/* uq/libs/mcmc/inc/uqFinDimLinearSpace.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_VECTOR_SPACE_H__
#define __UQ_VECTOR_SPACE_H__

#include <uqEnvironment.h>
#include <uqMiscellaneous.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <uqDefines.h>

template <class V, class M>
class uqVectorSpaceClass
{
public:
          uqVectorSpaceClass();
          uqVectorSpaceClass(const uqEnvironmentClass& env, // See template specialization
                             const char*               prefix,
                                   unsigned int        dimValue = 0);
         ~uqVectorSpaceClass();

  const   Epetra_Map&               map            ()                         const;
          unsigned int              dim            ()                         const;
  const   std::string&              componentName  (unsigned int componentId) const;
  const   std::vector<std::string>& componentsNames()                         const;

  const   V&                        zeroVector     ()                         const;
          V*                        newVector      ()                         const; // See template specialization
          V*                        newVector      (const V& v)               const;
          M*                        newMatrix      ()                         const; // See template specialization
          M*                        newDiagMatrix  (const V& v)               const;
          M*                        newDiagMatrix  (double diagValue)         const; // See template specialization

          void                      printCompsNames(std::ostream& os, bool printHorizontally) const;
          void                      print          (std::ostream& os) const;

protected:
        //void                      constructMap         ();
          void                      defineMyOptions      (po::options_description& optionsDesc) const;
          void                      getMyOptionValues    (po::options_description& optionsDesc);

          void                      readComponentsNamesFromSpecFile(std::string& specFileName);
          void                      setComponentName     (unsigned int componentId,
                                                          const std::string& name);
  const   uqEnvironmentClass&       m_env;
          std::string               m_prefix;
          unsigned int              m_dim;
  const   Epetra_Map*               m_map;
          V*                        m_zeroVector;

          po::options_description*  m_optionsDesc;
          std::string               m_option_help;
          std::string               m_option_dim;
          std::string               m_option_specFile;

  mutable std::vector<std::string>* m_componentsNames; // FIXME: will need to be a parallel vector in case of a very large number of components
};

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(
  const uqEnvironmentClass& env,
  const char*               prefix,
        unsigned int        dimValue)
  :
  m_env            (env),
  m_prefix         ((std::string)(prefix) + "space_"),
  m_dim            (dimValue),
  m_map            (NULL),
  m_zeroVector     (NULL),
  m_optionsDesc    (new po::options_description("Vector space options") ),
  m_option_help    (m_prefix + "help"    ),
  m_option_dim     (m_prefix + "dim"     ),
  m_option_specFile(m_prefix + "specFile"),
  m_componentsNames(new std::vector<std::string>(0))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqVectorSpaceClass<V,M>::constructor()"
              << std::endl;
  }

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting values of options with prefix '" << m_prefix
                                   << "', state of uqVectorSpaceClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  m_map        = new Epetra_Map(m_dim,0,m_env.comm());
  m_zeroVector = new V(m_env,*m_map);

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqVectorSpaceClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqVectorSpaceClass<V,M>::~uqVectorSpaceClass()
{
  //std::cout << "Entering uqVectorSpaceClass<V,M>::destructor()"
  //          << std::endl;

  if (m_componentsNames != NULL) delete m_componentsNames;
  if (m_zeroVector      != NULL) delete m_zeroVector;
  if (m_map             != NULL) delete m_map;

  //std::cout << "Leaving uqVectorSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
void
uqVectorSpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                                   "produce help message for vector space"              )
    (m_option_dim.c_str(),      po::value<unsigned int>()->default_value(0),  "Space dimension"                                    )
    (m_option_specFile.c_str(), po::value<std::string >()->default_value(""), "File with the specification of all components names")
  ;

  return;
}

template <class V, class M>
void
uqVectorSpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_dim.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap[m_option_dim.c_str()].as<unsigned int>();
  }

  // Read vector space spec file only if 0 dimension was passed to constructor
  if (m_componentsNames->size() == 0) {
    std::string specFileName("");
    if (m_env.allOptionsMap().count(m_option_specFile.c_str())) {
      const po::variables_map& tmpMap = m_env.allOptionsMap();
      specFileName = tmpMap[m_option_specFile.c_str()].as<std::string>();
      readComponentsNamesFromSpecFile(specFileName);
    }
  }

  return;
}

template <class V, class M>
const Epetra_Map&
uqVectorSpaceClass<V,M>::map() const
{
  UQ_FATAL_TEST_MACRO(m_map == NULL,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::map()",
                      "m_map is still NULL");
  return *m_map;
}

template<class V, class M>
const V&
uqVectorSpaceClass<V,M>::zeroVector() const
{
  UQ_FATAL_TEST_MACRO(m_zeroVector == NULL,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::zeroVector()",
                      "m_zeroVector is still NULL");
  return *m_zeroVector;
}

template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
V*
uqVectorSpaceClass<V,M>::newVector(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new V(v);
}

template <class V, class M>
M*
uqVectorSpaceClass<V,M>::newDiagMatrix(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new M(v);
}

template <class V, class M>
void
uqVectorSpaceClass<V,M>::readComponentsNamesFromSpecFile(std::string& specFileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(specFileName.c_str());

  // Determine number of lines
  unsigned int numLines = std::count(std::istreambuf_iterator<char>(ifs),
                                     std::istreambuf_iterator<char>(),
                                     '\n');

  // Determine number of components
  int iRC;
  ifs.seekg(0,std::ios_base::beg);
  unsigned int lineId = 0;
  unsigned int numComponents = 0;
  std::string tempString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    iRC = uqMiscReadStringAndDoubleFromFile(ifs,tempString,NULL);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqVectorSpaceClass<V,M>::constructor()",
                        "failed reading during the determination of the number of components");
    //std::cout << "lineId = "           << lineId
    //          << ", numComponents = " << numComponents
    //          << ", tempString = "     << tempString
    //          << std::endl;
    if (tempString[0] != '#') numComponents++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_dim != numComponents) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of components (%d) in space spec file does not match dimension (%d) passed in the main input file",numComponents,m_dim);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqVectorSpaceClass<V,M>::constructor()",
                        errorExplanation);
  }

  std::cout << "Space spec file '"     << specFileName
            << "' has "                << numLines
            << " lines and specifies " << numComponents
            << " components."
            << std::endl;
  m_componentsNames->resize(numComponents,"");

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int compId = 0;
  std::string  compName("");
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //std::cout << "Beginning read of line (in space spec file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, compName, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqVectorSpaceClass<V,M>::constructor()",
                        "failed reading a component name during the components names reading loop");

    lineId++;
    if (compName[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }

    // Check 'compId' before setting one more component
    if (compId >= m_componentsNames->size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"compId (%d) got too large during reading of space spec file",compId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqVectorSpaceClass<V,M>::constructor()",
                          errorExplanation);
    }

    std::cout << "Just read, for compId = " << compId
              << ": compName = "            << compName
              << std::endl;
    setComponentName(compId,
                     compName);
    compId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(compId != m_componentsNames->size(),
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "the number of components just read is nonconsistent");
  return;
}

template <class V, class M>
void
uqVectorSpaceClass<V,M>::setComponentName(
  unsigned int       componentId,
  const std::string& name)
{
  UQ_FATAL_TEST_MACRO((componentId > m_componentsNames->size()),
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::setComponentName()",
                      "componentId is too big");

  (*m_componentsNames)[componentId] = name;

  return;
}

template <class V, class M>
const std::string&
uqVectorSpaceClass<V,M>::componentName(unsigned int componentId) const
{
  UQ_FATAL_TEST_MACRO((componentId > m_componentsNames->size()),
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::componentName()",
                      "componentId is too big");

  return (*m_componentsNames)[componentId];
}

template<class V, class M>
void
uqVectorSpaceClass<V,M>::printCompsNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->componentName(i) << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->componentName(i) << "'"
         << std::endl;
    }
  }

  return;
}

template <class V, class M>
void
uqVectorSpaceClass<V,M>::print(std::ostream& os) const
{
  os << m_option_dim << " = " << m_dim
     << std::endl;

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqVectorSpaceClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_VECTOR_SPACE_H__

