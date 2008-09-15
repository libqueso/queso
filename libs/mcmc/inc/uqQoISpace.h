/* uq/libs/mcmc/inc/uqQoISpace.h
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

#ifndef __UQ_QOI_SPACE_H__
#define __UQ_QOI_SPACE_H__

#include <uqFinDimLinearSpace.h>
#include <uqMiscellaneous.h>
#include <uqQoI.h>
#include <vector>

template <class V, class M>
class uqQoISpaceClass : public uqFinDimLinearSpaceClass<V,M>
{
public:
           uqQoISpaceClass();
           uqQoISpaceClass(const uqEnvironmentClass& env,
                           const char*               prefix);
  virtual ~uqQoISpaceClass();

          unsigned int      dim                 () const;
          int               setQoI              (unsigned int       qoiId,
                                                 const std::string& name);
          const uqQoIClass& qoi                 (unsigned int qoiId) const;
  virtual void              print               (std::ostream& os) const;
          void              printQoINames       (std::ostream& os, bool printHorizontally) const; // See template specialization

protected:
          void              defineMyOptions     (po::options_description& optionsDesc) const;
          void              getMyOptionValues   (po::options_description& optionsDesc);
          void              readQoIsFromSpecFile(std::string& specFileName);
          void              resetValues         ();

  po::options_description* m_optionsDesc;
  std::string              m_option_help;
  std::string              m_option_dim;
  std::string              m_option_specificationFile;

  std::vector<uqQoIClass*> m_qois; // FIXME: will need to be a parallel vector in case of a very large number of qois
  uqQoIClass               m_dummyQoI;

  using uqFinDimLinearSpaceClass<V,M>::m_env;
  using uqFinDimLinearSpaceClass<V,M>::m_dim;
};

template <class V, class M>
uqQoISpaceClass<V,M>::uqQoISpaceClass()
  :
  uqFinDimLinearSpaceClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqQoISpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqQoISpaceClass<V,M>::uqQoISpaceClass(
  const uqEnvironmentClass& env,
  const char*               prefix)
  :
  uqFinDimLinearSpaceClass<V,M>(env,prefix),
  m_optionsDesc                (new po::options_description("QoI space options")),
  m_option_help                (uqFinDimLinearSpaceClass<V,M>::m_prefix + "qoiSpace_help"             ),
  m_option_dim                 (uqFinDimLinearSpaceClass<V,M>::m_prefix + "qoiSpace_dim"              ),
  m_option_specificationFile   (uqFinDimLinearSpaceClass<V,M>::m_prefix + "qoiSpace_specificationFile"),
  m_qois                       (0),//,NULL),
  m_dummyQoI                   ("NonExistentYet")
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqQoISpaceClass<V,M>::constructor()"
              << std::endl;
  }


  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  // Now that 'm_dim' has been set, construct Trilinos map
  this->constructMap();

  if (m_env.rank() == 0) std::cout << "After getting values of options with prefix '" << uqFinDimLinearSpaceClass<V,M>::m_prefix
                                   << "', state of uqQoISpaceClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqQoISpaceClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqQoISpaceClass<V,M>::~uqQoISpaceClass()
{
  //std::cout << "Entering uqQoISpaceClass<V,M>::destructor()"
  //          << std::endl;

  for (unsigned int i = 0; i < m_qois.size(); ++i) {
    if (m_qois[i]) delete m_qois[i];
  }

  if (m_optionsDesc) delete m_optionsDesc;

  //std::cout << "Leaving uqQoISpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
unsigned int
uqQoISpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
void
uqQoISpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                                            "produce help message for UQ qoi space"                 )
    (m_option_dim.c_str(),               po::value<unsigned int>()->default_value(0),  "Space dimension"                                       )
    (m_option_specificationFile.c_str(), po::value<std::string >()->default_value(""), "File with the specification of all qois to be computed")
  ;

  return;
}

template <class V, class M>
void
uqQoISpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_dim.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap[m_option_dim.c_str()].as<unsigned int>();
  }

  // Read qoi specification file only if 0 dimension was passed to constructor
  if (m_qois.size() == 0) {
    std::string specFileName("");
    if (m_env.allOptionsMap().count(m_option_specificationFile.c_str())) {
      const po::variables_map& tmpMap = m_env.allOptionsMap();
      specFileName = tmpMap[m_option_specificationFile.c_str()].as<std::string>();
      readQoIsFromSpecFile(specFileName);
    }
  }

  return;
}

template <class V, class M>
void
uqQoISpaceClass<V,M>::readQoIsFromSpecFile(std::string& specFileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(specFileName.c_str());

  // Determine number of lines
  unsigned int numLines = std::count(std::istreambuf_iterator<char>(ifs),
                                     std::istreambuf_iterator<char>(),
                                     '\n');

  // Determine number of qois
  int iRC;
  ifs.seekg(0,std::ios_base::beg);
  unsigned int lineId = 0;
  unsigned int numQoIs = 0;
  std::string tempString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    iRC = uqMiscReadStringAndDoubleFromFile(ifs,tempString,NULL);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqQoISpaceClass<V,M>::constructor()",
                        "failed reading during the determination of the number of qois");
    //std::cout << "lineId = "           << lineId
    //          << ", numQoIs = " << numQoIs
    //          << ", tempString = "     << tempString
    //          << std::endl;
    if (tempString[0] != '#') numQoIs++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqQoISpaceClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_dim != numQoIs) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of qois (%d) in qoi specification file does not match dimension (%d) passed in the main input file",numQoIs,m_dim);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqQoISpaceClass<V,M>::constructor()",
                        errorExplanation);
  }

  std::cout << "QoI specification file '" << specFileName
            << "' has "                          << numLines
            << " lines and specifies "           << numQoIs
            << " qois."
            << std::endl;
  m_qois.resize(numQoIs,NULL);

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int qoiId = 0;
  std::string  qoiName            ("");
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //std::cout << "Beginning read of line (in qoi specification file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, qoiName, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqQoISpaceClass<V,M>::constructor()",
                        "failed reading a qoi name during the qois reading loop");

    lineId++;
    if (qoiName[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }
    //UQ_FATAL_TEST_MACRO(endOfLineAchieved,
    //                    m_env.rank(),
    //                    "uqQoISpaceClass<V,M>::constructor()",
    //                    "failed to provide information beyond qoi name during the qois reading loop");

    // Check 'qoiId' before setting one more qoi
    if (qoiId >= m_qois.size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"qoiId (%d) got too large during reading of qoi specification file",qoiId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqQoISpaceClass<V,M>::constructor()",
                          errorExplanation);
    }

    std::cout << "Just read, for qoiId = " << qoiId
              << ": qoiName = "            << qoiName
              << std::endl;
    setQoI(qoiId,
           qoiName);
    qoiId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqQoISpaceClass<V,M>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(qoiId != m_qois.size(),
                      m_env.rank(),
                      "uqQoISpaceClass<V,M>::constructor()",
                      "the number of qois just read is nonconsistent");
  return;
}

template <class V, class M>
int
uqQoISpaceClass<V,M>::setQoI(
  unsigned int       qoiId,
  const std::string& name)
{
  UQ_TEST_MACRO((qoiId > m_qois.size()),
                m_env.rank(),
                "uqQoISpaceClass<V,M>::setQoI()",
                "qoiId is too big",
                UQ_INVALID_QOI_SPEC_RC);

  if (m_qois[qoiId] == NULL) {
    m_qois[qoiId] = new uqQoIClass(name);
  }
  else {
    m_qois[qoiId]->setName(name);
  }

  // These values cannot be trusted anymore
  // They need to be updated
  // They will be updated the next time they are requested
  resetValues();

  return 0;
}

template <class V, class M>
const uqQoIClass&
uqQoISpaceClass<V,M>::qoi(unsigned int qoiId) const
{
  if (qoiId > m_qois.size()) return m_dummyQoI;
  if (m_qois[qoiId] == NULL) return m_dummyQoI;
  return *(m_qois[qoiId]);
}

template <class V, class M>
void
uqQoISpaceClass<V,M>::resetValues()
{
}

template<class V, class M>
void
uqQoISpaceClass<V,M>::printQoINames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->qoi(i).name() << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->qoi(i).name() << "'"
         << std::endl;
    }
  }

  return;
}

template <class V, class M>
void
uqQoISpaceClass<V,M>::print(std::ostream& os) const
{
  os << uqFinDimLinearSpaceClass<V,M>::m_prefix << "dim = " << m_dim
     << "\nQoIs are:"
     << std::endl;
  for (unsigned int i = 0; i < this->dim(); ++i) {
    os << i << " ";
    if (m_qois[i]) {
      os << *(m_qois[i]);
    }
    else {
      os << "NULL";
    }
    os << std::endl;
  }

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqQoISpaceClass<V,M>& space)
{
  space.print(os);

  return os;
}
#endif // __UQ_QOI_SPACE_H__

