/* uq/libs/mcmc/inc/uqVectorRV.h
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

#ifndef __UQ_VECTOR_RV_H__
#define __UQ_VECTOR_RV_H__

#include <uqScalarRV.h>
#include <uqVectorSpace.h>
#include <uqVectorProbDensity.h>
#include <uqVectorRealizer.h>
#include <uqDefaultPrior.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorRVClass {
public:
  uqBaseVectorRVClass(const char*                    prefix,
                      const uqVectorSpaceClass<V,M>& imageSpace);
  virtual ~uqBaseVectorRVClass();

  const   uqEnvironmentClass&                env                       ()       const;
  const   uqVectorSpaceClass          <V,M>& imageSpace                ()       const;
  const   uqBaseVectorProbDensityClass<V,M>& probDensity               ()       const;
  const   uqBaseVectorRealizerClass   <V,M>& realizer                  ()       const;
  virtual void                               realization               (V& vec) const = 0;

          void                               setProbDensity            (const uqBaseVectorProbDensityClass<V,M>& probDensity);
          void                               setRealizer               (const uqBaseVectorRealizerClass   <V,M>& realizer   );

  const   uqBaseScalarRVClass&               component                 (unsigned int componentId) const;
  const   V&                                 minValues                 () const;
  const   V&                                 maxValues                 () const;
  const   V&                                 expectValues              () const;
  const   V&                                 stdDevValues              () const;
          void                               setComponent              (unsigned int componentId,
                                                                        double       minValue    = -INFINITY,
                                                                        double       maxValue    = INFINITY,
                                                                        double       expectValue = 0.,
                                                                        double       stdDevValue = INFINITY);

          bool                               outOfBounds               (const V& v) const;
          uqBaseVectorSequenceClass<V>&      chain                     ();

  virtual void                               print                     (std::ostream& os) const;

protected:
          void                               defineMyOptions           (po::options_description& optionsDesc) const;
          void                               getMyOptionValues         (po::options_description& optionsDesc);
          void                               readComponentsFromSpecFile(std::string& specFileName);
          void                               resetValues               ();
          void                               createMinValues           () const; // See template specialization
          void                               createMaxValues           () const; // See template specialization
          void                               createExpectValues        () const; // See template specialization
          void                               createStdDevValues        () const; // See template specialization

  const   uqEnvironmentClass&                m_env;
          std::string                        m_prefix;
  const   uqVectorSpaceClass          <V,M>& m_imageSpace;
  const   uqBaseVectorProbDensityClass<V,M>* m_probDensity;
  const   uqBaseVectorRealizerClass   <V,M>* m_realizer;

          po::options_description*           m_optionsDesc;
          std::string                        m_option_help;
          std::string                        m_option_specFile;

          std::vector<uqBaseScalarRVClass*>  m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
          uqBaseScalarRVClass                m_dummyComponent;
          mutable V*                         m_minValues;
          mutable V*                         m_maxValues;
          mutable V*                         m_expectValues;
          mutable V*                         m_stdDevValues;

          bool                               m_chainUse2;
          uqSequenceOfVectorsClass<V>        m_chain1;
          uqArrayOfSequencesClass<V>         m_chain2;
};

template<class V, class M>
uqBaseVectorRVClass<V,M>::uqBaseVectorRVClass(
  const char*                              prefix,
  const uqVectorSpaceClass          <V,M>& imageSpace)
  //const uqBaseVectorProbDensityClass<V,M>* probDensity,
  //const uqBaseVectorRealizerClass   <V,M>* realizer)
  :
  m_env            (imageSpace.env()),
  m_prefix         ((std::string)(prefix)+"rv_"),
  m_imageSpace     (imageSpace),
  m_probDensity    (NULL),//(probDensity),
  m_realizer       (NULL),//(realizer)
  m_optionsDesc    (new po::options_description("Vector random variable options")),
  m_option_help    (m_prefix + "help"),
  m_option_specFile(m_prefix + "specFile"),
  m_components     (m_imageSpace.dim(),NULL),
  m_dummyComponent (),
  m_minValues      (NULL),
  m_maxValues      (NULL),
  m_expectValues   (NULL),
  m_stdDevValues   (NULL),
  m_chainUse2      (false),
  m_chain1         (0,m_imageSpace.zeroVector()),
  m_chain2         (0,m_imageSpace.zeroVector())
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqBaseVectorRVClass<V,M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRVClass<V,M>::~uqBaseVectorRVClass()
{
  m_chain1.clear();
  m_chain2.clear();
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                                  "produce help message for vector random variable" )
    (m_option_specFile.c_str(), po::value<std::string>()->default_value(""), "File with the specification of all RV components")
  ;

  return;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRVClass<V,M>::getMyOptionValues()"
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  //if (m_env.allOptionsMap().count(m_option_chain_use2.c_str())) {
  //  m_chainUse2 = m_env.allOptionsMap()[m_option_chain_use2.c_str()].as<bool>();
  //}

  // Read RV components specification file only if 0 dimension was passed to constructor
  //if (m_components.size() == 0) { GAMBIARRA
    std::string specFileName("");
    if (m_env.allOptionsMap().count(m_option_specFile.c_str())) {
      const po::variables_map& tmpMap = m_env.allOptionsMap();
      specFileName = tmpMap[m_option_specFile.c_str()].as<std::string>();
      if (specFileName == ".") {
        if (m_env.rank() == 0) {
          std::cout << "In uqBaseVectorRVClass<V,M>::getMyOptionValues()"
                    << ": spec file '" << specFileName << "' is interpreted as 'no spec was speficied'"
                    << std::endl;
        }
      }
      else {
        readComponentsFromSpecFile(specFileName);
      }
    }
    else {
      if (m_env.rank() == 0) {
        std::cout << "In uqBaseVectorRVClass<V,M>::getMyOptionValues()"
                  << ": no spec file was specified"
                  << std::endl;
      }
    }
    //} // GAMBIARRA

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRVClass<V,M>::getMyOptionValues()"
              << std::endl;
  }

  return;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::readComponentsFromSpecFile(std::string& specFileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(specFileName.c_str());
  if (ifs.is_open() == false) {
    if (m_env.rank() == 0) {
      std::cout << "In uqBaseVectorRVClass<V,M>::readComponentsFromSpecFile()"
                << ", WARNING: spec file '" << specFileName
                << "' was not found"
                << std::endl;
    }
    return;
  }
  else {
    if (m_env.rank() == 0) {
      std::cout << "In uqBaseVectorRVClass<V,M>::readComponentsFromSpecFile()"
                << ": about to read file '" << specFileName
                << "'"
                << std::endl;
    }
  }

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
                        "uqBaseVectorRVClass<V,M>::constructor()",
                        "failed reading during the determination of the number of components");
    //std::cout << "lineId = "           << lineId
    //          << ", numObservables = " << numObservables
    //          << ", tempString = "     << tempString
    //          << std::endl;
    if (tempString[0] != '#') numComponents++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_imageSpace.dim() != numComponents) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of components (%d) in RV components specification file '%s' does not match dimension (%d) of the image space",numComponents,specFileName.c_str(),m_imageSpace.dim());
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqBaseVectorRVClass<V,M>::constructor()",
                        errorExplanation);
  }

  std::cout << "RV components specification file '" << specFileName
            << "' has "                             << numLines
            << " lines and specifies "              << numComponents
            << " components."
            << std::endl;

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int componentId = 0;
  std::string  componentExplan  ("");
  std::string  minValueString   ("");
  std::string  maxValueString   ("");
  std::string  expectValueString("");
  std::string  stdDevValueString("");
  double       minValue;
  double       maxValue;
  double       expectValue;
  double       stdDevValue;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //std::cout << "Beginning read of line (in RV components specification file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, componentExplan, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqBaseVectorRVClass<V,M>::constructor()",
                        "failed reading a component explanation during the components reading loop");

    lineId++;
    if (componentExplan[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }
    UQ_FATAL_TEST_MACRO(endOfLineAchieved,
                        m_env.rank(),
                        "uqBaseVectorRVClass<V,M>::constructor()",
                        "failed to provide information beyond component explanation during the components reading loop");

    // Check 'componentId' before setting one more component
    if (componentId >= m_components.size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"componentId (%d) got too large during reading of RV components specification file",componentId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqBaseVectorRVClass<V,M>::constructor()",
                          errorExplanation);
    }

    minValue    = -INFINITY;
    maxValue    = INFINITY;
    expectValue = nan("");
    stdDevValue = INFINITY;

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, minValueString, &minValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorRVClass<V,M>::constructor()",
                          "failed reading a minimal value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, maxValueString, &maxValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorRVClass<V,M>::constructor()",
                          "failed reading a maximum value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, expectValueString, &expectValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorRVClass<V,M>::constructor()",
                          "failed reading an expectation value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, stdDevValueString, &stdDevValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorRVClass<V,M>::constructor()",
                          "failed reading a std dev value during the components reading loop");
    }

    if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
    if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
      std::cout << "In uqBaseVectorRVClass<V,M>::constructor()"
                << ", just read, for componentId = " << componentId
                << ": componentExplan = "            << componentExplan
                << ", minValue = "                   << minValue
                << ", maxValue = "                   << maxValue
                << ", expectValue = "                << expectValue
                << ", stdDevValue = "                << stdDevValue
                << std::endl;
    }
    setComponent(componentId,
                 minValue,
                 maxValue,
                 expectValue,
                 stdDevValue);
    componentId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(componentId != m_components.size(),
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::constructor()",
                      "the number of components just read is nonconsistent");

  return;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::setComponent(
  unsigned int componentId,
  double       minValue,
  double       maxValue,
  double       expectValue,
  double       stdDevValue)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRVClass<V,M>::setComponent()"
              << ", componentId = " << componentId
              << ", minValue = "    << minValue
              << ", maxValue = "    << maxValue
              << ", expectValue = " << expectValue
              << ", stdDevValue = " << stdDevValue
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((componentId > m_components.size()),
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::setComponent()",
                      "componentId is too big");

  if (m_components[componentId] == NULL) {
    m_components[componentId] = new uqBaseScalarRVClass(minValue,
                                                        maxValue,
                                                        expectValue,
                                                        stdDevValue);
  }
  else {
    m_components[componentId]->setMinValue   (minValue);
    m_components[componentId]->setMaxValue   (maxValue);
    m_components[componentId]->setExpectValue(expectValue);
    m_components[componentId]->setStdDevValue(stdDevValue);
  }

  // These values cannot be trusted anymore
  // They need to be updated
  // They will be updated the next time they are requested
  resetValues();

  return;
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorRVClass<V,M>::imageSpace() const
{
  return m_imageSpace;
}

template<class V, class M>
const uqBaseVectorProbDensityClass<V,M>&
uqBaseVectorRVClass<V,M>::probDensity() const
{
  UQ_FATAL_TEST_MACRO(m_probDensity == NULL,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::probDensity()",
                      "m_probDensity is NULL");

  return *m_probDensity;
}

template<class V, class M>
const uqBaseVectorRealizerClass<V,M>&
uqBaseVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::realizer()",
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V, class M>
void
uqBaseVectorRVClass<V,M>::setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity)
{
  m_probDensity = &probDensity;
  return;
}

template<class V, class M>
void
uqBaseVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::resetValues()
{
  if (m_stdDevValues) delete m_stdDevValues;
  if (m_expectValues) delete m_expectValues;
  if (m_maxValues   ) delete m_maxValues;
  if (m_minValues   ) delete m_minValues;
  m_stdDevValues = NULL;
  m_expectValues = NULL;
  m_maxValues    = NULL;
  m_minValues    = NULL;
}

template <class V, class M>
const uqBaseScalarRVClass&
uqBaseVectorRVClass<V,M>::component(unsigned int componentId) const
{
  if (componentId > m_components.size()) return m_dummyComponent;
  if (m_components[componentId] == NULL) return m_dummyComponent;
  return *(m_components[componentId]);
}

template <class V, class M>
const V&
uqBaseVectorRVClass<V,M>::minValues() const
{
  if (m_minValues == NULL) this->createMinValues();
  return *m_minValues;
}

template <class V, class M>
const V&
uqBaseVectorRVClass<V,M>::maxValues() const
{
  if (m_maxValues == NULL) this->createMaxValues();
  return *m_maxValues;
}

template <class V, class M>
const V&
uqBaseVectorRVClass<V,M>::expectValues() const
{
  if (m_expectValues == NULL) this->createExpectValues();
  return *m_expectValues;
}

template <class V, class M>
const V&
uqBaseVectorRVClass<V,M>::stdDevValues() const
{
  if (m_stdDevValues == NULL) this->createStdDevValues();
  return *m_stdDevValues;
}

template <class V, class M>
bool
uqBaseVectorRVClass<V,M>::outOfBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->minValues()) ||
          v.atLeastOneComponentBiggerThan (this->maxValues()));
}
#if 0
template<class V, class M>
void
uqBaseVectorRVClass<V,M>::realization(V& vec) const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::realization()",
                      "m_realizer is NULL");

  m_realizer->nextSample(vec);

  return;
}
#endif

template<class V,class M>
uqBaseVectorSequenceClass<V>&
uqBaseVectorRVClass<V,M>::chain()
{
  if (m_chainUse2) return m_chain2;
  return m_chain1;
}

template <class V, class M>
const uqEnvironmentClass&
uqBaseVectorRVClass<V,M>::env() const
{
  return m_env;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::print(std::ostream& os) const
{
  os << "\nComponentsx are:"
     << std::endl;
  for (unsigned int i = 0; i < m_components.size(); ++i) {
    os << i << " ";
    if (m_components[i]) {
      os << *(m_components[i]);
    }
    else {
      os << "NULL";
    }
    os << std::endl;
  }
  return;
}

template<class V, class M>
std::ostream& operator<<(std::ostream& os, const uqBaseVectorRVClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

//*****************************************************
// Generic class
//*****************************************************
template<class V, class M>
class uqGenericVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGenericVectorRVClass(const char*                              prefix,
                         const uqVectorSpaceClass          <V,M>& imageSpace,
                         const uqBaseVectorProbDensityClass<V,M>* probDensity,
                         const uqBaseVectorRealizerClass   <V,M>* realizer);
  virtual ~uqGenericVectorRVClass();

private:
  void realization(V& vec) const;

  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_probDensity;
  using uqBaseVectorRVClass<V,M>::m_realizer;
};

template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                              prefix,
  const uqVectorSpaceClass          <V,M>& imageSpace,
  const uqBaseVectorProbDensityClass<V,M>* probDensity,
  const uqBaseVectorRealizerClass   <V,M>* realizer)
  :
  uqBaseVectorRVClass<V,M>(prefix,imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGenericVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_probDensity = probDensity;
  m_realizer    = realizer;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGenericVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGenericVectorRVClass<V,M>::~uqGenericVectorRVClass()
{
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::realization(V& vec) const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.rank(),
                      "uqGenericVectorRVClass<V,M>::realization()",
                      "m_realizer is NULL");

  m_realizer->nextSample(vec);

  return;
}

//*****************************************************
// Gaussian class
//*****************************************************
template<class V, class M>
class uqGaussianVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGaussianVectorRVClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& imageSpace,
                          const M*                       covMatrix);
  virtual ~uqGaussianVectorRVClass();

private:
  void realization(V& vec) const;

  const M* m_covMatrix;
  V*                                      m_paramPriorMus;
  V*                                      m_paramPriorSigmas;
  uqDefault_M2lPriorRoutine_DataType<V,M> m_m2lPriorRoutine_Data;

  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_probDensity;
  using uqBaseVectorRVClass<V,M>::m_realizer;
};

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  const M*                       covMatrix)
  :
  uqBaseVectorRVClass<V,M>(prefix,imageSpace),
  m_covMatrix             (covMatrix)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if (m_covMatrix == NULL) {
#if 0
    V tmpVec(m_imageSpace.zeroVector());
    for (unsigned int i = 0; i < m_imageSpace.dim(); ++i) {
      sigma = m_components[i]->stdDevValue();
      tmpVec[i] = sigma*sigma;
    }
    m_covMatrix = m_imageSpace.newDiagMatrix(tmpVec);
#endif

    m_paramPriorMus    = new V(this->expectValues());
    m_paramPriorSigmas = new V(this->stdDevValues());
    m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
    m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;
    m_probDensity = new uqRoutineVectorProbDensityClass<V,M>(m_imageSpace,
                                                             uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                             (void *) &m_m2lPriorRoutine_Data,
                                                             true); // the routine computes [-2.*ln(Likelihood)]
    if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
      std::cout << "In uqGaussianVectorRVClass<V,M>::constructor()"
                << ", prefix = " << m_prefix
                << ": priorMus = " << *m_paramPriorMus
                << ", priorSigmas = " << *m_paramPriorSigmas
                << std::endl;
    }
  }
  else {
    //m_probDensity = uqGaussianProbDensity(env,m_covMatrix);
  }

  m_realizer = NULL;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGaussianVectorRVClass<V,M>::~uqGaussianVectorRVClass()
{
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::realization(V& vec) const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::realization()",
                      "m_realizer is NULL");

  m_realizer->nextSample(vec);

  return;
}

#endif // __UQ_VECTOR_RV_H__
