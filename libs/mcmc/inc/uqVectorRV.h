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

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqVectorRVClass {
public:
  uqVectorRVClass(const uqEnvironmentClass&                env,
                  const char*                              prefix,
                  const uqVectorSpaceClass          <V,M>& imageSpace,
                  const uqBaseVectorProbDensityClass<V,M>* probDensity,
                  const uqBaseVectorRealizerClass   <V,M>* realizer);
  virtual ~uqVectorRVClass();

  const   uqVectorSpaceClass          <V,M>& imageSpace                ()       const;
  const   uqBaseVectorProbDensityClass<V,M>& probDensity               ()       const;
  const   uqBaseVectorRealizerClass   <V,M>& realizer                  ()       const;
          void                               realization               (V& vec) const;

          void                               setProbDensity            (const uqBaseVectorProbDensityClass<V,M>& probDensity);
          void                               setRealizer               (const uqBaseVectorRealizerClass   <V,M>& realizer   );

  const   uqBaseScalarRVClass&               component                 (unsigned int componentId) const;
  const   V&                                 minValues                 () const;
  const   V&                                 maxValues                 () const;
  const   V&                                 expectValues              () const;
  const   V&                                 stdDevValues              () const;

          bool                               outOfBounds               (const V& v) const;

  virtual void                               print                     (std::ostream& os) const;

protected:
          void                               defineMyOptions           (po::options_description& optionsDesc) const;
          void                               getMyOptionValues         (po::options_description& optionsDesc);
          void                               readComponentsFromSpecFile(std::string& specFileName);
          void                               setComponent              (unsigned int componentId,
                                                                        double       minValue    = -INFINITY,
                                                                        double       maxValue    = INFINITY,
                                                                        double       expectValue = 0.,
                                                                        double       stdDevValue = INFINITY);
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
          std::string                        m_option_specificationFile;

          std::vector<uqBaseScalarRVClass*>  m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
          uqBaseScalarRVClass                m_dummyComponent;
          mutable V*                         m_minValues;
          mutable V*                         m_maxValues;
          mutable V*                         m_expectValues;
          mutable V*                         m_stdDevValues;
};

template<class V, class M>
uqVectorRVClass<V,M>::uqVectorRVClass(
  const uqEnvironmentClass&                env,
  const char*                              prefix,
  const uqVectorSpaceClass          <V,M>& imageSpace,
  const uqBaseVectorProbDensityClass<V,M>* probDensity,
  const uqBaseVectorRealizerClass   <V,M>* realizer)
  :
  m_env                     (env),
  m_prefix                  (prefix),
  m_imageSpace              (imageSpace),
  m_probDensity             (probDensity),
  m_realizer                (realizer),
  m_optionsDesc             (new po::options_description("Vector random variable options")),
  m_option_help             (m_prefix + "help"),
  m_option_specificationFile(m_prefix + "specificationFile"),
  m_components              (0),//,NULL),
  m_dummyComponent          (".",0.),
  m_minValues               (NULL),
  m_maxValues               (NULL),
  m_expectValues            (NULL),
  m_stdDevValues            (NULL)
{
  if (m_probDensity == NULL) {
#if 0
    m_paramPriorMus    = new P_V(m_paramSpace->expectValues());
    m_paramPriorSigmas = new P_V(m_paramSpace->stdDevValues());
    m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
    m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;

    m_priorParamDensity = new uqRoutineProbDensity_Class<P_V,P_M>(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                  (void *) &m_m2lPriorRoutine_Data,
                                                                  true); // the routine computes [-2.*ln(Likelihood)]
#endif
  }

}

template<class V, class M>
uqVectorRVClass<V,M>::~uqVectorRVClass()
{
}

template <class V, class M>
void
uqVectorRVClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                                           "produce help message for vector random variable"             )
    (m_option_specificationFile.c_str(), po::value<std::string>()->default_value(""), "File with the specification of all components to be inferred")
  ;

  return;
}

template <class V, class M>
void
uqVectorRVClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqVectorRVClass<V,M>::getMyOptionValues()"
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  // Read component specification file only if 0 dimension was passed to constructor
  if (m_components.size() == 0) {
    std::string specFileName("");
    if (m_env.allOptionsMap().count(m_option_specificationFile.c_str())) {
      const po::variables_map& tmpMap = m_env.allOptionsMap();
      specFileName = tmpMap[m_option_specificationFile.c_str()].as<std::string>();
      readComponentsFromSpecFile(specFileName);
    }
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqVectorRVClass<V,M>::getMyOptionValues()"
              << std::endl;
  }

  return;
}

template <class V, class M>
void
uqVectorRVClass<V,M>::readComponentsFromSpecFile(std::string& specFileName)
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
                        "uqVectorRVClass<V,M>::constructor()",
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
                      "uqVectorRVClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_imageSpace.dim() != numComponents) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of components (%d) in component specification file does not match dimension (%d) of the image space",numComponents,m_imageSpace.dim());
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqVectorRVClass<V,M>::constructor()",
                        errorExplanation);
  }

  std::cout << "Component specification file '" << specFileName
            << "' has "                         << numLines
            << " lines and specifies "          << numComponents
            << " components."
            << std::endl;
  m_components.resize(numComponents,NULL);

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
    //std::cout << "Beginning read of line (in component specification file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, componentExplan, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqVectorRVClass<V,M>::constructor()",
                        "failed reading a component explanation during the components reading loop");

    lineId++;
    if (componentExplan[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }
    UQ_FATAL_TEST_MACRO(endOfLineAchieved,
                        m_env.rank(),
                        "uqVectorRVClass<V,M>::constructor()",
                        "failed to provide information beyond component explanation during the components reading loop");

    // Check 'componentId' before setting one more component
    if (componentId >= m_components.size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"componentId (%d) got too large during reading of component specification file",componentId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqVectorRVClass<V,M>::constructor()",
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
                          "uqVectorRVClass<V,M>::constructor()",
                          "failed reading a minimal value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, maxValueString, &maxValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqVectorRVClass<V,M>::constructor()",
                          "failed reading a maximum value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, expectValueString, &expectValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqVectorRVClass<V,M>::constructor()",
                          "failed reading an expectation value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, stdDevValueString, &stdDevValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqVectorRVClass<V,M>::constructor()",
                          "failed reading a std dev value during the components reading loop");
    }

    if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
    //std::cout << "Just read, for componentId = " << componentId
    //          << ": componentExplan = " << componentExplan
    //          << ", minValue = "        << minValue
    //          << ", maxValue = "        << maxValue
    //          << ", expectValue = "     << expectValue
    //          << ", stdDevValue = "     << stdDevValue
    //          << std::endl;
    setComponent(componentId,
                 minValue,
                 maxValue,
                 expectValue,
                 stdDevValue);
    componentId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqVectorRVClass<V,M>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(componentId != m_components.size(),
                      m_env.rank(),
                      "uqVectorRVClass<V,M>::constructor()",
                      "the number of components just read is nonconsistent");

  return;
}

template <class V, class M>
void
uqVectorRVClass<V,M>::setComponent(
  unsigned int componentId,
  double       minValue,
  double       maxValue,
  double       expectValue,
  double       stdDevValue)
{
  UQ_FATAL_TEST_MACRO((componentId > m_components.size()),
                      m_env.rank(),
                      "uqVectorRVClass<V,M>::setComponent()",
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
uqVectorRVClass<V,M>::imageSpace() const
{
  return m_imageSpace;
}

template<class V, class M>
const uqBaseVectorProbDensityClass<V,M>&
uqVectorRVClass<V,M>::probDensity() const
{
  UQ_FATAL_TEST_MACRO(m_probDensity == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqVectorRVClass<V,M>::probDensity()",
                      "m_probDensity is NULL");

  return *m_probDensity;
}

template<class V, class M>
const uqBaseVectorRealizerClass<V,M>&
uqVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqVectorRVClass<V,M>::realizer()",
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V, class M>
void
uqVectorRVClass<V,M>::setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity)
{
  m_probDensity = &probDensity;
  return;
}

template<class V, class M>
void
uqVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}

template <class V, class M>
void
uqVectorRVClass<V,M>::resetValues()
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
uqVectorRVClass<V,M>::component(unsigned int componentId) const
{
  if (componentId > m_components.size()) return m_dummyComponent;
  if (m_components[componentId] == NULL) return m_dummyComponent;
  return *(m_components[componentId]);
}

template <class V, class M>
const V&
uqVectorRVClass<V,M>::minValues() const
{
  if (m_minValues == NULL) this->createMinValues();
  return *m_minValues;
}

template <class V, class M>
const V&
uqVectorRVClass<V,M>::maxValues() const
{
  if (m_maxValues == NULL) this->createMaxValues();
  return *m_maxValues;
}

template <class V, class M>
const V&
uqVectorRVClass<V,M>::expectValues() const
{
  if (m_expectValues == NULL) this->createExpectValues();
  return *m_expectValues;
}

template <class V, class M>
const V&
uqVectorRVClass<V,M>::stdDevValues() const
{
  if (m_stdDevValues == NULL) this->createStdDevValues();
  return *m_stdDevValues;
}

template <class V, class M>
bool
uqVectorRVClass<V,M>::outOfBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->minValues()) ||
          v.atLeastOneComponentBiggerThan (this->maxValues()));
}

template<class V, class M>
void
uqVectorRVClass<V,M>::realization(V& vec) const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqVectorRVClass<V,M>::realization()",
                      "m_realizer is NULL");

  m_realizer->nextSample(vec);

  return;
}

template <class V, class M>
void
uqVectorRVClass<V,M>::print(std::ostream& os) const
{
  return;
}

template<class V, class M>
std::ostream& operator<<(std::ostream& os, const uqVectorRVClass<V,M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_VECTOR_RV_H__
