/* uq/libs/mcmc/inc/uqVectorProbDensity.h
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

#ifndef __UQ_VECTOR_PROB_DENSITY_H__
#define __UQ_VECTOR_PROB_DENSITY_H__

#include <uqEnvironment.h>
#include <math.h>
#include <uqDefaultPrior.h>

//*****************************************************
// Classes to accomodate a probability density.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorProbDensityClass {
public:
           uqBaseVectorProbDensityClass(const char*                    prefix,
                                        const uqVectorSpaceClass<V,M>& domainSpace,
                                        const V*                       minValues,
                                        const V*                       maxValues);
  virtual ~uqBaseVectorProbDensityClass();

  const   uqVectorSpaceClass<V,M>&              domainSpace               ()                     const;
  virtual double                                actualDensity             (const V& paramValues) const = 0;
  virtual double                                minus2LnDensity           (const V& paramValues) const = 0;

  //const   uqBaseScalarProbDensityClass<double>& component                 (unsigned int componentId) const;
  const   V&                                    minValues                 () const;
  const   V&                                    maxValues                 () const;
          bool                                  outOfBounds               (const V& v) const;

#if 0
  const   V&                                    expectValues              () const;
  const   V&                                    stdDevValues              () const;
          void                                  setComponent              (unsigned int componentId,
                                                                           double       minValue    = -INFINITY,
                                                                           double       maxValue    = INFINITY,
                                                                           double       expectValue = 0.,
                                                                           double       stdDevValue = INFINITY);
#endif

protected:
#if 0
          void                                  defineMyOptions           (po::options_description& optionsDesc) const;
          void                                  getMyOptionValues         (po::options_description& optionsDesc);
          void                                  readComponentsFromSpecFile(std::string& specFileName);
#endif
          void                                  resetValues               ();
          //void                                  createMinValues           () const; // See template specialization
          //void                                  createMaxValues           () const; // See template specialization
          //void                                  createExpectValues        () const; // See template specialization
          //void                                  createStdDevValues        () const; // See template specialization

  const   uqEnvironmentClass&                                m_env;
          std::string                                        m_prefix;
  const   uqVectorSpaceClass<V,M>&                           m_domainSpace;

  //po::options_description*                           m_optionsDesc;
  //std::string                                        m_option_help;
  //std::string                                        m_option_specFile;

  //std::vector<uqBaseScalarProbDensityClass<double>*> m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
  //uqBaseScalarProbDensityClass<double>               m_dummyComponent;
          mutable const V*                                         m_minValues;
          mutable const V*                                         m_maxValues;
          //mutable V*                                         m_expectValues;
          //mutable V*                                         m_stdDevValues;
};

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::uqBaseVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V*                       minValues,
  const V*                       maxValues)
  :
  m_env        (domainSpace.env()),
  m_prefix     ((std::string)(prefix)+"pd_"),
  m_domainSpace(domainSpace),
  m_minValues  (minValues),
  m_maxValues  (maxValues)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  //defineMyOptions                (*m_optionsDesc);
  //m_env.scanInputFileForMyOptions(*m_optionsDesc);
  //getMyOptionValues              (*m_optionsDesc);

  //if (m_env.rank() == 0) std::cout << "In uqBaseVectorProbDensityClass<V,M>::constructor()"
  //                                 << ": after getting values of options, state of object is:"
  //                                 << "\n" << *this
  //                                 << std::endl;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorProbDensityClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::~uqBaseVectorProbDensityClass()
{
}
#if 0
template <class V, class M>
void
uqBaseVectorProbDensityClass<V,M>::defineMyOptions(
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
uqBaseVectorProbDensityClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::getMyOptionValues()"
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  //if (m_env.allOptionsMap().count(m_option_chain_use2.c_str())) {
  //  m_chainUse2 = m_env.allOptionsMap()[m_option_chain_use2.c_str()].as<bool>();
  //}

  // Read RV components specification file
  std::string specFileName("");
  if (m_env.allOptionsMap().count(m_option_specFile.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    specFileName = tmpMap[m_option_specFile.c_str()].as<std::string>();
    if (specFileName == ".") {
      if (m_env.rank() == 0) {
        std::cout << "In uqBaseVectorProbDensityClass<V,M>::getMyOptionValues()"
                  << ": an spec file with name '" << specFileName << "' is interpreted as 'no spec file was speficied'"
                  << std::endl;
      }
    }
    else {
      readComponentsFromSpecFile(specFileName);
    }
  }
  else {
    if (m_env.rank() == 0) {
      std::cout << "In uqBaseVectorProbDensityClass<V,M>::getMyOptionValues()"
                << ": no spec file was specified"
                << std::endl;
    }
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorProbDensityClass<V,M>::getMyOptionValues()"
              << std::endl;
  }

  return;
}

template <class V, class M>
void
uqBaseVectorProbDensityClass<V,M>::readComponentsFromSpecFile(std::string& specFileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(specFileName.c_str());
  if (ifs.is_open() == false) {
    if (m_env.rank() == 0) {
      std::cout << "In uqBaseVectorProbDensityClass<V,M>::readComponentsFromSpecFile()"
                << ", WARNING: spec file '" << specFileName
                << "' was specified but it was not found"
                << std::endl;
    }
    return;
  }
  else {
    if (m_env.rank() == 0) {
      std::cout << "In uqBaseVectorProbDensityClass<V,M>::readComponentsFromSpecFile()"
                << ": about to read spec file '" << specFileName
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
                        "uqBaseVectorProbDensityClass<V,M>::constructor()",
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
                      "uqBaseVectorProbDensityClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_imageSpace.dim() != numComponents) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of components (%d) in RV components specification file '%s' does not match dimension (%d) of the image space",numComponents,specFileName.c_str(),m_imageSpace.dim());
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqBaseVectorProbDensityClass<V,M>::constructor()",
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
                        "uqBaseVectorProbDensityClass<V,M>::constructor()",
                        "failed reading a component explanation during the components reading loop");

    lineId++;
    if (componentExplan[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }
    UQ_FATAL_TEST_MACRO(endOfLineAchieved,
                        m_env.rank(),
                        "uqBaseVectorProbDensityClass<V,M>::constructor()",
                        "failed to provide information beyond component explanation during the components reading loop");

    // Check 'componentId' before setting one more component
    if (componentId >= m_components.size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"componentId (%d) got too large during reading of RV components specification file",componentId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqBaseVectorProbDensityClass<V,M>::constructor()",
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
                          "uqBaseVectorProbDensityClass<V,M>::constructor()",
                          "failed reading a minimal value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, maxValueString, &maxValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorProbDensityClass<V,M>::constructor()",
                          "failed reading a maximum value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, expectValueString, &expectValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorProbDensityClass<V,M>::constructor()",
                          "failed reading an expectation value during the components reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, stdDevValueString, &stdDevValue, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqBaseVectorProbDensityClass<V,M>::constructor()",
                          "failed reading a std dev value during the components reading loop");
    }

    if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
    if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
      std::cout << "In uqBaseVectorProbDensityClass<V,M>::constructor()"
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
                      "uqBaseVectorProbDensityClass<V,M>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(componentId != m_components.size(),
                      m_env.rank(),
                      "uqBaseVectorProbDensityClass<V,M>::constructor()",
                      "the number of components just read is nonconsistent");

  return;
}

template <class V, class M>
void
uqBaseVectorProbDensityClass<V,M>::setComponent(
  unsigned int componentId,
  double       minValue,
  double       maxValue,
  double       expectValue,
  double       stdDevValue)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::setComponent()"
              << ", componentId = " << componentId
              << ", minValue = "    << minValue
              << ", maxValue = "    << maxValue
              << ", expectValue = " << expectValue
              << ", stdDevValue = " << stdDevValue
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((componentId > m_components.size()),
                      m_env.rank(),
                      "uqBaseVectorProbDensityClass<V,M>::setComponent()",
                      "componentId is too big");

  if (m_components[componentId] == NULL) {
    m_components[componentId] = new uqBaseScalarProbDensityClass(minValue,
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

template <class V, class M>
void
uqBaseVectorProbDensityClass<V,M>::resetValues()
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
const uqBaseScalarProbDensityClass&
uqBaseVectorProbDensityClass<V,M>::component(unsigned int componentId) const
{
  if (componentId > m_components.size()) return m_dummyComponent;
  if (m_components[componentId] == NULL) return m_dummyComponent;
  return *(m_components[componentId]);
}
#endif

template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::minValues() const
{
  //if (m_minValues == NULL) this->createMinValues();
  return *m_minValues;
}

template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::maxValues() const
{
  //if (m_maxValues == NULL) this->createMaxValues();
  return *m_maxValues;
}
#if 0
template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::expectValues() const
{
  if (m_expectValues == NULL) this->createExpectValues();
  return *m_expectValues;
}

template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::stdDevValues() const
{
  if (m_stdDevValues == NULL) this->createStdDevValues();
  return *m_stdDevValues;
}
#endif
template <class V, class M>
bool
uqBaseVectorProbDensityClass<V,M>::outOfBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->minValues()) ||
          v.atLeastOneComponentBiggerThan (this->maxValues()));
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorProbDensityClass<V,M>::domainSpace() const
{
  return m_domainSpace;
}

//*****************************************************
// Generic probability density class
//*****************************************************
template<class V, class M>
class uqGenericVectorProbDensityClass : public uqBaseVectorProbDensityClass<V,M> {
public:
  uqGenericVectorProbDensityClass(const char*                    prefix,
                                  const uqVectorSpaceClass<V,M>& domainSpace,
                                  const V*                       minValues,
                                  const V*                       maxValues,
                                  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                                  const void* routineDataPtr,
                                  bool routineComputesMinus2LogOfDensity);
 ~uqGenericVectorProbDensityClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;

  bool m_routineComputesMinus2LogOfDensity;

  using uqBaseVectorProbDensityClass<V,M>::m_env;
  using uqBaseVectorProbDensityClass<V,M>::m_prefix;
  using uqBaseVectorProbDensityClass<V,M>::m_domainSpace;
};

template<class V, class M>
uqGenericVectorProbDensityClass<V,M>::uqGenericVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V*                       minValues,
  const V*                       maxValues,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr,
  bool        routineComputesMinus2LogOfDensity)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace,minValues,maxValues),
  m_routinePtr                       (routinePtr),
  m_routineDataPtr                   (routineDataPtr),
  m_routineComputesMinus2LogOfDensity(routineComputesMinus2LogOfDensity)
{
}

template<class V, class M>
uqGenericVectorProbDensityClass<V,M>::~uqGenericVectorProbDensityClass()
{
}

template<class V, class M>
double
uqGenericVectorProbDensityClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfDensity == false) {
    value = -2.*log(value);
  }

  return value;
}

template<class V, class M>
double
uqGenericVectorProbDensityClass<V,M>::actualDensity(const V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfDensity) {
    value = exp(-.5*value);
  }

  return value;
}

//*****************************************************
// Bayesian probability density class
//*****************************************************
template<class V, class M>
class uqBayesianVectorProbDensityClass : public uqBaseVectorProbDensityClass<V,M> {
public:
  uqBayesianVectorProbDensityClass(const char*                              prefix,
                                   const uqBaseVectorProbDensityClass<V,M>* priorDensity,
                                   const uqBaseVectorProbDensityClass<V,M>* likelihoodFunction); 
 ~uqBayesianVectorProbDensityClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  const uqBaseVectorProbDensityClass<V,M>* m_priorDensity;
  const uqBaseVectorProbDensityClass<V,M>* m_likelihoodFunction;

  using uqBaseVectorProbDensityClass<V,M>::m_env;
  using uqBaseVectorProbDensityClass<V,M>::m_prefix;
  using uqBaseVectorProbDensityClass<V,M>::m_domainSpace;
};

template<class V,class M>
uqBayesianVectorProbDensityClass<V,M>::uqBayesianVectorProbDensityClass(
  const char*                              prefix,
  const uqBaseVectorProbDensityClass<V,M>* priorDensity,
  const uqBaseVectorProbDensityClass<V,M>* likelihoodFunction)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"bay").c_str(),priorDensity->domainSpace(),NULL,NULL),
  m_priorDensity      (priorDensity),
  m_likelihoodFunction(likelihoodFunction)
{
}

template<class V,class M>
uqBayesianVectorProbDensityClass<V,M>::~uqBayesianVectorProbDensityClass()
{
}

template<class V, class M>
double
uqBayesianVectorProbDensityClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value1 = m_priorDensity->minus2LnDensity(paramValues);
  double value2 = m_likelihoodFunction->minus2LnDensity(paramValues);

  //if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
  //  std::cout << "In uqBayesianVectorProbDensityClass<P_V,P_M>::minus2LnDensity()"
  //            << ", -2ln(prior) = " << value1
  //            << ", -2ln(like) = "  << value2
  //            << std::endl;
  //}

  return value1+value2;
}

template<class V, class M>
double
uqBayesianVectorProbDensityClass<V,M>::actualDensity(const V& paramValues) const
{
  double value1 = m_priorDensity->actualDensity(paramValues);
  double value2 = m_likelihoodFunction->actualDensity(paramValues);

  return value1*value2;
}

//*****************************************************
// Guassian probability density class
//*****************************************************
template<class V, class M>
class uqGaussianVectorProbDensityClass : public uqBaseVectorProbDensityClass<V,M> {
public:
  uqGaussianVectorProbDensityClass(const char*                    prefix,
                                   const uqVectorSpaceClass<V,M>& domainSpace,
                                   const V*                       minValues,
                                   const V*                       maxValues,
                                   const M*                       covMatrix,
                                   const V*                       expectValues,  // GAMBIARRA
                                   const V*                       stdDevValues); // GAMBIARRA
 ~uqGaussianVectorProbDensityClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  const M* m_covMatrix;
  V*                                      m_paramPriorMus;
  V*                                      m_paramPriorSigmas;
  uqDefault_M2lPriorRoutine_DataType<V,M> m_m2lPriorRoutine_Data;
  const uqBaseVectorProbDensityClass<V,M>* m_probDensity;

  using uqBaseVectorProbDensityClass<V,M>::m_env;
  using uqBaseVectorProbDensityClass<V,M>::m_prefix;
  using uqBaseVectorProbDensityClass<V,M>::m_domainSpace;
};

template<class V,class M>
uqGaussianVectorProbDensityClass<V,M>::uqGaussianVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V*                       minValues,
  const V*                       maxValues,
  const M*                       covMatrix,
  const V*                       expectValues,
  const V*                       stdDevValues)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,minValues,maxValues),
  m_covMatrix                      (covMatrix)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorProbDensityClass<V,M>::constructor()"
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

    m_paramPriorMus    = new V(*expectValues);
    m_paramPriorSigmas = new V(*stdDevValues);
    m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
    m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;
    m_probDensity = new uqGenericVectorProbDensityClass<V,M>(m_prefix.c_str(),
                                                             m_domainSpace,
                                                             minValues,
                                                             maxValues,
                                                             uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                             (void *) &m_m2lPriorRoutine_Data,
                                                             true); // the routine computes [-2.*ln(Likelihood)]
    if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
      std::cout << "In uqGaussianVectorProbDensityClass<V,M>::constructor()"
                << ", prefix = "      << m_prefix
                << ": priorMus = "    << *m_paramPriorMus
                << ", priorSigmas = " << *m_paramPriorSigmas
                << std::endl;
    }
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorProbDensityClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorProbDensityClass<V,M>::~uqGaussianVectorProbDensityClass()
{
  delete m_probDensity;
  delete m_paramPriorSigmas;
  delete m_paramPriorMus;
}

template<class V, class M>
double
uqGaussianVectorProbDensityClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  return m_probDensity->minus2LnDensity(paramValues);
}

template<class V, class M>
double
uqGaussianVectorProbDensityClass<V,M>::actualDensity(const V& paramValues) const
{
  return m_probDensity->actualDensity(paramValues);
}

#endif // __UQ_VECTOR_PROB_DENSITY_H__
