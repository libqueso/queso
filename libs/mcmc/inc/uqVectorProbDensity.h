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
                                        const V&                       domainMinValues,
                                        const V&                       domainMaxValues,
                                        const V&                       domainExpectedValues,
                                        const V&                       domainVarianceValues);
           uqBaseVectorProbDensityClass(const char*                    prefix,
                                        const uqVectorSpaceClass<V,M>& domainSpace,
                                        const V&                       domainMinValues,
                                        const V&                       domainMaxValues,
                                        const V&                       domainExpectedValues);
           uqBaseVectorProbDensityClass(const char*                    prefix,
                                        const uqVectorSpaceClass<V,M>& domainSpace,
                                        const V&                       domainMinValues,
                                        const V&                       domainMaxValues);
           uqBaseVectorProbDensityClass(const char*                    prefix,
                                        const uqVectorSpaceClass<V,M>& domainSpace);
  virtual ~uqBaseVectorProbDensityClass();

  const   uqVectorSpaceClass<V,M>&              domainSpace               ()                     const;
  virtual double                                actualDensity             (const V& paramValues) const = 0;
  virtual double                                minus2LnDensity           (const V& paramValues) const = 0;

//const   uqBaseScalarProbDensityClass<double>& component                 (unsigned int componentId) const;
          bool                                  outOfDomainBounds         (const V& v) const;
  const   V&                                    domainMinValues           () const;
  const   V&                                    domainMaxValues           () const;
  const   V&                                    domainExpectedValues      () const;
  const   V&                                    domainVarianceValues      () const;

protected:

  const   uqEnvironmentClass&      m_env;
          std::string              m_prefix;
  const   uqVectorSpaceClass<V,M>& m_domainSpace;
          V*                       m_domainMinValues;
          V*                       m_domainMaxValues;
          V*                       m_domainExpectedValues;
          V*                       m_domainVarianceValues;

//std::vector<uqBaseScalarProbDensityClass<double>*> m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
//uqBaseScalarProbDensityClass<double>               m_dummyComponent;
};

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::uqBaseVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  m_env                 (domainSpace.env()),
  m_prefix              ((std::string)(prefix)+"pd_"),
  m_domainSpace         (domainSpace),
  m_domainMinValues     (new V(domainMinValues   )),
  m_domainMaxValues     (new V(domainMaxValues   )),
  m_domainExpectedValues(new V(domainExpectedValues)),
  m_domainVarianceValues(new V(domainVarianceValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorProbDensityClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::uqBaseVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues)
  :
  m_env                 (domainSpace.env()),
  m_prefix              ((std::string)(prefix)+"pd_"),
  m_domainSpace         (domainSpace),
  m_domainMinValues     (new V(domainMinValues         )),
  m_domainMaxValues     (new V(domainMaxValues         )),
  m_domainExpectedValues(new V(domainExpectedValues    )),
  m_domainVarianceValues(domainSpace.newVector(INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorProbDensityClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::uqBaseVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues)
  :
  m_env                 (domainSpace.env()),
  m_prefix              ((std::string)(prefix)+"pd_"),
  m_domainSpace         (domainSpace),
  m_domainMinValues     (new V(domainMinValues         )),
  m_domainMaxValues     (new V(domainMaxValues         )),
  m_domainExpectedValues(domainSpace.newVector(      0.)),
  m_domainVarianceValues(domainSpace.newVector(INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorProbDensityClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::uqBaseVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace)
  :
  m_env                 (domainSpace.env()),
  m_prefix              ((std::string)(prefix)+"pd_"),
  m_domainSpace         (domainSpace),
  m_domainMinValues     (domainSpace.newVector(-INFINITY)),
  m_domainMaxValues     (domainSpace.newVector( INFINITY)),
  m_domainExpectedValues(domainSpace.newVector(       0.)),
  m_domainVarianceValues(domainSpace.newVector( INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorProbDensityClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::~uqBaseVectorProbDensityClass()
{
  delete m_domainMinValues;
  delete m_domainMaxValues;
}
#if 0
template <class V, class M>
void
uqBaseVectorProbDensityClass<V,M>::setComponent(
  unsigned int componentId,
  double       minValue,
  double       maxValue,
  double       expectedValue,
  double       varianceValues)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorProbDensityClass<V,M>::setComponent()"
              << ", componentId = "    << componentId
              << ", minValue = "       << minValue
              << ", maxValue = "       << maxValue
              << ", expectedValue = "  << expectedValue
              << ", varianceValues = " << varianceValues
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((componentId > m_components.size()),
                      m_env.rank(),
                      "uqBaseVectorProbDensityClass<V,M>::setComponent()",
                      "componentId is too big");

  if (m_components[componentId] == NULL) {
    m_components[componentId] = new uqBaseScalarProbDensityClass(minValue,
                                                                 maxValue,
                                                                 expectedValue,
                                                                 varianceValues);
  }
  else {
    m_components[componentId]->setMinValue   (minValue);
    m_components[componentId]->setMaxValue   (maxValue);
    m_components[componentId]->setExpectValue(expectedValue);
    m_components[componentId]->setStdDevValue(varianceValues);
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
  if (m_varianceValuess) delete m_varianceValuess;
  if (m_expectedValues ) delete m_expectedValues;
  if (m_domainMaxValues) delete m_domainMaxValues;
  if (m_domainMinValues) delete m_domainMinValues;
  m_varianceValuess = NULL;
  m_expectedValues  = NULL;
  m_domainMaxValues = NULL;
  m_domainMinValues = NULL;
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
uqBaseVectorProbDensityClass<V,M>::domainMinValues() const
{
  return *m_domainMinValues;
}

template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::domainMaxValues() const
{
  return *m_domainMaxValues;
}

template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::domainExpectedValues() const
{
  return *m_domainExpectedValues;
}

template <class V, class M>
const V&
uqBaseVectorProbDensityClass<V,M>::domainVarianceValues() const
{
  return *m_domainVarianceValues;
}

template <class V, class M>
bool
uqBaseVectorProbDensityClass<V,M>::outOfDomainBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->domainMinValues()) ||
          v.atLeastOneComponentBiggerThan (this->domainMaxValues()));
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
                                  const V&                       domainMinValues,
                                  const V&                       domainMaxValues,
                                  const V&                       domainExpectedValues,
                                  const V&                       domainVarianceValues,
                                  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                                  const void* routineDataPtr,
                                  bool routineComputesMinus2LogOfDensity);
  uqGenericVectorProbDensityClass(const char*                    prefix,
                                  const uqVectorSpaceClass<V,M>& domainSpace,
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
  using uqBaseVectorProbDensityClass<V,M>::m_domainMinValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainMaxValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainVarianceValues;
};

template<class V, class M>
uqGenericVectorProbDensityClass<V,M>::uqGenericVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr,
  bool        routineComputesMinus2LogOfDensity)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace,domainMinValues,domainMaxValues,domainExpectedValues,domainVarianceValues),
  m_routinePtr                       (routinePtr),
  m_routineDataPtr                   (routineDataPtr),
  m_routineComputesMinus2LogOfDensity(routineComputesMinus2LogOfDensity)
{
}

template<class V, class M>
uqGenericVectorProbDensityClass<V,M>::uqGenericVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr,
  bool        routineComputesMinus2LogOfDensity)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace),
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
  using uqBaseVectorProbDensityClass<V,M>::m_domainMinValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainMaxValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainVarianceValues;
};

template<class V,class M>
uqBayesianVectorProbDensityClass<V,M>::uqBayesianVectorProbDensityClass(
  const char*                              prefix,
  const uqBaseVectorProbDensityClass<V,M>* priorDensity,
  const uqBaseVectorProbDensityClass<V,M>* likelihoodFunction)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"bay").c_str(),priorDensity->domainSpace()),
  m_priorDensity      (priorDensity),
  m_likelihoodFunction(likelihoodFunction)
{
#if 0
  uqMiscSelectMax(m_priorDensity->domainMinValues(),
                  m_likelihoodFunction->domainMinValues(),
                  *m_domainMinValues);
  uqMiscSelectMin(m_priorDensity->domainMaxValues(),
                  m_likelihoodFunction->domainMaxValues(),
                  *m_domainMaxValues);
#else
  for (unsigned int i = 0; i < m_domainMinValues->size(); ++i) {
    (*m_domainMinValues)[i] = std::max(m_priorDensity->domainMinValues()[i],
                                       m_likelihoodFunction->domainMinValues()[i]);
  }
  for (unsigned int i = 0; i < m_domainMaxValues->size(); ++i) {
    (*m_domainMaxValues)[i] = std::min(m_priorDensity->domainMaxValues()[i],
                                       m_likelihoodFunction->domainMaxValues()[i]);
  }
#endif
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
                                   const V&                       domainMinValues,
                                   const V&                       domainMaxValues,
                                   const V&                       domainExpectedValues,
                                   const V&                       domainVarianceValues);
  uqGaussianVectorProbDensityClass(const char*                    prefix,
                                   const uqVectorSpaceClass<V,M>& domainSpace,
                                   const V&                       domainMinValues,
                                   const V&                       domainMaxValues,
                                   const V&                       domainExpectedValues,
                                   const M&                       covMatrix);
 ~uqGaussianVectorProbDensityClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  const M*                                 m_covMatrix;
  uqDefault_M2lPriorRoutine_DataType<V,M>  m_m2lPriorRoutine_Data;
  const uqBaseVectorProbDensityClass<V,M>* m_probDensity;

  using uqBaseVectorProbDensityClass<V,M>::m_env;
  using uqBaseVectorProbDensityClass<V,M>::m_prefix;
  using uqBaseVectorProbDensityClass<V,M>::m_domainSpace;
  using uqBaseVectorProbDensityClass<V,M>::m_domainMinValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainMaxValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorProbDensityClass<V,M>::m_domainVarianceValues;

  void commonConstructor();
};

template<class V,class M>
uqGaussianVectorProbDensityClass<V,M>::uqGaussianVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues,domainExpectedValues,domainVarianceValues),
  m_covMatrix                      (m_domainSpace.newDiagMatrix(domainVarianceValues*domainVarianceValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorProbDensityClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  commonConstructor();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorProbDensityClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorProbDensityClass<V,M>::uqGaussianVectorProbDensityClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  uqBaseVectorProbDensityClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues,domainExpectedValues),
  m_covMatrix                      (new M(covMatrix))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorProbDensityClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  commonConstructor();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorProbDensityClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
void
uqGaussianVectorProbDensityClass<V,M>::commonConstructor()
{
#if 0
  V tmpVec(m_domainSpace.zeroVector());
  for (unsigned int i = 0; i < m_domainSpace.dim(); ++i) {
    sigma = m_components[i]->varianceValues();
    tmpVec[i] = sigma*sigma;
  }
  m_covMatrix = m_domainSpace.newDiagMatrix(tmpVec);
#endif

  m_m2lPriorRoutine_Data.paramPriorMus       = m_domainExpectedValues;
  m_m2lPriorRoutine_Data.paramPriorVariances = m_domainVarianceValues;
  m_probDensity = new uqGenericVectorProbDensityClass<V,M>(m_prefix.c_str(),
                                                           m_domainSpace,
                                                           *m_domainMinValues,
                                                           *m_domainMaxValues,
                                                           *m_domainExpectedValues,
                                                           *m_domainVarianceValues,
                                                           uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                           (void *) &m_m2lPriorRoutine_Data,
                                                           true); // the routine computes [-2.*ln(Likelihood)]
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqGaussianVectorProbDensityClass<V,M>::constructor()"
              << ", prefix = "      << m_prefix
              << ": priorMus = "    << *m_m2lPriorRoutine_Data.paramPriorMus
              << ", priorSigmas = " << *m_m2lPriorRoutine_Data.paramPriorVariances
              << std::endl;
  }

  return;
}

template<class V,class M>
uqGaussianVectorProbDensityClass<V,M>::~uqGaussianVectorProbDensityClass()
{
  delete m_covMatrix;
  delete m_probDensity;
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
