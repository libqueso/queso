/* uq/libs/queso/inc/uqVectorPdf.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
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
class uqBaseVectorPdfClass {
public:
           uqBaseVectorPdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace,
                                const V&                       domainMinValues,
                                const V&                       domainMaxValues,
                                const V&                       domainExpectedValues,
                                const V&                       domainVarianceValues);
           uqBaseVectorPdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace,
                                const V&                       domainMinValues,
                                const V&                       domainMaxValues,
                                const V&                       domainExpectedValues);
           uqBaseVectorPdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace,
                                const V&                       domainMinValues,
                                const V&                       domainMaxValues);
           uqBaseVectorPdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace);
  virtual ~uqBaseVectorPdfClass();

  const   uqVectorSpaceClass<V,M>&              domainSpace               ()                     const;
  virtual double                                actualDensity             (const V& paramValues) const = 0;
  virtual double                                minus2LnDensity           (const V& paramValues) const = 0;

//const   uqBaseScalarPdfClass<double>& component                 (unsigned int componentId) const;
          bool                                  outOfDomainBounds         (const V& v) const;
  const   V&                                    domainMinValues           () const;
  const   V&                                    domainMaxValues           () const;
  const   V&                                    domainExpectedValues      () const;
  const   V&                                    domainVarianceValues      () const;

protected:

  const   uqBaseEnvironmentClass&      m_env;
          std::string              m_prefix;
  const   uqVectorSpaceClass<V,M>& m_domainSpace;

          V*                       m_domainMinValues;
          V*                       m_domainMaxValues;
          V*                       m_domainExpectedValues;
          V*                       m_domainVarianceValues;

//std::vector<uqBaseScalarPdfClass<double>*> m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
//uqBaseScalarPdfClass<double>               m_dummyComponent;
};

template<class V, class M>
uqBaseVectorPdfClass<V,M>::uqBaseVectorPdfClass(
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
    std::cout << "Entering uqBaseVectorPdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorPdfClass<V,M>::uqBaseVectorPdfClass(
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
  //m_domainVarianceValues(domainSpace.newVector(INFINITY)) QUESTION: Ask Ernesto why this won't compile
  m_domainVarianceValues(new V(domainExpectedValues)) // FIXME: These variance values are bogus
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorPdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorPdfClass<V,M>::uqBaseVectorPdfClass(
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
    std::cout << "Entering uqBaseVectorPdfClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorPdfClass<V,M>::uqBaseVectorPdfClass(
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
    std::cout << "Entering uqBaseVectorPdfClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorPdfClass<V,M>::~uqBaseVectorPdfClass()
{
  delete m_domainMinValues;
  delete m_domainMaxValues;
}

#if 0
template <class V, class M>
const uqBaseScalarPdfClass&
uqBaseVectorPdfClass<V,M>::component(unsigned int componentId) const
{
  if (componentId > m_components.size()) return m_dummyComponent;
  if (m_components[componentId] == NULL) return m_dummyComponent;
  return *(m_components[componentId]);
}
#endif
#if 1
template <class V, class M>
const V&
uqBaseVectorPdfClass<V,M>::domainMinValues() const
{
  return *m_domainMinValues;
}

template <class V, class M>
const V&
uqBaseVectorPdfClass<V,M>::domainMaxValues() const
{
  return *m_domainMaxValues;
}
#endif
template <class V, class M>
const V&
uqBaseVectorPdfClass<V,M>::domainExpectedValues() const
{
  return *m_domainExpectedValues;
}

template <class V, class M>
const V&
uqBaseVectorPdfClass<V,M>::domainVarianceValues() const
{
  return *m_domainVarianceValues;
}

template <class V, class M>
bool
uqBaseVectorPdfClass<V,M>::outOfDomainBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->domainMinValues()) ||
          v.atLeastOneComponentBiggerThan (this->domainMaxValues()));
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorPdfClass<V,M>::domainSpace() const
{
  return m_domainSpace;
}

//*****************************************************
// Generic probability density class
//*****************************************************
template<class V, class M>
class uqGenericVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqGenericVectorPdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          const V&                       domainMinValues,
                          const V&                       domainMaxValues,
                          const V&                       domainExpectedValues,
                          const V&                       domainVarianceValues,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                          const void* routineDataPtr,
                          bool routineComputesMinus2LogOfDensity);
  uqGenericVectorPdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                          const void* routineDataPtr,
                          bool routineComputesMinus2LogOfDensity);
 ~uqGenericVectorPdfClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;

  bool m_routineComputesMinus2LogOfDensity;

  using uqBaseVectorPdfClass<V,M>::m_env;
  using uqBaseVectorPdfClass<V,M>::m_prefix;
  using uqBaseVectorPdfClass<V,M>::m_domainSpace;
  using uqBaseVectorPdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorPdfClass<V,M>::m_domainMaxValues;
  using uqBaseVectorPdfClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorPdfClass<V,M>::m_domainVarianceValues;
};

template<class V, class M>
uqGenericVectorPdfClass<V,M>::uqGenericVectorPdfClass(
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
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),
                                    domainSpace,
                                    domainMinValues,
                                    domainMaxValues,
                                    domainExpectedValues,
                                    domainVarianceValues),
  m_routinePtr                       (routinePtr),
  m_routineDataPtr                   (routineDataPtr),
  m_routineComputesMinus2LogOfDensity(routineComputesMinus2LogOfDensity)
{
}

template<class V, class M>
uqGenericVectorPdfClass<V,M>::uqGenericVectorPdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr,
  bool        routineComputesMinus2LogOfDensity)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace),
  m_routinePtr                       (routinePtr),
  m_routineDataPtr                   (routineDataPtr),
  m_routineComputesMinus2LogOfDensity(routineComputesMinus2LogOfDensity)
{
}

template<class V, class M>
uqGenericVectorPdfClass<V,M>::~uqGenericVectorPdfClass()
{
}

template<class V, class M>
double
uqGenericVectorPdfClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfDensity == false) {
    value = -2.*log(value);
  }

  return value;
}

template<class V, class M>
double
uqGenericVectorPdfClass<V,M>::actualDensity(const V& paramValues) const
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
class uqBayesianVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqBayesianVectorPdfClass(const char*                              prefix,
                           const uqBaseVectorPdfClass<V,M>* priorDensity,
                           const uqBaseVectorPdfClass<V,M>* likelihoodFunction); 
 ~uqBayesianVectorPdfClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  const uqBaseVectorPdfClass<V,M>* m_priorDensity;
  const uqBaseVectorPdfClass<V,M>* m_likelihoodFunction;

  using uqBaseVectorPdfClass<V,M>::m_env;
  using uqBaseVectorPdfClass<V,M>::m_prefix;
  using uqBaseVectorPdfClass<V,M>::m_domainSpace;
  using uqBaseVectorPdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorPdfClass<V,M>::m_domainMaxValues;
  using uqBaseVectorPdfClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorPdfClass<V,M>::m_domainVarianceValues;
};

template<class V,class M>
uqBayesianVectorPdfClass<V,M>::uqBayesianVectorPdfClass(
  const char*                              prefix,
  const uqBaseVectorPdfClass<V,M>* priorDensity,
  const uqBaseVectorPdfClass<V,M>* likelihoodFunction)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"bay").c_str(),priorDensity->domainSpace()),
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
uqBayesianVectorPdfClass<V,M>::~uqBayesianVectorPdfClass()
{
}

template<class V, class M>
double
uqBayesianVectorPdfClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value1 = m_priorDensity->minus2LnDensity(paramValues);
  double value2 = m_likelihoodFunction->minus2LnDensity(paramValues);

  //if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
  //  std::cout << "In uqBayesianVectorPdfClass<P_V,P_M>::minus2LnDensity()"
  //            << ", -2ln(prior) = " << value1
  //            << ", -2ln(like) = "  << value2
  //            << std::endl;
  //}

  return value1+value2;
}

template<class V, class M>
double
uqBayesianVectorPdfClass<V,M>::actualDensity(const V& paramValues) const
{
  double value1 = m_priorDensity->actualDensity(paramValues);
  double value2 = m_likelihoodFunction->actualDensity(paramValues);

  return value1*value2;
}

//*****************************************************
// Gaussian probability density class
//*****************************************************
template<class V, class M>
class uqGaussianVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqGaussianVectorPdfClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& domainSpace,
                           const V&                       domainMinValues,
                           const V&                       domainMaxValues,
                           const V&                       domainExpectedValues,
                           const V&                       domainVarianceValues);
  uqGaussianVectorPdfClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& domainSpace,
                           const V&                       domainMinValues,
                           const V&                       domainMaxValues,
                           const V&                       domainExpectedValues,
                           const M&                       covMatrix);
 ~uqGaussianVectorPdfClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  bool     m_diagonalCovMatrix;
  const M* m_covMatrix;

  using uqBaseVectorPdfClass<V,M>::m_env;
  using uqBaseVectorPdfClass<V,M>::m_prefix;
  using uqBaseVectorPdfClass<V,M>::m_domainSpace;
  using uqBaseVectorPdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorPdfClass<V,M>::m_domainMaxValues;
  using uqBaseVectorPdfClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorPdfClass<V,M>::m_domainVarianceValues;

};

template<class V,class M>
uqGaussianVectorPdfClass<V,M>::uqGaussianVectorPdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues,domainExpectedValues,domainVarianceValues),
  m_diagonalCovMatrix(true),
  m_covMatrix                      (m_domainSpace.newDiagMatrix(domainVarianceValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorPdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqGaussianVectorPdfClass<V,M>::constructor()"
              << ", prefix = "      << m_prefix
              << ": Mus = "    << this->domainExpectedValues()
	      << ", Variances = " << this->domainVarianceValues()
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorPdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorPdfClass<V,M>::uqGaussianVectorPdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues,domainExpectedValues),
  m_diagonalCovMatrix(false),
  m_covMatrix                      (new M(covMatrix))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorPdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqGaussianVectorPdfClass<V,M>::constructor()"
              << ", prefix = "      << m_prefix
              << ": Mus = "    << this->domainExpectedValues()
	      << ", Covariance Matrix = " << covMatrix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorPdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorPdfClass<V,M>::~uqGaussianVectorPdfClass()
{
  delete m_covMatrix;
}

template<class V, class M>
double
uqGaussianVectorPdfClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  if( m_diagonalCovMatrix ){
    V diffVec(paramValues - this->domainExpectedValues());
    return ((diffVec*diffVec)/this->domainVarianceValues()).sumOfComponents();
  } else{
    V diffVec(paramValues - this->domainExpectedValues());
    V tmpVec = this->m_covMatrix->invertMultiply(diffVec);
    return (diffVec*tmpVec).sumOfComponents();
  }
}

template<class V, class M>
double
uqGaussianVectorPdfClass<V,M>::actualDensity(const V& paramValues) const
{
  return exp(-0.5*this->minus2LnDensity(paramValues));
}

//*****************************************************
// Uniform probability density class
//*****************************************************
template<class V, class M>
class uqUniformVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqUniformVectorPdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          const V&                       domainMinValues,
                          const V&                       domainMaxValues);
 ~uqUniformVectorPdfClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  using uqBaseVectorPdfClass<V,M>::m_env;
  using uqBaseVectorPdfClass<V,M>::m_prefix;
  using uqBaseVectorPdfClass<V,M>::m_domainSpace;
  using uqBaseVectorPdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorPdfClass<V,M>::m_domainMaxValues;
  using uqBaseVectorPdfClass<V,M>::m_domainExpectedValues;
  using uqBaseVectorPdfClass<V,M>::m_domainVarianceValues;
};

template<class V,class M>
uqUniformVectorPdfClass<V,M>::uqUniformVectorPdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),
                                    domainSpace,
                                    domainMinValues,
                                    domainMaxValues)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqUniformVectorPdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqUniformVectorPdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqUniformVectorPdfClass<V,M>::~uqUniformVectorPdfClass()
{
}

template<class V, class M>
double
uqUniformVectorPdfClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  return 0.;
}

template<class V, class M>
double
uqUniformVectorPdfClass<V,M>::actualDensity(const V& paramValues) const
{
  return 1.;
}

#endif // __UQ_VECTOR_PROB_DENSITY_H__
