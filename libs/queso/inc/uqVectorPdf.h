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
#include <uqScalarFunction.h>

//*****************************************************
// Classes to accomodate a probability density.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorPdfClass : public uqBaseScalarFunctionClass<V,M> {
public:
           uqBaseVectorPdfClass(const char*                  prefix,
                                const uqVectorSetClass<V,M>& domainSet,
                                const V&                     domainExpVector,
                                const V&                     domainVarVector);
           uqBaseVectorPdfClass(const char*                  prefix,
                                const uqVectorSetClass<V,M>& domainSet,
                                const V&                     domainExpVector);
           uqBaseVectorPdfClass(const char*                  prefix,
                                const uqVectorSetClass<V,M>& domainSet);
  virtual ~uqBaseVectorPdfClass();

  virtual       double                 actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;
  virtual       double                 minus2LnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;

//const   uqBaseScalarPdfClass<double>& component           (unsigned int componentId) const;
          const V&                      domainExpVector() const;
          const V&                      domainVarVector() const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;

  V* m_domainExpVector;
  V* m_domainVarVector;

//std::vector<uqBaseScalarPdfClass<double>*> m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
//uqBaseScalarPdfClass<double>               m_dummyComponent;
};

template<class V, class M>
uqBaseVectorPdfClass<V,M>::uqBaseVectorPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     domainExpVector,
  const V&                     domainVarVector)
  :
  uqBaseScalarFunctionClass<V,M>(((std::string)(prefix)+"pd_").c_str(), domainSet),
  m_domainExpVector(new V(domainExpVector)),
  m_domainVarVector(new V(domainVarVector))
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
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     domainExpVector)
  :
  uqBaseScalarFunctionClass<V,M>(((std::string)(prefix)+"pd_").c_str(), domainSet),
  m_domainExpVector(new V(domainExpVector    )),
  m_domainVarVector(domainSet.vectorSpace().newVector(INFINITY))
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
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet)
  :
  uqBaseScalarFunctionClass<V,M>(((std::string)(prefix)+"pd_").c_str(), domainSet),
  m_domainExpVector(domainSet.vectorSpace().newVector(       0.)),
  m_domainVarVector(domainSet.vectorSpace().newVector( INFINITY))
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
uqBaseVectorPdfClass<V,M>::~uqBaseVectorPdfClass()
{
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

template <class V, class M>
const V&
uqBaseVectorPdfClass<V,M>::domainExpVector() const
{
  return *m_domainExpVector;
}

template <class V, class M>
const V&
uqBaseVectorPdfClass<V,M>::domainVarVector() const
{
  return *m_domainVarVector;
}

//*****************************************************
// Generic probability density class
//*****************************************************
template<class V, class M>
class uqGenericVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqGenericVectorPdfClass(const char*                           prefix,
                          const uqBaseScalarFunctionClass<V,M>& scalarFunction,
                          const V&                              domainExpVector,
                          const V&                              domainVarVector);
  uqGenericVectorPdfClass(const char*                           prefix,
                          const uqBaseScalarFunctionClass<V,M>& scalarFunction);
 ~uqGenericVectorPdfClass();

  double actualValue     (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double minus2LnValue   (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseVectorPdfClass<V,M>::m_domainExpVector;
  using uqBaseVectorPdfClass<V,M>::m_domainVarVector;

  const uqBaseScalarFunctionClass<V,M>& m_scalarFunction;
};

template<class V, class M>
uqGenericVectorPdfClass<V,M>::uqGenericVectorPdfClass(
  const char*                           prefix,
  const uqBaseScalarFunctionClass<V,M>& scalarFunction,
  const V&                              domainExpVector,
  const V&                              domainVarVector)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),
                            scalarFunction.domainSet(),
                            domainExpVector,
                            domainVarVector),
  m_scalarFunction(scalarFunction)
{
}

template<class V, class M>
uqGenericVectorPdfClass<V,M>::uqGenericVectorPdfClass(
  const char*                           prefix,
  const uqBaseScalarFunctionClass<V,M>& scalarFunction)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),scalarFunction.domainSet()),
  m_scalarFunction(scalarFunction)
{
}

template<class V, class M>
uqGenericVectorPdfClass<V,M>::~uqGenericVectorPdfClass()
{
}

template<class V, class M>
double
uqGenericVectorPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  return m_scalarFunction.actualValue(domainVector);
}

template<class V, class M>
double
uqGenericVectorPdfClass<V,M>::minus2LnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  return m_scalarFunction.minus2LnValue(domainVector);
}

//*****************************************************
// Bayesian probability density class
//*****************************************************
template<class V, class M>
class uqBayesianVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqBayesianVectorPdfClass(const char*                           prefix,
                           const uqBaseVectorPdfClass     <V,M>& priorDensity,
                           const uqBaseScalarFunctionClass<V,M>& likelihoodFunction,
                           const uqVectorSetClass         <V,M>& intersectionDomain); 
 ~uqBayesianVectorPdfClass();

  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double minus2LnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseVectorPdfClass<V,M>::m_domainExpVector;
  using uqBaseVectorPdfClass<V,M>::m_domainVarVector;

  const uqBaseVectorPdfClass     <V,M>& m_priorDensity;
  const uqBaseScalarFunctionClass<V,M>& m_likelihoodFunction;

  mutable V m_tmpVector1;
  mutable V m_tmpVector2;
  mutable M* m_tmpMatrix;
};

template<class V,class M>
uqBayesianVectorPdfClass<V,M>::uqBayesianVectorPdfClass(
  const char*                           prefix,
  const uqBaseVectorPdfClass     <V,M>& priorDensity,
  const uqBaseScalarFunctionClass<V,M>& likelihoodFunction,
  const uqVectorSetClass         <V,M>& intersectionDomain)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"bay").c_str(),intersectionDomain),
  m_priorDensity           (priorDensity),
  m_likelihoodFunction     (likelihoodFunction),
  m_tmpVector1             (m_domainSet.vectorSpace().zeroVector()),
  m_tmpVector2             (m_domainSet.vectorSpace().zeroVector()),
  m_tmpMatrix              (m_domainSet.vectorSpace().newMatrix())
{
}

template<class V,class M>
uqBayesianVectorPdfClass<V,M>::~uqBayesianVectorPdfClass()
{
  delete m_tmpMatrix;
}

template<class V, class M>
double
uqBayesianVectorPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  V* gradV = NULL;
  if (gradVector) gradV = &m_tmpVector1;

  M* hessianM = NULL;
  if (hessianMatrix) hessianM = m_tmpMatrix;

  V* hessianE = NULL;
  if (hessianEffect) hessianE = &m_tmpVector2;

  double value1 = m_priorDensity.actualValue      (domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);
  double value2 = m_likelihoodFunction.actualValue(domainVector,domainDirection,gradV     ,hessianM     ,hessianE     );

  return value1*value2;
}

template<class V, class M>
double
uqBayesianVectorPdfClass<V,M>::minus2LnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  V* gradV = NULL;
  if (gradVector) gradV = &m_tmpVector1;

  M* hessianM = NULL;
  if (hessianMatrix) hessianM = m_tmpMatrix;

  V* hessianE = NULL;
  if (hessianEffect) hessianE = &m_tmpVector2;

  double value1 = m_priorDensity.minus2LnValue      (domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);
  double value2 = m_likelihoodFunction.minus2LnValue(domainVector,domainDirection,gradV     ,hessianM     ,hessianE     );
  if (gradVector   ) *gradVector    += *gradV;
  if (hessianMatrix) *hessianMatrix += *hessianM;
  if (hessianEffect) *hessianEffect += *hessianE;

  //if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
  //  std::cout << "In uqBayesianVectorPdfClass<P_V,P_M>::minus2LnValue()"
  //            << ", -2ln(prior) = " << value1
  //            << ", -2ln(like) = "  << value2
  //            << std::endl;
  //}

  return value1+value2;
}

//*****************************************************
// Gaussian probability density class
//*****************************************************
template<class V, class M>
class uqGaussianVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqGaussianVectorPdfClass(const char*                  prefix,
                           const uqVectorSetClass<V,M>& domainSet,
                           const V&                     domainExpVector,
                           const V&                     domainVarVector);
  uqGaussianVectorPdfClass(const char*                  prefix,
                           const uqVectorSetClass<V,M>& domainSet,
                           const V&                     domainExpVector,
                           const M&                     covMatrix);
 ~uqGaussianVectorPdfClass();

  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double minus2LnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  void   updateExpVector  (const V& newExpVector);
  void   updateCovMatrix  (const M& newCovMatrix);

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseVectorPdfClass<V,M>::m_domainExpVector;
  using uqBaseVectorPdfClass<V,M>::m_domainVarVector;

  bool     m_diagonalCovMatrix;
  const M* m_covMatrix;
};

template<class V,class M>
uqGaussianVectorPdfClass<V,M>::uqGaussianVectorPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     domainExpVector,
  const V&                     domainVarVector)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet,domainExpVector,domainVarVector),
  m_diagonalCovMatrix(true),
  m_covMatrix                      (m_domainSet.newDiagMatrix(domainVarVector))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorPdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqGaussianVectorPdfClass<V,M>::constructor()"
              << ", prefix = "      << m_prefix
              << ": Mus = "    << this->domainExpVector()
	      << ", Variances = " << this->domainVarVector()
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
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     domainExpVector,
  const M&                     covMatrix)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet,domainExpVector),
  m_diagonalCovMatrix      (false),
  m_covMatrix              (new M(covMatrix))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorPdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqGaussianVectorPdfClass<V,M>::constructor()"
              << ", prefix = "      << m_prefix
              << ": Mus = "    << this->domainExpVector()
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
uqGaussianVectorPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  return exp(-0.5*this->minus2LnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}

template<class V, class M>
double
uqGaussianVectorPdfClass<V,M>::minus2LnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (m_diagonalCovMatrix) {
    V diffVec(domainVector - this->domainExpVector());
    return ((diffVec*diffVec)/this->domainVarVector()).sumOfComponents();
  }
  else {
    V diffVec(domainVector - this->domainExpVector());
    V tmpVec = this->m_covMatrix->invertMultiply(diffVec);
    return (diffVec*tmpVec).sumOfComponents();
  }
}

template<class V, class M>
void
uqGaussianVectorPdfClass<V,M>::updateExpVector(const V& newExpVector)
{
  // delete old expected values (alloced at construction or last call to this function)
  delete m_domainExpVector;
  m_domainExpVector = new V(newExpVector);
  return;
}

template<class V, class M>
void
uqGaussianVectorPdfClass<V,M>::updateCovMatrix(const M& newCovMatrix)
{
  // delete old expected values (alloced at construction or last call to this function)
  delete m_covMatrix;
  m_covMatrix = new M(newCovMatrix);
  return;
}

//*****************************************************
// Uniform probability density class
//*****************************************************
template<class V, class M>
class uqUniformVectorPdfClass : public uqBaseVectorPdfClass<V,M> {
public:
  uqUniformVectorPdfClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& domainSet);
 ~uqUniformVectorPdfClass();

  double actualValue     (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double minus2LnValue   (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseVectorPdfClass<V,M>::m_domainExpVector;
  using uqBaseVectorPdfClass<V,M>::m_domainVarVector;
};

template<class V,class M>
uqUniformVectorPdfClass<V,M>::uqUniformVectorPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet)
  :
  uqBaseVectorPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),
                            domainSet)
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
uqUniformVectorPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  return 1.;
}

template<class V, class M>
double
uqUniformVectorPdfClass<V,M>::minus2LnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  return 0.;
}
#endif // __UQ_VECTOR_PROB_DENSITY_H__
