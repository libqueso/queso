/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_VECTOR_PROB_DENSITY_H__
#define __UQ_VECTOR_PROB_DENSITY_H__

#include <uqEnvironment.h>
#include <math.h>
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqBaseVectorPdfClass<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [1]"
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqBaseVectorPdfClass<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [2]"
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqBaseVectorPdfClass<V,M>::constructor() [3]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqBaseVectorPdfClass<V,M>::constructor() [3]"
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqBayesianVectorPdfClass<V,M>::actualValue()"
                           << ": domainVector = " << domainVector
                           << std::endl;
  }

  V* gradVLike = NULL;
  if (gradVector) gradVLike = &m_tmpVector1;

  M* hessianMLike = NULL;
  if (hessianMatrix) hessianMLike = m_tmpMatrix;

  V* hessianELike = NULL;
  if (hessianEffect) hessianELike = &m_tmpVector2;

  double value1 = m_priorDensity.actualValue      (domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);
  double value2 = m_likelihoodFunction.actualValue(domainVector,domainDirection,gradVLike ,hessianMLike ,hessianELike );

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.fullRank(),
                      "uqBayesianVectorPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = value1*value2;

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqBayesianVectorPdfClass<V,M>::actualValue()"
                           << ": domainVector = " << domainVector
                           << ", returnValue = "  << returnValue
                           << std::endl;
  }

  return returnValue;
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                           << ": domainVector = " << domainVector
                           << std::endl;
  }

  V* gradVLike = NULL;
  if (gradVector) gradVLike = &m_tmpVector1;

  M* hessianMLike = NULL;
  if (hessianMatrix) hessianMLike = m_tmpMatrix;

  V* hessianELike = NULL;
  if (hessianEffect) hessianELike = &m_tmpVector2;

  double value1 = m_priorDensity.minus2LnValue      (domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "In uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                           << ", domainVector = " << domainVector
                           << ": about to call likelihood()"
                           << std::endl;
  }

  double value2 = m_likelihoodFunction.minus2LnValue(domainVector,domainDirection,gradVLike, hessianMLike, hessianELike );

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "In uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                           << ", domainVector = " << domainVector
                           << ": value1 = "       << value1
                           << ", value2 = "       << value2
                           << std::endl;
    if (gradVector) {
      *m_env.subDisplayOutputFile() << "In uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                             << ", domainVector = " << domainVector
                             << ": gradVector = "   << *gradVector
                             << ", gradVLike = "    << *gradVLike
                             << std::endl;
    }
    if (hessianMatrix) {
      *m_env.subDisplayOutputFile() << "In uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                             << ", domainVector = "  << domainVector
                             << ": hessianMatrix = " << *hessianMatrix
                             << ", hessianMLike = "  << *hessianMLike
                             << std::endl;
    }
    if (hessianEffect) {
      *m_env.subDisplayOutputFile() << "In uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                             << ", domainVector = "  << domainVector
                             << ": hessianEffect = " << *hessianEffect
                             << ", hessianELike = "  << *hessianELike
                             << std::endl;
    }
  }

  if (gradVector   ) *gradVector    += *gradVLike;
  if (hessianMatrix) *hessianMatrix += *hessianMLike;
  if (hessianEffect) *hessianEffect += *hessianELike;

  double returnValue = value1+value2;

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqBayesianVectorPdfClass<V,M>::minus2LnValue()"
                           << ": domainVector = " << domainVector
                           << ", returnValue = "  << returnValue
                           << std::endl;
  }

  return returnValue;
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

  double   actualValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double   minus2LnValue  (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  void     updateExpVector(const V& newExpVector);
  void     updateCovMatrix(const M& newCovMatrix);
  const M& covMatrix      () const;

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
  m_covMatrix        (m_domainSet.vectorSpace().newDiagMatrix(domainVarVector))
{
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqGaussianVectorPdfClass<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "In uqGaussianVectorPdfClass<V,M>::constructor()"
                         //<< ", prefix = "     << m_prefix
                           << ": meanVector = " << this->domainExpVector()
	                   << ", Variances = "  << this->domainVarVector()
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqGaussianVectorPdfClass<V,M>::constructor() [1]"
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqGaussianVectorPdfClass<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "In uqGaussianVectorPdfClass<V,M>::constructor()"
                         //<< ", prefix = "            << m_prefix
                           << ": meanVector = "        << this->domainExpVector()
	                   << ", Covariance Matrix = " << covMatrix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqGaussianVectorPdfClass<V,M>::constructor() [2]"
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqGaussianVectorPdfClass<V,M>::actualValue()"
                           << ", meanVector = "   << *m_domainExpVector
	                   << ", covMatrix = "    << *m_covMatrix
                           << ": domainVector = " << domainVector
                           << std::endl;
  }

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.fullRank(),
                      "uqGaussianVectorPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = std::exp(-0.5*this->minus2LnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqGaussianVectorPdfClass<V,M>::actualValue()"
                           << ", meanVector = "   << *m_domainExpVector
	                   << ", covMatrix = "    << *m_covMatrix
                           << ": domainVector = " << domainVector
                           << ", returnValue = "  << returnValue
                           << std::endl;
  }

  return returnValue;
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqGaussianVectorPdfClass<V,M>::minus2LnValue()"
                           << ", meanVector = "   << *m_domainExpVector
	                   << ", covMatrix = "    << *m_covMatrix
                           << ": domainVector = " << domainVector
                           << std::endl;
  }

  double returnValue = 0.;

  if (m_diagonalCovMatrix) {
    V diffVec(domainVector - this->domainExpVector());
    returnValue = ((diffVec*diffVec)/this->domainVarVector()).sumOfComponents();
  }
  else {
    V diffVec(domainVector - this->domainExpVector());
    V tmpVec = this->m_covMatrix->invertMultiply(diffVec);
    returnValue = (diffVec*tmpVec).sumOfComponents();
  }

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.fullRank(),
                      "uqGaussianVectorPdfClass<V,M>::minus2LnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqGaussianVectorPdfClass<V,M>::minus2LnValue()"
                           << ", meanVector = "   << *m_domainExpVector
	                   << ", covMatrix = "    << *m_covMatrix
                           << ": domainVector = " << domainVector
                           << ", returnValue = "  << returnValue
                           << std::endl;
  }

  return returnValue;
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

template<class V, class M>
const M&
uqGaussianVectorPdfClass<V,M>::covMatrix() const
{
  return *m_covMatrix;
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
  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Entering uqUniformVectorPdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayOutputFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayOutputFile() << "Leaving uqUniformVectorPdfClass<V,M>::constructor()"
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
