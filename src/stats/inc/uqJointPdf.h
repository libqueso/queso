//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_JOINT_PROB_DENSITY_H__
#define __UQ_JOINT_PROB_DENSITY_H__

#include <uqEnvironment.h>
#include <math.h>
#include <uqScalarFunction.h>
#include <boost/math/special_functions.hpp> // for Boost isnan. Note parantheses are important in function call.

//*****************************************************
// Classes to accomodate a probability density.
//*****************************************************

//*****************************************************
// Base class [PDF-00]
//*****************************************************
template<class V, class M>
class uqBaseJointPdfClass : public uqBaseScalarFunctionClass<V,M> {
public:
           uqBaseJointPdfClass(const char*                  prefix,
                               const uqVectorSetClass<V,M>& domainSet);
  virtual ~uqBaseJointPdfClass();

  virtual double actualValue                    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;
  virtual double lnValue                        (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;
  virtual void   setNormalizationStyle          (unsigned int value) const;
          void   setLogOfNormalizationFactor    (double value) const;
  virtual double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const = 0;
//const uqBaseScalarPdfClass<double>& component(unsigned int componentId) const;

protected:
          double commonComputeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;

  mutable unsigned int m_normalizationStyle;
  mutable double       m_logOfNormalizationFactor;
//std::vector<uqBaseScalarPdfClass<double>*> m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
//uqBaseScalarPdfClass<double>               m_dummyComponent;
};

template<class V, class M>
uqBaseJointPdfClass<V,M>::uqBaseJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet)
  :
  uqBaseScalarFunctionClass<V,M>(((std::string)(prefix)+"pd_").c_str(), domainSet),
  m_normalizationStyle(0),
  m_logOfNormalizationFactor(0.)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBaseJointPdfClass<V,M>::constructor() [3]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBaseJointPdfClass<V,M>::constructor() [3]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqBaseJointPdfClass<V,M>::~uqBaseJointPdfClass()
{
}

template<class V,class M>
void
uqBaseJointPdfClass<V,M>::setNormalizationStyle(unsigned int value) const
{
  m_normalizationStyle = value;
  return;
}

template<class V,class M>
void
uqBaseJointPdfClass<V,M>::setLogOfNormalizationFactor(double value) const
{
  m_logOfNormalizationFactor = value;
  return;
}

template<class V,class M>
double
uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  double volume = m_domainSet.volume();
  if (((boost::math::isnan)(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                )) {
    // Do nothing
  }
  else {
    const uqBoxSubsetClass<V,M>* boxSubset = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&m_domainSet);
    if (boxSubset == NULL) {
      // Do nothing
    }
    else {
      V tmpVec(m_domainSet.vectorSpace().zeroVector());
      double sum = 0.;
      for (unsigned int i = 0; i < numSamples; ++i) {
        tmpVec.cwSetUniform(boxSubset->minValues(),boxSubset->maxValues());
        sum += this->actualValue(tmpVec,NULL,NULL,NULL,NULL);
      }
      double avgValue = sum/((double) numSamples);
      value = -( log(avgValue) + log(volume) );
      if (updateFactorInternally) {
        m_logOfNormalizationFactor = value;
      }
    }
  }

  return value;
}

#if 0
template <class V, class M>
const uqBaseScalarPdfClass&
uqBaseJointPdfClass<V,M>::component(unsigned int componentId) const
{
  if (componentId > m_components.size()) return m_dummyComponent;
  if (m_components[componentId] == NULL) return m_dummyComponent;
  return *(m_components[componentId]);
}
#endif

//*****************************************************
// Generic probability density class [PDF-01]
//*****************************************************
template<class V, class M>
class uqGenericJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqGenericJointPdfClass(const char*                           prefix,
                         const uqBaseScalarFunctionClass<V,M>& scalarFunction);
 ~uqGenericJointPdfClass();

  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  const uqBaseScalarFunctionClass<V,M>& m_scalarFunction;
};

template<class V, class M>
uqGenericJointPdfClass<V,M>::uqGenericJointPdfClass(
  const char*                           prefix,
  const uqBaseScalarFunctionClass<V,M>& scalarFunction)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),scalarFunction.domainSet()),
  m_scalarFunction(scalarFunction)
{
}

template<class V, class M>
uqGenericJointPdfClass<V,M>::~uqGenericJointPdfClass()
{
}

template<class V, class M>
double
uqGenericJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  return ((exp(m_logOfNormalizationFactor))*m_scalarFunction.actualValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect)); // [PDF-01]
}

template<class V, class M>
double
uqGenericJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  return (m_logOfNormalizationFactor + m_scalarFunction.lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect)); // [PDF-01]
}

template<class V, class M>
double
uqGenericJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGenericJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGenericJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// Bayesian probability density class [PDF-02]
//*****************************************************
template<class V, class M>
class uqBayesianJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqBayesianJointPdfClass(const char*                           prefix,
                          const uqBaseJointPdfClass      <V,M>& priorDensity,
                          const uqBaseScalarFunctionClass<V,M>& likelihoodFunction,
                                double                          likelihoodExponent,
                          const uqVectorSetClass         <V,M>& intersectionDomain); 
 ~uqBayesianJointPdfClass();

  double actualValue              (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue                  (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  double lastComputedLogPrior     () const;
  double lastComputedLogLikelihood() const;
  void   setNormalizationStyle    (unsigned int value) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  const uqBaseJointPdfClass      <V,M>& m_priorDensity;
  const uqBaseScalarFunctionClass<V,M>& m_likelihoodFunction;
        double                          m_likelihoodExponent;
  mutable double                        m_lastComputedLogPrior;
  mutable double                        m_lastComputedLogLikelihood;

  mutable V  m_tmpVector1;
  mutable V  m_tmpVector2;
  mutable M* m_tmpMatrix;
};

template<class V,class M>
uqBayesianJointPdfClass<V,M>::uqBayesianJointPdfClass(
  const char*                           prefix,
  const uqBaseJointPdfClass     <V,M>&  priorDensity,
  const uqBaseScalarFunctionClass<V,M>& likelihoodFunction,
        double                          likelihoodExponent,
  const uqVectorSetClass         <V,M>& intersectionDomain)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"bay").c_str(),intersectionDomain),
  m_priorDensity             (priorDensity),
  m_likelihoodFunction       (likelihoodFunction),
  m_likelihoodExponent       (likelihoodExponent),
  m_lastComputedLogPrior     (0.),
  m_lastComputedLogLikelihood(0.),
  m_tmpVector1               (m_domainSet.vectorSpace().zeroVector()),
  m_tmpVector2               (m_domainSet.vectorSpace().zeroVector()),
  m_tmpMatrix                (m_domainSet.vectorSpace().newMatrix())
{
}

template<class V,class M>
uqBayesianJointPdfClass<V,M>::~uqBayesianJointPdfClass()
{
  delete m_tmpMatrix;
}

template<class V,class M>
void
uqBayesianJointPdfClass<V,M>::setNormalizationStyle(unsigned int value) const
{
  m_priorDensity.setNormalizationStyle(value);
  return;
}

template<class V,class M>
double
uqBayesianJointPdfClass<V,M>::lastComputedLogPrior() const
{
  return m_lastComputedLogPrior;
}

template<class V,class M>
double
uqBayesianJointPdfClass<V,M>::lastComputedLogLikelihood() const
{
  return m_lastComputedLogLikelihood;
}

template<class V, class M>
double
uqBayesianJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBayesianJointPdfClass<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqBayesianJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  V* gradVLike = NULL;
  if (gradVector) gradVLike = &m_tmpVector1;

  M* hessianMLike = NULL;
  if (hessianMatrix) hessianMLike = m_tmpMatrix;

  V* hessianELike = NULL;
  if (hessianEffect) hessianELike = &m_tmpVector2;

  double value1 = m_priorDensity.actualValue (domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);
  double value2 = 1.;
  if (m_likelihoodExponent != 0.) {
    value2 = m_likelihoodFunction.actualValue(domainVector,domainDirection,gradVLike ,hessianMLike ,hessianELike );
  }

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqBayesianJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = value1;
  if (m_likelihoodExponent == 0.) {
    // Do nothing
  }
  else if (m_likelihoodExponent == 1.) {
    returnValue *= value2;
  }
  else {
    returnValue *= pow(value2,m_likelihoodExponent);
  }
  returnValue *= exp(m_logOfNormalizationFactor); // [PDF-02] ???

  m_lastComputedLogPrior      = log(value1);
  m_lastComputedLogLikelihood = m_likelihoodExponent*log(value2);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBayesianJointPdfClass<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqBayesianJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBayesianJointPdfClass<V,M>::lnValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  V* gradVLike = NULL;
  if (gradVector) gradVLike = &m_tmpVector1;

  M* hessianMLike = NULL;
  if (hessianMatrix) hessianMLike = m_tmpMatrix;

  V* hessianELike = NULL;
  if (hessianEffect) hessianELike = &m_tmpVector2;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": about to call prior()"
                            << std::endl;
  }

  double value1 = m_priorDensity.lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": lnPrior = " << value1
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": about to call likelihood()"
                            << std::endl;
  }

  double value2 = 0.;
  if (m_likelihoodExponent != 0.) {
    value2 = m_likelihoodFunction.lnValue(domainVector,domainDirection,gradVLike, hessianMLike, hessianELike );
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": value1 = "       << value1
                            << ", value2 = "       << value2
                            << std::endl;
    if (gradVector) {
      *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                              << ", domainVector = " << domainVector
                              << ": gradVector = "   << *gradVector
                              << ", gradVLike = "    << *gradVLike
                              << std::endl;
    }
    if (hessianMatrix) {
      *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                              << ", domainVector = "  << domainVector
                              << ": hessianMatrix = " << *hessianMatrix
                              << ", hessianMLike = "  << *hessianMLike
                              << std::endl;
    }
    if (hessianEffect) {
      *m_env.subDisplayFile() << "In uqBayesianJointPdfClass<V,M>::lnValue()"
                              << ", domainVector = "  << domainVector
                              << ": hessianEffect = " << *hessianEffect
                              << ", hessianELike = "  << *hessianELike
                              << std::endl;
    }
  }

  if (gradVector   ) *gradVector    += *gradVLike;
  if (hessianMatrix) *hessianMatrix += *hessianMLike;
  if (hessianEffect) *hessianEffect += *hessianELike;

  double returnValue = value1;
  if (m_likelihoodExponent == 0.) {
    // Do nothing
  }
  else if (m_likelihoodExponent == 1.) {
    returnValue += value2;
  }
  else {
    returnValue += value2*m_likelihoodExponent;
  } // prudenci 2010/03/05
  returnValue += m_logOfNormalizationFactor; // [PDF-02] ???

  m_lastComputedLogPrior      = value1;
  m_lastComputedLogLikelihood = m_likelihoodExponent*value2;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBayesianJointPdfClass<V,M>::lnValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqBayesianJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  double volume = m_domainSet.volume();
  if (((boost::math::isnan)(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                )) {
    // Do nothing
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.worldRank(),
                        "uqBayesianJointPdfClass<V,M>::lnValue()",
                        "incomplete code for computeLogOfNormalizationFactor()");
  }

  return value;
}

//*****************************************************
// Gaussian probability density class [PDF-03]
//*****************************************************
template<class V, class M>
class uqGaussianJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqGaussianJointPdfClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& domainSet,
                          const V&                     lawExpVector,
                          const V&                     lawVarVector);
  uqGaussianJointPdfClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& domainSet,
                          const V&                     lawExpVector,
                          const M&                     lawCovMatrix);
 ~uqGaussianJointPdfClass();

  double   actualValue       (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double   lnValue           (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double   computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  void     updateLawExpVector(const V& newLawExpVector);
  void     updateLawCovMatrix(const M& newLawCovMatrix);
  const M& lawCovMatrix      () const;

  const V& lawExpVector() const;
  const V& lawVarVector() const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_normalizationStyle;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;
  V*       m_lawExpVector;
  V*       m_lawVarVector;
  bool     m_diagonalCovMatrix;
  const M* m_lawCovMatrix;
};

template<class V,class M>
uqGaussianJointPdfClass<V,M>::uqGaussianJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet),
  m_lawExpVector     (new V(lawExpVector)),
  m_lawVarVector     (new V(lawVarVector)),
  m_diagonalCovMatrix(true),
  m_lawCovMatrix     (m_domainSet.vectorSpace().newDiagMatrix(lawVarVector))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGaussianJointPdfClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "In uqGaussianJointPdfClass<V,M>::constructor()"
                          //<< ", prefix = "     << m_prefix
                            << ": meanVector = " << this->lawExpVector()
	                    << ", Variances = "  << this->lawVarVector()
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianJointPdfClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqGaussianJointPdfClass<V,M>::uqGaussianJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet),
  m_lawExpVector     (new V(lawExpVector)),
  m_lawVarVector     (domainSet.vectorSpace().newVector(INFINITY)), // FIX ME
  m_diagonalCovMatrix(false),
  m_lawCovMatrix     (new M(lawCovMatrix))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGaussianJointPdfClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "In uqGaussianJointPdfClass<V,M>::constructor()"
                          //<< ", prefix = "            << m_prefix
                            << ": meanVector = "        << this->lawExpVector()
	                    << ", Covariance Matrix = " << lawCovMatrix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianJointPdfClass<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqGaussianJointPdfClass<V,M>::~uqGaussianJointPdfClass()
{
  delete m_lawCovMatrix;
  delete m_lawVarVector;
  delete m_lawExpVector;
}

template <class V, class M>
const V&
uqGaussianJointPdfClass<V,M>::lawExpVector() const
{
  return *m_lawExpVector;
}

template <class V, class M>
const V&
uqGaussianJointPdfClass<V,M>::lawVarVector() const
{
  return *m_lawVarVector;
}

template<class V, class M>
double
uqGaussianJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering uqGaussianJointPdfClass<V,M>::actualValue()"
                            << ", meanVector = "   << *m_lawExpVector
	                    << ", lawCovMatrix = " << *m_lawCovMatrix
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqGaussianJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqGaussianJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = 0.;

  if (this->m_domainSet.contains(domainVector) == false) { // prudenci 2011-Oct-04
    returnValue = 0.;
  }
  else {
    returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
  }
  //returnValue *= exp(m_logOfNormalizationFactor); // No need, because 'lnValue()' is called right above [PDF-03]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianJointPdfClass<V,M>::actualValue()"
                            << ", meanVector = "   << *m_lawExpVector
	                    << ", lawCovMatrix = " << *m_lawCovMatrix
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqGaussianJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering uqGaussianJointPdfClass<V,M>::lnValue()"
                            << ", meanVector = "   << *m_lawExpVector
	                    << ", lawCovMatrix = " << *m_lawCovMatrix
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqGaussianJointPdfClass<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  if (domainDirection) {}; // just to remove compiler warning

  double returnValue = 0.;

  double lnDeterminant = 0.;
  if (this->m_domainSet.contains(domainVector) == false) { // prudenci 2011-Oct-04
    returnValue = -INFINITY;
  }
  else {
    V diffVec(domainVector - this->lawExpVector());
    if (m_diagonalCovMatrix) {
      returnValue = ((diffVec*diffVec)/this->lawVarVector()).sumOfComponents();
      if (m_normalizationStyle == 0) {
        unsigned int iMax = this->lawVarVector().sizeLocal();
        for (unsigned int i = 0; i < iMax; ++i) {
          lnDeterminant += log(this->lawVarVector()[i]);
        }
      }
    }
    else {
      V tmpVec = this->m_lawCovMatrix->invertMultiply(diffVec);
      returnValue = (diffVec*tmpVec).sumOfComponents();
      if (m_normalizationStyle == 0) {
        lnDeterminant = this->m_lawCovMatrix->lnDeterminant();
      }
    }
    if (m_normalizationStyle == 0) {
      returnValue += log(2*M_PI);   // normalization of pdf
      returnValue += lnDeterminant; // normalization of pdf
    }
    returnValue *= -0.5;
  }
  returnValue += m_logOfNormalizationFactor; // [PDF-03]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianJointPdfClass<V,M>::lnValue()"
                            << ", m_normalizationStyle = " << m_normalizationStyle
                            << ", meanVector = "           << *m_lawExpVector
	                    << ", lawCovMatrix = "         << *m_lawCovMatrix
                            << ": domainVector = "         << domainVector
                            << ", returnValue = "          << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqGaussianJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGaussianJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

template<class V, class M>
void
uqGaussianJointPdfClass<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // delete old expected values (alloced at construction or last call to this function)
  delete m_lawExpVector;
  m_lawExpVector = new V(newLawExpVector);
  return;
}

template<class V, class M>
void
uqGaussianJointPdfClass<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  // delete old expected values (alloced at construction or last call to this function)
  delete m_lawCovMatrix;
  m_lawCovMatrix = new M(newLawCovMatrix);
  return;
}

template<class V, class M>
const M&
uqGaussianJointPdfClass<V,M>::lawCovMatrix() const
{
  return *m_lawCovMatrix;
}

//*****************************************************
// Uniform probability density class [PDF-04]
//*****************************************************
template<class V, class M>
class uqUniformJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqUniformJointPdfClass(const char*                  prefix,
                         const uqVectorSetClass<V,M>& domainSet);
 ~uqUniformJointPdfClass();

  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_normalizationStyle;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;
};

template<class V,class M>
uqUniformJointPdfClass<V,M>::uqUniformJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),
                            domainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqUniformJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqUniformJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqUniformJointPdfClass<V,M>::~uqUniformJointPdfClass()
{
}

template<class V, class M>
double
uqUniformJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqUniformJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainDirection) {}; // just to remove compiler warning

  double volume = m_domainSet.volume();
  if (((boost::math::isnan)(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                ) ||
      (m_normalizationStyle != 0   )) {
    volume = 1.;
  }

  return 1./volume; // No need to multiply by exp(m_logOfNormalizationFactor) [PDF-04]
}

template<class V, class M>
double
uqUniformJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainVector[0]) {}; // just to remove compiler warning
  if (domainDirection) {}; // just to remove compiler warning

  double volume = m_domainSet.volume();
  if (((boost::math::isnan)(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                ) ||
      (m_normalizationStyle != 0   )) {
    volume = 1.;
  }

  return log(volume); // No need to add m_logOfNormalizationFactor [PDF-04]
}

template<class V, class M>
double
uqUniformJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqUniformJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqUniformJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// Beta probability density class [PDF-05]
//*****************************************************
template<class V, class M>
class uqBetaJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqBetaJointPdfClass(const char*                  prefix,
                      const uqVectorSetClass<V,M>& domainSet,
                      const V&                     alpha,
                      const V&                     beta);
 ~uqBetaJointPdfClass();

  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_normalizationStyle;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  V m_alpha;
  V m_beta;
};

template<class V,class M>
uqBetaJointPdfClass<V,M>::uqBetaJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),domainSet),
  m_alpha(alpha),
  m_beta (beta)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBetaJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBetaJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqBetaJointPdfClass<V,M>::~uqBetaJointPdfClass()
{
}

template<class V, class M>
double
uqBetaJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqBetaJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqBetaJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  // No need to multiply by exp(m_logOfNormalizationFactor) because 'lnValue()' is called [PDF-05]
  return exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}

template<class V, class M>
double
uqBetaJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqBetaJointPdfClass<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double aux = 0.;
  double result = 0.;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
    if (m_normalizationStyle == 0) {
      aux = log(gsl_ran_beta_pdf(domainVector[i],m_alpha[i],m_beta[i]));
    }
    else {
      aux = (m_alpha[i]-1.)*log(domainVector[i]) + (m_beta[i]-1.)*log(1.-domainVector[i]);
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqBetaJointPdfClass<V,M>::lnValue()"
                              << ", m_normalizationStyle = "      << m_normalizationStyle
                              << ": domainVector[" << i << "] = " << domainVector[i]
                              << ", m_alpha[" << i << "] = "      << m_alpha[i]
                              << ", m_beta[" << i << "] = "       << m_beta[i]
                              << ", log(pdf)= "                   << aux
                              << std::endl;
    }
    result += aux;
  }
  result += m_logOfNormalizationFactor; // [PDF-05]

  return result;
}

template<class V, class M>
double
uqBetaJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqBetaJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqBetaJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// Gamma probability density class [PDF-06]
//*****************************************************
template<class V, class M>
class uqGammaJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqGammaJointPdfClass(const char*                  prefix,
                       const uqVectorSetClass<V,M>& domainSet,
                       const V&                     a,
                       const V&                     b);
 ~uqGammaJointPdfClass();

  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_normalizationStyle;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  V m_a;
  V m_b;
};

template<class V,class M>
uqGammaJointPdfClass<V,M>::uqGammaJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     a,
  const V&                     b)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),domainSet),
  m_a(a),
  m_b(b)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGammaJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGammaJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqGammaJointPdfClass<V,M>::~uqGammaJointPdfClass()
{
}

template<class V, class M>
double
uqGammaJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqGammaJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqGammaJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  // No need to multiply by exp(m_logOfNormalizationFactor) because 'lnValue()' is called [PDF-06]
  return exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}

template<class V, class M>
double
uqGammaJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqGammaJointPdfClass<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double aux = 0.;
  double result = 0.;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
    if (m_normalizationStyle == 0) {
      aux = log(gsl_ran_gamma_pdf(domainVector[i],m_a[i],m_b[i]));
    }
    else {
      aux = (m_a[i]-1.)*log(domainVector[i]) - domainVector[i]/m_b[i];
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqGammaJointPdfClass<V,M>::lnValue()"
                              << ", m_normalizationStyle = "      << m_normalizationStyle
                              << ": domainVector[" << i << "] = " << domainVector[i]
                              << ", m_a[" << i << "] = "          << m_a[i]
                              << ", m_b[" << i << "] = "          << m_b[i]
                              << ", log(pdf)= "                   << aux
                              << std::endl;
    }
    result += aux;
  }
  result += m_logOfNormalizationFactor; // [PDF-06]

  return result;
}

template<class V, class M>
double
uqGammaJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGammaJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGammaJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// InverseGamma probability density class [PDF-07]
//*****************************************************
template<class V, class M>
class uqInverseGammaJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqInverseGammaJointPdfClass(const char*                  prefix,
                              const uqVectorSetClass<V,M>& domainSet,
                              const V&                     alpha,
                              const V&                     beta);
 ~uqInverseGammaJointPdfClass();

  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_normalizationStyle;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  V m_alpha;
  V m_beta;
};

template<class V,class M>
uqInverseGammaJointPdfClass<V,M>::uqInverseGammaJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     alpha,
  const V&                     beta)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),domainSet),
  m_alpha(alpha),
  m_beta (beta)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqInverseGammaJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqInverseGammaJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqInverseGammaJointPdfClass<V,M>::~uqInverseGammaJointPdfClass()
{
}

template<class V, class M>
double
uqInverseGammaJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqInverseGammaJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqInverseGammaJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  // No need to multiply by exp(m_logOfNormalizationFactor) because 'lnValue()' is called [PDF-07]
  return exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}

template<class V, class M>
double
uqInverseGammaJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqInverseGammaJointPdfClass<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double result = 0.;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
    result -= (m_alpha[i]+1.)*log(domainVector[i]);
    result -= m_beta[i]/domainVector[i];
    if (m_normalizationStyle == 0) {
      // Code needs to be done yet
    }
  }
  result += m_logOfNormalizationFactor; // [PDF-07]

  return result;
}

template<class V, class M>
double
uqInverseGammaJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqInverseGammaJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqInverseGammaJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// Powered probability density class [PDF-08]
//*****************************************************
template<class V, class M>
class uqPoweredJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqPoweredJointPdfClass(const char*                     prefix,
                         const uqBaseJointPdfClass<V,M>& srcDensity,
                               double                    exponent);
 ~uqPoweredJointPdfClass();

  double actualValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue              (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  void   setNormalizationStyle(unsigned int value) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  const uqBaseJointPdfClass<V,M>& m_srcDensity;
        double                    m_exponent;
};

template<class V,class M>
uqPoweredJointPdfClass<V,M>::uqPoweredJointPdfClass(
  const char*                     prefix,
  const uqBaseJointPdfClass<V,M>& srcDensity,
        double                    exponent)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"pow").c_str(),srcDensity.domainSet()),
  m_srcDensity            (srcDensity),
  m_exponent              (exponent)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqPoweredJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In uqPoweredJointPdfClass<V,M>::constructor()"
                          //<< ", prefix = "     << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqPoweredJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqPoweredJointPdfClass<V,M>::~uqPoweredJointPdfClass()
{
}

template<class V,class M>
void
uqPoweredJointPdfClass<V,M>::setNormalizationStyle(unsigned int value) const
{
  m_srcDensity.setNormalizationStyle(value);
  return;
}

template<class V, class M>
double
uqPoweredJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqPoweredJointPdfClass<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqPoweredJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  double value = m_srcDensity.actualValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqPoweredJointPdfClass<V,M>::actualValue()",
                      "incomplete code for domainDirection, gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = pow(value,m_exponent);
  returnValue *= exp(m_logOfNormalizationFactor); // [PDF-08] ???

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqPoweredJointPdfClass<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqPoweredJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqPoweredJointPdfClass<V,M>::lnValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  double value = m_srcDensity.lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqPoweredJointPdfClass<V,M>::lnValue()",
                      "incomplete code for domainDirection, gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = m_exponent*value;
  returnValue += m_logOfNormalizationFactor; // [PDF-08] ???

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqPoweredJointPdfClass<V,M>::lnValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqPoweredJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  double volume = m_domainSet.volume();
  if (((boost::math::isnan)(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                )) {
    // Do nothing
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.worldRank(),
                        "uqPoweredJointPdfClass<V,M>::lnValue()",
                        "incomplete code for computeLogOfNormalizationFactor()");
  }

  return value;
}

//*****************************************************
// Wigner probability density class [PDF-09]
//*****************************************************
template<class V, class M>
class uqWignerJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqWignerJointPdfClass(const char*                  prefix,
                        const uqVectorSetClass<V,M>& domainSet,
                        const V&                     centerPos,
                        double                       radius);
 ~uqWignerJointPdfClass();

  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;
  V*     m_centerPos;
  double m_radius;
};

template<class V,class M>
uqWignerJointPdfClass<V,M>::uqWignerJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     centerPos,
  double                       radius)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"uni").c_str(),
			   domainSet),
  m_centerPos(new V(centerPos)),
  m_radius   (radius)    
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqWignerJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_radius <= 0.,
                      m_env.worldRank(),
                      "uqWignerJointPdfClass<V,M>::constructor()",
                      "invalid radius");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqWignerJointPdfClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqWignerJointPdfClass<V,M>::~uqWignerJointPdfClass()
{
  delete m_centerPos;
}

template<class V, class M>
double
uqWignerJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqWignerJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  double returnValue = 0.;
  double distanceRatio = (domainVector - *m_centerPos).norm2()/m_radius;
  if (distanceRatio < 1.) {
    returnValue = 2.*m_radius*m_radius*sqrt(1. - distanceRatio*distanceRatio)/M_PI;
  }
  returnValue *= exp(m_logOfNormalizationFactor); // [PDF-09]

  return returnValue;
}

template<class V, class M>
double
uqWignerJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  // No need to add m_logOfNormalizationFactor because 'actualValue()' is called [PDF-09]
  return log(this->actualValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}

template<class V, class M>
double
uqWignerJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqWignerJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqWignerJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// LogNormal probability density class [PDF-10]
//*****************************************************
template<class V, class M>
class uqLogNormalJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqLogNormalJointPdfClass(const char*                  prefix,
                           const uqVectorSetClass<V,M>& domainSet,
                           const V&                     lawExpVector,
                           const V&                     lawVarVector);
 ~uqLogNormalJointPdfClass();

  double   actualValue (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double   lnValue     (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double   computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

  const V& lawExpVector() const;
  const V& lawVarVector() const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_normalizationStyle;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;
  V*   m_lawExpVector;
  V*   m_lawVarVector;
  bool m_diagonalCovMatrix;
};

template<class V,class M>
uqLogNormalJointPdfClass<V,M>::uqLogNormalJointPdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet),
  m_lawExpVector     (new V(lawExpVector)),
  m_lawVarVector     (new V(lawVarVector)),
  m_diagonalCovMatrix(true)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalJointPdfClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "In uqLogNormalJointPdfClass<V,M>::constructor()"
                          //<< ", prefix = "     << m_prefix
                            << ": meanVector = " << this->lawExpVector()
	                    << ", Variances = "  << this->lawVarVector()
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalJointPdfClass<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqLogNormalJointPdfClass<V,M>::~uqLogNormalJointPdfClass()
{
  delete m_lawVarVector;
  delete m_lawExpVector;
}

template <class V, class M>
const V&
uqLogNormalJointPdfClass<V,M>::lawExpVector() const
{
  return *m_lawExpVector;
}

template <class V, class M>
const V&
uqLogNormalJointPdfClass<V,M>::lawVarVector() const
{
  return *m_lawVarVector;
}

template<class V, class M>
double
uqLogNormalJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalJointPdfClass<V,M>::actualValue()"
                            << ", meanVector = "               << *m_lawExpVector
                            << ": domainVector = "             << domainVector
                            << ", domainVector.sizeLocal() = " << domainVector.sizeLocal()
                            << ", this->m_domainSet.vectorSpace().dimLocal() = " << this->m_domainSet.vectorSpace().dimLocal()
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqLogNormalJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqLogNormalJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = 0.;

  V zeroVector(domainVector);
  zeroVector.cwSet(0.);
  if (domainVector.atLeastOneComponentSmallerOrEqualThan(zeroVector)) {
    returnValue = 0.;
  }
  else if (this->m_domainSet.contains(domainVector) == false) { // prudenci 2011-Oct-04
    returnValue = 0.;
  }
  else {
    returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
  }
  //returnValue *= exp(m_logOfNormalizationFactor); // No need, because 'lnValue()' is called right above // [PDF-10]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalJointPdfClass<V,M>::actualValue()"
                            << ", meanVector = "   << *m_lawExpVector
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqLogNormalJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalJointPdfClass<V,M>::lnValue()"
                            << ", meanVector = "   << *m_lawExpVector
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO((gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqLogNormalJointPdfClass<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  if (domainDirection) {}; // just to remove compiler warning

  double returnValue = 0.;

  V zeroVector(domainVector);
  zeroVector.cwSet(0.);
  if (domainVector.atLeastOneComponentSmallerOrEqualThan(zeroVector)) {
    returnValue = -INFINITY;
  }
  else if (this->m_domainSet.contains(domainVector) == false) { // prudenci 2011-Oct-04
    returnValue = -INFINITY;
  }
  else {
    if (m_diagonalCovMatrix) {
      V diffVec(zeroVector);
      for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
        diffVec[i] = std::log(domainVector[i]) - this->lawExpVector()[i];
      }
      returnValue = ((diffVec*diffVec)/this->lawVarVector()).sumOfComponents();
      returnValue *= -0.5;
      if (m_normalizationStyle == 0) {
        for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
          returnValue -= std::log(domainVector[i] * std::sqrt(2. * M_PI * this->lawVarVector()[i])); // Contribution of 1/(x\sqrt{2\pi\sigma^2})
        }
      }
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqLogNormalJointPdfClass<V,M>::lnValue()",
                          "situation with a non-diagonal covariance matrix makes no sense");
    }
    returnValue += m_logOfNormalizationFactor; // [PDF-10]
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalJointPdfClass<V,M>::lnValue()"
                            << ", meanVector = "   << *m_lawExpVector
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqLogNormalJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqLogNormalJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = uqBaseJointPdfClass<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqLogNormalJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

//*****************************************************
// Concatenated probability density class [PDF-11]
//*****************************************************
template<class V, class M>
class uqConcatenatedJointPdfClass : public uqBaseJointPdfClass<V,M> {
public:
  uqConcatenatedJointPdfClass(const char*                     prefix,
                              const uqBaseJointPdfClass<V,M>& density1,
                              const uqBaseJointPdfClass<V,M>& density2,
                              const uqVectorSetClass   <V,M>& concatenatedDomain); 
  uqConcatenatedJointPdfClass(const char*                                          prefix,
                              const std::vector<const uqBaseJointPdfClass<V,M>* >& densities,
                              const uqVectorSetClass<V,M>&                         concatenatedDomain); 
 ~uqConcatenatedJointPdfClass();

  double actualValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue              (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  void   setNormalizationStyle(unsigned int value) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;
  using uqBaseJointPdfClass<V,M>::m_logOfNormalizationFactor;

  std::vector<const uqBaseJointPdfClass<V,M>* > m_densities;
};

template<class V,class M>
uqConcatenatedJointPdfClass<V,M>::uqConcatenatedJointPdfClass(
  const char*                     prefix,
  const uqBaseJointPdfClass<V,M>& density1,
  const uqBaseJointPdfClass<V,M>& density2,
  const uqVectorSetClass   <V,M>& concatenatedDomain)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"concat").c_str(),concatenatedDomain),
  m_densities             (2,(const uqBaseJointPdfClass<V,M>*) NULL)
{
  m_densities[0] = &density1;
  m_densities[1] = &density2;

  unsigned int size1 = m_densities[0]->domainSet().vectorSpace().dimLocal();
  unsigned int size2 = m_densities[1]->domainSet().vectorSpace().dimLocal();
  unsigned int size  = concatenatedDomain.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO((size1+size2) != size,
                      m_env.worldRank(),
                      "uqConcatenatedJointPdfClass<V,M>::constructor(1)",
                      "incompatible dimensions");
}

template<class V,class M>
uqConcatenatedJointPdfClass<V,M>::uqConcatenatedJointPdfClass(
  const char*                                          prefix,
  const std::vector<const uqBaseJointPdfClass<V,M>* >& densities,
  const uqVectorSetClass<V,M>&                         concatenatedDomain)
  :
  uqBaseJointPdfClass<V,M>(((std::string)(prefix)+"concat").c_str(),concatenatedDomain),
  m_densities             (densities.size(),(const uqBaseJointPdfClass<V,M>*) NULL)
{
  unsigned int sumSizes = 0;
  for (unsigned i = 0; i < m_densities.size(); ++i) {
    m_densities[i] = densities[i];
    sumSizes += m_densities[i]->domainSet().vectorSpace().dimLocal();
  }

  unsigned int size  = concatenatedDomain.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(sumSizes != size,
                      m_env.worldRank(),
                      "uqConcatenatedJointPdfClass<V,M>::constructor(2)",
                      "incompatible dimensions");
}

template<class V,class M>
uqConcatenatedJointPdfClass<V,M>::~uqConcatenatedJointPdfClass()
{
}

template<class V,class M>
void
uqConcatenatedJointPdfClass<V,M>::setNormalizationStyle(unsigned int value) const
{
  for (unsigned i = 0; i < m_densities.size(); ++i) {
    m_densities[i]->setNormalizationStyle(value);
  }
  return;
}

template<class V, class M>
double
uqConcatenatedJointPdfClass<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedJointPdfClass<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqConcatenatedJointPdfClass<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqConcatenatedJointPdfClass<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  std::vector<V*> vecs(m_densities.size(),(V*) NULL);
  std::vector<double> values(m_densities.size(),0.);
  double returnValue = 1.;
  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vecs[i] = new V(m_densities[i]->domainSet().vectorSpace().zeroVector());
    domainVector.cwExtract(cummulativeSize,*(vecs[i]));
    values[i] = m_densities[i]->actualValue(*(vecs[i]),NULL,NULL,NULL,NULL);
    returnValue *= values[i];
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqConcatenatedJointPdfClass<V,M>::actualValue()"
                              << ", *(vecs[" << i << "]) = "       << *(vecs[i])
                              << ": values[" << i << "] = "        << values[i]
                              << ", temporary cumulative value = " << returnValue
                              << std::endl;
    }
    cummulativeSize += vecs[i]->sizeLocal();
    delete vecs[i];
  }
  //returnValue *= exp(m_logOfNormalizationFactor); // No need, because each PDF should be already normalized [PDF-11]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedJointPdfClass<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqConcatenatedJointPdfClass<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedJointPdfClass<V,M>::lnValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "uqConcatenatedJointPdfClass<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  std::vector<V*> vecs(m_densities.size(),(V*) NULL);
  std::vector<double> values(m_densities.size(),0.);
  double returnValue = 0.;
  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vecs[i] = new V(m_densities[i]->domainSet().vectorSpace().zeroVector());
    domainVector.cwExtract(cummulativeSize,*(vecs[i]));
    values[i] = m_densities[i]->lnValue(*(vecs[i]),NULL,NULL,NULL,NULL);
    returnValue += values[i];
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqConcatenatedJointPdfClass<V,M>::lnValue()"
                              << ", *(vecs[" << i << "]) = "       << *(vecs[i])
                              << ": values[" << i << "] = "        << values[i]
                              << ", temporary cumulative value = " << returnValue
                              << std::endl;
    }
    cummulativeSize += vecs[i]->sizeLocal();
    delete vecs[i];
  }
  //returnValue += m_logOfNormalizationFactor; // No need, because each PDF should be already normalized [PDF-11]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedJointPdfClass<V,M>::lnValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
uqConcatenatedJointPdfClass<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqConcatenatedJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  double volume = m_domainSet.volume();
  if (((boost::math::isnan)(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                )) {
    // Do nothing
  }
  else {
    for (unsigned int i = 0; i < m_densities.size(); ++i) {
      m_densities[i]->computeLogOfNormalizationFactor(numSamples, updateFactorInternally);
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqConcatenatedJointPdfClass<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

#endif // __UQ_JOINT_PROB_DENSITY_H__
