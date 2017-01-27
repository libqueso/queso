//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <queso/BayesianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V,class M>
BayesianJointPdf<V,M>::BayesianJointPdf(
  const char*                           prefix,
  const BaseJointPdf     <V,M>&  priorDensity,
  const BaseScalarFunction<V,M>& likelihoodFunction,
        double                          likelihoodExponent,
  const VectorSet         <V,M>& intersectionDomain)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"bay").c_str(),intersectionDomain),
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
// Destructor ---------------------------------------
template<class V,class M>
BayesianJointPdf<V,M>::~BayesianJointPdf()
{
  delete m_tmpMatrix;
}
// Math methods -------------------------------------
template<class V,class M>
void
BayesianJointPdf<V,M>::setNormalizationStyle(unsigned int value) const
{
  m_priorDensity.setNormalizationStyle(value);
  return;
}
// --------------------------------------------------
template<class V,class M>
double
BayesianJointPdf<V,M>::lastComputedLogPrior() const
{
  return m_lastComputedLogPrior;
}
// --------------------------------------------------
template<class V,class M>
double
BayesianJointPdf<V,M>::lastComputedLogLikelihood() const
{
  return m_lastComputedLogLikelihood;
}
// --------------------------------------------------
template<class V, class M>
double
BayesianJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BayesianJointPdf<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");

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

  queso_require_msg(!(gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

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
    *m_env.subDisplayFile() << "Leaving BayesianJointPdf<V,M>::actualValue()"
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}

template<class V, class M>
double
BayesianJointPdf<V,M>::lnValue(const V & domainVector) const
{
  double value1 = m_priorDensity.lnValue(domainVector);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In BayesianJointPdf<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": lnPrior = " << value1
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In BayesianJointPdf<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": about to call likelihood()"
                            << std::endl;
  }

  double value2 = 0.;
  if (m_likelihoodExponent != 0.) {
    value2 = m_likelihoodFunction.lnValue(domainVector);
  }

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

  return returnValue;
}

template<class V, class M>
double
BayesianJointPdf<V,M>::lnValue(const V & domainVector, V & gradVector) const
{
  double value1 = m_priorDensity.lnValue(domainVector, gradVector);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In BayesianJointPdf<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": lnPrior = " << value1
                            << std::endl;
  }

  double value2 = 0.;
  if (m_likelihoodExponent != 0.) {
    value2 = m_likelihoodFunction.lnValue(domainVector, m_tmpVector1);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "In BayesianJointPdf<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": value1 = "       << value1
                            << ", value2 = "       << value2
                            << std::endl;
    *m_env.subDisplayFile() << "In BayesianJointPdf<V,M>::lnValue()"
                            << ", domainVector = " << domainVector
                            << ": gradVector = "   << gradVector
                            << ", gradVLike = "    << m_tmpVector1
                            << std::endl;
  }

  gradVector += m_tmpVector1;

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

  return returnValue;
}

// --------------------------------------------------
template<class V, class M>
double
BayesianJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  double volume = m_domainSet.volume();
  if ((queso_isnan(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                )) {
    // Do nothing
  }
  else {
    queso_error_msg("incomplete code for computeLogOfNormalizationFactor()");
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::BayesianJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
