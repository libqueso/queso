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

#include <queso/InverseGammaJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
InverseGammaJointPdf<V,M>::InverseGammaJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  const V&                     alpha,
  const V&                     beta)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"uni").c_str(),domainSet),
  m_alpha(alpha),
  m_beta (beta)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering InverseGammaJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving InverseGammaJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
InverseGammaJointPdf<V,M>::~InverseGammaJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
InverseGammaJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");

  queso_require_msg(!(domainDirection || gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  // No need to multiply by exp(m_logOfNormalizationFactor) because 'lnValue()' is called [PDF-07]
  return exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}
//--------------------------------------------------
template<class V, class M>
double
InverseGammaJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  queso_require_msg(!(domainDirection || gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

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
//--------------------------------------------------
template<class V, class M>
void
InverseGammaJointPdf<V,M>::distributionMean(V& meanVector) const
{
  queso_assert_equal_to(m_alpha.sizeLocal(), m_beta.sizeLocal());
  queso_assert_equal_to(m_alpha.sizeLocal(), meanVector.sizeLocal());

  for (unsigned int i = 0; i < m_alpha.sizeLocal(); ++i) {
    queso_assert_greater(m_alpha[i], 1);
    meanVector[i] = m_beta[i] / (m_alpha[i] - 1);
  }
}
//--------------------------------------------------
template<class V, class M>
void
InverseGammaJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  queso_assert_equal_to(m_alpha.sizeLocal(), m_beta.sizeLocal());
  queso_assert_equal_to(m_alpha.sizeLocal(), covMatrix.numCols());
  queso_assert_equal_to (covMatrix.numCols(), covMatrix.numRowsGlobal());

  covMatrix.zeroLower();
  covMatrix.zeroUpper();

  for (unsigned int i = 0; i < m_alpha.sizeLocal(); ++i) {
    queso_assert_greater(m_alpha[i], 2);
    covMatrix(i,i) = m_beta[i]*m_beta[i] / (m_alpha[i] - 1) /
                     (m_alpha[i] - 1) / (m_alpha[i] - 2);
  }
}
//--------------------------------------------------
template<class V, class M>
double
InverseGammaJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering InverseGammaJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving InverseGammaJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::InverseGammaJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
