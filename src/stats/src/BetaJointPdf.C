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

#include <queso/BetaJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BasicPdfsBase.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
BetaJointPdf<V,M>::BetaJointPdf(
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
    *m_env.subDisplayFile() << "Entering BetaJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BetaJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
BetaJointPdf<V,M>::~BetaJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
BetaJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");
  queso_require_msg(!(domainDirection || hessianMatrix || hessianEffect), "incomplete code for hessianMatrix and hessianEffect calculations");

  // No need to multiply by exp(m_logOfNormalizationFactor) because 'lnValue()' is called
  double value = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));

  if (gradVector) {
    // DM: This transformation may have issues near the boundary.  I haven't
    // fully explored this yet.
    //
    // This evaluation is also normalised (if m_normalizationStyle is zero) and
    // so the gradient in physical space contains the normalisation constant.
    (*gradVector) *= value;
  }

  return value;
}
//--------------------------------------------------
template<class V, class M>
double
BetaJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  queso_require_msg(!(domainDirection || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double aux = 0.;
  double result = 0.;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
    if (m_normalizationStyle == 0) {
      aux = log(m_env.basicPdfs()->betaPdfActualValue(domainVector[i],m_alpha[i],m_beta[i]));
    }
    else {
      aux = (m_alpha[i]-1.)*log(domainVector[i]) + (m_beta[i]-1.)*log(1.-domainVector[i]);
    }

    // The log of the beta pdf is
    // f(x) = (alpha - 1) * log(x) + (beta - 1) * log(1 - x)
    // Therefore
    // df/dx (x) = ((alpha - 1) / x) + ((1 - beta) / (1 - x))
    if (gradVector) {
      // We're computing grad of log of p which is p' / p
      (*gradVector)[i] = (m_alpha[i] - 1.0) / domainVector[i] +
        (1.0 - m_beta[i]) / (1.0 - domainVector[i]);
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In BetaJointPdf<V,M>::lnValue()"
                              << ", m_normalizationStyle = "      << m_normalizationStyle
                              << ": domainVector[" << i << "] = " << domainVector[i]
                              << ", m_alpha[" << i << "] = "      << m_alpha[i]
                              << ", m_beta[" << i << "] = "       << m_beta[i]
                              << ", log(pdf)= "                   << aux
                              << std::endl;
    }
    result += aux;
  }
  result += m_logOfNormalizationFactor;

  return result;
}
//--------------------------------------------------
template<class V, class M>
void
BetaJointPdf<V,M>::distributionMean(V& meanVector) const
{
  unsigned int n_params = meanVector.sizeLocal();
  queso_assert_equal_to (n_params, m_alpha.sizeLocal());

  for (unsigned int i = 0; i < n_params; ++i) {
    meanVector[i] = m_alpha[i] / (m_alpha[i] + m_beta[i]);
  }
}
//--------------------------------------------------
template<class V, class M>
void
BetaJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  unsigned int n_params = m_alpha.sizeLocal();
  queso_assert_equal_to (n_params, m_beta.sizeLocal());
  queso_assert_equal_to (n_params, covMatrix.numCols());
  queso_assert_equal_to (covMatrix.numCols(), covMatrix.numRowsGlobal());

  covMatrix.zeroLower();
  covMatrix.zeroUpper();

  for (unsigned int i = 0; i < n_params; ++i) {
    covMatrix(i,i) = (m_alpha[i] * m_beta[i]) /
      ((m_alpha[i] + m_beta[i]) *
       (m_alpha[i] + m_beta[i]) *
       (m_alpha[i] + m_beta[i] + 1));
  }
}
//--------------------------------------------------
template<class V, class M>
double
BetaJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering BetaJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving BetaJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::BetaJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
