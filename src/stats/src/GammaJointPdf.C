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

#include <queso/GammaJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BasicPdfsBase.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
GammaJointPdf<V,M>::GammaJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  const V&                     a,
  const V&                     b)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"uni").c_str(),domainSet),
  m_a(a),
  m_b(b)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GammaJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GammaJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
GammaJointPdf<V,M>::~GammaJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
GammaJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");

  queso_require_msg(!(domainDirection || gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  // No need to multiply by exp(m_logOfNormalizationFactor) because 'lnValue()' is called [PDF-06]
  return exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}
//--------------------------------------------------
template<class V, class M>
double
GammaJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  queso_require_msg(!(domainDirection || gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double aux = 0.;
  double result = 0.;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
    if (m_normalizationStyle == 0) {
      aux = log(m_env.basicPdfs()->gammaPdfActualValue(domainVector[i],m_a[i],m_b[i]));
    }
    else {
      aux = (m_a[i]-1.)*log(domainVector[i]) - domainVector[i]/m_b[i];
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GammaJointPdf<V,M>::lnValue()"
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
//--------------------------------------------------
template<class V, class M>
void
GammaJointPdf<V,M>::distributionMean(V& meanVector) const
{
  queso_assert_equal_to(m_a.sizeLocal(), m_b.sizeLocal());
  queso_assert_equal_to(m_a.sizeLocal(), meanVector.sizeLocal());

  for (unsigned int i = 0; i < m_a.sizeLocal(); ++i) {
    meanVector[i] = m_a[i] * m_b[i];
  }
}
//--------------------------------------------------
template<class V, class M>
void
GammaJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  queso_assert_equal_to(m_a.sizeLocal(), m_b.sizeLocal());
  queso_assert_equal_to(m_a.sizeLocal(), covMatrix.numRowsGlobal());
  queso_assert_equal_to(covMatrix.numCols(), covMatrix.numRowsGlobal());

  covMatrix.zeroLower();
  covMatrix.zeroUpper();

  for (unsigned int i = 0; i < m_a.sizeLocal(); ++i) {
    covMatrix(i,i) = m_a[i] * m_b[i] * m_b[i];
  }
}
//--------------------------------------------------
template<class V, class M>
double
GammaJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GammaJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GammaJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::GammaJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
