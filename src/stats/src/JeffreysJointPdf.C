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

#include <queso/JeffreysJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
JeffreysJointPdf<V,M>::JeffreysJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"jef").c_str(),
                            domainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering JeffreysJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving JeffreysJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
JeffreysJointPdf<V,M>::~JeffreysJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
JeffreysJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal()) {
    queso_error_msg("There is an invalid input while computing JeffreysJointPdf<V,M>::actualValue()");
  }

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainDirection) {}; // just to remove compiler warning

  double pdf = 1.0;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i){
    if (domainVector[i] < 0.0 ) {
      queso_error_msg("The domain for Jeffreys prior should be greater than zero.");
    }
    else if ((domainVector[i] == -INFINITY         ) ||
        (domainVector[i]      ==  INFINITY         ) ||
        (m_normalizationStyle != 0   )) {//TODO: R: not sure what this is doing?
      pdf = 0.0;
    }
  else {
    pdf = pdf * (1.0 / domainVector[i]);
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << " return pdf " << std::endl;
  }
  return pdf; // No need to multiply by exp(m_logOfNormalizationFactor) [PDF-04]
}

//--------------------------------------------------

template<class V, class M>
double
JeffreysJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.0;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainVector[0]) {}; // just to remove compiler warning
  if (domainDirection) {}; // just to remove compiler warning

  double pdf = 1.0;
  double result = 0.;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i){
    if (domainVector[i] < 0.0) {
      queso_error_msg("The domain for Jeffreys prior should be greater than zero.");
    }
    else if ((domainVector[i] == -INFINITY         ) ||
        (domainVector[i]      ==  INFINITY         ) ||
        (m_normalizationStyle != 0   )) {//TODO: R: not sure what this is doing?
      pdf = 0.0;
      result = -INFINITY; //TODO: what do we do here?
    }
  else {
    pdf = pdf * (1.0 / domainVector[i]);
    result = std::log(pdf);
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << " return log(pdf) " << std::endl;
  }
  return result; // No need to add m_logOfNormalizationFactor [PDF-04]
}
//--------------------------------------------------
template<class V, class M>
void
JeffreysJointPdf<V,M>::distributionMean(V& meanVector) const
{
  // FIXME: This is an improper prior; the "mean" makes little sense.
  // At least we can return something within the domain set.
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Warning: JeffreysJointPdf<V,M>::distributionMean() makes little sense"
                            << std::endl;
  }
  m_domainSet.centroid(meanVector);
}
//--------------------------------------------------
template<class V, class M>
void
JeffreysJointPdf<V,M>::distributionVariance (M & covMatrix) const
{
  // There's no way this is anything like well-defined
  queso_not_implemented();
}
//--------------------------------------------------
template<class V, class M>
double
JeffreysJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering JeffreysJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving JeffreysJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << " return normalization factor " << std::endl;
  }
  return value;
}

}  // End namespace QUESO

template class QUESO::JeffreysJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
