//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/LogNormalJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
LogNormalJointPdf<V,M>::LogNormalJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet),
  m_lawExpVector     (new V(lawExpVector)),
  m_lawVarVector     (new V(lawVarVector)),
  m_diagonalCovMatrix(true)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering LogNormalJointPdf<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "In LogNormalJointPdf<V,M>::constructor()"
                          //<< ", prefix = "     << m_prefix
                            << ": meanVector = " << this->lawExpVector()
                      << ", Variances = "  << this->lawVarVector()
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving LogNormalJointPdf<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
LogNormalJointPdf<V,M>::~LogNormalJointPdf()
{
  delete m_lawVarVector;
  delete m_lawExpVector;
}
// Math methods-------------------------------------
template <class V, class M>
const V&
LogNormalJointPdf<V,M>::lawExpVector() const
{
  return *m_lawExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
LogNormalJointPdf<V,M>::lawVarVector() const
{
  return *m_lawVarVector;
}
//--------------------------------------------------
template<class V, class M>
double
LogNormalJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering LogNormalJointPdf<V,M>::actualValue()"
                            << ", meanVector = "               << *m_lawExpVector
                            << ": domainVector = "             << domainVector
                            << ", domainVector.sizeLocal() = " << domainVector.sizeLocal()
                            << ", this->m_domainSet.vectorSpace().dimLocal() = " << this->m_domainSet.vectorSpace().dimLocal()
                            << std::endl;
  }

  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");

  queso_require_msg(!(gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

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
    *m_env.subDisplayFile() << "Leaving LogNormalJointPdf<V,M>::actualValue()"
                            << ", meanVector = "   << *m_lawExpVector
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}
//--------------------------------------------------
template<class V, class M>
double
LogNormalJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering LogNormalJointPdf<V,M>::lnValue()"
                            << ", meanVector = "   << *m_lawExpVector
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  queso_require_msg(!(gradVector || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

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
      queso_error_msg("situation with a non-diagonal covariance matrix makes no sense");
    }
    returnValue += m_logOfNormalizationFactor; // [PDF-10]
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Leaving LogNormalJointPdf<V,M>::lnValue()"
                            << ", meanVector = "   << *m_lawExpVector
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}
//--------------------------------------------------
template<class V, class M>
double
LogNormalJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering LogNormalJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving LogNormalJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::LogNormalJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
