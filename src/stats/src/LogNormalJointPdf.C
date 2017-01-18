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

  queso_require_msg(!(hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = 0.;

  V zeroVector(domainVector);
  zeroVector.cwSet(0.);
  if (domainVector.atLeastOneComponentSmallerOrEqualThan(zeroVector)) {
    // What should the gradient be here?
    returnValue = 0.;
  }
  else if (this->m_domainSet.contains(domainVector) == false) {
    // What should the gradient be here?
    returnValue = 0.;
  }
  else {
    // Already normalised
    returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));

    if (gradVector) {
      (*gradVector) *= returnValue;
    }
  }

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

  queso_require_msg(!(domainDirection || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = 0.;

  V zeroVector(domainVector);
  zeroVector.cwSet(0.);
  if (domainVector.atLeastOneComponentSmallerOrEqualThan(zeroVector)) {
    // What should the gradient be here?
    returnValue = -INFINITY;
  }
  else if (this->m_domainSet.contains(domainVector) == false) {
    // What should the gradient be here?
    returnValue = -INFINITY;
  }
  else {
    if (m_diagonalCovMatrix) {
      V diffVec(zeroVector);
      for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
        diffVec[i] = std::log(domainVector[i]) - this->lawExpVector()[i];

        // Compute the gradient of log of the PDF
        // The log of a log normal pdf is:
        // f(x) = -log(x \sigma sqrt(2 \pi)) - ((log(x) - \mu)^2 / (2 \sigma^2))
        // Therefore
        // \frac{df}{dx}(x) = -1/x - (log(x) - \mu) / (x \sigma^2)
        if (gradVector) {
          (*gradVector)[i] = -(1.0 / domainVector[i]) -
            diffVec[i] / (domainVector[i] * this->lawVarVector()[i]);
        }
      }
      returnValue = ((diffVec*diffVec)/this->lawVarVector()).sumOfComponents();
      returnValue *= -0.5;

      for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
        returnValue -= std::log(domainVector[i]);
        if (m_normalizationStyle == 0) {
          returnValue -= std::log(std::sqrt(2. * M_PI * this->lawVarVector()[i])); // Contribution of 1/(x\sqrt{2\pi\sigma^2})
        }
      }
    }
    else {
      queso_error_msg("situation with a non-diagonal covariance matrix makes no sense");
    }
    returnValue += m_logOfNormalizationFactor;
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
void
LogNormalJointPdf<V,M>::distributionMean(V& meanVector) const
{
  // FIXME - this is the mean of a non-truncated lognormal
  // distribution, and doesn't take into account a limited domainSet.

  if (m_diagonalCovMatrix) {
    unsigned int n_params = meanVector.sizeLocal();
    queso_assert_equal_to (n_params, this->lawExpVector().sizeLocal());

    for (unsigned int i = 0; i < n_params; ++i) {
      meanVector[i] = std::exp(this->lawExpVector()[i] + 0.5*this->lawVarVector()[i]);
    }
  }
  else {
    queso_error_msg("situation with a non-diagonal covariance matrix makes no sense");
  }
}
//--------------------------------------------------
template<class V, class M>
void
LogNormalJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  // FIXME - this is the variance of a non-truncated lognormal
  // distribution, and doesn't take into account a limited domainSet.

  if (m_diagonalCovMatrix) {
    unsigned int n_params = this->lawExpVector().sizeLocal();
    queso_assert_equal_to (n_params, this->lawVarVector().sizeLocal());
    queso_assert_equal_to (n_params, covMatrix.numCols());
    queso_assert_equal_to (covMatrix.numCols(), covMatrix.numRowsGlobal());

    covMatrix.zeroLower();
    covMatrix.zeroUpper();

    for (unsigned int i = 0; i < n_params; ++i) {
      covMatrix(i,i) = (std::exp(this->lawVarVector()[i]) - 1) *
                       std::exp(2*this->lawExpVector()[i] + this->lawVarVector()[i]);
    }
  }
  else {
    queso_error_msg("situation with a non-diagonal covariance matrix makes no sense");
  }
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
