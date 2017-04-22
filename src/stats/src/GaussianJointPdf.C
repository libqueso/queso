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

#include <queso/GaussianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
GaussianJointPdf<V,M>::GaussianJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet),
  m_lawExpVector     (new V(lawExpVector)),
  m_lawVarVector     (new V(lawVarVector)),
  m_diagonalCovMatrix(true),
  m_lawCovMatrix     (m_domainSet.vectorSpace().newDiagMatrix(lawVarVector))
{

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GaussianJointPdf<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "In GaussianJointPdf<V,M>::constructor()"
                          //<< ", prefix = "     << m_prefix
                            << ": meanVector = " << this->lawExpVector()
                      << ", Variances = "  << this->lawVarVector()
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GaussianJointPdf<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V,class M>
GaussianJointPdf<V,M>::GaussianJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"gau").c_str(),domainSet),
  m_lawExpVector     (new V(lawExpVector)),
  m_lawVarVector     (domainSet.vectorSpace().newVector(INFINITY)), // FIX ME
  m_diagonalCovMatrix(false),
  m_lawCovMatrix     (new M(lawCovMatrix))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GaussianJointPdf<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "In GaussianJointPdf<V,M>::constructor()"
                          //<< ", prefix = "            << m_prefix
                            << ": meanVector = "        << this->lawExpVector()
                      << ", Covariance Matrix = " << lawCovMatrix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GaussianJointPdf<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
GaussianJointPdf<V,M>::~GaussianJointPdf()
{
  delete m_lawCovMatrix;
  delete m_lawVarVector;
  delete m_lawExpVector;
}
// Math methods-------------------------------------
template <class V, class M>
const V&
GaussianJointPdf<V,M>::lawExpVector() const
{
  return *m_lawExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
GaussianJointPdf<V,M>::lawVarVector() const
{
  return *m_lawVarVector;
}

template <class V, class M>
void
GaussianJointPdf<V, M>::print(std::ostream & os) const
{
  // Print m_env?

  os << "Start printing GaussianJointPdf<V, M>" << std::endl;
  os << "m_prefix:" << std::endl;
  os << this->m_prefix << std::endl;
  os << "m_domainSet:" << std::endl;
  os << this->m_domainSet << std::endl;
  os << "m_normalizationStyle:" << std::endl;
  os << this->m_normalizationStyle << std::endl;
  os << "m_logOfNormalizationFactor:" << std::endl;
  os << this->m_logOfNormalizationFactor << std::endl;
  os << "Mean:" << std::endl;
  os << this->lawExpVector() << std::endl;
  os << "Variance vector:" << std::endl;
  os << this->lawVarVector() << std::endl;
  os << "Covariance matrix:" << std::endl;
  os << this->lawCovMatrix() << std::endl;
  os << "Diagonal covariance?" << std::endl;
  os << this->m_diagonalCovMatrix << std::endl;
  os << "End printing GaussianJointPdf<V, M>" << std::endl;
}

//--------------------------------------------------
template<class V, class M>
double
GaussianJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering GaussianJointPdf<V,M>::actualValue()"
                            << ", meanVector = "   << *m_lawExpVector
                      << ", lawCovMatrix = " << *m_lawCovMatrix
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");

  queso_require_msg(!(hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = 0.;

  if (this->m_domainSet.contains(domainVector) == false) {
    // What should the gradient be here?
    returnValue = 0.;
  }
  else {
    // Already normalised (so the gradient will have the normalisation constant
    // in it)
    returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));

    if (gradVector) {
      (*gradVector) *= returnValue;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Leaving GaussianJointPdf<V,M>::actualValue()"
                            << ", meanVector = "   << *m_lawExpVector
                      << ", lawCovMatrix = " << *m_lawCovMatrix
                            << ": domainVector = " << domainVector
                            << ", returnValue = "  << returnValue
                            << std::endl;
  }

  return returnValue;
}
//--------------------------------------------------
template<class V, class M>
double
GaussianJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 55)) {
    *m_env.subDisplayFile() << "Entering GaussianJointPdf<V,M>::lnValue()"
                            << ", meanVector = "   << *m_lawExpVector
                      << ", lawCovMatrix = " << *m_lawCovMatrix
                            << ": domainVector = " << domainVector
                            << std::endl;
  }

  queso_require_msg(!(domainDirection || hessianMatrix || hessianEffect), "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  double returnValue = 0.;

  double lnDeterminant = 0.;
  if (this->m_domainSet.contains(domainVector) == false) {
    // What should the gradient be here?
    returnValue = -INFINITY;
  }
  else {
    V diffVec(domainVector - this->lawExpVector());
    if (m_diagonalCovMatrix) {
      returnValue = ((diffVec*diffVec)/this->lawVarVector()).sumOfComponents();

      // Compute the gradient of log of the pdf.
      // The log of a Gaussian pdf is:
      // f(x) = - 1/2 (x - \mu)^T \Sigma^{-1} (x - \mu)
      // Therefore
      // \frac{df}{dx}(x) = - (x - \mu)^T \Sigma^{-1}
      //                  = - (\Sigma^{-1}^T (x - \mu))^T
      //                  = - (\Sigma^{-1} (x - \mu))^T
      //                  = - \Sigma^{-1} (x - \mu)  (row/column vector doesn't matter)
      //
      // So if \Sigma is diagonal we have a component-wise product of two
      // vectors (x - \mu) and the diagonal elements of \Sigma^{-1}
      if (gradVector) {
        (*gradVector) = diffVec;  // Copy
        (*gradVector) /= this->lawVarVector();
        (*gradVector) *= -1.0;
      }

      if (m_normalizationStyle == 0) {
        unsigned int iMax = this->lawVarVector().sizeLocal();
        for (unsigned int i = 0; i < iMax; ++i) {
          lnDeterminant += std::log(this->lawVarVector()[i]);
        }
      }
    }
    else {
      V tmpVec = this->m_lawCovMatrix->invertMultiply(diffVec);
      returnValue = (diffVec*tmpVec).sumOfComponents();

      // Compute the gradient of log of the pdf.
      // The log of a Gaussian pdf is:
      // f(x) = - 1/2 (x - \mu)^T \Sigma^{-1} (x - \mu)
      // Therefore
      // \frac{df}{dx}(x) = - (x - \mu)^T \Sigma^{-1}
      //                  = - (\Sigma^{-1}^T (x - \mu))^T
      //                  = - (\Sigma^{-1} (x - \mu))^T
      //                  = - \Sigma^{-1} (x - \mu)  (row/column vector doesn't matter)
      if (gradVector) {
        (*gradVector) = tmpVec;
        (*gradVector) *= -1.0;
      }

      if (m_normalizationStyle == 0) {
        lnDeterminant = this->m_lawCovMatrix->lnDeterminant();
      }
    }
    if (m_normalizationStyle == 0) {
      returnValue += ((double) this->lawVarVector().sizeLocal()) * std::log(2*M_PI);   // normalization of pdf
      returnValue += lnDeterminant; // normalization of pdf
    }
    returnValue *= -0.5;
  }
  returnValue += m_logOfNormalizationFactor;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "Leaving GaussianJointPdf<V,M>::lnValue()"
                            << ", m_normalizationStyle = " << m_normalizationStyle
                            << ", m_diagonalCovMatrix = " << m_diagonalCovMatrix
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << ", lnDeterminant = " << lnDeterminant
                            << ", meanVector = "           << *m_lawExpVector
                            << ", lawCovMatrix = "         << *m_lawCovMatrix
                            << ": domainVector = "         << domainVector
                            << ", returnValue = "          << returnValue
                            << std::endl;
  }

  return returnValue;
}
//--------------------------------------------------
template<class V, class M>
void
GaussianJointPdf<V,M>::distributionMean(V& meanVector) const
{
  meanVector = this->lawExpVector();
}
//--------------------------------------------------
template<class V, class M>
void
GaussianJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  queso_assert_equal_to (covMatrix.numCols(), covMatrix.numRowsGlobal());

  if (m_diagonalCovMatrix) {
    covMatrix.zeroLower();
    covMatrix.zeroUpper();

    unsigned int n_comp = this->lawVarVector().sizeLocal();
    queso_assert_equal_to (n_comp, covMatrix.numCols());

    for (unsigned int i = 0; i < n_comp; ++i) {
      covMatrix(i,i) = this->lawVarVector()[i];
    }
  } else {
    covMatrix = *this->m_lawCovMatrix;
  }
}
//--------------------------------------------------
template<class V, class M>
double
GaussianJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GaussianJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GaussianJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}
//--------------------------------------------------
template<class V, class M>
void
GaussianJointPdf<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lawExpVector;
  m_lawExpVector = new V(newLawExpVector);
  return;
}

template<class V, class M>
void
GaussianJointPdf<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lawCovMatrix;
  m_lawCovMatrix = new M(newLawCovMatrix);
  return;
}

template<class V, class M>
const M&
GaussianJointPdf<V,M>::lawCovMatrix() const
{
  return *m_lawCovMatrix;
}

}  // End namespace QUESO

template class QUESO::GaussianJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
