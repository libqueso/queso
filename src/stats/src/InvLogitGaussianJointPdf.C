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

#include <queso/InvLogitGaussianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
InvLogitGaussianJointPdf<V,M>::InvLogitGaussianJointPdf(
  const char*                  prefix,
  const BoxSubset<V,M>& domainBoxSubset,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      domainBoxSubset),
  m_lawExpVector(new V(lawExpVector)),
  m_lawVarVector(new V(lawVarVector)),
  m_diagonalCovMatrix(true),
  m_lawCovMatrix(m_domainSet.vectorSpace().newDiagMatrix(lawVarVector)),
  m_domainBoxSubset(domainBoxSubset)
{
}


template<class V,class M>
InvLogitGaussianJointPdf<V,M>::InvLogitGaussianJointPdf(
  const char*                  prefix,
  const BoxSubset<V,M>& domainBoxSubset,
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      domainBoxSubset),
  m_lawExpVector(new V(lawExpVector)),
  m_lawVarVector(domainBoxSubset.vectorSpace().newVector(INFINITY)), // FIX ME
  m_diagonalCovMatrix(false),
  m_lawCovMatrix(new M(lawCovMatrix)),
  m_domainBoxSubset(domainBoxSubset)
{
}

template<class V,class M>
InvLogitGaussianJointPdf<V,M>::~InvLogitGaussianJointPdf()
{
  delete m_lawCovMatrix;
  delete m_lawVarVector;
  delete m_lawExpVector;
}

template <class V, class M>
const V&
InvLogitGaussianJointPdf<V,M>::lawExpVector() const
{
  return *m_lawExpVector;
}

template <class V, class M>
const V&
InvLogitGaussianJointPdf<V,M>::lawVarVector() const
{
  return *m_lawVarVector;
}

template <class V, class M>
void
InvLogitGaussianJointPdf<V, M>::print(std::ostream & os) const
{
  // Print m_env?

  os << "Start printing InvLogitGaussianJointPdf<V, M>" << std::endl;
  os << "m_prefix:" << std::endl;
  os << this->m_prefix << std::endl;
  os << "m_domainSet:" << std::endl;
  os << this->m_domainBoxSubset << std::endl;
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
  os << "End printing InvLogitGaussianJointPdf<V, M>" << std::endl;
}

template<class V, class M>
double
InvLogitGaussianJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  double returnValue;

  returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));

  return returnValue;
}

template<class V, class M>
double
InvLogitGaussianJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  double returnValue;
  double lnDeterminant = 0.0;
  V transformedDomainVector(domainVector);

  V min_domain_bounds(this->m_domainBoxSubset.minValues());
  V max_domain_bounds(this->m_domainBoxSubset.maxValues());

  double lnjacobian = 0.0;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); i++) {
    double min_val = min_domain_bounds[i];
    double max_val = max_domain_bounds[i];

    if (queso_isfinite(min_val) &&
        queso_isfinite(max_val)) {

      if (domainVector[i] == min_val || domainVector[i] == max_val) {
        // Exit early if we can
        return -INFINITY;
      }

        // Left- and right-hand sides are finite.  Do full transform.
        transformedDomainVector[i] = std::log(domainVector[i] - min_val) -
            std::log(max_val - domainVector[i]);

        lnjacobian += std::log(max_val - min_val) -
          std::log(domainVector[i] - min_val) -
          std::log(max_val - domainVector[i]);
    }
    else if (queso_isfinite(min_val) &&
             !queso_isfinite(max_val)) {

      if (domainVector[i] == min_val) {
        // Exit early if we can
        return -INFINITY;
      }

      // Left-hand side finite, but right-hand side is not.
      // Do only left-hand transform.
      transformedDomainVector[i] = std::log(domainVector[i] - min_val);

      lnjacobian += -std::log(domainVector[i] - min_val);
    }
    else if (!queso_isfinite(min_val) &&
             queso_isfinite(max_val)) {

      if (domainVector[i] == max_val) {
        // Exit early if we can
        return -INFINITY;
      }

      // Right-hand side is finite, but left-hand side is not.
      // Do only right-hand transform.
      transformedDomainVector[i] = -std::log(max_val - domainVector[i]);

      lnjacobian += -std::log(max_val - domainVector[i]);
    }
    else {
      // No transform.
      transformedDomainVector[i] = domainVector[i];
    }
  }

  V diffVec(transformedDomainVector - this->lawExpVector());
  if (m_diagonalCovMatrix) {
    returnValue = ((diffVec * diffVec) /
        this->lawVarVector()).sumOfComponents();
    if (m_normalizationStyle == 0) {
      unsigned int iMax = this->lawVarVector().sizeLocal();
      for (unsigned int i = 0; i < iMax; ++i) {
        lnDeterminant += log(this->lawVarVector()[i]);
      }
    }
  }
  else {
    V tmpVec = this->m_lawCovMatrix->invertMultiply(diffVec);
    returnValue = (diffVec * tmpVec).sumOfComponents();
    if (m_normalizationStyle == 0) {
      lnDeterminant = this->m_lawCovMatrix->lnDeterminant();
    }
  }
  if (m_normalizationStyle == 0) {
    returnValue += ((double) this->lawVarVector().sizeLocal()) * log(2 * M_PI);
    returnValue += lnDeterminant;
  }
  returnValue *= -0.5;
  returnValue += m_logOfNormalizationFactor;
  returnValue += lnjacobian;

  return returnValue;
}

template<class V, class M>
void
InvLogitGaussianJointPdf<V,M>::distributionMean(V& meanVector) const
{
  // AFAIK there's no simple closed form here, and taking the inverse
  // transformation of the mean in the transformed space probably
  // isn't accurate enough in cases where the mean is too near the
  // bounds.
  queso_not_implemented();
}

//---------------------------------------------------
template<class V,class M>
void
InvLogitGaussianJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  // AFAIK there's no simple closed form here, and taking the inverse
  // transformation of the variance in the transformed space probably
  // isn't accurate enough in cases where the mean is too near the
  // bounds.
  queso_not_implemented();
}


template<class V, class M>
double
InvLogitGaussianJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
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

template<class V, class M>
void
InvLogitGaussianJointPdf<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lawExpVector;
  m_lawExpVector = new V(newLawExpVector);
}

template<class V, class M>
void
InvLogitGaussianJointPdf<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lawCovMatrix;
  m_lawCovMatrix = new M(newLawCovMatrix);
}

template<class V, class M>
const M&
InvLogitGaussianJointPdf<V,M>::lawCovMatrix() const
{
  return *m_lawCovMatrix;
}

}  // End namespace QUESO

template class QUESO::InvLogitGaussianJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
