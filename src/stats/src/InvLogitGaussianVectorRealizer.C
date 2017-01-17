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

#include <cmath>
#include <limits>
#include <queso/InvLogitGaussianVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/math_macros.h>

namespace QUESO {

template<class V, class M>
InvLogitGaussianVectorRealizer<V, M>::InvLogitGaussianVectorRealizer(
    const char * prefix,
    const BoxSubset<V, M> & unifiedImageBoxSubset,
    const V & lawExpVector,
    const M & lowerCholLawCovMatrix)
  : BaseVectorRealizer<V, M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      unifiedImageBoxSubset, std::numeric_limits<unsigned int>::max()),
    m_unifiedLawExpVector(new V(lawExpVector)),
    m_unifiedLawVarVector(
        unifiedImageBoxSubset.vectorSpace().newVector(INFINITY)),  // FIX ME
    m_lowerCholLawCovMatrix(new M(lowerCholLawCovMatrix)),
    m_matU(NULL),
    m_vecSsqrt(NULL),
    m_matVt(NULL),
    m_unifiedImageBoxSubset(unifiedImageBoxSubset)
{
  *m_unifiedLawExpVector = lawExpVector;
}

template<class V, class M>
InvLogitGaussianVectorRealizer<V, M>::InvLogitGaussianVectorRealizer(
    const char * prefix,
    const BoxSubset<V, M> & unifiedImageBoxSubset,
    const V & lawExpVector,
    const M & matU,
    const V & vecSsqrt,
    const M & matVt)
  : BaseVectorRealizer<V, M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      unifiedImageBoxSubset, std::numeric_limits<unsigned int>::max()),
    m_unifiedLawExpVector(new V(lawExpVector)),
    m_unifiedLawVarVector(
        unifiedImageBoxSubset.vectorSpace().newVector( INFINITY)), // FIX ME
    m_lowerCholLawCovMatrix(NULL),
    m_matU(new M(matU)),
    m_vecSsqrt(new V(vecSsqrt)),
    m_matVt(new M(matVt)),
    m_unifiedImageBoxSubset(unifiedImageBoxSubset)
{
  *m_unifiedLawExpVector = lawExpVector; // ????
}

template<class V, class M>
InvLogitGaussianVectorRealizer<V, M>::~InvLogitGaussianVectorRealizer()
{
  delete m_matVt;
  delete m_vecSsqrt;
  delete m_matU;
  delete m_lowerCholLawCovMatrix;
  delete m_unifiedLawVarVector;
  delete m_unifiedLawExpVector;
}

template <class V, class M>
const V &
InvLogitGaussianVectorRealizer<V, M>::unifiedLawExpVector() const
{
  return *m_unifiedLawExpVector;
}

template <class V, class M>
const V &
InvLogitGaussianVectorRealizer<V, M>::unifiedLawVarVector() const
{
  return *m_unifiedLawVarVector;
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::realization(V & nextValues) const
{
  V iidGaussianVector(m_unifiedImageSet.vectorSpace().zeroVector());

  iidGaussianVector.cwSetGaussian(0.0, 1.0);

  if (m_lowerCholLawCovMatrix) {
    nextValues = (*m_unifiedLawExpVector) +
      (*m_lowerCholLawCovMatrix) * iidGaussianVector;
  }
  else if (m_matU && m_vecSsqrt && m_matVt) {
    nextValues = (*m_unifiedLawExpVector) +
      (*m_matU) * ((*m_vecSsqrt) * ((*m_matVt) * iidGaussianVector));
  }
  else {
    queso_error_msg("GaussianVectorRealizer<V,M>::realization() inconsistent internal state");
  }

  V min_domain_bounds(this->m_unifiedImageBoxSubset.minValues());
  V max_domain_bounds(this->m_unifiedImageBoxSubset.maxValues());

  for (unsigned int i = 0; i < nextValues.sizeLocal(); i++) {
    double temp = std::exp(nextValues[i]);
    double min_val = min_domain_bounds[i];
    double max_val = max_domain_bounds[i];

    if (queso_isfinite(min_val) &&
        queso_isfinite(max_val)) {
        // Left- and right-hand sides are finite.  Do full transform.
        nextValues[i] = (max_val * temp + min_val) / (1.0 + temp);
      }
    else if (queso_isfinite(min_val) &&
             !queso_isfinite(max_val)) {
      // Left-hand side finite, but right-hand side is not.
      // Do only left-hand transform.
      nextValues[i] = temp + min_val;
    }
    else if (!queso_isfinite(min_val) &&
             queso_isfinite(max_val)) {
      // Right-hand side is finite, but left-hand side is not.
      // Do only right-hand transform.
      nextValues[i] = (max_val * temp - 1.0) / temp;
    }
  }
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::updateLawExpVector(
    const V & newLawExpVector)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_unifiedLawExpVector;

  m_unifiedLawExpVector = new V(newLawExpVector);
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::updateLowerCholLawCovMatrix(
    const M & newLowerCholLawCovMatrix)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lowerCholLawCovMatrix;
  delete m_matU;
  delete m_vecSsqrt;
  delete m_matVt;

  m_lowerCholLawCovMatrix = new M(newLowerCholLawCovMatrix);
  m_matU                  = NULL;
  m_vecSsqrt              = NULL;
  m_matVt                 = NULL;
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::updateLowerCholLawCovMatrix(
  const M & matU,
  const V & vecSsqrt,
  const M & matVt)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_lowerCholLawCovMatrix;
  delete m_matU;
  delete m_vecSsqrt;
  delete m_matVt;

  m_lowerCholLawCovMatrix = NULL;
  m_matU                  = new M(matU);
  m_vecSsqrt              = new V(vecSsqrt);
  m_matVt                 = new M(matVt);
}

}  // End namespace QUESO

template class QUESO::InvLogitGaussianVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
