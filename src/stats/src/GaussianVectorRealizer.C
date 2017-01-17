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

#include <limits>
#include <queso/GaussianVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
GaussianVectorRealizer<V,M>::GaussianVectorRealizer(const char* prefix,
                  const VectorSet<V,M>& unifiedImageSet,
                  const V& lawExpVector,
                  const M& lowerCholLawCovMatrix)
  :
  BaseVectorRealizer<V,M>( ((std::string)(prefix)+"gau").c_str(), unifiedImageSet, std::numeric_limits<unsigned int>::max()), // 2011/Oct/02 - Correction thanks to Corey
  m_unifiedLawExpVector  (new V(lawExpVector)),
  m_unifiedLawVarVector  (unifiedImageSet.vectorSpace().newVector( INFINITY)), // FIX ME
  m_lowerCholLawCovMatrix(new M(lowerCholLawCovMatrix)),
  m_matU                 (NULL),
  m_vecSsqrt             (NULL),
  m_matVt                (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorRealizer<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  *m_unifiedLawExpVector = lawExpVector; // ????

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorRealizer<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class V, class M>
GaussianVectorRealizer<V,M>::GaussianVectorRealizer(const char* prefix,
                  const VectorSet<V,M>& unifiedImageSet,
                  const V& lawExpVector,
                  const M& matU,
                  const V& vecSsqrt,
                  const M& matVt)
  :
  BaseVectorRealizer<V,M>( ((std::string)(prefix)+"gau").c_str(), unifiedImageSet, std::numeric_limits<unsigned int>::max()), // 2011/Oct/02 - Correction thanks to Corey
  m_unifiedLawExpVector  (new V(lawExpVector)),
  m_unifiedLawVarVector  (unifiedImageSet.vectorSpace().newVector( INFINITY)), // FIX ME
  m_lowerCholLawCovMatrix(NULL),
  m_matU                 (new M(matU)),
  m_vecSsqrt             (new V(vecSsqrt)),
  m_matVt                (new M(matVt))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorRealizer<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  *m_unifiedLawExpVector = lawExpVector; // ????

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorRealizer<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
GaussianVectorRealizer<V,M>::~GaussianVectorRealizer()
{
  delete m_matVt;
  delete m_vecSsqrt;
  delete m_matU;
  delete m_lowerCholLawCovMatrix;
  delete m_unifiedLawVarVector;
  delete m_unifiedLawExpVector;
}
// Realization-related methods----------------------
template <class V, class M>
const V&
GaussianVectorRealizer<V,M>::unifiedLawExpVector() const
{
  return *m_unifiedLawExpVector;
}
//--------------------------------------------------
template <class V, class M>
const V&
GaussianVectorRealizer<V,M>::unifiedLawVarVector() const
{
  return *m_unifiedLawVarVector;
}
//--------------------------------------------------
template<class V, class M>
void
GaussianVectorRealizer<V,M>::realization(V& nextValues) const
{
  V iidGaussianVector(m_unifiedImageSet.vectorSpace().zeroVector());

  bool outOfSupport = true;
  do {
    iidGaussianVector.cwSetGaussian(0.0, 1.0);

    if (m_lowerCholLawCovMatrix) {
      nextValues = (*m_unifiedLawExpVector) + (*m_lowerCholLawCovMatrix)*iidGaussianVector;
    }
    else if (m_matU && m_vecSsqrt && m_matVt) {
      nextValues = (*m_unifiedLawExpVector) + (*m_matU)*( (*m_vecSsqrt) * ((*m_matVt)*iidGaussianVector) );
    }
    else {
      queso_error_msg("inconsistent internal state");
    }

    outOfSupport = !(this->m_unifiedImageSet.contains(nextValues));
  } while (outOfSupport); // prudenci 2011-Oct-04

  return;
}
//--------------------------------------------------
template<class V, class M>
void
GaussianVectorRealizer<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // delete old expected values (allocated at construction or last call to this function)
  delete m_unifiedLawExpVector;

  m_unifiedLawExpVector = new V(newLawExpVector);

  return;
}
//--------------------------------------------------
template<class V, class M>
void
GaussianVectorRealizer<V,M>::updateLowerCholLawCovMatrix(const M& newLowerCholLawCovMatrix)
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

  return;
}
//--------------------------------------------------
template<class V, class M>
void
GaussianVectorRealizer<V,M>::updateLowerCholLawCovMatrix(
  const M& matU,
  const V& vecSsqrt,
  const M& matVt)
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

  return;
}

}  // End namespace QUESO

template class QUESO::GaussianVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
