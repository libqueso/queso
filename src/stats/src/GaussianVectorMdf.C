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

#include <queso/GaussianVectorMdf.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
GaussianVectorMdf<V,M>::GaussianVectorMdf(
  const char*                    prefix,
  const VectorSet<V,M>& domainSet,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  BaseVectorMdf<V,M>(prefix,domainSet),
  m_covMatrix              (m_domainSet.newDiagMatrix(domainVarianceValues*domainVarianceValues))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorMdf<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorMdf<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Constructor -------------------------------------
template<class V,class M>
GaussianVectorMdf<V,M>::GaussianVectorMdf(
  const char*                    prefix,
  const VectorSet<V,M>& domainSet,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  BaseVectorMdf<V,M>(prefix,domainSet),
  m_covMatrix              (new M(covMatrix))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorMdf<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorMdf<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor -------------------------------------
template<class V,class M>
GaussianVectorMdf<V,M>::~GaussianVectorMdf()
{
  delete m_covMatrix;
}
// Math method -------------------------------------
template<class V, class M>
void
GaussianVectorMdf<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  queso_error_msg("incomplete code");
  return;
}
// I/O method --------------------------------------
template <class V, class M>
void
GaussianVectorMdf<V,M>::print(std::ostream& os) const
{
  return;
}
// Protected member function------------------------
template<class V,class M>
void
GaussianVectorMdf<V,M>::commonConstructor()
{
  queso_error_msg("incomplete code");
  return;
}

}  // End namespace QUESO
