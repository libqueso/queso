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

#include <queso/VectorCdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V, class M>
BaseVectorCdf<V,M>::BaseVectorCdf(
  const char*                  prefix,
  const VectorSet<V,M>& pdfSupport)
  :
  m_env       (pdfSupport.env()),
  m_prefix    ((std::string)(prefix)+"Cdf_"),
  m_pdfSupport(pdfSupport)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering BaseVectorCdf<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving BaseVectorCdf<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BaseVectorCdf<V,M>::~BaseVectorCdf()
{
}
// Math methods--------------------------------------
template<class V, class M>
const VectorSet<V,M>&
BaseVectorCdf<V,M>::pdfSupport() const
{
  return m_pdfSupport;
}
// I/O methods---------------------------------------
template<class V, class M>
void
BaseVectorCdf<V,M>::subWriteContents(
  const std::string&            /* varNamePrefix */,
  const std::string&            /* fileName */,
  const std::string&            /* fileType */,
  const std::set<unsigned int>& /* allowedSubEnvIds */) const
{
  std::cerr << "WARNING: BaseVectorCdf<V,M>::subWriteContents() being used..."
            << std::endl;
  return;
}

//---------------------------------------------------
// Method outside either class definition------------
//---------------------------------------------------
//! It calculated the maximum horizontal distances between two vector CDFs.
//
template <class V, class M>
void
horizontalDistances(const BaseVectorCdf<V,M>& cdf1,
                    const BaseVectorCdf<V,M>& cdf2,
                    const V& epsilonVec,
                    V&       distances)
{
  for (unsigned int i = 0; i < cdf1.pdfSupport().vectorSpace().dimLocal(); ++i) {
    distances[i] = horizontalDistance(cdf1.cdf(i),
                                      cdf2.cdf(i),
                                      epsilonVec[i]);
  }

  return;
}

}  // End namespace QUESO
