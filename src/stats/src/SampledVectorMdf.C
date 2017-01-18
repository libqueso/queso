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

#include <queso/SampledVectorMdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V,class M>
SampledVectorMdf<V,M>::SampledVectorMdf(
  const char*                          prefix,
  const ArrayOfOneDGrids <V,M>& oneDGrids,
  const ArrayOfOneDTables<V,M>& mdfValues)
  :
  BaseVectorMdf<V,M>(prefix,oneDGrids.rowSpace()),
  m_oneDGrids(oneDGrids),
  m_mdfValues(mdfValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SampledVectorMdf<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SampledVectorMdf<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
SampledVectorMdf<V,M>::~SampledVectorMdf()
{
}
// Math methods--------------------------------------
template<class V, class M>
void
SampledVectorMdf<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  queso_error_msg("incomplete code");
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
SampledVectorMdf<V,M>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_oneDGrids;

  // Print *mdf* values *at* grid points
  os << m_mdfValues;

  return;
}

}  // End namespace QUESO
