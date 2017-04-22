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

#include <queso/GenericVectorMdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V, class M>
GenericVectorMdf<V,M>::GenericVectorMdf(
  const char*                    prefix,
  const VectorSet<V,M>& domainSet,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
  const void* routineDataPtr)
  :
  BaseVectorMdf<V,M>(prefix,domainSet),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}
// Destructor ---------------------------------------
template<class V, class M>
GenericVectorMdf<V,M>::~GenericVectorMdf()
{
}
// Math methods--------------------------------------
template<class V, class M>
void
GenericVectorMdf<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  m_routinePtr(paramValues, m_routineDataPtr, mdfVec);
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
GenericVectorMdf<V,M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO
