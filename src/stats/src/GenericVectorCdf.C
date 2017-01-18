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

#include <queso/GenericVectorCdf.h>

namespace QUESO {

// Constructor ----------------------------------------
template<class V, class M>
GenericVectorCdf<V,M>::GenericVectorCdf(
  const char*                    prefix,
  const VectorSet<V,M>& pdfSupport,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
  const void* routineDataPtr)
  :
  BaseVectorCdf<V,M>(prefix,pdfSupport),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}
// Destructor ---------------------------------------
template<class V, class M>
GenericVectorCdf<V,M>::~GenericVectorCdf()
{
}
// Math method --------------------------------------
template<class V, class M>
void
GenericVectorCdf<V,M>::values(
  const V& paramValues,
        V& cdfVec) const
{
  m_routinePtr(paramValues, m_routineDataPtr, cdfVec);
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
GenericVectorCdf<V,M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO
