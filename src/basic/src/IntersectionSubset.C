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

#include <queso/VectorSpace.h>
#include <queso/VectorSet.h>
#include <queso/IntersectionSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default, shaped constructor
template<class V, class M>
IntersectionSubset<V,M>::IntersectionSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    double volume,
    const VectorSet<V,M>& set1,
    const VectorSet<V,M>& set2)
  : VectorSubset<V,M>(prefix,vectorSpace,volume),
    m_set1(set1),
    m_set2(set2)
{
}

// Destructor
template<class V, class M>
IntersectionSubset<V,M>::~IntersectionSubset()
{
}

// Mathematical methods
template<class V, class M>
bool IntersectionSubset<V,M>::contains(const V& vec) const
{
  return (m_set1.contains(vec) && m_set2.contains(vec));
}

template<class V, class M>
void IntersectionSubset<V,M>::centroid(V& vec) const
{
  // No general way to do this?
  queso_not_implemented();
}

template<class V, class M>
void IntersectionSubset<V,M>::moments(M& vec) const
{
  // No general way to do this?
  queso_not_implemented();
}

// I/O methods
template <class V, class M>
void IntersectionSubset<V,M>::print(std::ostream& os) const
{
  os << "In IntersectionSubset<V,M>::print()"
     << ": m_set1 = " << m_set1
     << ", m_set2 = " << m_set2
     << std::endl;

  return;
}

}  // End namespace QUESO


template class QUESO::IntersectionSubset<QUESO::GslVector, QUESO::GslMatrix>;
