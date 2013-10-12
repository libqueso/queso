//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <queso/VectorSpace.h>
#include <queso/DiscreteSubset.h>

namespace QUESO {

// Default constructor
template<class V, class M>
DiscreteSubset<V,M>::DiscreteSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    const std::vector<V*>& elements)
  : VectorSubset<V,M>(prefix, vectorSpace, 0.),
    m_elements(elements.size(),NULL)
{
  m_volume = 0.;
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "DiscreteSubset<V,M>::contains()",
                      "incomplete code");
}

// Destructor
template<class V, class M>
DiscreteSubset<V,M>::~DiscreteSubset()
{
}

// Mathematical methods
template<class V, class M>
bool DiscreteSubset<V,M>::contains(const V& vec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "DiscreteSubset<V,M>::contains()",
                      "incomplete code");

  return false;
}

// I/O methods
template <class V, class M>
void DiscreteSubset<V,M>::print(std::ostream& os) const
{
  os << "In DiscreteSubset<V,M>::print()"
     << ": nothing to print"
     << std::endl;

  return;
}

}  // End namespace QUESO
