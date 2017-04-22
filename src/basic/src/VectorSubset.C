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

#include <queso/VectorSubset.h>
#include <queso/VectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Shaped constructor
template <class V, class M>
VectorSubset<V,M>::VectorSubset(const char* prefix,
    const VectorSpace<V,M>& vectorSpace,
    double volume)
  : VectorSet<V,M>(vectorSpace.env(), prefix, volume),
    m_vectorSpace(&vectorSpace)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering VectorSubset<V,M>::constructor()"
              << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving VectorSubset<V,M>::constructor()"
              << std::endl;
  }
}

// Destructor
template <class V, class M>
VectorSubset<V,M>::~VectorSubset()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering VectorSubset<V,M>::destructor()"
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving VectorSubset<V,M>::destructor()"
                            << std::endl;
  }
}

// Math methods
template <class V, class M>
const VectorSpace<V,M>& VectorSubset<V,M>::vectorSpace() const
{
  return *m_vectorSpace;
}

// I/O methods
template <class V, class M>
void VectorSubset<V,M>::print(std::ostream& os) const
{
  os << "In VectorSubset<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::VectorSubset<QUESO::GslVector, QUESO::GslMatrix>;
