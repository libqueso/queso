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

#include <queso/Environment.h>
#include <queso/Defines.h>
#include <queso/VectorSet.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TeuchosVector.h>
#include <queso/TeuchosMatrix.h>

namespace QUESO {

// Default constructor
template <class V, class M>
VectorSet<V,M>::VectorSet()
  : m_env(*(new EmptyEnvironment()))
{
}

// Shaped constructor
template <class V, class M>
VectorSet<V,M>::VectorSet(const BaseEnvironment& env, const char * prefix,
    double volume)
  : m_env(env),
    m_prefix(prefix),
    m_volume(volume)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering VectorSet<V,M>::constructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving VectorSet<V,M>::constructor()"
                           << std::endl;
  }
}

// Destructor
template <class V, class M>
VectorSet<V,M>::~VectorSet()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering VectorSet<V,M>::destructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving VectorSet<V,M>::destructor()"
                           << std::endl;
  }
}

// Environment methods
template <class V, class M>
const BaseEnvironment& VectorSet<V,M>::env() const
{
  return m_env;
}

template <class V, class M>
const std::string& VectorSet<V,M>::prefix() const
{
  return m_prefix;
}

// Mathematical methods
template <class V, class M>
double VectorSet<V,M>::volume() const
{
  return m_volume;
}

// I/O methods
template <class V, class M>
void VectorSet<V,M>::print(std::ostream& os) const
{
  os << "In VectorSet<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::VectorSet<QUESO::GslVector, QUESO::GslMatrix>;
#ifdef QUESO_HAS_TRILINOS
template class QUESO::VectorSet<QUESO::TeuchosVector, QUESO::TeuchosMatrix>;
#endif
