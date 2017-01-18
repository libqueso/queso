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

#include <queso/UniformOneDGrid.h>

namespace QUESO {

// Constructor-------------------------------------------
template<class T>
UniformOneDGrid<T>::UniformOneDGrid(
  const BaseEnvironment& env,
  const char*               prefix,
        unsigned int        size,
        T                   minPosition,
        T                   maxPosition)
  :
  BaseOneDGrid<T>(env,prefix),
  m_size       (size),
  m_minPosition(minPosition),
  m_maxPosition(maxPosition)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering UniformOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving UniformOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor--------------------------------------------
template<class T>
UniformOneDGrid<T>::~UniformOneDGrid()
{
}
// Math methods------------------------------------------
template<class T>
unsigned int
UniformOneDGrid<T>::size() const
{
  return m_size;
}
//-------------------------------------------------------
template<class T>
T
UniformOneDGrid<T>::operator[](unsigned int i) const
{
  queso_require_less_msg(i, m_size, "too large i");

  T ratio = ((T) i)/(((T)m_size)-1.); // IMPORTANT: Yes, '-1.'
  T position = (1.-ratio)*m_minPosition + ratio*m_maxPosition;
  return position;
}
//-------------------------------------------------------
template<class T>
unsigned int
UniformOneDGrid<T>::findIntervalId(const T& paramValue) const
{
  queso_require_msg(!((paramValue < m_minPosition) || (m_maxPosition < paramValue)), "paramValue is out of domain");

  T ratio = (paramValue - m_minPosition)/(m_maxPosition - m_minPosition);
  unsigned int i = (unsigned int) (ratio*(m_size-1.));
  if ((i > 0                  ) &&
      ((*this)[i] > paramValue)) {
    i--;
  }

  return i;
}

}  // End namespace QUESO

template class QUESO::UniformOneDGrid<double>;
