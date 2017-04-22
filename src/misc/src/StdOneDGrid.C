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

#include <queso/StdOneDGrid.h>

namespace QUESO {

template<class T>
StdOneDGrid<T>::StdOneDGrid(
  const BaseEnvironment& env,
  const char*                   prefix,
  const std::vector<T>&         points)
  :
  BaseOneDGrid<T>(env,prefix),
  m_points              (points)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering StdOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving StdOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
StdOneDGrid<T>::~StdOneDGrid()
{
}

template<class T>
unsigned int
StdOneDGrid<T>::size() const
{
  return m_points.size();
}

template<class T>
T
StdOneDGrid<T>::operator[](unsigned int i) const
{
  queso_require_less_msg(i, m_points.size(), "too large i");

  return m_points[i];
}

template<class T>
unsigned int
StdOneDGrid<T>::findIntervalId(const T& paramValue) const
{
  queso_require_msg(!((paramValue < m_points[0]) || (m_points[m_points.size()-1] < paramValue)), "paramValue is out of domain");

  unsigned int iMax = m_points.size();
  unsigned int i = 1; // Yes, '1'
  for (i = 1; i < iMax; ++i) { // Yes, '1'
    if (paramValue < m_points[i]) {
      i--;
      break;
    }
  }

  return i;
}

}  // End namespace QUESO
