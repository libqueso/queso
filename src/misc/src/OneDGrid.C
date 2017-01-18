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

#include <queso/OneDGrid.h>

namespace QUESO {

template<class T>
BaseOneDGrid<T>::BaseOneDGrid(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix)+"grid")
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering BaseOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving BaseOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
BaseOneDGrid<T>::~BaseOneDGrid()
{
}

template <class T>
void
BaseOneDGrid<T>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_prefix << "_sub" << m_env.subIdString() << " = zeros(" << this->size()
     << ","                                                      << 1
     << ");"
     << std::endl;
  os << m_prefix << "_sub" << m_env.subIdString() << " = [";
  for (unsigned int j = 0; j < this->size(); ++j) {
    os << (*this)[j] << " ";
  }
  os << "];"
     << std::endl;

  return;
}

}  // End namespace QUESO

template class QUESO::BaseOneDGrid<double>;
