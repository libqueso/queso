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

#include <queso/ModelValidation.h>

namespace QUESO {

// Constructor -------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
ModelValidation<P_V,P_M,Q_V,Q_M>::ModelValidation(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix) + ""),
  m_cycle (NULL)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering ModelValidation<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving ModelValidation<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  return;
}
// Destructor---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
ModelValidation<P_V,P_M,Q_V,Q_M>::~ModelValidation()
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering ModeValidation::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_cycle) delete m_cycle;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving ModeValidation::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Misc methods-------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const BaseEnvironment&
ModelValidation<P_V,P_M,Q_V,Q_M>::env() const
{
  return m_env;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const ValidationCycle<P_V,P_M,Q_V,Q_M>&
ModelValidation<P_V,P_M,Q_V,Q_M>::cycle() const
{
  return *m_cycle;
}

}  // End namespace QUESO
