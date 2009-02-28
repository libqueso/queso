/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MODEL_VALIDATION_H__
#define __UQ_MODEL_VALIDATION_H__

#include <uqValidationCycle.h>

template <class P_V,class P_M,class Q_V,class Q_M>
class uqModelValidationClass
{
public:
  uqModelValidationClass(const uqBaseEnvironmentClass& env,
                         const char*                   prefix);
 ~uqModelValidationClass();

  virtual void run() = 0;

  const uqBaseEnvironmentClass&                  env  () const;
  const uqValidationCycleClass<P_V,P_M,Q_V,Q_M>& cycle() const;

protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;

  uqValidationCycleClass<P_V,P_M,Q_V,Q_M>* m_cycle;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::uqModelValidationClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix) + ""),
  m_cycle (NULL)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Entering uqModelValidationClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Leaving uqModelValidationClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::~uqModelValidationClass()
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Entering uqModeValidation::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_cycle) delete m_cycle;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Leaving uqModeValidation::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqBaseEnvironmentClass&
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::env() const
{
  return m_env;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqValidationCycleClass<P_V,P_M,Q_V,Q_M>&
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::cycle() const
{
  return *m_cycle;
}

#endif // __UQ_MODEL_VALIDATION_H__
