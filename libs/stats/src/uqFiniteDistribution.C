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

#include <uqFiniteDistribution.h>

uqFiniteDistributionClass::uqFiniteDistributionClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
  const std::vector<double>&    inpWeights)
  :
  m_env    (env),
  m_prefix ((std::string)(prefix)+"fd_"),
  m_weights(inpWeights.size(),0.)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqFiniteDistributionClass::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  for (unsigned int i = 0; i < m_weights.size(); ++i) {
    m_weights[i] = inpWeights[i];
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqFiniteDistributionClass::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

uqFiniteDistributionClass::~uqFiniteDistributionClass()
{
  m_weights.clear();
}

const uqBaseEnvironmentClass&
uqFiniteDistributionClass::env() const
{
  return m_env;
}

const std::vector<double>&
uqFiniteDistributionClass::weights() const
{
  return m_weights;
}

unsigned int
uqFiniteDistributionClass::sample() const
{
  unsigned int result = 0;

  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqFiniteDistributionClass::sample()",
                      "incomplete code");

  return result;
}
