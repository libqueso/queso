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

  unsigned int numOfZeroWeights = 0;
  double sumCheck = 0.;
  for (unsigned int i = 0; i < m_weights.size(); ++i) {
    m_weights[i] = inpWeights[i];
    sumCheck += m_weights[i];
    m_map.insert(std::map<double,unsigned int>::value_type(sumCheck,i));
    if (inpWeights[i] == 0.) numOfZeroWeights++;
  }

  std::cout << "\n  numOfZeroWeights = " << numOfZeroWeights
            << "\n  m_map.size() = "     << m_map.size()
            << "\n  m_weights.size() = " << m_weights.size()
            << std::endl;

  UQ_FATAL_TEST_MACRO((m_map.size() != m_weights.size()),
                      m_env.fullRank(),
                      "uqFiniteDistributionClass::constructor()",
                      "map and inpWeights have different sizes");

  UQ_FATAL_TEST_MACRO(((1 - sumCheck) > 1.e-12) || ((sumCheck - 1) > 1.e-12),
                      m_env.fullRank(),
                      "uqFiniteDistributionClass::constructor()",
                      "weights sum is too far from 1.");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqFiniteDistributionClass::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

uqFiniteDistributionClass::~uqFiniteDistributionClass()
{
  m_map.empty();
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

  double aux = gsl_rng_uniform(m_env.rng());
  if (aux == 0.) {
    result = 0;
  }
  else if (aux == 1.) {
    result = m_map.size()-1;
  }
  else {
    result = m_map.upper_bound(aux)->second-1;
  }

  UQ_FATAL_TEST_MACRO((result >= m_map.size()),
                      m_env.fullRank(),
                      "uqFiniteDistributionClass::sample()",
                      "invalid result");

  return result;
}
