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

#include <uqStatisticalForwardProblemOptions.h>
#include <uqMiscellaneous.h>

uqStatisticalForwardProblemOptionsClass::uqStatisticalForwardProblemOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_prefix                                   ((std::string)(prefix) + "mc_"),
  m_env                                      (env),
  m_optionsDesc                              (new po::options_description("Statistical Forward Problem options")),
  m_option_help                              (m_prefix + "help"                                                 )
{
}

uqStatisticalForwardProblemOptionsClass::~uqStatisticalForwardProblemOptionsClass()
{
  if (m_optionsDesc                    ) delete m_optionsDesc;
} 

void
uqStatisticalForwardProblemOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqStatisticalForwardProblemOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
};

void
uqStatisticalForwardProblemOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                                             "produce help message for Bayesian Markov chain distr. calculator")
  ;

  return;
}

void
uqStatisticalForwardProblemOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  return;
}

void
uqStatisticalForwardProblemOptionsClass::print(std::ostream& os) const
{
  return;
}

std::ostream& operator<<(std::ostream& os, const uqStatisticalForwardProblemOptionsClass& obj)
{
  obj.print(os);

  return os;
}
