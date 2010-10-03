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

#ifndef __UQ_ENVIRONMENT_OPTIONS_H__
#define __UQ_ENVIRONMENT_OPTIONS_H__

#include <uqEnvironment.h>

#define UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE "."
#define UQ_ENV_FILENAME_FOR_NO_INPUT_FILE  "."

#define UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV    1
#define UQ_ENV_SUB_SCREEN_WRITE_ODV        0
#define UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV   UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV   0
#define UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV ""
#define UQ_ENV_DISPLAY_VERBOSITY_ODV       0
#define UQ_ENV_SYNC_VERBOSITY_ODV          0
#define UQ_ENV_SEED_ODV                    0
#define UQ_ENV_IDENTIFYING_STRING_ODV      ""
#define UQ_ENV_NUM_DEBUG_PARAMS_ODV        0
#define UQ_ENV_DEBUG_PARAM_ODV             0.

class uqEnvironmentOptionsClass
{
public:
  uqEnvironmentOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
 ~uqEnvironmentOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  unsigned int             m_numSubEnvironments;
  std::string              m_subDisplayFileName;
  bool                     m_subDisplayAllowAll;
  std::set<unsigned int>   m_subDisplayAllowedSet;
  unsigned int             m_displayVerbosity;
  unsigned int             m_syncVerbosity;
  int                      m_seed;
  std::string              m_identifyingString;
  unsigned int             m_numDebugParams;
  std::vector<double>      m_debugParams;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  std::string              m_prefix;
  po::options_description* m_optionsDesc;

  std::string              m_option_help;

  std::string              m_option_numSubEnvironments;
  std::string              m_option_subDisplayFileName;
  std::string              m_option_subDisplayAllowAll;
  std::string              m_option_subDisplayAllowedSet;
  std::string              m_option_displayVerbosity;
  std::string              m_option_syncVerbosity;
  std::string              m_option_seed;
  std::string              m_option_identifyingString;
};

std::ostream& operator<<(std::ostream& os, const uqEnvironmentOptionsClass& obj);

#endif // __UQ_ENVIRONMENT_CLASS_H__
