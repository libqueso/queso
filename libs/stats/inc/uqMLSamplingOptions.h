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

#ifndef __UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H__
#define __UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H__

#include <uqEnvironment.h>
#include <uqMLSamplingLevelOptions.h>

#define UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV   UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV ""

class uqMLSamplingOptionsClass
{
public:
  uqMLSamplingOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
 ~uqMLSamplingOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  std::string            m_prefix;

  std::string            m_restartInputFileName;
  std::string            m_dataOutputFileName;
  std::set<unsigned int> m_dataOutputAllowedSet;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;
  std::string                   m_option_restartInputFileName;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
};

std::ostream& operator<<(std::ostream& os, const uqMLSamplingOptionsClass& obj);
#endif // __UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H__
