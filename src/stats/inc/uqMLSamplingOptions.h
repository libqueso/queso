//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H__
#define __UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H__

#include <uqEnvironment.h>
#include <uqMLSamplingLevelOptions.h>

#define UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV "m"
#define UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV      100
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
  std::string            m_restartInputFileType;
  unsigned int           m_restartChainSize;
  std::string            m_dataOutputFileName;
  std::set<unsigned int> m_dataOutputAllowedSet;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;
  std::string                   m_option_restartInputFileName;
  std::string                   m_option_restartInputFileType;
  std::string                   m_option_restartChainSize;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
};

std::ostream& operator<<(std::ostream& os, const uqMLSamplingOptionsClass& obj);
#endif // __UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H__
