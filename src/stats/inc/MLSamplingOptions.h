//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H
#define UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H

#include <queso/Environment.h>
#include <queso/MLSamplingLevelOptions.h>

#define UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE "."

// _ODV = option default value

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY

#define UQ_ML_SAMPLING_RESTART_OUTPUT_LEVEL_PERIOD_ODV         0
#define UQ_ML_SAMPLING_RESTART_OUTPUT_BASE_NAME_FOR_FILES_ODV  UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_RESTART_OUTPUT_FILE_TYPE_ODV            UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_RESTART_INPUT_BASE_NAME_FOR_FILES_ODV   UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT

#else

#define UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV             UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV             UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV                  100

#endif

#define UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV               UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_DATA_OUTPUT_ALLOW_ALL_ODV               0
#define UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV             ""

namespace QUESO {

class MLSamplingOptions
{
public:
  MLSamplingOptions(const BaseEnvironment& env, const char* prefix);
 ~MLSamplingOptions();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  std::string            m_prefix;

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  unsigned int           m_restartOutput_levelPeriod;
  std::string            m_restartOutput_baseNameForFiles;
  std::string            m_restartOutput_fileType;
  std::string            m_restartInput_baseNameForFiles;
  std::string            m_restartInput_fileType;
#else
  std::string            m_restartInputFileName;
  std::string            m_restartInputFileType;
  unsigned int           m_restartChainSize;
#endif
  std::string            m_dataOutputFileName;
  bool                   m_dataOutputAllowAll;
  std::set<unsigned int> m_dataOutputAllowedSet;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const BaseEnvironment& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  std::string                   m_option_restartOutput_levelPeriod;
  std::string                   m_option_restartOutput_baseNameForFiles;
  std::string                   m_option_restartOutput_fileType;
  std::string                   m_option_restartInput_baseNameForFiles;
  std::string                   m_option_restartInput_fileType;
#else
  std::string                   m_option_restartInputFileName;
  std::string                   m_option_restartInputFileType;
  std::string                   m_option_restartChainSize;
#endif
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowAll;
  std::string                   m_option_dataOutputAllowedSet;
};

std::ostream& operator<<(std::ostream& os, const MLSamplingOptions& obj);

}  // End namespace QUESO

#endif // UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H
