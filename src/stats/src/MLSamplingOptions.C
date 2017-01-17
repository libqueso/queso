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

#include <queso/MLSamplingOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

MLSamplingOptions::MLSamplingOptions(const BaseEnvironment& env, const char* prefix)
  :
    m_prefix                               ((std::string)(prefix) + "ml_"                        ),
    m_help            (UQ_ML_SAMPLING_HELP),
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
    m_restartOutput_levelPeriod            (UQ_ML_SAMPLING_RESTART_OUTPUT_LEVEL_PERIOD_ODV       ),
    m_restartOutput_baseNameForFiles       (UQ_ML_SAMPLING_RESTART_OUTPUT_BASE_NAME_FOR_FILES_ODV),
    m_restartOutput_fileType               (UQ_ML_SAMPLING_RESTART_OUTPUT_FILE_TYPE_ODV          ),
    m_restartInput_baseNameForFiles        (UQ_ML_SAMPLING_RESTART_INPUT_BASE_NAME_FOR_FILES_ODV ),
    m_restartInput_fileType                (UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV           ),
#else
    m_restartInputFileName                 (UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV),
    m_restartInputFileType                 (UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV),
    m_restartChainSize                     (UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV     ),
#endif
    m_dataOutputFileName                   (UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  ),
  //m_dataOutputAllowedSet                 (),
    m_env                                  (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(env.optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help                          (m_prefix + "help"                          ),
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
    m_option_restartOutput_levelPeriod     (m_prefix + "restartOutput_levelPeriod"     ),
    m_option_restartOutput_baseNameForFiles(m_prefix + "restartOutput_baseNameForFiles"),
    m_option_restartOutput_fileType        (m_prefix + "restartOutput_fileType"        ),
    m_option_restartInput_baseNameForFiles (m_prefix + "restartInput_baseNameForFiles" ),
    m_option_restartInput_fileType         (m_prefix + "restartInput_fileType"         ),
#else
    m_option_restartInputFileName          (m_prefix + "restartInputFileName"),
    m_option_restartInputFileType          (m_prefix + "restartInputFileType"),
    m_option_restartChainSize              (m_prefix + "restartChainSize"    ),
#endif
    m_option_dataOutputFileName            (m_prefix + "dataOutputFileName"  ),
    m_option_dataOutputAllowedSet          (m_prefix + "dataOutputAllowedSet")
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string >(m_option_help,                           UQ_ML_SAMPLING_HELP,                                   "produce help msg for ML sampling options"      );
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  m_parser->registerOption<unsigned int>(m_option_restartOutput_levelPeriod,      UQ_ML_SAMPLING_RESTART_OUTPUT_LEVEL_PERIOD_ODV,        "restartOutput_levelPeriod"                     );
  m_parser->registerOption<std::string >(m_option_restartOutput_baseNameForFiles, UQ_ML_SAMPLING_RESTART_OUTPUT_BASE_NAME_FOR_FILES_ODV, "restartOutput_baseNameForFiles"                );
  m_parser->registerOption<std::string >(m_option_restartOutput_fileType,         UQ_ML_SAMPLING_RESTART_OUTPUT_FILE_TYPE_ODV,           "restartOutput_fileType"                        );
  m_parser->registerOption<std::string >(m_option_restartInput_baseNameForFiles,  UQ_ML_SAMPLING_RESTART_INPUT_BASE_NAME_FOR_FILES_ODV,  "restartInput_baseNameForFiles"                 );
  m_parser->registerOption<std::string >(m_option_restartInput_fileType,          UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV,            "restartInput_fileType"                         );
#else
  m_parser->registerOption<std::string >(m_option_restartInputFileName,           UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV,            "name of restart input file"                    );
  m_parser->registerOption<std::string >(m_option_restartInputFileType,           UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV,            "type of restart input file"                    );
  m_parser->registerOption<unsigned int>(m_option_restartChainSize,               UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV,                 "size of restart chain"                         );
#endif
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,             UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  ,            "name of generic output file"                   );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,           UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV,            "subEnvs that will write to generic output file");

  m_parser->scanInputFile();

  m_parser->getOption<std::string >(m_option_help,                           m_help);
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  m_parser->getOption<unsigned int>(m_option_restartOutput_levelPeriod,      m_restartOutput_levelPeriod);
  m_parser->getOption<std::string >(m_option_restartOutput_baseNameForFiles, m_restartOutput_baseNameForFiles);
  m_parser->getOption<std::string >(m_option_restartOutput_fileType,         m_restartOutput_fileType);
  m_parser->getOption<std::string >(m_option_restartInput_baseNameForFiles,  m_restartInput_baseNameForFiles);
  m_parser->getOption<std::string >(m_option_restartInput_fileType,          m_restartInput_fileType);
#else
  m_parser->getOption<std::string >(m_option_restartInputFileName,           m_restartInputFileName);
  m_parser->getOption<std::string >(m_option_restartInputFileType,           m_restartInputFileType);
  m_parser->getOption<unsigned int>(m_option_restartChainSize,               m_restartChainSize);
#endif
  m_parser->getOption<std::string >(m_option_dataOutputFileName,             m_dataOutputFileName);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,           m_dataOutputAllowedSet);
#else
  m_help = m_env.input()(m_option_help, UQ_ML_SAMPLING_HELP);
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  m_restartOutput_levelPeriod = m_env.input()(m_option_restartOutput_levelPeriod, UQ_ML_SAMPLING_RESTART_OUTPUT_LEVEL_PERIOD_ODV);
  m_restartOutput_baseNameForFiles = m_env.input()(m_option_restartOutput_baseNameForFiles, UQ_ML_SAMPLING_RESTART_OUTPUT_BASE_NAME_FOR_FILES_ODV);
  m_restartOutput_fileType = m_env.input()(m_option_restartOutput_fileType, UQ_ML_SAMPLING_RESTART_OUTPUT_FILE_TYPE_ODV);
  m_restartInput_baseNameForFiles = m_env.input()(m_option_restartInput_baseNameForFiles, UQ_ML_SAMPLING_RESTART_INPUT_BASE_NAME_FOR_FILES_ODV);
  m_restartInput_fileType = m_env.input()(m_option_restartInput_fileType, UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV);
#else
  m_restartInputFileName = m_env.input()(m_option_restartInputFileName, UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV);
  m_restartInputFileType = m_env.input()(m_option_restartInputFileType, UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV);
  m_restartChainSize = m_env.input()(m_option_restartChainSize, UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV);
#endif
  m_dataOutputFileName = m_env.input()(m_option_dataOutputFileName, UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  );

  // UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by
  // default
  unsigned int size = m_env.input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env.input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions(&env);
}

MLSamplingOptions::~MLSamplingOptions()
{
}

void
MLSamplingOptions::checkOptions(const BaseEnvironment * env)
{
  // DM: I'm printing here because I'm not sure where this object is created
  //     in the statistical inverse problem
  if (m_help != "") {
    if (env->subDisplayFile()) {
      *(env->subDisplayFile()) << (*this) << std::endl;
    }
  }

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  if ((m_restartOutput_levelPeriod > 0)) queso_require_not_equal_to_msg(m_restartOutput_baseNameForFiles, std::string("."), std::string("Option 'restartOutput_levelPeriod' is > 0, but 'restartOutput_baseNameForFiles' is not specified..."));
#endif
}

void
MLSamplingOptions::print(std::ostream& os) const
{
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  os <<         m_option_restartOutput_levelPeriod      << " = " << m_restartOutput_levelPeriod
     << "\n" << m_option_restartOutput_baseNameForFiles << " = " << m_restartOutput_baseNameForFiles
     << "\n" << m_option_restartOutput_fileType         << " = " << m_restartOutput_fileType
     << "\n" << m_option_restartInput_baseNameForFiles  << " = " << m_restartInput_baseNameForFiles
     << "\n" << m_option_restartInput_fileType          << " = " << m_restartInput_fileType
#else
  os <<         m_option_restartInputFileName           << " = " << m_restartInputFileName
     << "\n" << m_option_restartInputFileType           << " = " << m_restartInputFileType
     << "\n" << m_option_restartChainSize               << " = " << m_restartChainSize
#endif
     << "\n" << m_option_dataOutputFileName             << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet           << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n";

  return;
}

std::ostream& operator<<(std::ostream& os, const MLSamplingOptions& obj)
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  obj.print(os);
  return os;
}

}  // End namespace QUESO
