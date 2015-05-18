//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <boost/program_options.hpp>

#include <queso/MLSamplingOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

MLSamplingOptions::MLSamplingOptions(const BaseEnvironment& env, const char* prefix)
  :
    m_prefix                               ((std::string)(prefix) + "ml_"                        ),
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
    m_parser(new BoostInputOptionsParser(env.optionsInputFileName())),
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
  m_parser->registerOption(m_option_help,                                                                                                                            "produce help msg for ML sampling options"      );
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
}

MLSamplingOptions::~MLSamplingOptions()
{
}

// void
// MLSamplingOptions::defineOptions()
// {
//   (*m_optionsDescription).add_options()
//     (m_option_help.c_str(),                                                                                                                            "produce help msg for ML sampling options"      )
// #ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
//     (m_option_restartOutput_levelPeriod.c_str(),      boost::program_options::value<unsigned int>()->default_value(UQ_ML_SAMPLING_RESTART_OUTPUT_LEVEL_PERIOD_ODV),        "restartOutput_levelPeriod"                     )
//     (m_option_restartOutput_baseNameForFiles.c_str(), boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_OUTPUT_BASE_NAME_FOR_FILES_ODV), "restartOutput_baseNameForFiles"                )
//     (m_option_restartOutput_fileType.c_str(),         boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_OUTPUT_FILE_TYPE_ODV),           "restartOutput_fileType"                        )
//     (m_option_restartInput_baseNameForFiles.c_str(),  boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_BASE_NAME_FOR_FILES_ODV),  "restartInput_baseNameForFiles"                 )
//     (m_option_restartInput_fileType.c_str(),          boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV),            "restartInput_fileType"                         )
// #else
//     (m_option_restartInputFileName.c_str(),           boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV),            "name of restart input file"                    )
//     (m_option_restartInputFileType.c_str(),           boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV),            "type of restart input file"                    )
//     (m_option_restartChainSize.c_str(),               boost::program_options::value<unsigned int>()->default_value(UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV),                 "size of restart chain"                         )
// #endif
//     (m_option_dataOutputFileName.c_str(),             boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  ),            "name of generic output file"                   )
//     (m_option_dataOutputAllowedSet.c_str(),           boost::program_options::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV),            "subEnvs that will write to generic output file")
//   ;
// }
//
// void
// MLSamplingOptions::getOptionValues()
// {
//   if ((*m_optionsMap).count(m_option_help.c_str())) {
//     if (m_env.subDisplayFile()) {
//       *m_env.subDisplayFile() << (*m_optionsDescription)
//                               << std::endl;
//     }
//   }
//
// #ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
//   if ((*m_optionsMap).count(m_option_restartOutput_levelPeriod.c_str())) {
//     m_restartOutput_levelPeriod = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartOutput_levelPeriod.c_str()]).as<unsigned int>();
//   }
//
//   if ((*m_optionsMap).count(m_option_restartOutput_baseNameForFiles.c_str())) {
//     m_restartOutput_baseNameForFiles = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartOutput_baseNameForFiles.c_str()]).as<std::string>();
//   }
//
//   if ((m_restartOutput_levelPeriod > 0)) queso_require_not_equal_to_msg(m_restartOutput_baseNameForFiles, ".", "Option 'restartOutput_levelPeriod' is > 0, but 'restartOutput_baseNameForFiles' is not specified...");
//
//   if ((*m_optionsMap).count(m_option_restartOutput_fileType.c_str())) {
//     m_restartOutput_fileType = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartOutput_fileType.c_str()]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_restartInput_baseNameForFiles.c_str())) {
//     m_restartInput_baseNameForFiles = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartInput_baseNameForFiles.c_str()]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_restartInput_fileType.c_str())) {
//     m_restartInput_fileType = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartInput_fileType.c_str()]).as<std::string>();
//   }
// #else
//   if ((*m_optionsMap).count(m_option_restartInputFileName.c_str())) {
//     m_restartInputFileName = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartInputFileName.c_str()]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_restartInputFileType.c_str())) {
//     m_restartInputFileType = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartInputFileType.c_str()]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_restartChainSize.c_str())) {
//     m_restartChainSize = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_restartChainSize.c_str()]).as<unsigned int>();
//   }
// #endif
//   if ((*m_optionsMap).count(m_option_dataOutputFileName.c_str())) {
//     m_dataOutputFileName = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_dataOutputFileName.c_str()]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_dataOutputAllowedSet.c_str())) {
//     m_dataOutputAllowedSet.clear();
//     std::vector<double> tmpAllow(0,0.);
//     std::string inputString = (*m_optionsMap)[m_option_dataOutputAllowedSet.c_str()].as<std::string>();
//     MiscReadDoublesFromString(inputString,tmpAllow);
//
//     if (tmpAllow.size() > 0) {
//       for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
//         m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
//       }
//     }
//   }
//
//   return;
// }

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
  obj.print(os);

  return os;
}

}  // End namespace QUESO
