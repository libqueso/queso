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

#include <uqMLSamplingOptions.h>
#include <uqMiscellaneous.h>

namespace QUESO {

uqMLSamplingOptionsClass::uqMLSamplingOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
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
  m_optionsDesc                          (new po::options_description("Multilevel sampling options")),
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
}

uqMLSamplingOptionsClass::~uqMLSamplingOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
} 

void
uqMLSamplingOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqMLSamplingOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
uqMLSamplingOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                                            "produce help msg for ML sampling options"      )
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
    (m_option_restartOutput_levelPeriod.c_str(),      po::value<unsigned int>()->default_value(UQ_ML_SAMPLING_RESTART_OUTPUT_LEVEL_PERIOD_ODV),        "restartOutput_levelPeriod"                     )
    (m_option_restartOutput_baseNameForFiles.c_str(), po::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_OUTPUT_BASE_NAME_FOR_FILES_ODV), "restartOutput_baseNameForFiles"                )
    (m_option_restartOutput_fileType.c_str(),         po::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_OUTPUT_FILE_TYPE_ODV),           "restartOutput_fileType"                        )
    (m_option_restartInput_baseNameForFiles.c_str(),  po::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_BASE_NAME_FOR_FILES_ODV),  "restartInput_baseNameForFiles"                 )
    (m_option_restartInput_fileType.c_str(),          po::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV),            "restartInput_fileType"                         )
#else
    (m_option_restartInputFileName.c_str(),           po::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV),            "name of restart input file"                    )
    (m_option_restartInputFileType.c_str(),           po::value<std::string >()->default_value(UQ_ML_SAMPLING_RESTART_INPUT_FILE_TYPE_ODV),            "type of restart input file"                    )
    (m_option_restartChainSize.c_str(),               po::value<unsigned int>()->default_value(UQ_ML_SAMPLING_RESTART_CHAIN_SIZE_ODV),                 "size of restart chain"                         )
#endif
    (m_option_dataOutputFileName.c_str(),             po::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  ),            "name of generic output file"                   )
    (m_option_dataOutputAllowedSet.c_str(),           po::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV),            "subEnvs that will write to generic output file")
  ;

  return;
}

void
uqMLSamplingOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  if (m_env.allOptionsMap().count(m_option_restartOutput_levelPeriod.c_str())) {
    m_restartOutput_levelPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartOutput_levelPeriod.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_restartOutput_baseNameForFiles.c_str())) {
    m_restartOutput_baseNameForFiles = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartOutput_baseNameForFiles.c_str()]).as<std::string>();
  }

  UQ_FATAL_TEST_MACRO((m_restartOutput_levelPeriod > 0) && (m_restartOutput_baseNameForFiles == "."),
                      m_env.worldRank(),
                      "uqMLSamplingOptionsClass::getMyOptionsValues()",
                      "Option 'restartOutput_levelPeriod' is > 0, but 'restartOutput_baseNameForFiles' is not specified...");

  if (m_env.allOptionsMap().count(m_option_restartOutput_fileType.c_str())) {
    m_restartOutput_fileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartOutput_fileType.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_restartInput_baseNameForFiles.c_str())) {
    m_restartInput_baseNameForFiles = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartInput_baseNameForFiles.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_restartInput_fileType.c_str())) {
    m_restartInput_fileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartInput_fileType.c_str()]).as<std::string>();
  }
#else
  if (m_env.allOptionsMap().count(m_option_restartInputFileName.c_str())) {
    m_restartInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartInputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_restartInputFileType.c_str())) {
    m_restartInputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartInputFileType.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_restartChainSize.c_str())) {
    m_restartChainSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_restartChainSize.c_str()]).as<unsigned int>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_dataOutputFileName.c_str())) {
    m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet.c_str())) {
    m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  return;
}

void
uqMLSamplingOptionsClass::print(std::ostream& os) const
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

std::ostream& operator<<(std::ostream& os, const uqMLSamplingOptionsClass& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO
