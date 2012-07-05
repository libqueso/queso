//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#include <uqEnvironmentOptions.h>
#include <uqMiscellaneous.h>

uqEnvOptionsValuesClass::uqEnvOptionsValuesClass()
  :
  m_numSubEnvironments  (UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
  m_subDisplayFileName  (UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),
  m_subDisplayAllowAll  (UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),
//m_subDisplayAllowedSet(),
  m_displayVerbosity    (UQ_ENV_DISPLAY_VERBOSITY_ODV),
  m_syncVerbosity       (UQ_ENV_SYNC_VERBOSITY_ODV),
  m_checkingLevel       (UQ_ENV_CHECKING_LEVEL_ODV),
  m_seed                (UQ_ENV_SEED_ODV),
  m_platformName        (UQ_ENV_PLATFORM_NAME_ODV),
  m_identifyingString   (UQ_ENV_IDENTIFYING_STRING_ODV),
  m_numDebugParams      (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams         (m_numDebugParams,0.)
{
}

uqEnvOptionsValuesClass::~uqEnvOptionsValuesClass()
{
}

uqEnvOptionsValuesClass::uqEnvOptionsValuesClass(const uqEnvOptionsValuesClass& src)
{
  this->copy(src);
}

uqEnvOptionsValuesClass&
uqEnvOptionsValuesClass::operator=(const uqEnvOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}

void
uqEnvOptionsValuesClass::copy(const uqEnvOptionsValuesClass& src)
{
  m_numSubEnvironments   = src.m_numSubEnvironments;
  m_subDisplayFileName   = src.m_subDisplayFileName;
  m_subDisplayAllowAll   = src.m_subDisplayAllowAll;
  m_subDisplayAllowedSet = src.m_subDisplayAllowedSet;
  m_displayVerbosity     = src.m_displayVerbosity;
  m_syncVerbosity        = src.m_syncVerbosity;
  m_checkingLevel        = src.m_checkingLevel;
  m_seed                 = src.m_seed;
  m_platformName         = src.m_platformName;
  m_identifyingString    = src.m_identifyingString;
  m_numDebugParams       = src.m_numDebugParams;
  m_debugParams          = src.m_debugParams;

  return;
}

uqEnvironmentOptionsClass::uqEnvironmentOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_ov                         (),
  m_env                        (env),
  m_prefix                     ((std::string)(prefix) + "env_"),
  m_optionsDesc                (new po::options_description("Environment options")),
  m_option_help                (m_prefix + "help"                ),
  m_option_numSubEnvironments  (m_prefix + "numSubEnvironments"  ),
  m_option_subDisplayFileName  (m_prefix + "subDisplayFileName"  ),
  m_option_subDisplayAllowAll  (m_prefix + "subDisplayAllowAll"  ),
  m_option_subDisplayAllowedSet(m_prefix + "subDisplayAllowedSet"),
  m_option_displayVerbosity    (m_prefix + "displayVerbosity"    ),
  m_option_syncVerbosity       (m_prefix + "syncVerbosity"       ),
  m_option_checkingLevel       (m_prefix + "checkingLevel"       ),
  m_option_seed                (m_prefix + "seed"                ),
  m_option_platformName        (m_prefix + "platformName"        ),
  m_option_identifyingString   (m_prefix + "identifyingString"   )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "uqEnvironmentOptionsClass::constructor(1)",
                      "this constructor is incompatible with the abscense of an options input file");
}

uqEnvironmentOptionsClass::uqEnvironmentOptionsClass(
  const uqBaseEnvironmentClass&  env,
  const char*                    prefix,
  const uqEnvOptionsValuesClass& alternativeOptionsValues)
  :
  m_ov                         (alternativeOptionsValues),
  m_env                        (env),
  m_prefix                     ((std::string)(prefix) + "env_"),
  m_optionsDesc                (NULL),
  m_option_help                (m_prefix + "help"                ),
  m_option_numSubEnvironments  (m_prefix + "numSubEnvironments"  ),
  m_option_subDisplayFileName  (m_prefix + "subDisplayFileName"  ),
  m_option_subDisplayAllowAll  (m_prefix + "subDisplayAllowAll"  ),
  m_option_subDisplayAllowedSet(m_prefix + "subDisplayAllowedSet"),
  m_option_displayVerbosity    (m_prefix + "displayVerbosity"    ),
  m_option_syncVerbosity       (m_prefix + "syncVerbosity"       ),
  m_option_checkingLevel       (m_prefix + "checkingLevel"       ),
  m_option_seed                (m_prefix + "seed"                ),
  m_option_platformName        (m_prefix + "platformName"        ),
  m_option_identifyingString   (m_prefix + "identifyingString"   )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "uqEnvironmentOptionsClass::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqEnvironmentOptionsClass::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

uqEnvironmentOptionsClass::~uqEnvironmentOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
}

void
uqEnvironmentOptionsClass::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "uqEnvironmentOptionsClass::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  // 'm_subDisplayOutputFile' is still not available at this moment. Use 'std::cout'
  //if (m_env.subScreenFile() != NULL) {
  //  *m_env.subScreenFile()
  if ((m_env.fullRank() == 0) && (m_env.displayVerbosity() >= 3)) {
    std::cout << "In uqEnvironmentOptionsClass::scanOptionsValues()"
              << ": after reading values of options with prefix '" << m_prefix
              << "', state of object is:"
              << "\n" << *this
              << std::endl;
  }

  return;
}

void
uqEnvironmentOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In uqEnvOptionsClass::defineMyOptions(), before add_options()" << std::endl;
#endif
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                               "produce help message for uq environment"  )
    (m_option_numSubEnvironments.c_str(),   po::value<unsigned int>()->default_value(UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),    "number of subEnvironments"                )
    (m_option_subDisplayFileName.c_str(),   po::value<std::string >()->default_value(UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),   "output filename for subscreen writing"    )
    (m_option_subDisplayAllowAll.c_str(),   po::value<bool        >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),   "Allow all subEnvs to write to output file")
    (m_option_subDisplayAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV), "subEnvs that will write to output file"   )
    (m_option_displayVerbosity.c_str(),     po::value<unsigned int>()->default_value(UQ_ENV_DISPLAY_VERBOSITY_ODV),       "set verbosity"                            )
    (m_option_syncVerbosity.c_str(),        po::value<unsigned int>()->default_value(UQ_ENV_SYNC_VERBOSITY_ODV),          "set sync verbosity"                       )
    (m_option_checkingLevel.c_str(),        po::value<unsigned int>()->default_value(UQ_ENV_CHECKING_LEVEL_ODV),          "set checking level"                       )
    (m_option_seed.c_str(),                 po::value<int         >()->default_value(UQ_ENV_SEED_ODV),                    "set seed"                                 )
    (m_option_platformName.c_str(),         po::value<std::string >()->default_value(UQ_ENV_PLATFORM_NAME_ODV),           "platform name"                            )
    (m_option_identifyingString.c_str(),    po::value<std::string >()->default_value(UQ_ENV_IDENTIFYING_STRING_ODV),      "identifying string"                       )
  //(m_option_numDebugParams.c_str(),       po::value<unsigned int>()->default_value(UQ_ENV_NUM_DEBUG_PARAMS_ODV),        "set number of debug parameters"           )
  ;
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In uqEnvOptionsClass::defineMyOptions(), after add_options()" << std::endl;
#endif

  return;
}

void
uqEnvironmentOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering uqEnvOptionsClass::getMyOptionsValues()" << std::endl;
#endif
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    // 'm_subDisplayOutputFile' is still not available at this moment. Use 'std::cout'
    if (m_env.fullRank() == 0) std::cout << optionsDesc
                                         << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_numSubEnvironments.c_str())) {
    m_ov.m_numSubEnvironments = m_env.allOptionsMap()[m_option_numSubEnvironments].as<unsigned int>();
  }
  if ((m_env.fullComm().NumProc()%m_ov.m_numSubEnvironments) != 0) {
    std::cerr << "In uqBaseEnvironmentClass::getMyOptionValues()"
              << ": m_env.fullComm().NumProc() = " << m_env.fullComm().NumProc()
              << ", m_numSubEnvironments = "       << m_ov.m_numSubEnvironments
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_env.fullComm().NumProc()%m_ov.m_numSubEnvironments) != 0,
                      m_env.worldRank(),
                      "uqBaseEnvironmentClass::getMyOptionValues()",
                      "total number of processors in environment must be multiple of the specified number of subEnvironments");

  if (m_env.allOptionsMap().count(m_option_subDisplayFileName.c_str())) {
    m_ov.m_subDisplayFileName = m_env.allOptionsMap()[m_option_subDisplayFileName].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_subDisplayAllowAll.c_str())) {
    m_ov.m_subDisplayAllowAll = m_env.allOptionsMap()[m_option_subDisplayAllowAll].as<bool>();
  }

  if (m_ov.m_subDisplayAllowAll) {
    m_ov.m_subDisplayAllowedSet.clear();
    // The line below is commented because 'm_subId' is not set at this point yet
    //m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }
  else if (m_env.allOptionsMap().count(m_option_subDisplayAllowedSet.c_str())) {
    m_ov.m_subDisplayAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_subDisplayAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);
    //if (m_subDisplayOutputFile) {
    //  *m_subDisplayOutputFile << "In uqEnvironmentOptionsClass::getMyOptionValues(): allow = ";
    //  for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
    //    *m_subDisplayOutputFile << " " << tmpAllow[i];
    //  }
    //  *m_subDisplayOutputFile << std::endl;
    //}

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_subDisplayAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_displayVerbosity.c_str())) {
    m_ov.m_displayVerbosity = m_env.allOptionsMap()[m_option_displayVerbosity].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_syncVerbosity.c_str())) {
    m_ov.m_syncVerbosity = m_env.allOptionsMap()[m_option_syncVerbosity].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_checkingLevel.c_str())) {
    m_ov.m_checkingLevel = m_env.allOptionsMap()[m_option_checkingLevel].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_seed.c_str())) {
    m_ov.m_seed = m_env.allOptionsMap()[m_option_seed].as<int>();
  }

  if (m_env.allOptionsMap().count(m_option_platformName.c_str())) {
    m_ov.m_platformName = m_env.allOptionsMap()[m_option_platformName].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_identifyingString.c_str())) {
    m_ov.m_identifyingString = m_env.allOptionsMap()[m_option_identifyingString].as<std::string>();
  }

  //if (m_env.allOptionsMap().count(m_option_numDebugParams.c_str())) {
  //  m_seed = m_env.allOptionsMap()[m_option_numDebugParams].as<unsigned int>();
  //}

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Leaving uqEnvOptionsClass::getMyOptionsValues()" << std::endl;
#endif

  return;
}

void
uqEnvironmentOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_numSubEnvironments   << " = " << m_ov.m_numSubEnvironments
     << "\n" << m_option_subDisplayFileName   << " = " << m_ov.m_subDisplayFileName
     << "\n" << m_option_subDisplayAllowAll   << " = " << m_ov.m_subDisplayAllowAll
     << "\n" << m_option_subDisplayAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_subDisplayAllowedSet.begin(); setIt != m_ov.m_subDisplayAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_displayVerbosity  << " = " << m_ov.m_displayVerbosity
     << "\n" << m_option_syncVerbosity     << " = " << m_ov.m_syncVerbosity
     << "\n" << m_option_checkingLevel     << " = " << m_ov.m_checkingLevel
     << "\n" << m_option_seed              << " = " << m_ov.m_seed
     << "\n" << m_option_platformName      << " = " << m_ov.m_platformName
     << "\n" << m_option_identifyingString << " = " << m_ov.m_identifyingString
   //<< "\n" << m_option_numDebugParams    << " = " << m_ov.m_numDebugParams
     << std::endl;
  return;
}

std::ostream& operator<<(std::ostream& os, const uqEnvironmentOptionsClass& obj)
{
  obj.print(os);

  return os;
}
