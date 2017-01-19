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

#include <queso/Defines.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/EnvironmentOptions.h>
#include <queso/Miscellaneous.h>

#include <queso/asserts.h>

namespace QUESO {

// --------------------------------------------------
// EnvOptionsValues --------------------------
// --------------------------------------------------

// Default constructor ------------------------------
EnvOptionsValues::EnvOptionsValues()
  :
    m_prefix("env_"),
    m_help(UQ_ENV_HELP),
    m_numSubEnvironments(UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
    m_subDisplayFileName(UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),
    m_subDisplayAllowAll(UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),
    m_subDisplayAllowInter0(UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV),
    //m_subDisplayAllowedSet (),
    m_displayVerbosity(UQ_ENV_DISPLAY_VERBOSITY_ODV),
    m_syncVerbosity(UQ_ENV_SYNC_VERBOSITY_ODV),
    m_checkingLevel(UQ_ENV_CHECKING_LEVEL_ODV),
    m_rngType(UQ_ENV_RNG_TYPE_ODV),
    m_seed(UQ_ENV_SEED_ODV),
    m_platformName(UQ_ENV_PLATFORM_NAME_ODV),
    m_identifyingString(UQ_ENV_IDENTIFYING_STRING_ODV),
    m_numDebugParams(UQ_ENV_NUM_DEBUG_PARAMS_ODV),
    m_debugParams(m_numDebugParams,0.),
    m_env(NULL),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help(m_prefix + "help"),
    m_option_numSubEnvironments(m_prefix + "numSubEnvironments"),
    m_option_subDisplayFileName(m_prefix + "subDisplayFileName"),
    m_option_subDisplayAllowAll(m_prefix + "subDisplayAllowAll"),
    m_option_subDisplayAllowInter0(m_prefix + "subDisplayAllowInter0"),
    m_option_subDisplayAllowedSet(m_prefix + "subDisplayAllowedSet"),
    m_option_displayVerbosity(m_prefix + "displayVerbosity"),
    m_option_syncVerbosity(m_prefix + "syncVerbosity"),
    m_option_checkingLevel(m_prefix + "checkingLevel"),
    m_option_rngType(m_prefix + "rngType"),
    m_option_seed(m_prefix + "seed"),
    m_option_platformName(m_prefix + "platformName"),
    m_option_identifyingString(m_prefix + "identifyingString")
{
}

EnvOptionsValues::EnvOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
    m_prefix((std::string) + "env_"),
    m_help(UQ_ENV_HELP),
    m_numSubEnvironments(UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
    m_subDisplayFileName(UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),
    m_subDisplayAllowAll(UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),
    m_subDisplayAllowInter0(UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV),
    //m_subDisplayAllowedSet (),
    m_displayVerbosity(UQ_ENV_DISPLAY_VERBOSITY_ODV),
    m_syncVerbosity(UQ_ENV_SYNC_VERBOSITY_ODV),
    m_checkingLevel(UQ_ENV_CHECKING_LEVEL_ODV),
    m_rngType(UQ_ENV_RNG_TYPE_ODV),
    m_seed(UQ_ENV_SEED_ODV),
    m_platformName(UQ_ENV_PLATFORM_NAME_ODV),
    m_identifyingString(UQ_ENV_IDENTIFYING_STRING_ODV),
    m_numDebugParams(UQ_ENV_NUM_DEBUG_PARAMS_ODV),
    m_debugParams(m_numDebugParams,0.),
    m_env(env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help(m_prefix + "help"),
    m_option_numSubEnvironments(m_prefix + "numSubEnvironments"),
    m_option_subDisplayFileName(m_prefix + "subDisplayFileName"),
    m_option_subDisplayAllowAll(m_prefix + "subDisplayAllowAll"),
    m_option_subDisplayAllowInter0(m_prefix + "subDisplayAllowInter0"),
    m_option_subDisplayAllowedSet(m_prefix + "subDisplayAllowedSet"),
    m_option_displayVerbosity(m_prefix + "displayVerbosity"),
    m_option_syncVerbosity(m_prefix + "syncVerbosity"),
    m_option_checkingLevel(m_prefix + "checkingLevel"),
    m_option_rngType(m_prefix + "rngType"),
    m_option_seed(m_prefix + "seed"),
    m_option_platformName(m_prefix + "platformName"),
    m_option_identifyingString(m_prefix + "identifyingString")
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  // Register all options with parser
  m_parser->registerOption<std::string >(m_option_help, UQ_ENV_HELP, "produce help message for environment");
  m_parser->registerOption<unsigned int>(m_option_numSubEnvironments, UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV, "number of subEnvironments");
  m_parser->registerOption<std::string >(m_option_subDisplayFileName, UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV, "output filename for subscreen writing");
  m_parser->registerOption<bool        >(m_option_subDisplayAllowAll, UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV, "Allow all processors to write to output file");
  m_parser->registerOption<bool        >(m_option_subDisplayAllowInter0, UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV, "Allow all inter0 nodes to write to output file");
  m_parser->registerOption<std::string >(m_option_subDisplayAllowedSet, UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV, "subEnvs that will write to output file");
  m_parser->registerOption<unsigned int>(m_option_displayVerbosity, UQ_ENV_DISPLAY_VERBOSITY_ODV, "set verbosity");
  m_parser->registerOption<unsigned int>(m_option_syncVerbosity, UQ_ENV_SYNC_VERBOSITY_ODV, "set sync verbosity");
  m_parser->registerOption<unsigned int>(m_option_checkingLevel, UQ_ENV_CHECKING_LEVEL_ODV, "set checking level");
  m_parser->registerOption<std::string >(m_option_rngType, UQ_ENV_RNG_TYPE_ODV, "set rngType");
  m_parser->registerOption<int         >(m_option_seed, UQ_ENV_SEED_ODV, "set seed");
  m_parser->registerOption<std::string >(m_option_platformName, UQ_ENV_PLATFORM_NAME_ODV, "platform name");
  m_parser->registerOption<std::string >(m_option_identifyingString, UQ_ENV_IDENTIFYING_STRING_ODV, "identifying string");

  // Read the input file
  m_parser->scanInputFile();

  m_parser->getOption<std::string >(m_option_help, m_help);
  m_parser->getOption<unsigned int>(m_option_numSubEnvironments, m_numSubEnvironments);
  m_parser->getOption<std::string>(m_option_subDisplayFileName, m_subDisplayFileName);
  m_parser->getOption<bool>(m_option_subDisplayAllowAll, m_subDisplayAllowAll);
  m_parser->getOption<bool>(m_option_subDisplayAllowInter0, m_subDisplayAllowInter0);
  m_parser->getOption<std::set<unsigned int> >(m_option_subDisplayAllowedSet, m_subDisplayAllowedSet);
  m_parser->getOption<unsigned int>(m_option_displayVerbosity, m_displayVerbosity);
  m_parser->getOption<unsigned int>(m_option_syncVerbosity, m_syncVerbosity);
  m_parser->getOption<unsigned int>(m_option_checkingLevel, m_checkingLevel);
  m_parser->getOption<std::string>(m_option_rngType, m_rngType);
  m_parser->getOption<int>(m_option_seed, m_seed);
  m_parser->getOption<std::string>(m_option_platformName, m_platformName);
  m_parser->getOption<std::string>(m_option_identifyingString, m_identifyingString);
#else
  m_help = m_env->input()(m_option_help, UQ_ENV_HELP);
  m_numSubEnvironments = m_env->input()(m_option_numSubEnvironments, UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV);
  m_subDisplayFileName = m_env->input()(m_option_subDisplayFileName, UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV);
  m_subDisplayAllowAll = m_env->input()(m_option_subDisplayAllowAll, UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV);
  m_subDisplayAllowInter0 = m_env->input()(m_option_subDisplayAllowInter0, UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV);

  // UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = m_env->input().vector_variable_size(m_option_subDisplayAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env->input()(m_option_subDisplayAllowedSet, i, i);
    m_subDisplayAllowedSet.insert(allowed);
  }

  m_displayVerbosity = m_env->input()(m_option_displayVerbosity, UQ_ENV_DISPLAY_VERBOSITY_ODV);
  m_syncVerbosity = m_env->input()(m_option_syncVerbosity, UQ_ENV_SYNC_VERBOSITY_ODV);
  m_checkingLevel = m_env->input()(m_option_checkingLevel, UQ_ENV_CHECKING_LEVEL_ODV);
  m_rngType = m_env->input()(m_option_rngType, UQ_ENV_RNG_TYPE_ODV);
  m_seed = m_env->input()(m_option_seed, UQ_ENV_SEED_ODV);
  m_platformName = m_env->input()(m_option_platformName, UQ_ENV_PLATFORM_NAME_ODV);
  m_identifyingString = m_env->input()(m_option_identifyingString, UQ_ENV_IDENTIFYING_STRING_ODV);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

void
EnvOptionsValues::checkOptions()
{
  // Clear the permitted set of ranks if the user specifies that all are
  // allowed to display or the Inter0 communicator is allowed to display
  if (m_subDisplayAllowAll) {
    // This will get filled by the Environment after the sub communicators are
    // created
    m_subDisplayAllowedSet.clear();
  }
  else if (m_subDisplayAllowInter0) {
    m_subDisplayAllowedSet.clear();
  }
}

// Copy constructor ---------------------------------
EnvOptionsValues::EnvOptionsValues(const EnvOptionsValues& src)
{
  this->copy(src);
}

// Destructor ---------------------------------------
EnvOptionsValues::~EnvOptionsValues()
{
}

// Set methods---------------------------------------
EnvOptionsValues&
EnvOptionsValues::operator=(const EnvOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}
// Private methods-----------------------------------
void
EnvOptionsValues::copy(const EnvOptionsValues& src)
{
  m_numSubEnvironments    = src.m_numSubEnvironments;
  m_subDisplayFileName    = src.m_subDisplayFileName;
  m_subDisplayAllowAll    = src.m_subDisplayAllowAll;
  m_subDisplayAllowInter0 = src.m_subDisplayAllowInter0;
  m_subDisplayAllowedSet  = src.m_subDisplayAllowedSet;
  m_displayVerbosity      = src.m_displayVerbosity;
  m_syncVerbosity         = src.m_syncVerbosity;
  m_checkingLevel         = src.m_checkingLevel;
  m_rngType               = src.m_rngType;
  m_seed                  = src.m_seed;
  m_platformName          = src.m_platformName;
  m_identifyingString     = src.m_identifyingString;
  m_numDebugParams        = src.m_numDebugParams;
  m_debugParams           = src.m_debugParams;

  return;
}

std::ostream& operator<<(std::ostream& os, const EnvOptionsValues & obj)
{
  // Print the parser
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#else
  obj.m_env->input().print(os);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  // Print the option names and current values
  os <<         obj.m_option_numSubEnvironments    << " = " << obj.m_numSubEnvironments
     << "\n" << obj.m_option_subDisplayFileName    << " = " << obj.m_subDisplayFileName
     << "\n" << obj.m_option_subDisplayAllowAll    << " = " << obj.m_subDisplayAllowAll
   //<< "\n" << obj.m_option_subDisplayAllowInter0 << " = " << obj.m_subDisplayAllowInter0
     << "\n" << obj.m_option_subDisplayAllowedSet  << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_subDisplayAllowedSet.begin(); setIt != obj.m_subDisplayAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << obj.m_option_displayVerbosity  << " = " << obj.m_displayVerbosity
     << "\n" << obj.m_option_syncVerbosity     << " = " << obj.m_syncVerbosity
     << "\n" << obj.m_option_checkingLevel     << " = " << obj.m_checkingLevel
     << "\n" << obj.m_option_rngType           << " = " << obj.m_rngType
     << "\n" << obj.m_option_seed              << " = " << obj.m_seed
     << "\n" << obj.m_option_platformName      << " = " << obj.m_platformName
     << "\n" << obj.m_option_identifyingString << " = " << obj.m_identifyingString
   //<< "\n" << obj.m_option_numDebugParams    << " = " << obj.m_numDebugParams
     << std::endl;
  return os;
}

// --------------------------------------------------
// EnvironmentOptions ------------------------
// --------------------------------------------------

// Default constructor ------------------------------
EnvironmentOptions::EnvironmentOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_ov                          (),
  m_env                         (env),
  m_prefix                      ((std::string)(prefix) + "env_"),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                 (new boost::program_options::options_description("Environment options")),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                 (m_prefix + "help"                 ),
  m_option_numSubEnvironments   (m_prefix + "numSubEnvironments"   ),
  m_option_subDisplayFileName   (m_prefix + "subDisplayFileName"   ),
  m_option_subDisplayAllowAll   (m_prefix + "subDisplayAllowAll"   ),
  m_option_subDisplayAllowInter0(m_prefix + "subDisplayAllowInter0"),
  m_option_subDisplayAllowedSet (m_prefix + "subDisplayAllowedSet" ),
  m_option_displayVerbosity     (m_prefix + "displayVerbosity"     ),
  m_option_syncVerbosity        (m_prefix + "syncVerbosity"        ),
  m_option_checkingLevel        (m_prefix + "checkingLevel"        ),
  m_option_rngType              (m_prefix + "rngType"              ),
  m_option_seed                 (m_prefix + "seed"                 ),
  m_option_platformName         (m_prefix + "platformName"         ),
  m_option_identifyingString    (m_prefix + "identifyingString"    )
{
  queso_deprecated();
  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the abscense of an options input file"));
}
// Constructor with alternative values --------------
EnvironmentOptions::EnvironmentOptions(
  const BaseEnvironment&  env,
  const char*                    prefix,
  const EnvOptionsValues& alternativeOptionsValues)
  :
  m_ov                          (alternativeOptionsValues),
  m_env                         (env),
  m_prefix                      ((std::string)(prefix) + "env_"),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                 (),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                 (m_prefix + "help"                 ),
  m_option_numSubEnvironments   (m_prefix + "numSubEnvironments"   ),
  m_option_subDisplayFileName   (m_prefix + "subDisplayFileName"   ),
  m_option_subDisplayAllowAll   (m_prefix + "subDisplayAllowAll"   ),
  m_option_subDisplayAllowInter0(m_prefix + "subDisplayAllowInter0"),
  m_option_subDisplayAllowedSet (m_prefix + "subDisplayAllowedSet" ),
  m_option_displayVerbosity     (m_prefix + "displayVerbosity"     ),
  m_option_syncVerbosity        (m_prefix + "syncVerbosity"        ),
  m_option_checkingLevel        (m_prefix + "checkingLevel"        ),
  m_option_rngType              (m_prefix + "rngType"              ),
  m_option_seed                 (m_prefix + "seed"                 ),
  m_option_platformName         (m_prefix + "platformName"         ),
  m_option_identifyingString    (m_prefix + "identifyingString"    )
{
  queso_deprecated();
  queso_require_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the existence of an options input file"));

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In EnvironmentOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}
// Destructor ---------------------------------------
EnvironmentOptions::~EnvironmentOptions()
{
  queso_deprecated();
}

// I/O methods---------------------------------------
void
EnvironmentOptions::scanOptionsValues()
{
  queso_deprecated();
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  // 'm_subDisplayOutputFile' is still not available at this moment. Use 'std::cout'
  //if (m_env.subScreenFile() != NULL) {
  //  *m_env.subScreenFile()
  if ((m_env.fullRank() == 0) && (m_env.displayVerbosity() >= 3)) {
    std::cout << "In EnvironmentOptions::scanOptionsValues()"
              << ": after reading values of options with prefix '" << m_prefix
              << "', state of object is:"
              << "\n" << *this
              << std::endl;
  }

  return;
}
// --------------------------------------------------
void
EnvironmentOptions::print(std::ostream& os) const
{
  queso_deprecated();

  os <<         m_option_numSubEnvironments    << " = " << m_ov.m_numSubEnvironments
     << "\n" << m_option_subDisplayFileName    << " = " << m_ov.m_subDisplayFileName
     << "\n" << m_option_subDisplayAllowAll    << " = " << m_ov.m_subDisplayAllowAll
   //<< "\n" << m_option_subDisplayAllowInter0 << " = " << m_ov.m_subDisplayAllowInter0
     << "\n" << m_option_subDisplayAllowedSet  << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_subDisplayAllowedSet.begin(); setIt != m_ov.m_subDisplayAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_displayVerbosity  << " = " << m_ov.m_displayVerbosity
     << "\n" << m_option_syncVerbosity     << " = " << m_ov.m_syncVerbosity
     << "\n" << m_option_checkingLevel     << " = " << m_ov.m_checkingLevel
     << "\n" << m_option_rngType           << " = " << m_ov.m_rngType
     << "\n" << m_option_seed              << " = " << m_ov.m_seed
     << "\n" << m_option_platformName      << " = " << m_ov.m_platformName
     << "\n" << m_option_identifyingString << " = " << m_ov.m_identifyingString
   //<< "\n" << m_option_numDebugParams    << " = " << m_ov.m_numDebugParams
     << std::endl;
  return;
}

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
// Private methods ----------------------------------
void
EnvironmentOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In EnvOptions::defineMyOptions(), before add_options()" << std::endl;
#endif
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                 "produce help message for  environment"       )
    (m_option_numSubEnvironments.c_str(),    boost::program_options::value<unsigned int>()->default_value(UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),     "number of subEnvironments"                     )
    (m_option_subDisplayFileName.c_str(),    boost::program_options::value<std::string >()->default_value(UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),    "output filename for subscreen writing"         )
    (m_option_subDisplayAllowAll.c_str(),    boost::program_options::value<bool        >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),    "Allow all processors to write to output file"  )
    (m_option_subDisplayAllowInter0.c_str(), boost::program_options::value<bool        >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV), "Allow all inter0 nodes to write to output file")
    (m_option_subDisplayAllowedSet.c_str(),  boost::program_options::value<std::string >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV),  "subEnvs that will write to output file"        )
    (m_option_displayVerbosity.c_str(),      boost::program_options::value<unsigned int>()->default_value(UQ_ENV_DISPLAY_VERBOSITY_ODV),        "set verbosity"                                 )
    (m_option_syncVerbosity.c_str(),         boost::program_options::value<unsigned int>()->default_value(UQ_ENV_SYNC_VERBOSITY_ODV),           "set sync verbosity"                            )
    (m_option_checkingLevel.c_str(),         boost::program_options::value<unsigned int>()->default_value(UQ_ENV_CHECKING_LEVEL_ODV),           "set checking level"                            )
    (m_option_rngType.c_str(),               boost::program_options::value<std::string >()->default_value(UQ_ENV_RNG_TYPE_ODV),                 "set rngType"                                   )
    (m_option_seed.c_str(),                  boost::program_options::value<int         >()->default_value(UQ_ENV_SEED_ODV),                     "set seed"                                      )
    (m_option_platformName.c_str(),          boost::program_options::value<std::string >()->default_value(UQ_ENV_PLATFORM_NAME_ODV),            "platform name"                                 )
    (m_option_identifyingString.c_str(),     boost::program_options::value<std::string >()->default_value(UQ_ENV_IDENTIFYING_STRING_ODV),       "identifying string"                            )
  //(m_option_numDebugParams.c_str(),        boost::program_options::value<unsigned int>()->default_value(UQ_ENV_NUM_DEBUG_PARAMS_ODV),         "set number of debug parameters"                )
  ;
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In EnvOptions::defineMyOptions(), after add_options()" << std::endl;
#endif

  return;
}
// --------------------------------------------------
void
EnvironmentOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering EnvOptions::getMyOptionsValues()" << std::endl;
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
    std::cerr << "In BaseEnvironment::getMyOptionValues()"
              << ": m_env.fullComm().NumProc() = " << m_env.fullComm().NumProc()
              << ", m_numSubEnvironments = "       << m_ov.m_numSubEnvironments
              << std::endl;
  }
  queso_require_equal_to_msg((m_env.fullComm().NumProc()%m_ov.m_numSubEnvironments), 0, "total number of processors in environment must be multiple of the specified number of subEnvironments");

  if (m_env.allOptionsMap().count(m_option_subDisplayFileName.c_str())) {
    m_ov.m_subDisplayFileName = m_env.allOptionsMap()[m_option_subDisplayFileName].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_subDisplayAllowAll.c_str())) {
    m_ov.m_subDisplayAllowAll = m_env.allOptionsMap()[m_option_subDisplayAllowAll].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_subDisplayAllowInter0.c_str())) {
    m_ov.m_subDisplayAllowInter0 = m_env.allOptionsMap()[m_option_subDisplayAllowInter0].as<bool>();
  }

  if (m_ov.m_subDisplayAllowAll) {
    m_ov.m_subDisplayAllowedSet.clear();
    // The line below is commented because 'm_subId' is not set at this point yet
    //m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }
  else if (m_ov.m_subDisplayAllowInter0) {
    m_ov.m_subDisplayAllowedSet.clear();
  }
  else if (m_env.allOptionsMap().count(m_option_subDisplayAllowedSet.c_str())) {
    m_ov.m_subDisplayAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_subDisplayAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);
    //if (m_subDisplayOutputFile) {
    //  *m_subDisplayOutputFile << "In EnvironmentOptions::getMyOptionValues(): allow = ";
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

  if (m_env.allOptionsMap().count(m_option_rngType.c_str())) {
    m_ov.m_rngType = m_env.allOptionsMap()[m_option_rngType].as<std::string>();
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
  //  m_numDebugParams = m_env.allOptionsMap()[m_option_numDebugParams].as<unsigned int>();
  //}

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Leaving EnvOptions::getMyOptionsValues()" << std::endl;
#endif

  return;
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

// Operator outside class definition ----------------
std::ostream& operator<<(std::ostream& os, const EnvironmentOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
