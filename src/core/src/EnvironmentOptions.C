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

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#define GETPOT_NAMESPACE QUESO // So we don't clash with other getpots
#include <queso/getpot.h>
#undef GETPOT_NAMESPACE
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

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
  m_env(NULL)
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  , m_parser(new BoostInputOptionsParser())
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
{
  this->set_defaults();
  this->set_prefix("");
}


EnvOptionsValues::EnvOptionsValues(const BaseEnvironment* env,
				   const char* prefix)
{
  this->set_defaults();
  this->parse(*env, prefix);
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
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#else
  obj.m_env->input().print(os);
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

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


void
EnvOptionsValues::set_defaults()
{
  m_help = UQ_ENV_HELP;
  m_numSubEnvironments = UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV;
  m_subDisplayFileName = UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV;
  m_subDisplayAllowAll = UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV;
  m_subDisplayAllowInter0 = UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV;
  //m_subDisplayAllowedSet  = ;
  m_displayVerbosity = UQ_ENV_DISPLAY_VERBOSITY_ODV;
  m_syncVerbosity = UQ_ENV_SYNC_VERBOSITY_ODV;
  m_checkingLevel = UQ_ENV_CHECKING_LEVEL_ODV;
  m_rngType = UQ_ENV_RNG_TYPE_ODV;
  m_seed = UQ_ENV_SEED_ODV;
  m_platformName = UQ_ENV_PLATFORM_NAME_ODV;
  m_identifyingString = UQ_ENV_IDENTIFYING_STRING_ODV;
  m_numDebugParams = UQ_ENV_NUM_DEBUG_PARAMS_ODV;
  m_debugParams.assign(m_numDebugParams, 0.);
}


void
EnvOptionsValues::set_prefix(const std::string& prefix)
{
  m_prefix = prefix + "env_";

  m_option_help = m_prefix + "help";
  m_option_numSubEnvironments = m_prefix + "numSubEnvironments";
  m_option_subDisplayFileName = m_prefix + "subDisplayFileName";
  m_option_subDisplayAllowAll = m_prefix + "subDisplayAllowAll";
  m_option_subDisplayAllowInter0 = m_prefix + "subDisplayAllowInter0";
  m_option_subDisplayAllowedSet = m_prefix + "subDisplayAllowedSet";
  m_option_displayVerbosity = m_prefix + "displayVerbosity";
  m_option_syncVerbosity = m_prefix + "syncVerbosity";
  m_option_checkingLevel = m_prefix + "checkingLevel";
  m_option_rngType = m_prefix + "rngType";
  m_option_seed = m_prefix + "seed";
  m_option_platformName = m_prefix + "platformName";
  m_option_identifyingString = m_prefix + "identifyingString";

}


void
EnvOptionsValues::parse(const BaseEnvironment& env, const std::string& prefix)
{
  m_env = &env;

  if (m_env->optionsInputFileName().empty()) {
    queso_error_msg("Missing input file is required");
  }

  this->set_prefix(prefix);

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  m_parser.reset(new BoostInputOptionsParser(env.optionsInputFileName()));

  // Register all options with parser
  m_parser->registerOption<std::string>
    (m_option_help, m_help, 
    "produce help message for environment");
  m_parser->registerOption<unsigned int>
    (m_option_numSubEnvironments, m_numSubEnvironments,
    "number of subEnvironments");
  m_parser->registerOption<std::string >
    (m_option_subDisplayFileName, m_subDisplayFileName,
    "output filename for subscreen writing");
  m_parser->registerOption<bool>
    (m_option_subDisplayAllowAll, m_subDisplayAllowAll,
    "Allow all processors to write to output file");
  m_parser->registerOption<bool>
    (m_option_subDisplayAllowInter0, m_subDisplayAllowInter0,
    "Allow all inter0 nodes to write to output file");

  // convert the current member set of integers to a string for
  // default handling, omitting trailing whitespace
  std::ostringstream sdas_oss;
  std::set<unsigned int>::const_iterator sdas_it = m_subDisplayAllowedSet.begin();
  for ( ; sdas_it != m_subDisplayAllowedSet.end(); ++sdas_it) {
    if (sdas_it != m_subDisplayAllowedSet.begin())
      sdas_oss << ' ';
    sdas_oss << *sdas_it;
  }
  m_parser->registerOption<std::string>
    (m_option_subDisplayAllowedSet, sdas_oss.str(),
    "subEnvs that will write to output file");

  m_parser->registerOption<unsigned int>
    (m_option_displayVerbosity, m_displayVerbosity,
    "set verbosity");
  m_parser->registerOption<unsigned int>
    (m_option_syncVerbosity, m_syncVerbosity,
    "set sync verbosity");
  m_parser->registerOption<unsigned int>
    (m_option_checkingLevel, m_checkingLevel,
    "set checking level");
  m_parser->registerOption<std::string>
    (m_option_rngType, m_rngType,
    "set rngType");
  m_parser->registerOption<int>
    (m_option_seed, m_seed,
    "set seed");
  m_parser->registerOption<std::string>
    (m_option_platformName, m_platformName,
    "platform name");
  m_parser->registerOption<std::string>
    (m_option_identifyingString, m_identifyingString,
    "identifying string");

  // Read the input file
  m_parser->scanInputFile();

  // Retrieve parsed options
  m_parser->getOption<std::string>(m_option_help, m_help);
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

  m_help = m_env->input()(m_option_help, m_help);

  m_numSubEnvironments =
    m_env->input()(m_option_numSubEnvironments, m_numSubEnvironments);

  m_subDisplayFileName =
    m_env->input()(m_option_subDisplayFileName, m_subDisplayFileName);

  m_subDisplayAllowAll =
    m_env->input()(m_option_subDisplayAllowAll, m_subDisplayAllowAll);

  m_subDisplayAllowInter0 =
    m_env->input()(m_option_subDisplayAllowInter0, m_subDisplayAllowInter0);

  // UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = m_env->input().vector_variable_size(m_option_subDisplayAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env->input()(m_option_subDisplayAllowedSet, i, i);
    m_subDisplayAllowedSet.insert(allowed);
  }

  m_displayVerbosity =
    m_env->input()(m_option_displayVerbosity, m_displayVerbosity);

  m_syncVerbosity = m_env->input()(m_option_syncVerbosity, m_syncVerbosity);

  m_checkingLevel = m_env->input()(m_option_checkingLevel, m_checkingLevel);

  m_rngType = m_env->input()(m_option_rngType, m_rngType);

  m_seed = m_env->input()(m_option_seed, m_seed);

  m_platformName = m_env->input()(m_option_platformName, m_platformName);

  m_identifyingString =
    m_env->input()(m_option_identifyingString, m_identifyingString);

#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

}  // End namespace QUESO
