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

#ifndef UQ_ENVIRONMENT_OPTIONS_H
#define UQ_ENVIRONMENT_OPTIONS_H

#include <string>
#include <set>
#include <vector>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/ScopedPtr.h>

#define UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE "."
#define UQ_ENV_FILENAME_FOR_NO_INPUT_FILE  "."

#define UQ_ENV_HELP                         ""
#define UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV     1
#define UQ_ENV_SUB_SCREEN_WRITE_ODV         0
#define UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV    UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV    0
#define UQ_ENV_SUB_DISPLAY_ALLOW_INTER0_ODV 0
#define UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV  ""
#define UQ_ENV_DISPLAY_VERBOSITY_ODV        0
#define UQ_ENV_SYNC_VERBOSITY_ODV           0
#define UQ_ENV_CHECKING_LEVEL_ODV           0
#define UQ_ENV_RNG_TYPE_ODV                 "gsl"
#define UQ_ENV_SEED_ODV                     0
#define UQ_ENV_IDENTIFYING_STRING_ODV       ""
#define UQ_ENV_PLATFORM_NAME_ODV            ""
#define UQ_ENV_NUM_DEBUG_PARAMS_ODV         0
#define UQ_ENV_DEBUG_PARAM_ODV              0.

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
// Forward declarations
namespace boost {
  namespace program_options {
    class options_description;
    }
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

namespace QUESO {

// Forward declarations
class BaseEnvironment;

/*! \file EnvironmentOptions.h
    \brief Class to allow options to be passed to a QUESO environment.
*/

/*! \class EnvOptionsValues
 *  \brief This class provides a suite options one can pass to a QUESO environment.
 *
 *  QUESO expects the user to provide an input file with environment options for the library variables.
 *  If no input file, a collection of default values is assigned to some of the variables. The class
 *  EnvOptionsValues is responsible for this task.
 */

class EnvOptionsValues
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  EnvOptionsValues();
  EnvOptionsValues(const BaseEnvironment * env, const char * prefix);

  //! Copy constructor
  EnvOptionsValues(const EnvOptionsValues& src);

  //! Destructor
  virtual ~EnvOptionsValues();
  //@}

  //! @name Set methods
  //@{
  //! Operator for copying the options of an environment.
  EnvOptionsValues & operator=(const EnvOptionsValues& rhs);
  //@}

  std::string m_prefix;

  //! @name Attributes
  //! If this string is non-empty, print the options object to the output file
  std::string m_help;

  //! Number of sub-environments (chains).  Each chain may have multiple
  //! processes.
  unsigned int m_numSubEnvironments;

  //! Output filename for sub-screen writing.
  std::string m_subDisplayFileName;

  //! Allows (or not) all sub-environments to write to output file.
  /*!
   * If this option is true, m_subDisplayAllowedSet is ignored.
   */
  bool m_subDisplayAllowAll;

  //! Allows (or not) all inter0 nodes to write to output file
  /*!
   * If this option is true, m_subDisplayAllowedSet is ignored.
   */
  bool m_subDisplayAllowInter0;

  //! Sub-environments that will write to output.
  /*!
   * This option is ignored if either m_subDisplayAllowAll or
   * m_subDisplayAllowInter0 are true.
   */
  std::set<unsigned int> m_subDisplayAllowedSet;

  //! Verbosity.
  unsigned int m_displayVerbosity;

  //! Synchronized verbosity.
  unsigned int m_syncVerbosity;

  //! Checking level
  unsigned int m_checkingLevel;

  //! Type of the random number generator.
  std::string m_rngType;

  //! Seed of the random number generator.
  /*!
   * If env_seed = -z, with z>=1, then each processor sets the seed to value
   * MPI_RANK + z.  It is crucial that \verb+env_seed+ takes a
   * \underline{negative} value, otherwise all chain samples are going to be
   * the same.
   */
  int m_seed;

  //! Platform name.
  std::string m_platformName;

  //! Identifying string.
  std::string m_identifyingString;

  //! Number of debug parameters.  Unused?
  unsigned int m_numDebugParams;

  //! Debug parameters.  Unused?
  std::vector<double> m_debugParams;
  //@}

private:
  const BaseEnvironment * m_env;

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  //! Input file option name for flagging helpful printing output
  std::string m_option_help;

  //! Input file option name for m_numSubEnvironments.
  std::string m_option_numSubEnvironments;

  //! Input file option name for m_subDisplayFileName
  std::string m_option_subDisplayFileName;

  //! Input file option name for m_subDisplayAllowAll
  std::string m_option_subDisplayAllowAll;

  //! Input file option name for m_subDisplayAllowInter0
  std::string m_option_subDisplayAllowInter0;

  //! Input file option name for m_subDisplayAllowedSet
  std::string m_option_subDisplayAllowedSet;

  //! Input file option name for m_displayVerbosity
  std::string m_option_displayVerbosity;

  //! Input file option name for m_syncVerbosity
  std::string m_option_syncVerbosity;

  //! Input file option name for m_checkingLevel
  std::string m_option_checkingLevel;

  //! Input file option name for m_rngType
  std::string m_option_rngType;

  //! Input file option name for m_seed
  std::string m_option_seed;

  //! Input file option name for m_platformName
  std::string m_option_platformName;

  //! Input file option name for m_identifyingString
  std::string m_option_identifyingString;

  //! Makes an exact copy of an existing EnvOptionsValues instance.
  void copy(const EnvOptionsValues& src);

  //! Sorts out any inter-option conflicts
  void checkOptions();

  //! Print values of the options chosen.
  friend std::ostream& operator<<(std::ostream& os,
      const EnvOptionsValues & obj);
};

/*! \class EnvironmentOptions
 *  \brief This class reads options one can pass to a QUESO environment through an input file.
 *
 *  QUESO expects the user to provide an input file with environment options for the library variables.
 *  This class reads the input options for QUESO environment variables. */

class EnvironmentOptions
{
public:
      //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Assigns the default suite of options to the environment.*/
  EnvironmentOptions(const BaseEnvironment& env, const char* prefix);

  //! Constructor with alternative options values.
  EnvironmentOptions(const BaseEnvironment& env, const char* prefix, const EnvOptionsValues& alternativeOptionsValues);

  //! Destructor
 ~EnvironmentOptions();
 //@}

 //! @name I/O methods
 //@{
 //! Scans option values from input file.
  void scanOptionsValues();

  //! Print values of the options chosen.
  void print            (std::ostream& os) const;
  //@}


  //! Instance of EnvOptionsValues, a class with default values for QUESO environment.
  EnvOptionsValues  m_ov;

private:
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  //! Define my environment options as the default options
  void   defineMyOptions  (boost::program_options::options_description& optionsDesc) const;

  //! Gets the option values of the environment.
  void   getMyOptionValues(boost::program_options::options_description& optionsDesc);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  //! Environment.
  const BaseEnvironment& m_env;

  //! Options prefix.
  std::string              m_prefix;

  //! Environment options description.
  /*! Uses boost::program_options::options_description. A set of option descriptions.
   * This provides convenient interface for adding new option method, and facilities
   * to search for options by name.*/
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<boost::program_options::options_description>::Type m_optionsDesc;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  std::string              m_option_help;

  //! My number of sub-environments.
  std::string              m_option_numSubEnvironments;

  //! My output filename for sub-screen writing.
  std::string              m_option_subDisplayFileName;

  //! Allows (or not) all sub-environments to write to output file.
  std::string              m_option_subDisplayAllowAll;

  //! Allows (or not) all inter0 nodes to write to output file
  std::string              m_option_subDisplayAllowInter0;

  //! Sub-environments that will write to output.
  std::string              m_option_subDisplayAllowedSet;

  //! Verbosity.
  std::string              m_option_displayVerbosity;

  //! Synchronized verbosity.
  std::string              m_option_syncVerbosity;

  //! Checking level
  std::string              m_option_checkingLevel;

  //! Type of the random number generator.
  std::string              m_option_rngType;

  //! Seed of the random number generator.
  /*! If env_seed = -z, with z>=1, then each processor sets the seed to value MPI_RANK + z.
   It is crucial that \verb+env_seed+ takes a \underline{negative} value, otherwise all chain samples are going to be the same.*/
  std::string              m_option_seed;

  //! Platform name.
  std::string              m_option_platformName;

  //! Identifying string.
  std::string              m_option_identifyingString;
};

//! Print values of the options chosen.
std::ostream& operator<<(std::ostream& os, const EnvironmentOptions& obj);

}  // End namespace QUESO

#endif // UQ_ENVIRONMENT_CLASS_H
