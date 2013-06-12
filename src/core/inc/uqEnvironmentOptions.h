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

#ifndef __UQ_ENVIRONMENT_OPTIONS_H__
#define __UQ_ENVIRONMENT_OPTIONS_H__

#include <uqEnvironment.h>

#define UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE "."
#define UQ_ENV_FILENAME_FOR_NO_INPUT_FILE  "."

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


/*! \file uqEnvironmentOptions.h
    \brief Class to allow options to be passed to a QUESO environment.
*/

/*! \class uqEnvironmentOptionsClass
 *  \brief This class reads options one can pass to a QUESO environment through an input file.
 * 
 *  QUESO expects the user to provide an input file with environment options for the library variables. 
 *  This class reads the input options for QUESO environment variables. */
 
class uqEnvironmentOptionsClass
{
public:
      //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor
  /*! Assigns the default suite of options to the environment.*/
  uqEnvironmentOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
  
  //! Constructor with alternative options values.
  uqEnvironmentOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix, const uqEnvOptionsValuesClass& alternativeOptionsValues);
  
  //! Destructor
 ~uqEnvironmentOptionsClass();
 //@}

 //! @name I/O methods
 //@{
 //! Scans option values from input file.
  void scanOptionsValues();
  
  //! Print values of the options chosen.
  void print            (std::ostream& os) const;
  //@}
  
  
  //! Instance of uqEnvOptionsValuesClass, a class with default values for QUESO environment.
  uqEnvOptionsValuesClass  m_ov;

private:
  //! Define my environment options as the default options
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  
  //! Gets the option values of the environment.
  void   getMyOptionValues(po::options_description& optionsDesc);

  //! Environment.
  const uqBaseEnvironmentClass& m_env;
  
  //! Options prefix.
  std::string              m_prefix;
  
  //! Environment options description.
  /*! Uses boost::program_options::options_description. A set of option descriptions. 
   * This provides convenient interface for adding new option method, and facilities 
   * to search for options by name.*/
  po::options_description* m_optionsDesc;

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
std::ostream& operator<<(std::ostream& os, const uqEnvironmentOptionsClass& obj);

#endif // __UQ_ENVIRONMENT_CLASS_H__
