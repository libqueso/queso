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

#ifndef UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H
#define UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H

#include <queso/Environment.h>
#include <queso/MLSamplingLevelOptions.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#define UQ_ML_SAMPLING_FILENAME_FOR_NO_FILE "."

// _ODV = option default value

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY

#define UQ_ML_SAMPLING_HELP                                    ""
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

/*! \file MLSamplingOptions.h
    \brief Classes to allow options to be passed to the Multilevel algorithm.
*/

/*! \class MLSamplingOptions
 *  \brief This class provides options for the Multilevel sequence generator if no input file is available.
 *
 *  Multilevel sequence generator expects options for its methods. This class provides default
 *  values for such options if no input file is available. */

class MLSamplingOptions
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Assigns the default suite of options to the Multilevel sequence generator.*/
  MLSamplingOptions(const BaseEnvironment& env, const char* prefix);

  //! Destructor
  virtual ~MLSamplingOptions();
  //@}

  //! @name I/O methods
  //@{
  //!  It prints the option values.
  void print            (std::ostream& os) const;
  //@}

  //! Class prefix. (ml)
  std::string            m_prefix;

  //! If non-empty string, options and values are printed to the output file
  std::string m_help;

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
  //! Period of restart output file (level).
  unsigned int           m_restartOutput_levelPeriod;

  //! Base name of restart output file.
  std::string            m_restartOutput_baseNameForFiles;

  //! Type of restart output file.
  std::string            m_restartOutput_fileType;

  //! Base name of restart input file.
  std::string            m_restartInput_baseNameForFiles;

  //! Type of restart input file
  std::string            m_restartInput_fileType;
#else
  //! Name of restart input file.
  std::string            m_restartInputFileName;

  //! Type of restart input file
  std::string            m_restartInputFileType;

  //! Size of restart chain
  unsigned int           m_restartChainSize;
#endif
  //! Name of generic output file
  std::string            m_dataOutputFileName;

  //! subEnvs that will write to generic output file
  bool                   m_dataOutputAllowAll;
  std::set<unsigned int> m_dataOutputAllowedSet;

private:
  const BaseEnvironment& m_env;

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  BoostInputOptionsParser * m_parser;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

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

  void checkOptions(const BaseEnvironment * env);

  friend std::ostream & operator<<(std::ostream & os,
      const MLSamplingOptions & obj);
};


}  // End namespace QUESO

#endif // UQ_MULTI_LEVEL_SAMPLING_OPTIONS_H
