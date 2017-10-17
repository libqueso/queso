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

#include <queso/config_queso.h>
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/Environment.h>

#ifndef UQ_SFP_OPTIONS_H
#define UQ_SFP_OPTIONS_H

#undef UQ_SFP_READS_SOLVER_OPTION

#define UQ_SFP_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_SFP_HELP        ""
#define UQ_SFP_COMPUTE_SOLUTION_ODV        1
#define UQ_SFP_COMPUTE_COVARIANCES_ODV     1
#define UQ_SFP_COMPUTE_CORRELATIONS_ODV    1
#define UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV   UQ_SFP_FILENAME_FOR_NO_FILE
#define UQ_SFP_DATA_OUTPUT_ALLOWED_SET_ODV ""
#ifdef UQ_SFP_READS_SOLVER_OPTION
#define UQ_SFP_SOLVER_ODV                  "mc" // Monte Carlo
#endif

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
namespace boost {
  namespace program_options {
    class options_description;
  }
}
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

namespace QUESO {

/*! \file StatisticalForwardProblemOptions.h
    \brief Classes to allow options to be passed to a Statistical Forward Problem.
*/

/*! \class SfpOptionsValues
 *  \brief This class provides options for a Statistical Forward Problem if no input file is available.
 *
 * In order to solve a Statistical Forward Problem (SFP), QUESO expects some options for its methods to be
 * fully defined. This class provides default values for such options if no input file is available. */

class SfpOptionsValues
{
public:
  //! Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Assigns the default suite of options to the Statistical Forward Problem.*/
  SfpOptionsValues            ();
  SfpOptionsValues(const BaseEnvironment * env, const char * prefix);

  //! Copy constructor.
  /*! It assigns the same options values from  \c src to \c this.*/
  SfpOptionsValues            (const SfpOptionsValues& src);

  //! Set parameter option names to begin with prefix
  void set_prefix(const std::string& prefix);

  //! Set default values for parameter options
  void set_defaults();

  //! Given prefix, read the input file for parameters named "prefix"+*
  void parse(const BaseEnvironment& env, const std::string& prefix);

  //! Destructor
  virtual ~SfpOptionsValues            ();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator; it copies \c rhs to \c this.
  SfpOptionsValues& operator= (const SfpOptionsValues& rhs);
  //@}

  std::string                   m_prefix;

  //! If non-empty string, options and values are printed to the output file
  std::string m_help;

  bool                   m_computeSolution;
  bool                   m_computeCovariances;
  bool                   m_computeCorrelations;
  std::string            m_dataOutputFileName;
  std::set<unsigned int> m_dataOutputAllowedSet;
#ifdef UQ_SFP_READS_SOLVER_OPTION
  std::string            m_solverString;
#endif

  //McOptionsValues m_mcOptionsValues;

private:
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  // The input options as strings so we can parse the input file later
  std::string                   m_option_help;
  std::string                   m_option_computeSolution;
  std::string                   m_option_computeCovariances;
  std::string                   m_option_computeCorrelations;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
#ifdef UQ_SFP_READS_SOLVER_OPTION
  std::string                   m_option_solver;
#endif

  //! Copies the option values from \c src to \c this.
  void copy(const SfpOptionsValues& src);

  void checkOptions();

  friend std::ostream & operator<<(std::ostream & os,
      const SfpOptionsValues & obj);
};

}  // End namespace QUESO

#endif // UQ_SFP_OPTIONS_H
