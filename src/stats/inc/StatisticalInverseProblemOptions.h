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

#ifndef UQ_SIP_OPTIONS_H
#define UQ_SIP_OPTIONS_H

#undef UQ_SIP_READS_SOLVER_OPTION

#define UQ_SIP_HELP ""
#define UQ_SIP_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_SIP_COMPUTE_SOLUTION_ODV        1
#define UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV   UQ_SIP_FILENAME_FOR_NO_FILE
#define UQ_SIP_DATA_OUTPUT_ALLOWED_SET_ODV ""
#ifdef UQ_SIP_READS_SOLVER_OPTION
#define UQ_SIP_SOLVER_ODV                  "bayes_mc" // Bayesian formula + Metropolis-Hastings
#endif

#define UQ_SIP_SEEDWITHMAPESTIMATOR 0
#define UQ_SIP_USEOPTIMIZERMONITOR 1

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
namespace boost {
  namespace program_options {
    class options_description;
  }
}
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

namespace QUESO {

/*!
 * \file StatisticalInverseProblemOptions.h
 * \brief Classes to allow options to be passed to a Statistical Inverse Problem
 */

/*!
 * \class SipOptionsValues
 * \brief This class provides options for a Statistical Inverse Problem if no input file is available.
 *
 * In order to solve a Statistical Inverse Problem (SIP), QUESO expects some
 * options for its methods to be fully defined.  This class provides default
 * values for such options if no input file is available.
 */

class SipOptionsValues
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Assigns the default suite of options to the Statistical Inverse Problem.*/
  SipOptionsValues            ();
  SipOptionsValues(const BaseEnvironment * env, const char * prefix);

  //! Copy constructor.
  /*! It assigns the same options values from  \c src to \c this.*/
  SipOptionsValues            (const SipOptionsValues& src);

  //! Set parameter option names to begin with prefix
  void set_prefix(const std::string& prefix);

  //! Set default values for parameter options
  void set_defaults();

  //! Given prefix, read the input file for parameters named "prefix"+*
  void parse(const BaseEnvironment& env, const std::string& prefix);

  //! Destructor
  virtual ~SipOptionsValues            ();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator; it copies \c rhs to \c this.
  SipOptionsValues& operator= (const SipOptionsValues& rhs);
  //@}

  std::string m_prefix;

  //! If this string is non-empty, options are print to the output file
  std::string m_help;

  bool                   m_computeSolution;
  std::string            m_dataOutputFileName;
  std::set<unsigned int> m_dataOutputAllowedSet;
#ifdef UQ_SIP_READS_SOLVER_OPTION
  std::string            m_solverString;
#endif

  bool m_seedWithMAPEstimator;
  bool m_useOptimizerMonitor;

private:
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  // The input options as strings so we can parse the input file later
  std::string                   m_option_help;
  std::string                   m_option_computeSolution;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
#ifdef UQ_SIP_READS_SOLVER_OPTION
  std::string                   m_option_solver;
#endif
  std::string m_option_seedWithMAPEstimator;
  std::string m_option_useOptimizerMonitor;

  //! Copies the option values from \c src to \c this.
  void copy(const SipOptionsValues& src);

  void checkOptions();

  friend std::ostream & operator<<(std::ostream & os,
      const SipOptionsValues & obj);
};

}  // End namespace QUESO

#endif // UQ_SIP_OPTIONS_H
