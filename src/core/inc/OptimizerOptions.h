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

#ifndef UQ_OPT_OPTIONS_H
#define UQ_OPT_OPTIONS_H

#include <queso/SharedPtr.h>

#define UQ_OPT_HELP ""

// _ODV = option default value
#define UQ_OPT_MAX_ITERATIONS 100
#define UQ_OPT_TOLERANCE 1e-3
#define UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE 1e-4

namespace QUESO {

class BaseEnvironment;
class BoostInputOptionsParser;

/*!
 * \file OptimizerOptions.h
 * \brief Classes to allow options to be passed to an instance of Optimizer
 */

/*!
 * \class OptimizerOptions
 * \brief This class provides options for a Optimizer
 */

class OptimizerOptions
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.  All options have their default values.
  OptimizerOptions();

  //! A constructor that takes the environment for parsing input file options
  OptimizerOptions(const BaseEnvironment * env, const char * prefix = "");

  //! Destructor
  virtual ~OptimizerOptions();
  //@}

  SharedPtr<const OptimizerOptions>::Type clone() const;

  std::string m_prefix;

  //! If this string is non-empty, options are print to the output file
  std::string m_help;

  unsigned int m_maxIterations;
  double m_tolerance;
  double m_finiteDifferenceStepSize;

private:
  BoostInputOptionsParser * m_parser;

  // The input options as strings so we can parse the input file later
  std::string m_option_help;
  std::string m_option_maxIterations;
  std::string m_option_tolerance;
  std::string m_option_finiteDifferenceStepSize;

  void checkOptions();

  friend std::ostream & operator<<(std::ostream & os,
      const OptimizerOptions & obj);
};

}  // End namespace QUESO

#endif // UQ_OPT_OPTIONS_H
