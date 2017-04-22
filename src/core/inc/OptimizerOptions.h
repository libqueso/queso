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

#ifndef UQ_OPT_OPTIONS_H
#define UQ_OPT_OPTIONS_H

// C++
#include <string>

// QUESO
#include <queso/SharedPtr.h>

#define UQ_OPT_HELP ""

// _ODV = option default value
#define UQ_OPT_MAX_ITERATIONS 100
#define UQ_OPT_TOLERANCE 1e-3
#define UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE 1e-4
#define UQ_OPT_SOLVER_TYPE "bfgs2"
#define UQ_OPT_FSTEP_SIZE 0.1
#define UQ_OPT_FDFSTEP_SIZE 1.0
#define UQ_OPT_LINE_TOLERANCE 0.1

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

  //! Copy constructor
  OptimizerOptions(const OptimizerOptions & rhs);

  //! Destructor
  virtual ~OptimizerOptions();
  //@}

  std::string m_prefix;

  //! If this string is non-empty, options are print to the output file
  std::string m_help;

  //! The maximum number of iterations to do for optimization.  Default is 100.
  unsigned int m_maxIterations;

  //! The tolerance at which optimization stops.  Default is 1e-3.
  /*!
   * The tolerance is measured in terms of the 2-norm of the gradient of the
   * objective function for optimizers that use derivatives.
   *
   * For optimizers that do not use derivatives, the tolerance is measured in
   * terms of the size of the size of the simplex for the simplex optimizers.
   */
  double m_tolerance;

  //! The step size used to compute gradients from finite differencing.  Default is 1e-4.
  /*!
   * This is only done when analytical gradients are not provided (either by
   * the user or by QUESO).
   */
  double m_finiteDifferenceStepSize;

  //! The optimization algorithm to use.  Default is bfgs2.
  /*!
   *  Choices are:
   *    - fletcher_reeves_cg
   *    - polak_ribiere_cg
   *    - bfgs
   *    - bfgs2
   *    - steepest_descent
   *    - nelder_mead
   *    - nelder_mead2
   *    - nelder_mead2_rand
   */
  std::string m_solverType;

  //! The size of the initial trial steps for optimizing without gradients.  Default is 0.1
  /*!
   * The initial trial steps is actually a vector of size equal to the number
   * of parameters, and this vector is set to \c m_fstepSize in each
   * component.
   */
  double m_fstepSize;

  //! The size of the first step when optimizing with gradients.  Default is 1.0.
  double m_fdfstepSize;

  //! Accuracy to which to solve line minization to.  Default is 0.1.
  /*!
   * The exact meaning of this parameter depends on the method used.  It's not
   * used for algorithms that don't use derivatives.
   */
  double m_lineTolerance;

private:
  const BaseEnvironment * m_env;

  BoostInputOptionsParser * m_parser;

  // The input options as strings so we can parse the input file later
  std::string m_option_help;

  //! Option name for OptimizerOptions::m_maxIterations.  Default is m_prefix + "optimizer_maxIterations"
  std::string m_option_maxIterations;
  //! Option name for OptimizerOptions::m_tolerance.  Default is m_prefix + "optimizer_tolerance"
  std::string m_option_tolerance;
  //! Option name for OptimizerOptions::m_finiteDifferenceStepSize.  Default is m_prefix + "optimizer_finiteDifferenceStepSize"
  std::string m_option_finiteDifferenceStepSize;
  //! Option name for OptimizerOptions::m_solverType.  Default is m_prefix + "optimizer_solverType"
  std::string m_option_solverType;
  //! Option name for OptimizerOptions::m_fstepSize.  Default is m_prefix + "optimizer_fstepSize"
  std::string m_option_fstepSize;
  //! Option name for OptimizerOptions::m_fdfstepSize.  Default is m_prefix + "optimizer_fdfStepSize"
  std::string m_option_fdfstepSize;
  //! Option name for OptimizerOptions::m_lineTolerance.  Default is m_prefix + "optimizer_lineTolerance"
  std::string m_option_lineTolerance;

  void checkOptions();

  friend std::ostream & operator<<(std::ostream & os,
      const OptimizerOptions & obj);
};

}  // End namespace QUESO

#endif // UQ_OPT_OPTIONS_H
