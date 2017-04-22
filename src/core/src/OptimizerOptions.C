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
#include <queso/Environment.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#else
#include <queso/getpot.h>
#endif

#include <queso/OptimizerOptions.h>

namespace QUESO {

OptimizerOptions::OptimizerOptions()
  : m_prefix("ip_"),
    m_help(UQ_OPT_HELP),
    m_maxIterations(UQ_OPT_MAX_ITERATIONS),
    m_tolerance(UQ_OPT_TOLERANCE),
    m_finiteDifferenceStepSize(UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE),
    m_solverType(UQ_OPT_SOLVER_TYPE),
    m_fstepSize(UQ_OPT_FSTEP_SIZE),
    m_fdfstepSize(UQ_OPT_FDFSTEP_SIZE),
    m_lineTolerance(UQ_OPT_LINE_TOLERANCE),
    m_env(NULL),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(NULL),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help(m_prefix + "help"),
    m_option_maxIterations(m_prefix + "maxIterations"),
    m_option_tolerance(m_prefix + "tolerance"),
    m_option_finiteDifferenceStepSize(m_prefix + "finiteDifferenceStepSize"),
    m_option_solverType(m_prefix + "solverType"),
    m_option_fstepSize(m_prefix + "fstepSize"),
    m_option_fdfstepSize(m_prefix + "fdfStepSize"),
    m_option_lineTolerance(m_prefix + "lineTolerance")
{
}

OptimizerOptions::OptimizerOptions(const BaseEnvironment * env, const char *
    prefix)
  : m_prefix((std::string)(prefix) + "optimizer_"),
    m_help(UQ_OPT_HELP),
    m_maxIterations(UQ_OPT_MAX_ITERATIONS),
    m_tolerance(UQ_OPT_TOLERANCE),
    m_finiteDifferenceStepSize(UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE),
    m_solverType(UQ_OPT_SOLVER_TYPE),
    m_fstepSize(UQ_OPT_FSTEP_SIZE),
    m_fdfstepSize(UQ_OPT_FDFSTEP_SIZE),
    m_lineTolerance(UQ_OPT_LINE_TOLERANCE),
    m_env(env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help(m_prefix + "help"),
    m_option_maxIterations(m_prefix + "maxIterations"),
    m_option_tolerance(m_prefix + "tolerance"),
    m_option_finiteDifferenceStepSize(m_prefix + "finiteDifferenceStepSize"),
    m_option_solverType(m_prefix + "solverType"),
    m_option_fstepSize(m_prefix + "fstepSize"),
    m_option_fdfstepSize(m_prefix + "fdfStepSize"),
    m_option_lineTolerance(m_prefix + "lineTolerance")
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string>(m_option_help, UQ_OPT_HELP,
      "produce help message for statistical inverse problem");

  m_parser->registerOption<unsigned int>(m_option_maxIterations,
      UQ_OPT_MAX_ITERATIONS,
      "max number of optimizer iterations to do");
  m_parser->registerOption<double>(m_option_tolerance, UQ_OPT_TOLERANCE,
      "optimize until gradient is less than tolerance");
  m_parser->registerOption<double>(m_option_finiteDifferenceStepSize,
      UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE,
      "if no deriv is given, do finite difference with this step size");
  m_parser->registerOption<std::string>(m_option_solverType, UQ_OPT_SOLVER_TYPE,
      "which optimisation algorithm to use");
  m_parser->registerOption<double>(m_option_fstepSize, UQ_OPT_FSTEP_SIZE,
      "sets the step size used in gradient-free solvers");
  m_parser->registerOption<double>(m_option_fdfstepSize, UQ_OPT_FDFSTEP_SIZE,
      "sets the step size used in gradient-based solvers");
  m_parser->registerOption<double>(m_option_lineTolerance, UQ_OPT_LINE_TOLERANCE,
      "sets the line minimisation tolerance");

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help, m_help);
  m_parser->getOption<unsigned int>(m_option_maxIterations, m_maxIterations);
  m_parser->getOption<double>(m_option_tolerance, m_tolerance);
  m_parser->getOption<double>(m_option_finiteDifferenceStepSize,
      m_finiteDifferenceStepSize);
  m_parser->getOption<std::string>(m_option_solverType, m_solverType);
  m_parser->getOption<double>(m_option_fstepSize, m_fstepSize);
  m_parser->getOption<double>(m_option_fdfstepSize, m_fdfstepSize);
  m_parser->getOption<double>(m_option_lineTolerance, m_lineTolerance);
#else
  m_help = m_env->input()(m_option_help, UQ_OPT_HELP);
  m_maxIterations = m_env->input()(m_option_maxIterations, UQ_OPT_MAX_ITERATIONS);
  m_tolerance = m_env->input()(m_option_tolerance, UQ_OPT_TOLERANCE);
  m_finiteDifferenceStepSize = m_env->input()(m_option_finiteDifferenceStepSize, UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE);
  m_solverType = m_env->input()(m_option_solverType, UQ_OPT_SOLVER_TYPE);
  m_fstepSize = m_env->input()(m_option_fstepSize, UQ_OPT_FSTEP_SIZE);
  m_fdfstepSize = m_env->input()(m_option_fdfstepSize, UQ_OPT_FDFSTEP_SIZE);
  m_lineTolerance = m_env->input()(m_option_lineTolerance, UQ_OPT_LINE_TOLERANCE);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

OptimizerOptions::OptimizerOptions(const OptimizerOptions & rhs)
  :
    m_prefix(rhs.m_prefix),
    m_help(rhs.m_help),
    m_maxIterations(rhs.m_maxIterations),
    m_tolerance(rhs.m_tolerance),
    m_finiteDifferenceStepSize(rhs.m_finiteDifferenceStepSize),
    m_solverType(rhs.m_solverType),
    m_fstepSize(rhs.m_fstepSize),
    m_fdfstepSize(rhs.m_fdfstepSize),
    m_lineTolerance(rhs.m_lineTolerance),
    m_parser(rhs.m_parser),  // We'll never touch the input file in a copied object
    m_option_help(rhs.m_option_help),
    m_option_maxIterations(rhs.m_option_maxIterations),
    m_option_tolerance(rhs.m_option_tolerance),
    m_option_finiteDifferenceStepSize(rhs.m_option_finiteDifferenceStepSize),
    m_option_solverType(rhs.m_option_solverType),
    m_option_fstepSize(rhs.m_option_fstepSize),
    m_option_fdfstepSize(rhs.m_option_fdfstepSize),
    m_option_lineTolerance(rhs.m_option_lineTolerance)
{
}

OptimizerOptions::~OptimizerOptions()
{
}

void
OptimizerOptions::checkOptions()
{
  queso_require_greater_msg(m_tolerance, 0.0, "optimizer tolerance must be > 0");
  queso_require_greater_msg(m_finiteDifferenceStepSize, 0.0, "finite difference step must be > 0");
  queso_require_greater_msg(m_maxIterations, 0, "max iterations must be > 0");
  queso_require_greater_msg(m_fstepSize, 0.0, "fstepSize must be > 0");
  queso_require_greater_msg(m_fdfstepSize, 0.0, "fdfstepSize must be > 0");
  queso_require_greater_msg(m_lineTolerance, 0.0, "line tolerance must be > 0");
}

std::ostream &
operator<<(std::ostream& os, const OptimizerOptions & obj)
{
  os << "\n" << obj.m_option_maxIterations << " = " << obj.m_maxIterations
     << "\n" << obj.m_option_tolerance << " = " << obj.m_tolerance;
  os << "\n" << obj.m_option_finiteDifferenceStepSize << " = "
             << obj.m_finiteDifferenceStepSize;
  os << "\n" << obj.m_option_solverType << " = " << obj.m_solverType;
  os << "\n" << obj.m_option_fstepSize << " = " << obj.m_fstepSize;
  os << "\n" << obj.m_option_fdfstepSize << " = " << obj.m_fdfstepSize;
  os << "\n" << obj.m_option_lineTolerance << " = " << obj.m_lineTolerance;
  os << std::endl;
  return os;
}

}  // End namespace QUESO
