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

#include <queso/Optimizer.h>

namespace QUESO {

BaseOptimizer::BaseOptimizer()
  : m_maxIterations(UQ_OPT_MAX_ITERATIONS),
    m_tolerance(UQ_OPT_TOLERANCE),
    m_finiteDifferenceStepSize(UQ_OPT_FINITE_DIFFERENCE_STEP_SIZE),
    m_solverType(UQ_OPT_SOLVER_TYPE),
    m_fstepSize(UQ_OPT_FSTEP_SIZE),
    m_fdfstepSize(UQ_OPT_FDFSTEP_SIZE),
    m_lineTolerance(UQ_OPT_LINE_TOLERANCE)
{
  m_optionsObj.reset(new OptimizerOptions);
}

BaseOptimizer::BaseOptimizer(OptimizerOptions options)
{
  m_optionsObj.reset(new OptimizerOptions(options));
}

BaseOptimizer::~BaseOptimizer()
{
}

unsigned int
BaseOptimizer::getMaxIterations() const
{
  return this->m_optionsObj->m_maxIterations;
}

double
BaseOptimizer::getTolerance() const
{
  return this->m_optionsObj->m_tolerance;
}

double
BaseOptimizer::getFiniteDifferenceStepSize() const
{
  return this->m_optionsObj->m_finiteDifferenceStepSize;
}

std::string
BaseOptimizer::getSolverType() const
{
  return this->m_optionsObj->m_solverType;
}

double
BaseOptimizer::getFstepSize() const
{
  return this->m_optionsObj->m_fstepSize;
}

double
BaseOptimizer::getFdfstepSize() const
{
  return this->m_optionsObj->m_fdfstepSize;
}

double
BaseOptimizer::getLineTolerance() const
{
  return this->m_optionsObj->m_lineTolerance;
}

void
BaseOptimizer::setMaxIterations(unsigned int maxIterations)
{
  this->m_optionsObj->m_maxIterations = maxIterations;
}

void
BaseOptimizer::setTolerance(double tolerance)
{
  this->m_optionsObj->m_tolerance = tolerance;
}

void
BaseOptimizer::setFiniteDifferenceStepSize(double h)
{
  this->m_optionsObj->m_finiteDifferenceStepSize = h;
}

void
BaseOptimizer::setSolverType(std::string solverType)
{
  this->m_optionsObj->m_solverType = solverType;
}

void
BaseOptimizer::setFstepSize(double fstepSize)
{
  this->m_optionsObj->m_fstepSize = fstepSize;
}

void
BaseOptimizer::setFdfstepSize(double fdfstepSize)
{
  this->m_optionsObj->m_fdfstepSize = fdfstepSize;
}

void
BaseOptimizer::setLineTolerance(double lineTolerance)
{
  this->m_optionsObj->m_lineTolerance = lineTolerance;
}

}  // End namespace QUESO
