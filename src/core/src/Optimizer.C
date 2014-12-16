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

#include <queso/Optimizer.h>

namespace QUESO {

BaseOptimizer::BaseOptimizer()
  : m_maxIterations(100),
    m_tolerance(1e-3),
    m_finiteDifferenceStepSize(1e-4)
{
}

BaseOptimizer::~BaseOptimizer()
{
}

unsigned int
BaseOptimizer::getMaxIterations() const
{
  return this->m_maxIterations;
}

double
BaseOptimizer::getTolerance() const
{
  return this->m_tolerance;
}

double
BaseOptimizer::getFiniteDifferenceStepSize() const
{
  return this->m_finiteDifferenceStepSize;
}

void
BaseOptimizer::setMaxIterations(unsigned int maxIterations)
{
  this->m_maxIterations = maxIterations;
}

void
BaseOptimizer::setTolerance(double tolerance)
{
  this->m_tolerance = tolerance;
}

void
BaseOptimizer::setFiniteDifferenceStepSize(double h)
{
  this->m_finiteDifferenceStepSize = h;
}

}  // End namespace QUESO
