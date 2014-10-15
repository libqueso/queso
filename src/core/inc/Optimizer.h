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

#ifndef UQ_BASE_OPTIMIZER_H
#define UQ_BASE_OPTIMIZER_H

namespace QUESO {

/*!
 * \file Optimizer.h
 * \brief Class for handling optimization of scalar functions
 *
 * \class BaseOptimizer
 * \brief A base class for handling optimisation of scalar functions
 *
 * This class defines the interface every optimizer must adhere to.  All
 * optimizers should subclass this class.
 */

class Vector;

class BaseOptimizer {
public:
  //! Default constructor.
  BaseOptimizer();

  //! Destructor
  virtual ~BaseOptimizer();

  virtual const Vector * minimize(const Vector & initialPoint) = 0;

  //! Returns the maximum number of iterations the optimizer will do
  unsigned int getMaxIterations() const;

  //! Returns the tolerance used to test for an extremum in the optimizer
  double getTolerance() const;

  //! Returns the step size used in the finite difference formula
  double getFiniteDifferenceStepSize() const;

  //! Sets the maximum number of iterations to be used by the optimizer
  void setMaxIterations(unsigned int maxIterations);

  //! Sets the tolerance the optimizer will use to test for an extremum
  void setTolerance(double tolerance);

  //! Sets the step to use in the finite difference derivative
  void setFiniteDifferenceStepSize(double h);

protected:
  Vector * minimizer;

private:
  unsigned int m_maxIterations;
  double m_tolerance;
  double m_finiteDifferenceStepSize;
};

}  // End namespace QUESO

#endif // UQ_BASE_OPTIMIZER_H
