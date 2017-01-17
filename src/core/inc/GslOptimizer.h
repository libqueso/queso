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

#ifndef UQ_GSL_OPTIMIZER_H
#define UQ_GSL_OPTIMIZER_H

#include <queso/Optimizer.h>

namespace QUESO {

/*!
 * \file GslOptimizer.h
 * \brief Class for handling optimization of scalar functions
 *
 * \class GslOptimizer
 * \brief A base class for handling optimisation of scalar functions
 *
 * WRITE DOCS HERE
 */

class Vector;
class GslVector;
class GslMatrix;
class OptimizerMonitor;

template <class V, class M>
class BaseScalarFunction;

class GslOptimizer : public BaseOptimizer {
public:
  //! Constructs an object that will maximize a scalar function
  /*!
   * The function \c objectiveFunction is the function that will be maximized.
   */
  GslOptimizer(
      const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction);


  //! Constructs an object that will maximize a scalar function
  /*!
   * The function \c objectiveFunction is the function that will be maximized.
   * This constructor allows the passing of custom options to optimizer to
   * modify things like tolerance, maximum number of iterations, and finite
   * difference step size.
   */
  GslOptimizer(OptimizerOptions options,
               const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction);

  //! Destructor
  virtual ~GslOptimizer();

  //! Minimize the objective function, starting at \c m_initialPoint
  /*!
   * m_initialPoint is handled in the derived class
   */
  virtual void minimize(OptimizerMonitor* monitor = NULL);

  //! Returns the objective function
  const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction() const;

  enum SolverType { FLETCHER_REEVES_CG,
                    POLAK_RIBIERE_CG,
                    BFGS,
                    BFGS2,
                    STEEPEST_DESCENT,
                    NELDER_MEAD,
                    NELDER_MEAD2,
                    NELDER_MEAD2_RAND };

  //! Set the point at which the optimization starts
  void setInitialPoint(const GslVector & intialPoint);

  //! Return the point that minimizes the objective function
  /*!
   * This state will be filled with GSL_NAN if, for some reason, the
   * optimization failed
   */
  const GslVector & minimizer() const;

  void set_solver_type( SolverType solver );

  void set_solver_type( std::string& solver );

  SolverType string_to_enum( std::string& solver );

  //! Sets step size used in gradient-free solvers
  /*!
   * By default, the step size used will be a vector of 0.1.
   * Use this method to reset the step_size to the desired
   * values.
   */
  void set_step_size( const GslVector& step_size );

  //! Sets step size used in gradient-based solvers
  /*!
   * GSL doesn't document this parameter well, but it seems to be related
   * to the line search, so we default to 1.0 for full step.
   */
  void set_step_size( double step_size );

  //! Set GSL line minimization tolerance
  /*!
   *  Applicable only to gradient-based solvers. Default is 0.1, as
   *  recommended by GSL documentation. See GSL documentation
   *  for more details.
   */
  void set_line_tol( double tol );

  //! Gets the algorithm to use for minimisation
  virtual std::string getSolverType() const;

  //! Gets the step size to use in gradient-free solvers
  virtual double getFstepSize() const;

  //! Gets the step to use in gradient-based solvers
  virtual double getFdfstepSize() const;

  //! Gets the tolerance to use for line minimisation
  virtual double getLineTolerance() const;

  //! Sets the algorithm to use for minimisation
  virtual void setSolverType(std::string solverType);

  //! Sets the step size to use in gradient-free solvers
  virtual void setFstepSize(double fstepSize);

  //! Sets the step to use in gradient-based solvers
  virtual void setFdfstepSize(double fdfstepSize);

  //! Sets the tolerance to use for line minimisation
  virtual void setLineTolerance(double lineTolerance);

private:
  const BaseScalarFunction<GslVector, GslMatrix> & m_objectiveFunction;

  GslVector * m_initialPoint;
  GslVector * m_minimizer;

  SolverType m_solver_type;

  //! For use in gradient-free algorithms
  GslVector m_fstep_size;

  //! For use in gradient-based algorithms
  double m_fdfstep_size;

  //! Line minimization tolerance in gradient-based algorithms
  double m_line_tol;

  //! Helper function
  bool solver_needs_gradient(SolverType solver);

  void minimize_with_gradient( unsigned int dim, OptimizerMonitor* monitor );

  void minimize_no_gradient( unsigned int dim, OptimizerMonitor* monitor );

};

}  // End namespace QUESO

#endif // UQ_GSL_OPTIMIZER_H
