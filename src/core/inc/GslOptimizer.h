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

  //! Destructor
  virtual ~GslOptimizer();

  //! Minimize the objective function, starting at \c m_initialPoint
  /*!
   * m_initialPoint is handled in the derived class
   */
  virtual void minimize();

  //! Returns the objective function
  const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction() const;
  
  enum SolverType { FLETCHER_REEVES_CG,
                    POLAK_RIBIERE_CG,
                    BFGS,
                    BFGS2,
                    STEEPEST_DECENT,
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

private:
  const BaseScalarFunction<GslVector, GslMatrix> & m_objectiveFunction;
  
  GslVector * m_initialPoint;
  GslVector * m_minimizer;

  SolverType m_solver_type;

  //! Helper function
  bool solver_needs_gradient(SolverType solver);

  const GslVector* minimize_with_gradient( unsigned int dim, const GslVector& initial_guess );

};

}  // End namespace QUESO

#endif // UQ_GSL_OPTIMIZER_H
