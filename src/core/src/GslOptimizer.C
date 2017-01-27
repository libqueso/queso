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

#include <iostream>

#include <queso/Defines.h>
#include <queso/GslVector.h>
#include <queso/VectorSpace.h>
#include <queso/ScalarFunction.h>
#include <queso/GslOptimizer.h>
#include <queso/OptimizerMonitor.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

namespace QUESO {

// We need to extern "C" because gsl needs a pointer to a C function to
// minimize
extern "C" {
  // This evaluate -log posterior
  double c_evaluate(const gsl_vector * x, void * context) {

    GslOptimizer * optimizer = static_cast<GslOptimizer * >(context);

    GslVector state(
        optimizer->objectiveFunction().domainSet().vectorSpace().zeroVector());

    // DM: Doing this copy sucks, but whatever.  It'll do for now.
    for (unsigned int i = 0; i < state.sizeLocal(); i++) {
      state[i] = gsl_vector_get(x, i);
    }

    // Bail early if GSL tries to evaluate outside of the domain
    if (!optimizer->objectiveFunction().domainSet().contains(state)) {
      return GSL_NAN;
    }

    // Should cache derivative here so we don't recompute stuff
    double result = -optimizer->objectiveFunction().lnValue(state);

    return result;
  }

  // This evaluates the derivative of -log posterior
  void c_evaluate_derivative(const gsl_vector * x, void * context,
      gsl_vector * derivative) {
    GslOptimizer * optimizer = static_cast<GslOptimizer * >(context);

    GslVector state(
        optimizer->objectiveFunction().domainSet().vectorSpace().zeroVector());
    GslVector deriv(
        optimizer->objectiveFunction().domainSet().vectorSpace().zeroVector());

    // DM: Doing this copy sucks, but whatever.  It'll do for now.
    for (unsigned int i = 0; i < state.sizeLocal(); i++) {
      state[i] = gsl_vector_get(x, i);

      // We fill with GSL_NAN and use it as a flag to check later that the user
      // actually fills the derivative vector with stuff
      deriv[i] = GSL_NAN;
    }

    if (!optimizer->objectiveFunction().domainSet().contains(state)) {
      // Fill derivative with error codes if the point is outside of the
      // domain
      for (unsigned int i = 0; i < deriv.sizeLocal(); i++) {
        gsl_vector_set(derivative, i, GSL_NAN);
      }
    }
    else {
      // We should cache the return value here so we don't recompute stuff
      double fx = -optimizer->objectiveFunction().lnValue(state, deriv);

      // Decide whether or not we need to do a finite difference based on
      // whether the user actually filled deriv with values that are not
      // GSL_NAN
      //
      // We're currently doing this check every time this function gets called.
      // We could probably pull this logic out of here and put it somewhere
      // where it only happens once
      bool userComputedDerivative = true;
      for (unsigned int i = 0; i < deriv.sizeLocal(); i++) {
        // If the user missed out a derivative in any direction, fall back to
        // a finite difference
        if (gsl_isnan(deriv[i])) {
          userComputedDerivative = false;
          break;
        }
      }

      if (userComputedDerivative) {
        for (unsigned int i = 0; i < deriv.sizeLocal(); i++) {
          gsl_vector_set(derivative, i, -deriv[i]);  // We need the minus sign
        }
      }
      else {
        // Finite difference step-size
        double h = optimizer->getFiniteDifferenceStepSize();

        // User did not provide a derivative, so do a finite difference
        for (unsigned int i = 0; i < deriv.sizeLocal(); i++) {
          double tempState = state[i];
          state[i] += h;

          // User didn't provide a derivative, so we don't bother passing in
          // the derivative vector again
          double fxph = -optimizer->objectiveFunction().lnValue(state);

          // Reset the state back to what it was before
          state[i] = tempState;

          // Make sure we didn't do anything dumb and tell gsl if we did
          if (!gsl_isnan(fx) && !gsl_isnan(fxph)) {
            gsl_vector_set(derivative, i, (fxph - fx) / h);
          }
          else {
            gsl_vector_set(derivative, i, GSL_NAN);
          }
        }
      }
    }
  }

  // This evaluates -log posterior and the derivative of -log posterior
  void c_evaluate_with_derivative(const gsl_vector * x, void * context,
      double * f, gsl_vector * derivative) {
    // We don't need to call both of these
    *f = c_evaluate(x, context);
    c_evaluate_derivative(x, context, derivative);
  }
}  // End extern "C"

GslOptimizer::GslOptimizer(
    const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction)
  : BaseOptimizer(),
    m_objectiveFunction(objectiveFunction),
    m_initialPoint(new GslVector(objectiveFunction.domainSet().
          vectorSpace().zeroVector())),
    m_minimizer(new GslVector(this->m_objectiveFunction.domainSet().
        vectorSpace().zeroVector())),
    m_solver_type(BFGS2),
    m_fstep_size(this->m_objectiveFunction.domainSet().vectorSpace().zeroVector()),
    m_fdfstep_size(getFdfstepSize()),
    m_line_tol(getLineTolerance())
{
  // We initialize the minimizer to GSL_NAN just in case the optimization fails
  m_minimizer->cwSet(GSL_NAN);

  // Set to documented default value.
  m_fstep_size.cwSet(getFstepSize());

  // Set solver type to the one set in the options object
  setSolverType(getSolverType());
}

GslOptimizer::GslOptimizer(
    OptimizerOptions options,
    const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction)
  : BaseOptimizer(options),
    m_objectiveFunction(objectiveFunction),
    m_initialPoint(new GslVector(objectiveFunction.domainSet().
          vectorSpace().zeroVector())),
    m_minimizer(new GslVector(this->m_objectiveFunction.domainSet().
        vectorSpace().zeroVector())),
    m_solver_type(BFGS2),
    m_fstep_size(this->m_objectiveFunction.domainSet().vectorSpace().zeroVector()),
    m_fdfstep_size(getFdfstepSize()),
    m_line_tol(getLineTolerance())
{
  // We initialize the minimizer to GSL_NAN just in case the optimization fails
  m_minimizer->cwSet(GSL_NAN);

  // Set to documented default value.
  m_fstep_size.cwSet(getFstepSize());

  // Set solver type to the one set in the options object
  setSolverType(getSolverType());
}

GslOptimizer::~GslOptimizer()
{
  delete this->m_initialPoint;
}

void GslOptimizer::minimize(OptimizerMonitor* monitor) {

  // First check that initial guess is reasonable
  if (!this->m_objectiveFunction.domainSet().contains(*(this->m_initialPoint)))
    {
      if( m_objectiveFunction.domainSet().env().fullRank() == 0 )
        {
          std::cerr << "Minimization was given initial point outside of domain"
                    << std::endl;
        }
      queso_error();
    }

  unsigned int dim = this->m_objectiveFunction.domainSet().vectorSpace().
    zeroVector().sizeLocal();

  // We use m_solver_type here because we need the enum
  if( this->solver_needs_gradient(m_solver_type) )
    {
      this->minimize_with_gradient( dim, monitor );
    }
  else
    {
      this->minimize_no_gradient( dim, monitor );
    }

  return;
}

const BaseScalarFunction<GslVector, GslMatrix> &
GslOptimizer::objectiveFunction() const
{
  return this->m_objectiveFunction;
}

void
GslOptimizer::setInitialPoint(const GslVector & initialPoint)
{
  for (unsigned int i = 0; i < initialPoint.sizeLocal(); i++) {
    (*(this->m_initialPoint))[i] = initialPoint[i];
  }
}

const GslVector &
GslOptimizer::minimizer() const
{
  return *(this->m_minimizer);
}

void GslOptimizer::set_solver_type( SolverType solver )
{
  queso_deprecated();
  m_solver_type = solver;
}

bool GslOptimizer::solver_needs_gradient( SolverType solver )
{
  bool gradient_needed = false;

  switch(solver)
    {
    case(FLETCHER_REEVES_CG):
    case(POLAK_RIBIERE_CG):
    case(BFGS):
    case(BFGS2):
    case(STEEPEST_DESCENT):
      {
        gradient_needed = true;
        break;
      }
    case(NELDER_MEAD):
    case(NELDER_MEAD2):
    case(NELDER_MEAD2_RAND):
      {
        break;
      }
    default:
      {
        // Wat?!
        queso_error();
      }
    } // switch(solver)

  return gradient_needed;
}

void GslOptimizer::minimize_with_gradient( unsigned int dim, OptimizerMonitor* monitor )
{
  // Set initial point
  gsl_vector * x = gsl_vector_alloc(dim);
  for (unsigned int i = 0; i < dim; i++) {
    gsl_vector_set(x, i, (*m_initialPoint)[i]);
  }

  // Tell GSL which solver we're using
  const gsl_multimin_fdfminimizer_type* type = NULL;

  // We use m_solver_type here because we need the enum
  switch(m_solver_type)
    {
    case(FLETCHER_REEVES_CG):
      type = gsl_multimin_fdfminimizer_conjugate_fr;
      break;
    case(POLAK_RIBIERE_CG):
      type = gsl_multimin_fdfminimizer_conjugate_pr;
      break;
    case(BFGS):
      type = gsl_multimin_fdfminimizer_vector_bfgs;
      break;
    case(BFGS2):
      type = gsl_multimin_fdfminimizer_vector_bfgs2;
      break;
    case(STEEPEST_DESCENT):
      type = gsl_multimin_fdfminimizer_steepest_descent;
      break;
    case(NELDER_MEAD):
    case(NELDER_MEAD2):
    case(NELDER_MEAD2_RAND):
    default:
      // Wat?!
      queso_error();
    }

  // Init solver
  gsl_multimin_fdfminimizer * solver =
    gsl_multimin_fdfminimizer_alloc(type, dim);

  // Point GSL to the right functions
  gsl_multimin_function_fdf minusLogPosterior;
  minusLogPosterior.n = dim;
  minusLogPosterior.f = &c_evaluate;
  minusLogPosterior.df = &c_evaluate_derivative;
  minusLogPosterior.fdf = &c_evaluate_with_derivative;
  minusLogPosterior.params = (void *)(this);

  gsl_multimin_fdfminimizer_set(solver, &minusLogPosterior, x, getFdfstepSize(), getLineTolerance());

  int status;
  size_t iter = 0;

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(solver);

    if (status) {
      if( m_objectiveFunction.domainSet().env().fullRank() == 0 )
        {
          std::cerr << "Error while GSL does optimisation. "
                    << "See below for GSL error type." << std::endl;
          std::cerr << "Gsl error: " << gsl_strerror(status) << std::endl;
        }
      break;
    }

    status = gsl_multimin_test_gradient(solver->gradient, this->getTolerance());

    if(monitor)
      {
        gsl_vector* x = gsl_multimin_fdfminimizer_x(solver);
        std::vector<double> x_min(dim);
        for( unsigned int i = 0; i < dim; i++)
          x_min[i] = gsl_vector_get(x,i);

        double f = gsl_multimin_fdfminimizer_minimum(solver);

        gsl_vector* grad = gsl_multimin_fdfminimizer_gradient(solver);
        double grad_norm = gsl_blas_dnrm2(grad);

        monitor->append( x_min, f, grad_norm );
      }

  } while ((status == GSL_CONTINUE) && (iter < this->getMaxIterations()));

  for (unsigned int i = 0; i < dim; i++) {
    (*m_minimizer)[i] = gsl_vector_get(solver->x, i);
  }

  // We're being good human beings and cleaning up the memory we allocated
  gsl_multimin_fdfminimizer_free(solver);
  gsl_vector_free(x);

  return;
}

void GslOptimizer::minimize_no_gradient( unsigned int dim, OptimizerMonitor* monitor )
{
  // Set initial point
  gsl_vector* x = gsl_vector_alloc(dim);
  for (unsigned int i = 0; i < dim; i++) {
    gsl_vector_set(x, i, (*m_initialPoint)[i]);
  }

  // Tell GSL which solver we're using
  const gsl_multimin_fminimizer_type* type = NULL;

  switch(m_solver_type)
    {
    case(NELDER_MEAD):
      type = gsl_multimin_fminimizer_nmsimplex;
      break;
    case(NELDER_MEAD2):
      type = gsl_multimin_fminimizer_nmsimplex2;
      break;
    case(NELDER_MEAD2_RAND):
      type = gsl_multimin_fminimizer_nmsimplex2rand;
      break;
    case(FLETCHER_REEVES_CG):
    case(POLAK_RIBIERE_CG):
    case(BFGS):
    case(BFGS2):
    case(STEEPEST_DESCENT):
    default:
      // Wat?!
      queso_error();
    }

  // Init solver
  gsl_multimin_fminimizer* solver =
    gsl_multimin_fminimizer_alloc(type, dim);

  // Point GSL at the right functions
  gsl_multimin_function minusLogPosterior;
  minusLogPosterior.n = dim;
  minusLogPosterior.f = &c_evaluate;
  minusLogPosterior.params = (void *)(this);

  // Needed for these gradient free algorithms.
  gsl_vector* step_size = gsl_vector_alloc(dim);

  for(unsigned int i = 0; i < dim; i++) {
    gsl_vector_set(step_size, i, m_fstep_size[i]);
  }

  gsl_multimin_fminimizer_set(solver, &minusLogPosterior, x, step_size);

  int status;
  size_t iter = 0;
  double size = 0.0;

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(solver);

      if (status) {
        if( m_objectiveFunction.domainSet().env().fullRank() == 0 )
          {
            std::cerr << "Error while GSL does optimisation. "
                      << "See below for GSL error type." << std::endl;
            std::cerr << "Gsl error: " << gsl_strerror(status) << std::endl;
          }
        break;
      }

      size = gsl_multimin_fminimizer_size(solver);

      status = gsl_multimin_test_size (size, this->getTolerance());

      if(monitor)
      {
        gsl_vector* x = gsl_multimin_fminimizer_x(solver);
        std::vector<double> x_min(dim);
        for( unsigned int i = 0; i < dim; i++)
          x_min[i] = gsl_vector_get(x,i);

        double f = gsl_multimin_fminimizer_minimum(solver);

        monitor->append( x_min, f, size );
      }

    }

  while ((status == GSL_CONTINUE) && (iter < this->getMaxIterations()));

  for (unsigned int i = 0; i < dim; i++) {
    (*m_minimizer)[i] = gsl_vector_get(solver->x, i);
  }

  // We're being good human beings and cleaning up the memory we allocated
  gsl_vector_free(step_size);
  gsl_multimin_fminimizer_free(solver);
  gsl_vector_free(x);

  return;
}

void GslOptimizer::set_step_size( const GslVector& step_size )
{
  queso_deprecated();
  m_fstep_size = step_size;
}

void GslOptimizer::set_step_size( double step_size )
{
  queso_deprecated();
  m_fdfstep_size = step_size;
}

GslOptimizer::SolverType GslOptimizer::string_to_enum( std::string& solver )
{
  SolverType solver_type;

  if( solver == std::string("fletcher_reeves_cg") )
    solver_type = FLETCHER_REEVES_CG;
  else if( solver == std::string("polak_ribiere_cg") )
    solver_type = POLAK_RIBIERE_CG;
  else if( solver == std::string("bfgs") )
    solver_type = BFGS;
  else if( solver == std::string("bfgs2") )
    solver_type = BFGS2;
  else if( solver == std::string("steepest_decent") ) {
    queso_deprecated();
    solver_type = STEEPEST_DESCENT;
  }
  else if( solver == std::string("steepest_descent") )
    solver_type = STEEPEST_DESCENT;
  else if( solver == std::string("nelder_mead") )
    solver_type = NELDER_MEAD;
  else if( solver == std::string("nelder_mead2") )
    solver_type = NELDER_MEAD2;
  else if( solver == std::string("nelder_mead2_rand") )
    solver_type = NELDER_MEAD2_RAND;
  else
    {
      if( m_objectiveFunction.domainSet().env().fullRank() == 0 )
        {
          std::cerr << "Error: Invalid GslOptimizer solver name: " << solver << std::endl
                    << "       Valids choices are: fletcher_reeves_cg" << std::endl
                    << "                           polak_ribiere_cg" << std::endl
                    << "                           bfgs" << std::endl
                    << "                           bfgs2" << std::endl
                    << "                           steepest_descent" << std::endl
                    << "                           nelder_mead" << std::endl
                    << "                           nelder_mead2" << std::endl
                    << "                           nelder_mead2_rand" << std::endl;
        }
      queso_error();
    }

  return solver_type;
}

void GslOptimizer::set_solver_type( std::string& solver )
{
  queso_deprecated()
  this->set_solver_type( this->string_to_enum(solver) );
}

void
GslOptimizer::setSolverType(std::string solverType)
{
  this->m_optionsObj->m_solverType = solverType;
  this->set_solver_type(solverType);
}

void
GslOptimizer::setFstepSize(double fstepSize)
{
  this->m_optionsObj->m_fstepSize = fstepSize;

  GslVector fstepSizeVector(
      objectiveFunction().domainSet().vectorSpace().zeroVector());
  fstepSizeVector.cwSet(fstepSize);

  this->set_step_size(fstepSizeVector);
}

void
GslOptimizer::setFdfstepSize(double fdfstepSize)
{
  this->m_optionsObj->m_fdfstepSize = fdfstepSize;
  this->set_step_size(fdfstepSize);
}

void
GslOptimizer::setLineTolerance(double lineTolerance)
{
  this->m_optionsObj->m_lineTolerance = lineTolerance;
}

std::string
GslOptimizer::getSolverType() const
{
  return this->m_optionsObj->m_solverType;
}

double
GslOptimizer::getFstepSize() const
{
  return this->m_optionsObj->m_fstepSize;
}

double
GslOptimizer::getFdfstepSize() const
{
  return this->m_optionsObj->m_fdfstepSize;
}

double
GslOptimizer::getLineTolerance() const
{
  return this->m_optionsObj->m_lineTolerance;
}

}  // End namespace QUESO
