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

#include <iostream>

#include <gsl/gsl_multimin.h>
#include <queso/GslVector.h>
#include <queso/VectorSpace.h>
#include <queso/ScalarFunction.h>
#include <queso/GslOptimizer.h>

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

    // Should cache derivative here so we don't a) segfault in the user's code
    // and b) so we don't recompute stuff
    double result = -optimizer->objectiveFunction().lnValue(state, NULL, NULL,
        NULL, NULL);

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
      double fx = -optimizer->objectiveFunction().lnValue(state, NULL, &deriv,
          NULL, NULL);

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
          double fxph = -optimizer->objectiveFunction().lnValue(state, NULL,
              NULL, NULL, NULL);

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
    m_solver_type(BFGS2)
{
  // We initialize the minimizer to GSL_NAN just in case the optimization fails
  for (unsigned int i = 0; i < this->m_minimizer->sizeLocal(); i++) {
    (*(this->m_minimizer))[i] = GSL_NAN;
  }
}

GslOptimizer::~GslOptimizer()
{
  delete this->m_initialPoint;
}

const Vector *
GslOptimizer::minimize(const Vector & initialPoint) {

  unsigned int dim = this->m_objectiveFunction.domainSet().vectorSpace().
    zeroVector().sizeLocal();

  // DM: The dynamic cast is needed because the abstract base class does not
  // have a pure virtual operator[] method.  We can implement one and remove
  // this dynamic cast in the future.
  const GslVector & initial_guess = dynamic_cast<const GslVector &>(initialPoint);

  const GslVector* minimizer = this->minimize_with_gradient( dim, initial_guess );

  return minimizer;
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
    m_solver_type = solver;
  }

  bool GslOptimizer::solver_needs_gradient( SolverType solver )
  {
    bool gradient_needed = false;

    switch(solver)
      {
      case(FLETCHER_REEVES):
      case(CONJUGATE_GRADIENT):
      case(BFGS):
      case(BFGS2):
      case(STEEPEST_DECENT):
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

  const GslVector* GslOptimizer::minimize_with_gradient( unsigned int dim, const GslVector& initial_guess )
  {
    // Set initial point
    gsl_vector * x = gsl_vector_alloc(dim);
    for (unsigned int i = 0; i < dim; i++) {
      gsl_vector_set(x, i, initial_guess[i]);
    }

    const gsl_multimin_fdfminimizer_type * T =
      gsl_multimin_fdfminimizer_vector_bfgs2;

    gsl_multimin_fdfminimizer * s =
      gsl_multimin_fdfminimizer_alloc(T, dim);

    gsl_multimin_function_fdf minusLogPosterior;
    minusLogPosterior.n = dim;
    minusLogPosterior.f = &c_evaluate;
    minusLogPosterior.df = &c_evaluate_derivative;
    minusLogPosterior.fdf = &c_evaluate_with_derivative;
    minusLogPosterior.params = (void *)(&(this->m_objectiveFunction));


    /*!
     * \todo Allow the user to tweak these hard-coded values
     */
    gsl_multimin_fdfminimizer_set(s, &minusLogPosterior, x, 0.01, 0.1);

    int status;
    size_t iter = 0;

    do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if (status) {
        std::cerr << "Error while GSL does optimisation. "
                  << "See below for GSL error type." << std::endl;
        std::cerr << "Gsl error: " << gsl_strerror(status) << std::endl;
        break;
      }

      // TODO: Allow the user to tweak this hard-coded value
      status = gsl_multimin_test_gradient(s->gradient, this->getTolerance());

      /*!
       * \todo We shouldn't be hard-coding the max number of iterations
       */
    } while ((status == GSL_CONTINUE) && (iter < this->getMaxIterations())); 

    // Get the minimizer
    GslVector * minimizer = new GslVector(this->m_objectiveFunction.domainSet().
                                          vectorSpace().zeroVector());
    for (unsigned int i = 0; i < dim; i++) {
      (*minimizer)[i] = gsl_vector_get(s->x, i);
    }

    // We're being good human beings and cleaning up the memory we allocated
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

    return minimizer;
  }

}  // End namespace QUESO
