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

    BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction =
      *(static_cast<BaseScalarFunction<GslVector, GslMatrix> * >(context));

    GslVector state(objectiveFunction.domainSet().vectorSpace().zeroVector());

    // DM: Doing this copy sucks, but whatever.  It'll do for now.
    for (unsigned int i = 0; i < state.sizeLocal(); i++) {
      state[i] = gsl_vector_get(x, i);
    }

    // Bail early if GSL tries to evaluate outside of the domain
    if (!objectiveFunction.domainSet().contains(state)) {
      return GSL_NAN;
    }

    // Should cace derivative here so we don't a) segfault in the user's code
    // and b) so we don't recompute stuff
    double result = -objectiveFunction.lnValue(state, NULL, NULL, NULL, NULL);

    return result;
  }

  // This evaluates the derivative of -log posterior
  void c_evaluate_derivative(const gsl_vector * x, void * context,
      gsl_vector * derivative) {
    BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction =
      *(static_cast<BaseScalarFunction<GslVector, GslMatrix> * >(context));

    GslVector state(objectiveFunction.domainSet().vectorSpace().zeroVector());
    GslVector deriv(objectiveFunction.domainSet().vectorSpace().zeroVector());

    // DM: Doing this copy sucks, but whatever.  It'll do for now.
    for (unsigned int i = 0; i < state.sizeLocal(); i++) {
      state[i] = gsl_vector_get(x, i);
    }

    if (!objectiveFunction.domainSet().contains(state)) {
      // Fill derivative with error codes if the point is outside of the
      // domain
      for (unsigned int i = 0; i < deriv.sizeLocal(); i++) {
        gsl_vector_set(derivative, i, GSL_NAN);
      }
    }
    else {
      // We should cache the return value here so we don't recompute stuff
      double fx = -objectiveFunction.lnValue(state, NULL, &deriv, NULL, NULL);

      for (unsigned int i = 0; i < deriv.sizeLocal(); i++) {
        gsl_vector_set(derivative, i, -deriv[i]);  // We need the minus sign
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
    m_objectiveFunction(objectiveFunction)
{
}

GslOptimizer::~GslOptimizer()
{
}

const Vector *
GslOptimizer::minimize() {
  size_t iter = 0;
  int status;
  unsigned int dim = this->m_objectiveFunction.domainSet().vectorSpace().
    zeroVector().sizeLocal();

  const gsl_multimin_fdfminimizer_type * T =
    gsl_multimin_fdfminimizer_conjugate_fr;

  gsl_multimin_fdfminimizer * s =
    gsl_multimin_fdfminimizer_alloc(T, dim);

  gsl_multimin_function_fdf minusLogPosterior;
  minusLogPosterior.n = dim;
  minusLogPosterior.f = &c_evaluate;
  minusLogPosterior.df = &c_evaluate_derivative;
  minusLogPosterior.fdf = &c_evaluate_with_derivative;
  minusLogPosterior.params = (void *)(&(this->m_objectiveFunction));

  // Set initial point (should take this from the user instead);
  gsl_vector * x = gsl_vector_alloc(dim);
  for (unsigned int i = 0; i < dim; i++) {
    gsl_vector_set(x, i, 9.0);  // DM: Initial point should not be hard-coded
  }

  // What are these hard-coded values?
  gsl_multimin_fdfminimizer_set(s, &minusLogPosterior, x, 0.01, 0.1);

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status) {
      if (status == GSL_ENOPROG) {
        std::cerr << "Minimizer could not make progress" << std::endl;
      }
      std::cerr << "Gsl error: " << gsl_strerror(status) << std::endl;
      break;
    }

    // What is this hard-coded value?
    status = gsl_multimin_test_gradient(s->gradient, 1e-3);

    if (status == GSL_SUCCESS) {
      std::cerr << "Minimizer will do another iteration" << std::endl;
    }
  } while (status == GSL_CONTINUE && iter < 100);  // We shouldn't be hard-coding the max number of iterations

  // Get the minimizer
  GslVector * minimizer = new GslVector(this->m_objectiveFunction.domainSet().
      vectorSpace().zeroVector());
  for (unsigned int i = 0; i < dim; i++) {
    (*minimizer)[i] = gsl_vector_get(s->x, i);
  }

  // We're being good human beings and cleaning up the memory we allocated
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  // Return the minimizer.  That's all she wrote.
  return minimizer;
}

}  // End namespace QUESO
