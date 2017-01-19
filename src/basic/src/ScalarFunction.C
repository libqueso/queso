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

#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#else
#include <queso/getpot.h>
#endif

#include <cstdlib>

#define QUESO_BASESCALARFN_FD_STEPSIZE_ODV "1e-6"

namespace QUESO {

// Default constructor
template<class V, class M>
BaseScalarFunction<V, M>::BaseScalarFunction(const char * prefix,
    const VectorSet<V, M> & domainSet)
  : m_env(domainSet.env()),
    m_prefix((std::string)(prefix) + "func_"),
    m_domainSet(domainSet),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(m_env.optionsInputFileName())),
#endif
    m_fdStepSize(1, std::atof(QUESO_BASESCALARFN_FD_STEPSIZE_ODV))
{
  unsigned int dim = m_domainSet.vectorSpace().dimLocal();

  // Snarf fd step size from input file.
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string>(m_prefix + "fdStepSize",
                                        QUESO_BASESCALARFN_FD_STEPSIZE_ODV,
                                        "step size for finite difference");
  m_parser->scanInputFile();
  m_parser->getOption<std::vector<double> >(m_prefix + "fdStepSize",
                                            m_fdStepSize);

  // Check size of finite difference vector the user provided.
  queso_require_msg((m_fdStepSize.size() == 1) || (m_fdStepSize.size() == dim),
                    "Finite difference vector is not the correct size");

  // If the user provided a scalar for a multi-dimensional function...
  if (dim > 1 && m_fdStepSize.size() == 1) {
    // ...get the only element
    double stepSize = m_fdStepSize[0];

    // and use it to fill a vector of length dim
    m_fdStepSize.resize(dim, stepSize);
  }
#else
  unsigned int size = m_env.input().vector_variable_size(m_prefix + "fdStepSize");

  if (size == 0) {
    m_fdStepSize.resize(dim, std::atof(QUESO_BASESCALARFN_FD_STEPSIZE_ODV));
  }
  else if (size == 1) {
    double value = m_env.input()(m_prefix + "fdStepSize",
                                 std::atof(QUESO_BASESCALARFN_FD_STEPSIZE_ODV),
                                 0);

    m_fdStepSize.resize(dim, value);
  }
  else if (size == dim) {
    for (unsigned int i = 0; i < size; i++) {
      m_fdStepSize[i] = m_env.input()(m_prefix + "fdStepSize",
                                      std::atof(QUESO_BASESCALARFN_FD_STEPSIZE_ODV),
                                      i);
    }
  }
  else {
    // Either the user provides nothing, a scalar, or the whole vector.
    // Any other possiblities are not allowed so we error in this case.
    queso_error_msg("Finite difference vector must be a scalar or a vector of length parameter dimension");
  }
#endif

  // Check all the elements of the finite difference vector are positive
  for (unsigned int i = 0; i < dim; i++) {
    queso_require_greater_msg(m_fdStepSize[i],
                              0.0,
                              "Finite difference step sizes must be positive");
  }
}

// Destructor
template<class V, class M>
BaseScalarFunction<V, M>::~BaseScalarFunction()
{
}

// Math methods
template<class V, class M>
const VectorSet<V, M> & BaseScalarFunction<V, M>::domainSet() const
{
  return m_domainSet;
}

template <class V, class M>
double
BaseScalarFunction<V, M>::lnValue(const V & domainVector,
                                  const V * domainDirection,
                                  V * gradVector,
                                  M * hessianMatrix,
                                  V * hessianEffect) const
{
  std::string msg;

  msg += "Implementation of all lnValue methods is missing.  Please implement";
  msg += " at least lnValue(const V &).";

  queso_error_msg(msg);
}

template <class V, class M>
double
BaseScalarFunction<V, M>::lnValue(const V & domainVector) const
{
  return this->lnValue(domainVector, NULL, NULL, NULL, NULL);
}

template <class V, class M>
double
BaseScalarFunction<V, M>::lnValue(const V & domainVector, V & gradVector) const
{
  double value = this->lnValue(domainVector);

  // Create perturbed version of domainVector to use in finite difference
  V perturbedVector(domainVector);

  // Fill up gradVector with a finite difference approximation
  for (unsigned int i = 0; i < domainVector.sizeLocal(); ++i) {
    // Store the old value of the perturbed element so we can undo it later
    double tmp = perturbedVector[i];

    perturbedVector[i] += m_fdStepSize[i];
    gradVector[i] = (this->lnValue(perturbedVector) - value) / m_fdStepSize[i];

    // Restore the old value of the perturbedVector element
    perturbedVector[i] = tmp;
  }

  return value;
}

template <class V, class M>
double
BaseScalarFunction<V, M>::lnValue(const V & domainVector,
                                  V & gradVector,
                                  const V & domainDirection,
                                  V & hessianEffect) const
{
  std::string msg;

  msg += "QUESO asked for Hessian information from an lnValue method, but the";
  msg += " implementation of is missing.  Please implement";
  msg += " lnValue(const V &, V &, const V &, V &).";

  queso_error_msg(msg);
}

template <class V, class M>
void
BaseScalarFunction<V, M>::setFiniteDifferenceStepSize(double fdStepSize)
{
  queso_require_greater_msg(fdStepSize, 0.0,
                            "Must provide a finite difference step > 0");

  for (unsigned int i = 0; i < this->m_fdStepSize.size(); i++) {
    this->m_fdStepSize[i] = fdStepSize;
  }
}

template <class V, class M>
void
BaseScalarFunction<V, M>::setFiniteDifferenceStepSize(unsigned int i,
                                                      double fdStepSize)
{
  queso_require_greater_msg(fdStepSize, 0.0,
                            "Must provide a finite difference step > 0");

  queso_require_greater_equal_msg(i, 0, "Must provide a nonnegative index");

  unsigned int size = this->m_fdStepSize.size();
  queso_require_less_msg(i, size, "Must provide an index less than size of parameter dimension");

  this->m_fdStepSize[i] = fdStepSize;
}

}  // End namespace QUESO

template class QUESO::BaseScalarFunction<QUESO::GslVector, QUESO::GslMatrix>;
