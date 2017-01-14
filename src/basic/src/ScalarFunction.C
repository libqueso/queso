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

#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoostInputOptionsParser.h>

#define QUESO_BASESCALARFN_FD_STEPSIZE_ODV 1e-6

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
    m_fdStepSize(QUESO_BASESCALARFN_FD_STEPSIZE_ODV)
{
  // Snarf fd step size from input file.
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<double>(m_prefix + "fdStepSize",
                                   QUESO_BASESCALARFN_FD_STEPSIZE_ODV,
                                   "step size for finite difference");
  m_parser->scanInputFile();
  m_parser->getOption<double>(m_prefix + "fdStepSize", m_fdStepSize);
#else
  m_fdStepSize = m_env.input()(m_prefix + "fdStepSize",
                               QUESO_BASESCALARFN_FD_STEPSIZE_ODV);
#endif

  queso_require_greater_msg(m_fdStepSize,
                            0.0,
                            "Finite difference step size must be positive");
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

    perturbedVector[i] += m_fdStepSize;
    gradVector[i] = (this->lnValue(perturbedVector) - value) / m_fdStepSize;

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

}  // End namespace QUESO

template class QUESO::BaseScalarFunction<QUESO::GslVector, QUESO::GslMatrix>;
