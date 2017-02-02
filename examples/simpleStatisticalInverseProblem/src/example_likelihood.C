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

#include <example_likelihood.h>

template <class V, class M>
double
Likelihood<V, M>::lnValue(const V & paramValues) const
{
  // Just checking: the user, at the application level, expects
  // vector 'paramValues' to have size 2.
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 2,
                      QUESO::UQ_UNAVAILABLE_RANK,
                      "likelihoodRoutine()",
                      "paramValues vector does not have size 2");

  // Actual code
  //
  // This code exemplifies multiple Metropolis-Hastings solvers, each calling this likelihood
  // routine. In this simple example, only node 0 in each subenvironment does the job even
  // though there might be more than one node per sub-environment. In a more realistic
  // situation, if the user is asking for multiple nodes per subenvironment, then the model
  // code in the likelihood routines might really demand more than one node. Here we use
  // 'env.subRank()' only. A realistic application might want to use either 'env.subComm()'
  // or 'env.subComm().Comm()'.

  double result = 0.;
  const QUESO::BaseEnvironment& env = paramValues.env();
  if (env.subRank() == 0) {

    QUESO::GslVector diffVec(paramValues - *meanVector);

    result= scalarProduct(diffVec,covMatrix->invertMultiply(diffVec));
  }
  else {
    // Do nothing;
  }

  return -.5*result;
}

template class Likelihood<QUESO::GslVector, QUESO::GslMatrix>;
