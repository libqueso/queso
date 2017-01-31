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

#ifndef EX_LIKELIHOOD_H
#define EX_LIKELIHOOD_H

#include <queso/GslMatrix.h>
#include <queso/ScalarFunction.h>
#include <cmath>

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet)
  {}

  ~Likelihood() {};

  virtual double lnValue(const V & paramValues) const;
  virtual double actualValue(const V & paramValues,
                             const V * paramDirection,
                             V * gradVector,
                             M * hessianMatrix,
                             V * hessianEffect) const
  {
    return std::exp(this->lnValue(paramValues));
  }

  const QUESO::GslVector* meanVector;
  const QUESO::GslMatrix* covMatrix;
};

#endif
