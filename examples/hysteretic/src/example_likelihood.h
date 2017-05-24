/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __EX_LIKELIHOOD_H__
#define __EX_LIKELIHOOD_H__

#include <queso/GslMatrix.h>
#include <queso/ScalarFunction.h>
#include <cmath>

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet)
  {
  }

  virtual ~Likelihood() { }

  virtual double lnValue(const V & paramValues) const;
  virtual double actualValue(const V & paramValues,
                             const V * /* paramDirection */,
                             V *       /* gradVector */,
                             M *       /* hessianMatrix */,
                             V *       /* hessianEffect */) const
  {
    return std::exp(this->lnValue(paramValues));
  }

  std::vector<std::vector<double>* > floor;
              std::vector<double>    accel;
};

#endif
