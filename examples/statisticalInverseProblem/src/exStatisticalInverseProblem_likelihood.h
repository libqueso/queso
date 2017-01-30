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
 *-------------------------------------------------------------------------- */

#ifndef EX_STATISTICAL_INVERSE_PROBLEM_LIKELIHOOD_H
#define EX_STATISTICAL_INVERSE_PROBLEM_LIKELIHOOD_H

template<class P_V, class P_M>
class Likelihood : public QUESO::BaseScalarFunction<P_V, P_M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<P_V, P_M> domainSet)
    : QUESO::BaseScalarFunction<P_V, P_M>(prefix, domainSet)
  {
    // Do nothing
  }

  virtual ~Likelihood()
  {
    // Do nothing
  }

  virtual double lnValue(const P_V & paramValues) const
  {
    double result = 0.0;

    const P_V& paramMeans        = *paramMeans;
    const P_M& matrix            = *matrix;
    bool       applyMatrixInvert =  applyMatrixInvert;

    P_V diffVec(paramValues - paramMeans);
    if (applyMatrixInvert) {
      result = scalarProduct(diffVec, matrix.invertMultiply(diffVec));
    }
    else {
      result = scalarProduct(diffVec, matrix * diffVec);
    }

    return -.5*result;
  }

  const P_V* paramMeans;
  const P_M* matrix;
  bool       applyMatrixInvert;
}

#endif // EX_STATISTICAL_INVERSE_PROBLEM_LIKELIHOOD_H
