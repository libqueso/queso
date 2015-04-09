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

#include <example_likelihood.h>

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  double theta1  = paramValues[0];
  double theta2  = paramValues[1];
  double sigmaSq = paramValues[2];

  unsigned int numModes = ((likelihoodRoutine_DataType *) functionDataPtr)->numModes;

  double aux1 = theta1 + 2.*theta2;
  double aux2 = sqrt(theta1*theta1+4.*theta2*theta2);

  double w1 = 10.*sqrt(10.*(aux1+aux2));
  double w2 = 10.*sqrt(10.*(aux1-aux2));

  double sum1 = 0.;
  double sum2 = 0.;

  double aux = (w1 - 72.0470);
  sum1 += aux*aux;
  aux = (w1 - 71.8995);
  sum1 += aux*aux;
  aux = (w1 - 72.2801);
  sum1 += aux*aux;
  aux = (w1 - 71.9421);
  sum1 += aux*aux;
  aux = (w1 - 72.3578);
  sum1 += aux*aux;

  if (numModes == 1) {
    // Ok, do nothing
  }
  else if (numModes == 2) {
    aux = (w2 - 28.0292);
    sum2 += aux*aux;
    aux = (w2 - 27.3726);
    sum2 += aux*aux;
    aux = (w2 - 27.5388);
    sum2 += aux*aux;
    aux = (w2 - 27.0357);
    sum2 += aux*aux;
    aux = (w2 - 27.1588);
    sum2 += aux*aux;
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        paramValues.env().fullRank(),
                        "example_likelihood()",
                        "invalid 'numModes'");
  }

  double result = -0.5*((double) numModes)*5.*log(2.*M_PI*sigmaSq) - 0.5*(sum1+sum2)/sigmaSq;


  return result;
}
