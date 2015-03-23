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
#include <example_hyst.h>

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  const QUESO::BaseEnvironment& env = paramValues.env();
  UQ_FATAL_TEST_MACRO(paramValues.sizeLocal() != 15,
                      env.fullRank(),
                      "example_likelihood()",
                      "invalid parameter size");

  const std::vector<std::vector<double>* >&
    floor = ((likelihoodRoutine_DataType *) functionDataPtr)->floor;
  const std::vector<double>&
    accel = ((likelihoodRoutine_DataType *) functionDataPtr)->accel;

  unsigned int numFloors    = floor.size();
  unsigned int numTimeSteps = accel.size();

  UQ_FATAL_TEST_MACRO((numFloors != 4),
                      env.fullRank(),
                      "example_likelihood()",
                      "invalid 'numFloors'");

  UQ_FATAL_TEST_MACRO((numTimeSteps != 401),
                      env.fullRank(),
                      "example_likelihood()",
                      "invalid 'numTimeSteps'");

  for (unsigned int i = 0; i < numFloors; ++i) {
    UQ_FATAL_TEST_MACRO(floor[i]->size() != numTimeSteps,
                        env.fullRank(),
                        "example_likelihood()",
                        "invalid number of steps");
  }

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> floorSpace(env, "floor_", numFloors, NULL);

  double sigmaSq = paramValues[0];

  QUESO::GslVector kVec(floorSpace.zeroVector());
  kVec[0] = 2.20e+7 * exp(paramValues[1]);
  kVec[1] = 2.00e+7 * exp(paramValues[2]);
  kVec[2] = 1.70e+7 * exp(paramValues[3]);
  kVec[3] = 1.45e+7 * exp(paramValues[4]);

  QUESO::GslVector rVec(floorSpace.zeroVector());
  rVec[0] = 1.0e-1 * exp(paramValues[5]);
  rVec[1] = 1.0e-1 * exp(paramValues[6]);
  rVec[2] = 1.0e-1 * exp(paramValues[7]);
  rVec[3] = 1.0e-1 * exp(paramValues[8]);

  QUESO::GslVector uVec(floorSpace.zeroVector());
  uVec[0] = 8.0e-3 * exp(paramValues[9]);
  uVec[1] = 8.0e-3 * exp(paramValues[10]);
  uVec[2] = 7.0e-3 * exp(paramValues[11]);
  uVec[3] = 7.0e-3 * exp(paramValues[12]);

  double rho   = 7.959e-1 * exp(paramValues[13]);

  double gamma = 2.500e-3 * exp(paramValues[14]);

  std::vector<double> t(numTimeSteps,0.);
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> u     (floorSpace,numTimeSteps,""); // absolute displacement
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> ud    (floorSpace,numTimeSteps,""); // velocity
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> udd   (floorSpace,numTimeSteps,""); // acceleration
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> resfor(floorSpace,numTimeSteps,""); // restoring force
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> ru    (floorSpace,numTimeSteps,""); // relative displacement

  QUESO::GslVector massVec(floorSpace.zeroVector());
  massVec.cwSet(2.0e+4);

  hystereticModel(env,
                  massVec,
                  kVec,
                  rVec,
                  uVec,
                  rho,
                  gamma,
                  accel,
                  t, // output
                  u,
                  ud,
                  udd,
                  resfor,
                  ru);

  QUESO::GslVector auxVec(floorSpace.zeroVector());

  double sum = 0.;
  for (unsigned int i = 0; i < numFloors; ++i) {
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
    udd.getPositionValues(j,auxVec);
    sum += ((*floor[i])[j]-auxVec[i]-accel[j])*((*floor[i])[j]-auxVec[i]-accel[j]);
    }
  }

  double result = -0.5*((double) numFloors)*((double) numTimeSteps)*log(2.*M_PI*sigmaSq) - 0.5*sum/sigmaSq;
  return result;
}
