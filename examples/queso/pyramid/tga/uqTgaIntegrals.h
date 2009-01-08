/* uq/examples/queso/pyramid/uqTgaIntegrals.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_INTEGRALS_H__
#define __UQ_TGA_INTEGRALS_H__

#include <uqTgaDefines.h>
#include <uqTgaW.h>
#include <uqTgaLambda.h>
#include <uqDefines.h>

template<class P_V,class P_M>
void
uqTgaIntegrals(
  const P_V&                       paramValues,
  const P_V*                       paramDirection,
  const uqBase1D1DFunctionClass&   temperatureFunctionObj,
  double                           lowerIntegralTime,
  double                           upperIntegralTime,
  unsigned int                     reqNumIntervals,
  const uqTgaWClass     <P_V,P_M>& wObj,
  const uqTgaLambdaClass<P_V,P_M>& lambdaObj,
  const uqBase1D1DFunctionClass*   weightFunction,
  P_V*                             reducedGrad,
  P_M*                             reducedHessian,
  P_V*                             hessianEffect)
{
  UQ_FATAL_TEST_MACRO((hessianEffect != NULL) && (paramDirection == NULL),
                      paramValues.env().rank(),
                      "uqTgaIntegrals<P_V,P_M>()",
                      "hessianEffect is being requested but no paramDirection is supplied");

  // Set all possible return values to zero
  if (reducedGrad   ) *reducedGrad    *= 0.;
  if (reducedHessian) *reducedHessian *= 0.;

  unsigned int wSize = wObj.times().size();
  unsigned int lambdaSize = lambdaObj.times().size();

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering uqTgaIntegrals()"
              << ", wObj.times()[0] = "        << wObj.times()[0]
              << ", wObj.times()[max] = "      << wObj.times()[wSize-1]
              << ", lambdaObj.times()[0] = "   << lambdaObj.times()[0]
              << ", lambdaObj.times()[max] = " << lambdaObj.times()[lambdaSize-1]
              << ", lowerIntegralTime = "      << lowerIntegralTime
              << ", upperIntegralTime = "      << upperIntegralTime
              << ", reqNumIntervals = "        << reqNumIntervals
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO((wObj.times()[0] != lambdaObj.times()[0]),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "wObj.times()[0] and lambdaObj.times()[0] are different");

  UQ_FATAL_TEST_MACRO((wObj.times()[0] > lowerIntegralTime),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "wObj.times()[0] is too large");

  UQ_FATAL_TEST_MACRO((lambdaObj.times()[0] > lowerIntegralTime),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "lambdaObj.times()[0] is too large");

  UQ_FATAL_TEST_MACRO((wObj.times()[wSize-1] < upperIntegralTime),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "wObj.times()[last] is too small");

  UQ_FATAL_TEST_MACRO((lambdaObj.times()[lambdaSize-1] < upperIntegralTime),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "lambdaObj.times()[last] is too small");

  double A = paramValues[0];
  double E = paramValues[1];

  unsigned int currentWIntervalId      = 0;
  unsigned int currentLambdaIntervalId = 0;

  double w = 0.;
  P_V wGrad(paramValues);
  P_V* wGradPtr = &wGrad;
  double wDir = 0.;
  double* wDirPtr = &wDir;

  double lambda = 0.;
  P_V lambdaGrad(paramValues);
  P_V* lambdaGradPtr = &lambdaGrad;
  double lambdaDir = 0.;
  double* lambdaDirPtr = &lambdaDir;

  if (reducedHessian == NULL) {
    wGradPtr      = NULL;
    lambdaGradPtr = NULL;
  }

  if (hessianEffect == NULL) {
    lambdaDirPtr = NULL;
  }

  const uqDeltaSeq1D1DFunctionClass* deltaSeqFunction = NULL;
  if (weightFunction) {
    deltaSeqFunction = dynamic_cast< const uqDeltaSeq1D1DFunctionClass* >(weightFunction);
  }

  unsigned int numIntervals = reqNumIntervals;
  if (reqNumIntervals == 0) numIntervals = 1000;

  if (deltaSeqFunction == NULL) {
    double timeIntervalSize = (upperIntegralTime-lowerIntegralTime)/((double)numIntervals);
    double firstTime = lowerIntegralTime+.5*timeIntervalSize;
 
    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      std::cout << "In uqTgaIntegrals()"
                << ", continuous weight case"
                << ": beginning integration loop on time interval "
                << "["          << lowerIntegralTime
                << ", "         << upperIntegralTime
                << "], with = " << numIntervals
                << " subintervals"
                << "; wObj.times().size() = "      << wObj.times().size()
                << ", lambdaObj.times().size() = " << lambdaObj.times().size()
                << std::endl;
    }

    for (unsigned int i = 0; i < numIntervals; ++i) {
      double time = firstTime + ((double) i)*timeIntervalSize;
      double temp = temperatureFunctionObj.value(time);

      wObj.interpolate(time,
                       currentWIntervalId,
                       1.,
                       &w,
                       wGradPtr,
                       wDirPtr,
                       NULL);

      lambdaObj.interpolate(time,
                            currentLambdaIntervalId,
                            &lambda,
                            lambdaGradPtr,
                            lambdaDirPtr,
                            NULL);

      double expTerm = exp(-E/(R_CONSTANT*temp));
      if (reducedGrad) {
        (*reducedGrad)[0] +=                          lambda * w * expTerm;
        (*reducedGrad)[1] += -(A/(R_CONSTANT*temp)) * lambda * w * expTerm;
      }
      if (reducedHessian) {
        double wA      = (*wGradPtr     )[0];
        double wE      = (*wGradPtr     )[1];
        double lambdaA = (*lambdaGradPtr)[0];
        double lambdaE = (*lambdaGradPtr)[1];
        (*reducedHessian)(0,0) += (lambdaA*w + lambda*wA) * expTerm;
        (*reducedHessian)(0,1) += ( -lambda*w/(R_CONSTANT*temp) + lambdaE*w + lambda*wE ) * expTerm; // (d/dE)(dF/dA)
        (*reducedHessian)(1,0) += ( -lambda*w/(R_CONSTANT*temp) - A*(lambdaA*w + lambda*wA)/(R_CONSTANT*temp) ) * expTerm; // (d/dA)(dF/dE)
        (*reducedHessian)(1,1) += ( A*lambda*w/(R_CONSTANT*R_CONSTANT*temp*temp) - A*(lambdaE*w + lambda*wE)/(R_CONSTANT*temp) ) * expTerm;
      }
    }
    if (reducedGrad   ) (*reducedGrad   ) *= timeIntervalSize;
    if (reducedHessian) (*reducedHessian) *= timeIntervalSize;
  }
  else {
    // Determine number of milestone points (= number of milestone intervals + 1)
    unsigned int numDeltaTimes = deltaSeqFunction->domainValues().size();
    const std::vector<double>& deltaTimes = deltaSeqFunction->domainValues();

    unsigned int jMin = 0;
    unsigned int jMax = 0;
    for (jMin = 0; jMin < numDeltaTimes; ++jMin) {
      if (lowerIntegralTime < deltaTimes[jMin]) break;
    }
    for (jMax = jMin; jMax < numDeltaTimes; ++jMax) {
      if (upperIntegralTime <= deltaTimes[jMax]) break;
    }

    unsigned int numMilestones = jMax-jMin+2;
    std::vector<double> milestones(numMilestones,0.);
    milestones[0              ] = lowerIntegralTime;
    milestones[numMilestones-1] = upperIntegralTime;
    for (unsigned int j = 1; j < (numMilestones-1); ++j) {
      milestones[j] = deltaTimes[jMin+j-1];
    }

    // Determine the number of subintervals per milestone interval: needs to be at least 10
    double aux = ((double) numIntervals)/((double) (numMilestones-1));
    aux = std::max(aux,9.5);
    numIntervals = (unsigned int) (aux+0.5);

    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      std::cout << "In uqTgaIntegrals()"
                << ", delta seq weight case"
                << ": beginning integration loop on time interval "
                << "["          << lowerIntegralTime
                << ", "         << upperIntegralTime
                << "], with = " << numIntervals
                << " subintervals for each of the " << numMilestones-1
                << " milestone intervals"
                << "; wObj.times().size() = "      << wObj.times().size()
                << ", lambdaObj.times().size() = " << lambdaObj.times().size()
                << ", deltaTimes[0] = " << deltaTimes[0]
                << ", deltaTimes[max] = " << deltaTimes[numDeltaTimes-1]
                << ", jMin = " << jMin
                << ", jMax = " << jMax
                << std::endl;
    }

    for (unsigned int j = 1; j < numMilestones; ++j) { // Yes, from '1'
      double timeIntervalSize = (milestones[j] - milestones[j-1])/((double)numIntervals);
      double firstTime = milestones[j-1]+.5*timeIntervalSize;

      for (unsigned int i = 0; i < numIntervals; ++i) {
        double time = firstTime + ((double) i)*timeIntervalSize;
        double temp = temperatureFunctionObj.value(time);

        wObj.interpolate(time,
                         currentWIntervalId,
                         1.,
                         &w,
                         wGradPtr,
                         wDirPtr,
                         NULL);

        lambdaObj.interpolate(time,
                              currentLambdaIntervalId,
                              &lambda,
                              lambdaGradPtr,
                              lambdaDirPtr,
                              NULL);

        double expTerm = exp(-E/(R_CONSTANT*temp));
        if (reducedGrad) {
          (*reducedGrad)[0] +=                          lambda * w * expTerm * timeIntervalSize;
          (*reducedGrad)[1] += -(A/(R_CONSTANT*temp)) * lambda * w * expTerm * timeIntervalSize;
        }
        if (reducedHessian) {
          double wA      = (*wGradPtr     )[0];
          double wE      = (*wGradPtr     )[1];
          double lambdaA = (*lambdaGradPtr)[0];
          double lambdaE = (*lambdaGradPtr)[1];
          (*reducedHessian)(0,0) += (lambdaA*w + lambda*wA) * expTerm * timeIntervalSize;
          (*reducedHessian)(0,1) += ( -lambda*w/(R_CONSTANT*temp) + lambdaE*w + lambda*wE ) * expTerm * timeIntervalSize; // (d/dE)(dF/dA)
          (*reducedHessian)(1,0) += ( -lambda*w/(R_CONSTANT*temp) - A*(lambdaA*w + lambda*wA)/(R_CONSTANT*temp) ) * expTerm * timeIntervalSize; // (d/dA)(dF/dE)
          (*reducedHessian)(1,1) += ( A*lambda*w/(R_CONSTANT*R_CONSTANT*temp*temp) - A*(lambdaE*w + lambda*wE)/(R_CONSTANT*temp) ) * expTerm * timeIntervalSize;
        }
      } // i
    } // j
  }

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "In uqTgaIntegrals()"
              << ": finsihed integration loop with"
              << " currentWIntervalId = " << currentWIntervalId
              << ", currentLambdaIntervalId = " << currentLambdaIntervalId
              << std::endl;
  }

#if 0
  UQ_FATAL_TEST_MACRO((currentWIntervalId != (allWTimes.size()-2)),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "currentWIntervalId finished with an invalid value");

  UQ_FATAL_TEST_MACRO((currentLambdaIntervalId != (allLambdaTimes.size()-2)),
                      paramValues.env().rank(),
                      "uqTgaIntegrals()",
                      "currentLambdaIntervalId finished with an invalid value");
#endif

  return;
}
#endif // __UQ_TGA_INTEGRALS_H__
