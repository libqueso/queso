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
#include <uqTgaStorageW.h>
#include <uqTgaComputableW.h>
#include <uqTgaLambda.h>
#include <uqDefines.h>

template<class P_V,class P_M>
void
uqTgaLagrangianGradientWrtDesignParameters(
  const P_V&                     paramValues,
  const uqBase1D1DFunctionClass& temperatureFunctionObj,
  double                         upperIntegralLimit,
  const std::vector<double>&     allWTimes,
  const std::vector<double>&     allWs,
  const std::vector<double>&     allLambdaTimes,
  const std::vector<double>&     allLambdas,
  P_V&                           LagrangianGrad)
{
  std::cout << "Entering uqTgaLagrangianGradientWrtDesignParameters()"
            << ", allWTimes[0] = "        << allWTimes[0]
            << ", allWTimes[max] = "      << allWTimes[allWTimes.size()-1]
            << ", allLambdaTimes[0] = "   << allLambdaTimes[0]
            << ", allLambdaTimes[max] = " << allLambdaTimes[allLambdaTimes.size()-1]
            << ", upperIntegralLimit = "  << upperIntegralLimit
            << std::endl;

  UQ_FATAL_TEST_MACRO((allWTimes[0] != allLambdaTimes[0]),
                      paramValues.env().rank(),
                      "uqTgaLagrangianGradientWrtDesignParameters()",
                      "allWTimes[0] and allLambdaTimes[0] are different");

  UQ_FATAL_TEST_MACRO((allWTimes[allWTimes.size()-1] < upperIntegralLimit),
                      paramValues.env().rank(),
                      "uqTgaLagrangianGradientWrtDesignParameters()",
                      "allWTimes[last] is too small");

  UQ_FATAL_TEST_MACRO((allLambdaTimes[allLambdaTimes.size()-1] < upperIntegralLimit),
                      paramValues.env().rank(),
                      "uqTgaLagrangianGradientWrtDesignParameters()",
                      "allLambdaTimes[last] is too small");

  double A = paramValues[0];
  double E = paramValues[1];

  unsigned int currentWIntervalId      = 0;
  unsigned int currentLambdaIntervalId = 0;

  unsigned int numIntervals = 1000;
  double timeIntervalSize = (upperIntegralLimit-allWTimes[0])/((double)numIntervals);
  double firstTime = allWTimes[0]+.5*timeIntervalSize;

  std::cout << "In tgaAdjointEquation()"
            << ": beginning integration loop on time interval "
            << "["          << allWTimes[0]
            << ", "         << upperIntegralLimit
            << "], with = " << numIntervals
            << " subintervals"
            << "; allWTimes.size() = "      << allWTimes.size()
            << ", allLambdaTimes.size() = " << allLambdaTimes.size()
            << std::endl;

  LagrangianGrad *= 0.;

  for (unsigned int i = 0; i < numIntervals; ++i) {
    double time = firstTime + ((double) i)*timeIntervalSize;
    double temp = temperatureFunctionObj.value(time);

    while (allWTimes[currentWIntervalId+1] <= time) {
      currentWIntervalId++;
    }
    if ((allWTimes[currentWIntervalId] <= time ) &&
        (time < allWTimes[currentWIntervalId+1])) {
      // Ok
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          paramValues.env().rank(),
                          "uqTgaLagrangianGradientWrtDesignParameters()",
                          "invalid situation with allWTimes");
    }
    double wTimeDelta = allWTimes[currentWIntervalId+1] - allWTimes[currentWIntervalId];
    double wTimeRatio = (time - allWTimes[currentWIntervalId])/wTimeDelta;
    double wDelta = allWs[currentWIntervalId+1] - allWs[currentWIntervalId];
    double w      = allWs[currentWIntervalId] + wTimeRatio * wDelta;

    while (allLambdaTimes[currentLambdaIntervalId+1] < time) {
      currentLambdaIntervalId++;
    }
    if ((allLambdaTimes[currentLambdaIntervalId] <= time ) &&
        (time < allLambdaTimes[currentLambdaIntervalId+1])) {
      // Ok
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          paramValues.env().rank(),
                          "uqTgaLagrangianGradientWrtDesignParameters()",
                          "invalid situation with allLambdaTimes");
    }
    double lambdaTimeDelta = allLambdaTimes[currentLambdaIntervalId+1] - allLambdaTimes[currentLambdaIntervalId];
    double lambdaTimeRatio = (time - allLambdaTimes[currentLambdaIntervalId])/lambdaTimeDelta;
    double lambdaDelta = allLambdas[currentLambdaIntervalId+1] - allLambdas[currentLambdaIntervalId];
    double lambda      = allLambdas[currentLambdaIntervalId] + lambdaTimeRatio * lambdaDelta;

    LagrangianGrad[0] +=                          lambda * w * exp(-E/(R_CONSTANT*temp));
    LagrangianGrad[1] += -(A/(R_CONSTANT*temp)) * lambda * w * exp(-E/(R_CONSTANT*temp));
  }
  LagrangianGrad *= timeIntervalSize;

  std::cout << "In tgaAdjointEquation()"
            << ": finsihed integration loop with"
            << " currentWIntervalId = " << currentWIntervalId
            << ", currentLambdaIntervalId = " << currentLambdaIntervalId
            << std::endl;

#if 0
  UQ_FATAL_TEST_MACRO((currentWIntervalId != (allWTimes.size()-2)),
                      paramValues.env().rank(),
                      "uqTgaLagrangianGradientWrtDesignParameters()",
                      "currentWIntervalId finished with an invalid value");

  UQ_FATAL_TEST_MACRO((currentLambdaIntervalId != (allLambdaTimes.size()-2)),
                      paramValues.env().rank(),
                      "uqTgaLagrangianGradientWrtDesignParameters()",
                      "currentLambdaIntervalId finished with an invalid value");
#endif

  return;
}
#endif // __UQ_TGA_INTEGRALS_H__
