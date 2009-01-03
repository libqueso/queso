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
  const uqBase1D1DFunctionClass&   temperatureFunctionObj,
  double                           lowerIntegralTime,
  double                           upperIntegralTime,
  unsigned int                     reqNumIntervals,
  const uqTgaWClass     <P_V,P_M>& wObj,
  const uqTgaLambdaClass<P_V,P_M>& lambdaObj,
  P_V*                             LagrangianGrad,
  P_M*                             LagrangianHessian)
{
  unsigned int wSize = wObj.times().size();
  unsigned int lambdaSize = lambdaObj.times().size();

  if ((paramValues.env().verbosity() >= 0) && (paramValues.env().rank() == 0)) {
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

  unsigned int numIntervals = reqNumIntervals;
  if (reqNumIntervals == 0) numIntervals = 1000;
  double timeIntervalSize = (upperIntegralTime-lowerIntegralTime)/((double)numIntervals);
  double firstTime = lowerIntegralTime+.5*timeIntervalSize;

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "In uqTgaIntegrals()"
              << ": beginning integration loop on time interval "
              << "["          << lowerIntegralTime
              << ", "         << upperIntegralTime
              << "], with = " << numIntervals
              << " subintervals"
              << "; wObj.times().size() = "      << wObj.times().size()
              << ", lambdaObj.times().size() = " << lambdaObj.times().size()
              << std::endl;
  }

  if (LagrangianGrad   ) *LagrangianGrad    *= 0.;
  if (LagrangianHessian) *LagrangianHessian *= 0.;

  unsigned int currentWIntervalId      = 0;
  unsigned int currentLambdaIntervalId = 0;

  double w = 0.;
  P_V wGrad(paramValues);
  P_V* wGradPtr = &wGrad;

  double lambda = 0.;
  P_V lambdaGrad(paramValues);
  P_V* lambdaGradPtr = &lambdaGrad;

  if (LagrangianHessian == NULL) {
    wGradPtr      = NULL;
    lambdaGradPtr = NULL;
  }

  for (unsigned int i = 0; i < numIntervals; ++i) {
    double time = firstTime + ((double) i)*timeIntervalSize;
    double temp = temperatureFunctionObj.value(time);

    wObj.interpolate(time,
                     currentWIntervalId,
                     &w,
                     wGradPtr,
                     NULL);

    lambdaObj.interpolate(time,
                          currentLambdaIntervalId,
                          &lambda,
                          lambdaGradPtr,
                          NULL);

    double expTerm = exp(-E/(R_CONSTANT*temp));
    if (LagrangianGrad) {
      (*LagrangianGrad)[0] +=                          lambda * w * expTerm;
      (*LagrangianGrad)[1] += -(A/(R_CONSTANT*temp)) * lambda * w * expTerm;
    }
    if (LagrangianHessian) {
      double wA      = (*wGradPtr     )[0];
      double wE      = (*wGradPtr     )[1];
      double lambdaA = (*lambdaGradPtr)[0];
      double lambdaE = (*lambdaGradPtr)[1];
      (*LagrangianHessian)(0,0) += (lambdaA*w + lambda*wA) * expTerm;
      (*LagrangianHessian)(0,1) += ( -lambda*w/(R_CONSTANT*temp) + lambdaE*w + lambda*wE ) * expTerm; // (d/dE)(dF/dA)
      (*LagrangianHessian)(1,0) += ( -lambda*w/(R_CONSTANT*temp) - A*(lambdaA*w + lambda*wA)/(R_CONSTANT*temp) ) * expTerm; // (d/dA)(dF/dE)
      (*LagrangianHessian)(1,1) += ( A*lambda*w/(R_CONSTANT*R_CONSTANT*temp*temp) - A*(lambdaE*w + lambda*wE)/(R_CONSTANT*temp) ) * expTerm;
    }
  }
  if (LagrangianGrad   ) (*LagrangianGrad   ) *= timeIntervalSize;
  if (LagrangianHessian) (*LagrangianHessian) *= timeIntervalSize;

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
