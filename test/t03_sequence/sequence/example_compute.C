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

#include <example_compute.h>
#include <uqGslMatrix.h>
#include <uqVectorRV.h>

void compute(const uqFullEnvironmentClass& env) {
  // Step 1 of 9: Instantiate the parameter space
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 9: Instantiate the parameter domain
  uqGslVectorClass paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  uqGslVectorClass paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  uqBoxSubsetClass<uqGslVectorClass,uqGslMatrixClass>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 9: Instantiate the vector RV
  uqGslVectorClass meanVector(paramSpace.zeroVector());
  meanVector[0] = -1;
  meanVector[1] =  2;
  uqGslMatrixClass* covMatrix = paramSpace.newMatrix();
  (*covMatrix)(0,0) = 4.; (*covMatrix)(0,1) = 0.;
  (*covMatrix)(1,0) = 0.; (*covMatrix)(1,1) = 1.;
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>
    auxRv("", paramDomain,meanVector,*covMatrix);

  // Return
  delete covMatrix;

  return;
}
