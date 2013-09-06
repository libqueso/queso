//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id: example_compute.C 37704 2013-03-08 21:10:36Z karl $
//
//--------------------------------------------------------------------------

#include <example_compute.h>
#include <example_likelihood.h>
#include <uqGslMatrix.h>
#include <uqStatisticalInverseProblem.h>

void compute(const QUESO::uqFullEnvironmentClass& env) {
  // Step 1 of 5: Instantiate the parameter space
  QUESO::uqVectorSpaceClass<QUESO::uqGslVectorClass,QUESO::uqGslMatrixClass>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 5: Instantiate the parameter domain
  QUESO::uqGslVectorClass paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  QUESO::uqGslVectorClass paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  QUESO::uqBoxSubsetClass<QUESO::uqGslVectorClass,QUESO::uqGslMatrixClass>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 5: Instantiate the likelihood function object
  QUESO::uqGslVectorClass meanVector(paramSpace.zeroVector());
  meanVector[0] = -1;
  meanVector[1] =  2;
  
  QUESO::uqGslMatrixClass covMatrix(paramSpace.zeroVector());
  covMatrix(0,0) = 4.; covMatrix(0,1) = 0.;
  covMatrix(1,0) = 0.; covMatrix(1,1) = 1.;
  
  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.meanVector = &meanVector;
  likelihoodRoutine_Data.covMatrix  = &covMatrix;
  
  QUESO::uqGenericScalarFunctionClass<QUESO::uqGslVectorClass,QUESO::uqGslMatrixClass>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [ln(function)]

  // Step 4 of 5: Instantiate the inverse problem
  QUESO::uqUniformVectorRVClass<QUESO::uqGslVectorClass,QUESO::uqGslMatrixClass>
    priorRv("prior_", paramDomain);
  QUESO::uqGenericVectorRVClass<QUESO::uqGslVectorClass,QUESO::uqGslMatrixClass>
    postRv("post_", paramSpace);
  QUESO::uqStatisticalInverseProblemClass<QUESO::uqGslVectorClass,QUESO::uqGslMatrixClass>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  // Step 5 of 5: Solve the inverse problem
  QUESO::uqGslVectorClass paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.1;
  paramInitials[1] = -1.4;
  
  QUESO::uqGslMatrixClass proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = 8.; proposalCovMatrix(0,1) = 4.;
  proposalCovMatrix(1,0) = 4.; proposalCovMatrix(1,1) = 16.;
  
  ip.solveWithBayesMetropolisHastings(NULL,paramInitials, &proposalCovMatrix);
 
  return;
}
