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

#include <example_compute.h>
#include <example_likelihood.h>
#include <example_qoi.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>

void compute(const QUESO::FullEnvironment& env) {
  // Step 1 of 9: Instantiate the parameter space
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 9: Instantiate the parameter domain
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 9: Instantiate the likelihood function object
  QUESO::GslVector meanVector(paramSpace.zeroVector());
  meanVector[0] = -1;
  meanVector[1] =  2;
  QUESO::GslMatrix covMatrix(paramSpace.zeroVector());
  covMatrix(0,0) = 4.; covMatrix(0,1) = 0.;
  covMatrix(1,0) = 0.; covMatrix(1,1) = 1.;
  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.meanVector = &meanVector;
  likelihoodRoutine_Data.covMatrix  = &covMatrix;
  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [ln(function)]

  // Step 4 of 9: Instantiate the inverse problem
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_", paramDomain);
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_", paramSpace);
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  // Step 5 of 9: Solve the inverse problem
  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.1;
  paramInitials[1] = -1.4;
  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = 8.; proposalCovMatrix(0,1) = 4.;
  proposalCovMatrix(1,0) = 4.; proposalCovMatrix(1,1) = 16.;
  ip.solveWithBayesMetropolisHastings(NULL,paramInitials, &proposalCovMatrix);

  // Step 6 of 9: Instantiate the qoi space
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    qoiSpace(env, "qoi_", 1, NULL);

  // Step 7 of 9: Instantiate the qoi function object
  qoiRoutine_DataType qoiRoutine_Data;
  qoiRoutine_Data.coef1 = 1.;
  qoiRoutine_Data.coef2 = 1.;
  QUESO::GenericVectorFunction<QUESO::GslVector,QUESO::GslMatrix,
                               QUESO::GslVector,QUESO::GslMatrix>
    qoiFunctionObj("qoi_",
                   paramDomain,
                   qoiSpace,
                   qoiRoutine,
                   (void *) &qoiRoutine_Data);

  // Step 8 of 9: Instantiate the forward problem
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    qoiRv("qoi_", qoiSpace);
  QUESO::StatisticalForwardProblem<QUESO::GslVector,QUESO::GslMatrix,
                                   QUESO::GslVector,QUESO::GslMatrix>
    fp("", NULL, postRv, qoiFunctionObj, qoiRv);

  // Step 9 of 9: Solve the forward problem
  fp.solveWithMonteCarlo(NULL);

  return;
}
