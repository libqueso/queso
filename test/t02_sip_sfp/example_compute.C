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
// $Id$
//
//--------------------------------------------------------------------------

#include <example_compute.h>
#include <example_likelihood.h>
#include <example_qoi.h>
#include <uqGslMatrix.h>
#include <uqStatisticalInverseProblem.h>
#include <uqStatisticalForwardProblem.h>

void compute(const QUESO::FullEnvironmentClass& env) {
  // Step 1 of 9: Instantiate the parameter space
  QUESO::VectorSpaceClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 9: Instantiate the parameter domain
  QUESO::GslVectorClass paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  QUESO::GslVectorClass paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  QUESO::BoxSubsetClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 9: Instantiate the likelihood function object
  QUESO::GslVectorClass meanVector(paramSpace.zeroVector());
  meanVector[0] = -1;
  meanVector[1] =  2;
  QUESO::GslMatrixClass covMatrix(paramSpace.zeroVector());
  covMatrix(0,0) = 4.; covMatrix(0,1) = 0.;
  covMatrix(1,0) = 0.; covMatrix(1,1) = 1.;
  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.meanVector = &meanVector;
  likelihoodRoutine_Data.covMatrix  = &covMatrix;
  QUESO::GenericScalarFunctionClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [ln(function)]

  // Step 4 of 9: Instantiate the inverse problem
  QUESO::UniformVectorRVClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    priorRv("prior_", paramDomain);
  QUESO::GenericVectorRVClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    postRv("post_", paramSpace);
  QUESO::StatisticalInverseProblemClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  // Step 5 of 9: Solve the inverse problem
  QUESO::GslVectorClass paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.1;
  paramInitials[1] = -1.4;
  QUESO::GslMatrixClass proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = 8.; proposalCovMatrix(0,1) = 4.;
  proposalCovMatrix(1,0) = 4.; proposalCovMatrix(1,1) = 16.;
  ip.solveWithBayesMetropolisHastings(NULL,paramInitials, &proposalCovMatrix);

  // Step 6 of 9: Instantiate the qoi space
  QUESO::VectorSpaceClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    qoiSpace(env, "qoi_", 1, NULL);

  // Step 7 of 9: Instantiate the qoi function object
  qoiRoutine_DataType qoiRoutine_Data;
  qoiRoutine_Data.coef1 = 1.;
  qoiRoutine_Data.coef2 = 1.;
  QUESO::GenericVectorFunctionClass<QUESO::GslVectorClass,QUESO::GslMatrixClass,
                               QUESO::GslVectorClass,QUESO::GslMatrixClass>
    qoiFunctionObj("qoi_",
                   paramDomain,
                   qoiSpace,
                   qoiRoutine,
                   (void *) &qoiRoutine_Data);

  // Step 8 of 9: Instantiate the forward problem
  QUESO::GenericVectorRVClass<QUESO::GslVectorClass,QUESO::GslMatrixClass>
    qoiRv("qoi_", qoiSpace);
  QUESO::StatisticalForwardProblemClass<QUESO::GslVectorClass,QUESO::GslMatrixClass,
                                   QUESO::GslVectorClass,QUESO::GslMatrixClass>
    fp("", NULL, postRv, qoiFunctionObj, qoiRv);

  // Step 9 of 9: Solve the forward problem
  fp.solveWithMonteCarlo(NULL);

  return;
}
