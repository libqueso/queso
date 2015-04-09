//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <simple_sfp_example_compute.h>
#include <simple_sfp_example_qoi.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalForwardProblem.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>

void compute(const QUESO::FullEnvironment& env) {

  // Step 1 of 6: Instantiate the parameter space
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 6: Instantiate the parameter domain
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 6: Instantiate the qoi space
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    qoiSpace(env, "qoi_", 1, NULL);

  // Step 4 of 6: Instantiate the qoi function object
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

  // Step 5 of 6: Instantiate the forward problem
  // Parameters are Gaussian RV
  QUESO::GslVector meanVector( paramSpace.zeroVector() );
  meanVector[0] = -1;
  meanVector[1] = 2;

  QUESO::GslMatrix covMatrix = QUESO::GslMatrix(paramSpace.zeroVector());
  covMatrix(0,0) = 4.;
  covMatrix(0,1) = 0.;
  covMatrix(1,0) = 0.;
  covMatrix(1,1) = 1.;

  QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    paramRv("param_",paramDomain,meanVector,covMatrix);

  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    qoiRv("qoi_", qoiSpace);

  QUESO::StatisticalForwardProblem<QUESO::GslVector,QUESO::GslMatrix,
                                   QUESO::GslVector,QUESO::GslMatrix>
    fp("", NULL, paramRv, qoiFunctionObj, qoiRv);

  // Step 6 of 6: Solve the forward problem
  fp.solveWithMonteCarlo(NULL);

  return;
}
