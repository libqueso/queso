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

/*
 * Brief description of this file:
 *
 * The SIP consists on calibrating the magnitude 'g' of acceleration gravity
 * using measurements of the time that it takes for an object in free fall to
 * reach the ground from a given height and zero initial velocity. The
 * solution of the SIP is the posterior probability density function (PDF) of
 * 'g'.
 *
 * The SFP consists of calculating the maximum distance traveled by an object
 * in projectile motion. The posterior PDF of 'g' from the SIP might be used
 * as input to the SFP.
 *
 * The code consists of 7 files:
 * - 'gravity_main.C' (this file)
 * - 'gravity_likelihood.C' (necessary for the SIP)
 * - 'gravity_likelihood.h'
 * - 'gravity_qoi.C' (necessary for the SFP)
 * - 'gravity_qoi.h'
 */

#include <gravity_likelihood.h>
#include <gravity_qoi.h>

#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

#include <cmath>
#include <cstdlib>

int main(int argc, char* argv[])
{
#ifndef QUESO_HAS_MPI
  // Skip this test if we're not in parallel
  return 77;
#else

  MPI_Init(&argc, &argv);

  std::string inputFileName = argv[1];
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

  // Initialize QUESO environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);

  //================================================================
  // Statistical inverse problem (SIP): find posterior PDF for 'g'
  //================================================================

  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<> paramSpace(env, "param_", 1, NULL);

  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());

  paramMinValues[0] = 8.;
  paramMaxValues[0] = 11.;

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMinValues,
      paramMaxValues);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function
  // object to be used by QUESO.
  //------------------------------------------------------
  Likelihood<> lhood("like_", paramDomain);

  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  // Extra prefix before the default "rv_" prefix
  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  // No extra prefix before the default "ip_" prefix
  QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV
  //------------------------------------------------------
  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  priorRv.realizer().realization(paramInitials);

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = std::pow(std::abs(paramInitials[0]) / 20.0, 2.0);

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  //================================================================
  // Statistical forward problem (SFP): find the max distance
  // traveled by an object in projectile motion; input pdf for 'g'
  // is the solution of the SIP above.
  //================================================================

  //------------------------------------------------------
  // SFP Step 1 of 6: Instantiate the parameter *and* qoi spaces.
  // SFP input RV = FIP posterior RV, so SFP parameter space
  // has been already defined.
  //------------------------------------------------------
  QUESO::VectorSpace<> qoiSpace(env, "qoi_", 1, NULL);

  //------------------------------------------------------
  // SFP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------

  // Not necessary because input RV of the SFP = output RV of SIP.
  // Thus, the parameter domain has been already defined.

  //------------------------------------------------------
  // SFP Step 3 of 6: Instantiate the qoi object
  // to be used by QUESO.
  //------------------------------------------------------
  Qoi<> qoi("qoi_", paramDomain, qoiSpace);

  //------------------------------------------------------
  // SFP Step 4 of 6: Define the input RV
  //------------------------------------------------------

  // Not necessary because input RV of SFP = output RV of SIP
  // (postRv).

  //------------------------------------------------------
  // SFP Step 5 of 6: Instantiate the forward problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<> qoiRv("qoi_", qoiSpace);

  QUESO::StatisticalForwardProblem<> fp("", NULL, postRv, qoi, qoiRv);

  //------------------------------------------------------
  // SFP Step 6 of 6: Solve the forward problem
  //------------------------------------------------------
  fp.solveWithMonteCarlo(NULL);

  MPI_Finalize();

  return 0;
#endif  // QUESO_HAS_MPI
}
