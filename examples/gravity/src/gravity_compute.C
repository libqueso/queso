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
 * This file is divided in two parts:
 * - the first one handles the statistical inverse problem (SIP) for estimating
 *   the magnitude 'g' of gravity acceleration; and
 * - the second part handles the statistical forward problem (SFP) for
 *   predicting the maximum distance traveled by a projectile.
 *
 * The SIP definition requires a user defined likelihood function; refer to
 * files 'gravity_likelihood.h' and 'gravity_likelihood.C'. The SFP definition
 * requires a user defined qoi function; refer to files 'gravity_qoi.h' and
 * 'gravity_qoi.C'.
 */

#include <cmath>
#include <sys/time.h>

#include <queso/GslMatrix.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

#include <gravity_compute.h>
#include <gravity_likelihood.h>
#include <gravity_qoi.h>

void computeGravityAndTraveledDistance(const QUESO::FullEnvironment& env) {
  struct timeval timevalNow;

  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\nBeginning run of 'Gravity + Projectile motion' example at "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  // Just examples of possible calls
  if ((env.subDisplayFile()       ) &&
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Beginning run of 'Gravity + Projectile motion' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  env.fullComm().Barrier();
  env.subComm().Barrier();  // Just an example of a possible call

  //================================================================
  // Statistical inverse problem (SIP): find posterior PDF for 'g'
  //================================================================
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> Gravity estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

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
  gettimeofday(&timevalNow, NULL);
  std::cout << "Beginning 'SFP -> Projectile motion' at "
            << ctime(&timevalNow.tv_sec)
            << std::endl;

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
  std::cout << "Solving the SFP with Monte Carlo"
            << std::endl << std::endl;
  fp.solveWithMonteCarlo(NULL);

  gettimeofday(&timevalNow, NULL);
  if ((env.subDisplayFile()       ) &&
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Ending run of 'Gravity + Projectile motion' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  if (env.fullRank() == 0) {
    std::cout << "Ending run of 'Gravity + Projectile motion' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }
}
