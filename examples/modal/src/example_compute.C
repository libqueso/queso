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
#include <example_likelihood.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/ConcatenatedVectorRV.h>
#include <queso/InverseGammaVectorRV.h>
#include <queso/ConcatenationSubset.h>

#define APPLS_MODAL_USES_CONCATENATION

void compute(const QUESO::FullEnvironment& env, unsigned int numModes) {
  ////////////////////////////////////////////////////////
  // Step 1 of 5: Instantiate the parameter space
  ////////////////////////////////////////////////////////


#ifdef APPLS_MODAL_USES_CONCATENATION
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpaceA(env, "paramA_", 2, NULL);
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpaceB(env, "paramB_", 1, NULL);
#endif
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace(env, "param_", 3, NULL);

  ////////////////////////////////////////////////////////
  // Step 2 of 5: Instantiate the parameter domain
  ////////////////////////////////////////////////////////
#ifdef APPLS_MODAL_USES_CONCATENATION
  QUESO::GslVector paramMinsA(paramSpaceA.zeroVector());
  paramMinsA[0] = 0.;
  paramMinsA[1] = 0.;
  QUESO::GslVector paramMaxsA(paramSpaceA.zeroVector());
  paramMaxsA[0] = 3.;
  paramMaxsA[1] = 3.;
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
  
  paramDomainA("paramA_",paramSpaceA,paramMinsA,paramMaxsA);

  QUESO::GslVector paramMinsB(paramSpaceB.zeroVector());
  paramMinsB[0] = 0.;
  QUESO::GslVector paramMaxsB(paramSpaceB.zeroVector());
  paramMaxsB[0] = INFINITY;
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomainB("paramB_",paramSpaceB,paramMinsB,paramMaxsB);

  QUESO::ConcatenationSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("",paramSpace,paramDomainA,paramDomainB);
#else
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = 0.;
  paramMins[1] = 0.;
  paramMins[2] = 0.;
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 3.;
  paramMaxs[1] = 3.;
  paramMaxs[2] = .3;
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);
#endif

  ////////////////////////////////////////////////////////
  // Step 3 of 5: Instantiate the likelihood function object
  ////////////////////////////////////////////////////////
  likelihoodRoutine_DataType likelihoodRoutine_Data;
 
  likelihoodRoutine_Data.numModes = numModes;
  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [-2.*ln(function)]
  
    ////////////////////////////////////////////////////////
  // Step 4 of 5: Instantiate the inverse problem
  ////////////////////////////////////////////////////////
#ifdef APPLS_MODAL_USES_CONCATENATION
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRvA("priorA_", paramDomainA);
 
  QUESO::GslVector alpha(paramSpaceB.zeroVector());
  alpha[0] = 3.;
  QUESO::GslVector beta(paramSpaceB.zeroVector());
  if (numModes == 1) {
    beta[0] = 0.09709133373799;
  }
  else {
    beta[0] = 0.08335837191688;
  }
  QUESO::InverseGammaVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRvB("priorB_", paramDomainB,alpha,beta);

  QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_", priorRvA, priorRvB, paramDomain);
#else
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_", paramDomain);
#endif

  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  ////////////////////////////////////////////////////////
  // Step 5 of 5: Solve the inverse problem
  ////////////////////////////////////////////////////////
  ip.solveWithBayesMLSampling();

  ////////////////////////////////////////////////////////
  // Print some statistics
  ////////////////////////////////////////////////////////
  unsigned int numPosTotal = postRv.realizer().subPeriod();
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "numPosTotal = " << numPosTotal
                          << std::endl;
  }

  QUESO::GslVector auxVec(paramSpace.zeroVector());
  unsigned int numPosTheta1SmallerThan1dot5 = 0;
  for (unsigned int i = 0; i < numPosTotal; ++i) {
    postRv.realizer().realization(auxVec);
    if (auxVec[0] < 1.5) numPosTheta1SmallerThan1dot5++;
  }

  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "numPosTheta1SmallerThan1dot5 = " << numPosTheta1SmallerThan1dot5
                          << ", ratio = " << ((double) numPosTheta1SmallerThan1dot5)/((double) numPosTotal)
                          << std::endl;
  }

  return;
}
