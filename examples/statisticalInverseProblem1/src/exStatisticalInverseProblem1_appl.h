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

#ifndef __EX_STATISTICAL_INVERSE_PROBLEM_1_APPL_H__
#define __EX_STATISTICAL_INVERSE_PROBLEM_1_APPL_H__

#include <exStatisticalInverseProblem1_likelihood.h>
#include <uqStatisticalInverseProblem.h>
#include <uqAsciiTable.h>
#include <uqCovCond.h>

//********************************************************
// The driving routine: called by main()
//********************************************************
template<class P_V,class P_M>
void 
uqAppl(const uqBaseEnvironmentClass& env)
{
  if (env.fullRank() == 0) {
    std::cout << "Beginning run of 'exStatisticalInverseProblem1_example'\n"
              << std::endl;
  }

  //******************************************************
  // Step 1 of 5: Instantiate the parameter space
  //******************************************************
  uqVectorSpaceClass<P_V,P_M> paramSpace(env,"param_",4,NULL);

  //******************************************************
  // Step 2 of 5: Instantiate the parameter domain
  //******************************************************
  P_V paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  P_V paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",paramSpace,paramMins,paramMaxs);

  //******************************************************
  // Step 3 of 5: Instantiate the likelihood function object (data + routine), to be used by QUESO.
  //******************************************************
  P_V paramMeans(paramSpace.zeroVector());

  double condNumber = 100.0;
  P_V direction(paramSpace.zeroVector());
  direction.cwSet(1.);
  P_M* covMatrixInverse = paramSpace.newMatrix();
  P_M* covMatrix        = paramSpace.newMatrix();
  uqCovCond(condNumber,direction,*covMatrix,*covMatrixInverse);
  std::cout << "covMatrix = [" << *covMatrix
            << "];"
            << "\n"
            << "covMatrixInverse = [" << *covMatrixInverse
            << "];"
            << std::endl;

  likelihoodRoutine_DataType<P_V,P_M> likelihoodRoutine_Data;
  likelihoodRoutine_Data.paramMeans        = &paramMeans;
  likelihoodRoutine_Data.matrix            = covMatrixInverse;
  likelihoodRoutine_Data.applyMatrixInvert = false;

  uqGenericScalarFunctionClass<P_V,P_M> likelihoodFunctionObj("like_",
                                                              paramDomain,
                                                              likelihoodRoutine<P_V,P_M>,
                                                              (void *) &likelihoodRoutine_Data,
                                                              true); // the routine computes [-2.*ln(function)]

  //******************************************************
  // Step 4 of 5: Instantiate the inverse problem
  //******************************************************
  uqUniformVectorRVClass<P_V,P_M> priorRv("prior_", // Extra prefix before the default "rv_" prefix
                                          paramDomain);

  uqGenericVectorRVClass<P_V,P_M> postRv("post_", // Extra prefix before the default "rv_" prefix
                                         paramSpace);

  uqStatisticalInverseProblemClass<P_V,P_M> ip("", // No extra prefix before the default "ip_" prefix
                                               priorRv,
                                               likelihoodFunctionObj,
                                               postRv);

  //******************************************************
  // Step 5 of 5: Solve the inverse problem
  //******************************************************
  P_V paramInitials(paramSpace.zeroVector());
  P_M* proposalCovMatrix = postRv.imageSet().vectorSpace().newGaussianMatrix(NULL,&paramInitials);
  ip.solveWithBayesMarkovChain(paramInitials,
                               proposalCovMatrix);
  delete proposalCovMatrix;

  P_V tmpVec (paramSpace.zeroVector());
  P_V diffVec(paramSpace.zeroVector());
  const uqBaseVectorRealizerClass<P_V,P_M>& postRealizer = postRv.realizer();
  for (unsigned int i = 0; i < postRealizer.subPeriod(); ++i) {
    postRealizer.realization(tmpVec);
    diffVec = tmpVec - paramMeans;
    std::cout << "12345 " << scalarProduct(diffVec, *covMatrixInverse * diffVec)
              << std::endl;
  }

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************
  delete covMatrixInverse;
  delete covMatrix;

  if (env.fullRank() == 0) {
    std::cout << "Finishing run of 'exStatisticalInverseProblem1_example'"
              << std::endl;
  }

  return;
}
#endif // __EX_STATISTICAL_INVERSE_PROBLEM_1_APPL_H__
