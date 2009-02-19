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

#ifndef __UQ_MARKOV_CHAIN_1_APPL_H__
#define __UQ_MARKOV_CHAIN_1_APPL_H__

#include <uqMarkovChain1Likelihood.h>
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
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqMarkovChain1Example'\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //******************************************************
  // Read Ascii file with information on parameters
  //******************************************************
  uqAsciiTableClass<P_V,P_M> paramsTable(env,
                                         4,    // # of rows
                                         3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                         NULL, // All extra columns are of 'double' type
                                         "inputData/params.tab");

  const EpetraExt::DistArray<std::string>& paramNames = paramsTable.stringColumn(0);
  P_V                                      paramMins    (paramsTable.doubleColumn(1));
  P_V                                      paramMaxs    (paramsTable.doubleColumn(2));
  P_V                                      paramInitials(paramsTable.doubleColumn(3));

  uqVectorSpaceClass<P_V,P_M> paramSpace(env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         paramsTable.numRows(),
                                         &paramNames);

  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",
                                        paramSpace,
                                        paramMins,
                                        paramMaxs);

  //******************************************************
  // Instantiate a likelihood function object (data + routine), to be used by QUESO.
  //******************************************************
  double condNumber = 100.0;
  P_V direction(paramSpace.zeroVector());
  direction.cwSet(1.);
  P_M* covMatrixInverse = paramSpace.newMatrix();
  P_M* covMatrix        = paramSpace.newMatrix();
  uqCovCond(condNumber,direction,*covMatrix,*covMatrixInverse);

  std::cout << "covMatrix = [" << *covMatrix
            << "];"
            << std::endl;

  std::cout << "covMatrixInverse = [" << *covMatrixInverse
            << "];"
            << std::endl;

  likelihoodRoutine_DataType<P_V,P_M> likelihoodRoutine_Data;
  likelihoodRoutine_Data.paramMeans        = &paramInitials;
  likelihoodRoutine_Data.matrix            = covMatrixInverse;
  likelihoodRoutine_Data.applyMatrixInvert = false;

  uqGenericScalarFunctionClass<P_V,P_M> likelihoodFunctionObj("like_",
                                                              paramDomain,
                                                              likelihoodRoutine<P_V,P_M>,
                                                              (void *) &likelihoodRoutine_Data,
                                                              true); // the routine computes [-2.*ln(function)]

  //******************************************************
  // Instantiate inverse problem
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
  // Solve inverse problem, that is, set 'pdf' and 'realizer' of 'postRv'
  //******************************************************
  P_M* proposalCovMatrix = postRv.imageSet().vectorSpace().newGaussianMatrix(priorRv.pdf().domainVarVector(),
                                                                             paramInitials);
  ip.solveWithBayesMarkovChain(paramInitials,
                               proposalCovMatrix);
  delete proposalCovMatrix;

  P_V tmpVec (paramSpace.zeroVector());
  P_V meanVec(paramInitials);
  P_V diffVec(paramSpace.zeroVector());
  const uqBaseVectorRealizerClass<P_V,P_M>& postRealizer = postRv.realizer();
  for (unsigned int i = 0; i < postRealizer.period(); ++i) {
    postRealizer.realization(tmpVec);
    diffVec = tmpVec - meanVec;
    std::cout << "12345 " << scalarProduct(diffVec, *covMatrixInverse * diffVec)
              << std::endl;
  }

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************
  delete covMatrixInverse;
  delete covMatrix;

#if 0
  //******************************************************
  // Step 1 of 6: Define the finite dimensional linear spaces.
  //******************************************************
  uqVectorSpaceClass<P_V,P_M> paramSpace     (env,"calib_",4,NULL);

  P_V paramMins    (paramSpace.zeroVector());
  P_V paramMaxs    (paramSpace.zeroVector());
  P_V paramInitials(paramSpace.zeroVector());
  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",
                                        paramSpace,
                                        paramMins,
                                        paramMaxs);

  //******************************************************
  // Step 2 of 6: Define the prior prob. density function object: -2*ln[prior]
  //******************************************************
  uqDefault_M2lPriorRoutine_DataType<P_V,P_M> calib_M2lPriorRoutine_Data; // use default prior() routine
  P_V calib_ParamPriorMus   (paramSpace.priorMuValues   ());
  P_V calib_ParamPriorSigmas(paramSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uqM2lProbDensity_Class<P_V,P_M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                (void *) &calib_M2lPriorRoutine_Data); 
  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************
  double condNumber = 100.0;
  P_V direction(paramSpace.zeroVector());
  direction.cwSet(1.);
  P_M* covMatrixInverse = paramSpace.newMatrix();
  P_M* covMatrix        = paramSpace.newMatrix();
  uqCovCond(condNumber,direction,*covMatrix,*covMatrixInverse);

  likelihoodRoutine_DataType<P_V,P_M> likelihoodRoutine_Data;
  P_V calib_ParamInitials(paramSpace.initialValues());
  likelihoodRoutine_Data.paramMeans        = &calib_ParamInitials;
  likelihoodRoutine_Data.matrix            = covMatrixInverse;
  likelihoodRoutine_Data.applyMatrixInvert = false;

  uqMisfitLikelihoodFunction_Class<P_V,P_M likelihoodFunction_Obj(likelihoodRoutine<P_V,P_M>,
                                                                  (void *) &likelihoodRoutine_Data);

  uqGenericScalarFunctionClass<P_V,P_M> likelihoodFunctionObj("like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &likelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(function)]
  //******************************************************
  // Step 4 of 6: Define the Markov chain generator.
  //******************************************************
  uqDRAM_MarkovChainGeneratorClass<P_V,P_M> mcg(env,
                                                "calib_",
                                                paramSpace,
                                                calib_M2lPriorProbDensity_Obj,
                                                likelihoodFunction_Obj);

  //******************************************************
  // Step 5 of 6: Compute the proposal covariance matrix.
  //******************************************************
  P_M proposalCovMatrix(*covMatrix);
  proposalCovMatrix *= (2.4*2.4/(double) paramSpace.dim());

  //******************************************************
  // Step 6 of 6: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //******************************************************
  mcg.generateChains(&proposalCovMatrix,
                     covMatrix,
                     true);

#endif
  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqMarkovChain1Example'"
              << std::endl;
  }

  return;
}
#endif // __UQ_MARKOV_CHAIN_1_APPL_H__
