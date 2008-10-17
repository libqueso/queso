/* uq/examples/mcmc/normal/uqNormalEx.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_NORMAL_EX_H__
#define __UQ_NORMAL_EX_H__

#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>
#include <uqCovCond.h>

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class P_V,class P_M>
struct
calib_MisfitLikelihoodRoutine_DataType
{
  const P_V* paramInitials;
  const P_M* matrix;
  bool       applyMatrixInvert;
};

template<class P_V,class P_M,class L_V,class L_M>
void calib_MisfitLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr, L_V& resultValues)
{
  const P_V& paramInitials     = *((calib_MisfitLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->paramInitials;
  const P_M& matrix            = *((calib_MisfitLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->matrix;
  bool       applyMatrixInvert =  ((calib_MisfitLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->applyMatrixInvert;

  P_V diffVec(paramValues - paramInitials);
  if (applyMatrixInvert) {
    resultValues[0] = scalarProduct(diffVec, matrix.invertMultiply(diffVec));
  }
  else {
    resultValues[0] = scalarProduct(diffVec, matrix * diffVec);
  }

  return;
}

//********************************************************
// The MCMC driving routine: called by main()
//********************************************************
template<class P_V,class P_M,class L_V,class L_M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqNormalEx' example\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //******************************************************
  // Step 1 of 6: Define the finite dimensional linear spaces.
  //******************************************************
  uqParamSpaceClass<P_V,P_M>      calib_ParamSpace     (env,"calib_");
  uqObservableSpaceClass<L_V,L_M> calib_ObservableSpace(env,"calib_");

  //******************************************************
  // Step 2 of 6: Define the prior prob. density function object: -2*ln[prior]
  //******************************************************
  uqDefault_M2lPriorRoutine_DataType<P_V,P_M> calib_M2lPriorRoutine_Data; // use default prior() routine
  P_V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  P_V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uqM2lProbDensity_Class<P_V,P_M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                (void *) &calib_M2lPriorRoutine_Data); 
  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************
  double condNumber = 100.0;
  P_V direction(calib_ParamSpace.zeroVector());
  direction.cwSet(1.);
  P_M* precMatrix = calib_ParamSpace.newMatrix();
  P_M* covMatrix  = calib_ParamSpace.newMatrix();
  uqCovCond(condNumber,direction,*covMatrix,*precMatrix);

  calib_MisfitLikelihoodRoutine_DataType<P_V,P_M> calib_MisfitLikelihoodRoutine_Data;
  P_V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_MisfitLikelihoodRoutine_Data.paramInitials     = &calib_ParamInitials;
  calib_MisfitLikelihoodRoutine_Data.matrix            = precMatrix;
  calib_MisfitLikelihoodRoutine_Data.applyMatrixInvert = false;

  uqMisfitLikelihoodFunction_Class<P_V,P_M,L_V,L_M> calib_MisfitLikelihoodFunction_Obj(calib_MisfitLikelihoodRoutine<P_V,P_M,L_V,L_M>,
                                                                                       (void *) &calib_MisfitLikelihoodRoutine_Data);

  //******************************************************
  // Step 4 of 6: Define the Markov chain generator.
  //******************************************************
  uqDRAM_MarkovChainGeneratorClass<P_V,P_M,L_V,L_M> mcg(env,
                                                        "calib_",
                                                        calib_ParamSpace,
                                                        calib_ObservableSpace,
                                                        calib_M2lPriorProbDensity_Obj,
                                                        calib_MisfitLikelihoodFunction_Obj);

  //******************************************************
  // Step 5 of 6: Compute the proposal covariance matrix.
  //******************************************************
  P_M proposalCovMatrix(*covMatrix);
  proposalCovMatrix *= (2.4*2.4/(double) calib_ParamSpace.dim());

  //******************************************************
  // Step 6 of 6: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //******************************************************
  mcg.generateChains(&proposalCovMatrix,
                     covMatrix,
                     true);

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************
  delete precMatrix;
  delete covMatrix;

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqNormalEx' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_NORMAL_EX_H__
