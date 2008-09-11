/* uq/examples/mcmc/etc/uqEtcEx.h
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_ETC_EX_H__
#define __UQ_ETC_EX_H__

//#include <uqStateSpace.h>
#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> The application can use the default prior() routine provided by the MCMC tool.
//********************************************************
template<class P_V,class P_M>
struct
calib_M2lPriorRoutine_DataType
{
  const P_V* maybeYouNeedAVector;
  const P_M* maybeYouNeedAMatrix;
};

template<class P_V, class P_M>
double calib_M2lPriorRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  const P_V& v1 = *((calib_M2lPriorRoutine_DataType<P_V,P_M> *) functionDataPtr)->maybeYouNeedAVector;
  const P_M& m1 = *((calib_M2lPriorRoutine_DataType<P_V,P_M> *) functionDataPtr)->maybeYouNeedAMatrix;

  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_M2lPriorRoutine(), in uqEtcEx.h",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");

  return (m1 * v1).sumOfComponents();
}

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class S_V,class S_M>
struct
calib_CompleteLikelihoodRoutine_DataType
{
  const S_V* aVectorForInstance;
  const S_M* aMatrixForInstance;
};

template<class P_V,class P_M,class S_V,class S_M,class L_V,class L_M>
void calib_CompleteLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr, L_V& resultValues)
{
  const S_V& v1 = *((calib_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->aVectorForInstance;
  const S_M& m1 = *((calib_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->aMatrixForInstance;

  const S_V& v2 = *((calib_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->aVectorForInstance;
  const S_M& m2 = *((calib_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->aMatrixForInstance;

  L_V misfit(resultValues);
  misfit[0] = (m1 * v1).norm2Sq();
  misfit[1] = (m2 * v2).norm2Sq();

  L_V sigma(resultValues);
  sigma[0] = 1.;
  sigma[1] = .34;

  L_V minus2LogLikelihood(resultValues);
  minus2LogLikelihood = misfit/sigma;

  // The MCMC library will always expect the value [-2. * ln(Likelihood)], so this routine should use this line below
  resultValues = minus2LogLikelihood;

  // If, however, you really need to return the Likelihood(s), then use this line below
  //resultsValues[0] = exp(-.5*minus2LogLikelihood[0]);
  //resultsValues[1] = exp(-.5*minus2LogLikelihood[1]);

  return;
}

//********************************************************
// The MCMC driving routine: called by main()
//********************************************************
template<class P_V,class P_M,class S_V,class S_M,class L_V,class L_M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqEtcEx' example\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //*****************************************************
  // Step 1 of 6: Define the finite dimensional linear spaces.
  //*****************************************************
  uqParamSpaceClass<P_V,P_M>      calib_ParamSpace     (env,"calib_");
  uqObservableSpaceClass<L_V,L_M> calib_ObservableSpace(env,"calib_");

  //******************************************************
  // Step 2 of 6: Define the prior prob. density function object: -2*ln[prior]
  //******************************************************
#if 1 // You might want to use the default prior function provided by the MCMC tool
  uqDefault_M2lPriorRoutine_DataType<P_V,P_M> calib_M2lPriorRoutine_Data; // use default prior() routine
  P_V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  P_V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uq_M2lProbDensity_Class<P_V,P_M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                 (void *) &calib_M2lPriorRoutine_Data);
#else // Or you might need to use your own prior function
  calib_M2lPriorRoutine_DataType<P_V,P_M> calib_M2lPriorRoutine_Data;
  P_V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.maybeYouNeedAVector = &calib_ParamPriorSigmas;
  calib_M2lPriorRoutine_Data.maybeYouNeedAMatrix = NULL;

  uq_M2lProbDensity_Class<P_V,P_M> calib_M2lPriorProbDensity_Obj(calib_M2lPriorRoutine<P_V,P_M>,
                                                                 (void *) &calib_M2lPriorRoutine_Data);
#endif

  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************
  calib_CompleteLikelihoodRoutine_DataType<S_V,S_M> calib_CompleteLikelihoodRoutine_Data;
  P_V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_CompleteLikelihoodRoutine_Data.aVectorForInstance = &calib_ParamInitials;
  calib_CompleteLikelihoodRoutine_Data.aMatrixForInstance = NULL;

  uq_CompleteLikelihoodFunction_Class<P_V,P_M,L_V,L_M> calib_CompleteLikelihoodFunction_Obj(calib_CompleteLikelihoodRoutine<P_V,P_M,S_V,S_M,L_V,L_M>,
                                                                                            (void *) &calib_CompleteLikelihoodRoutine_Data,
                                                                                            true); // the routine computes [-2.*ln(Likelihood)]

  //******************************************************
  // Step 4 of 6: Define the Markov chain generator.
  //******************************************************
  uqDRAM_MarkovChainGeneratorClass<P_V,P_M,L_V,L_M> mcg(env,
                                                        "calib_",
                                                        calib_ParamSpace,
                                                        calib_ObservableSpace,
                                                        calib_M2lPriorProbDensity_Obj,
                                                        calib_CompleteLikelihoodFunction_Obj);

  //******************************************************
  // Step 5 of 6: Compute the proposal covariance matrix.
  //******************************************************
  // Proposal covariance matrix will be computed internally
  // by the Markov chain generator.
  P_M* proposalCovMatrix = NULL;

  //******************************************************
  // Step 6 of 6: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //******************************************************
  mcg.generateChains(proposalCovMatrix);

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqEtcEx' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_ETC_EX_H__
