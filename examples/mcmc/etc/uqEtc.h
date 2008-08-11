/* uq/examples/mcmc/etc/uqEtcEx.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
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

#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> The application can use the default prior() routine provided by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_M2lPriorRoutine_DataType
{
  const V* maybeYouNeedAVector;
  const M* maybeYouNeedAMatrix;
};

template<class V, class M>
double calib_M2lPriorRoutine(const V& paramValues, const void* functionDataPtr)
{
  const V& v1 = *((calib_M2lPriorRoutine_DataType<V,M> *) functionDataPtr)->maybeYouNeedAVector;
  const M& m1 = *((calib_M2lPriorRoutine_DataType<V,M> *) functionDataPtr)->maybeYouNeedAMatrix;

  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_M2lPriorRoutine(), in uqEtcEx.h",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");

  return (m1 * v1).sumOfComponents();
}

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_MisfitLikelihoodRoutine_DataType
{
  const V* aVectorForInstance;
  const M* aMatrixForInstance;
};

template<class V, class M>
void calib_MisfitLikelihoodRoutine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const V& v1 = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->aVectorForInstance;
  const M& m1 = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->aMatrixForInstance;

  resultValues[0] = (m1 * v1).norm2Sq();

  return;
}

//********************************************************
// The MCMC driving routine: called by main()
//********************************************************
template<class V, class M>
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
  uqParamSpaceClass<V,M>      calib_ParamSpace     (env,"calib");
  uqObservableSpaceClass<V,M> calib_ObservableSpace(env,"calib");

  //******************************************************
  // Step 2 of 6: Define the prior prob. density function object: -2*ln[prior]
  //******************************************************
#if 1 // You might need to use your own prior function
  calib_M2lPriorRoutine_DataType<V,M> calib_M2lPriorRoutine_Data;
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.maybeYouNeedAVector = &calib_ParamPriorSigmas;
  calib_M2lPriorRoutine_Data.maybeYouNeedAMatrix = NULL;

  uq_M2lProbDensity_Class<V,M> calib_M2lPriorProbDensity_Obj(calib_M2lPriorRoutine<V,M>,
                                                             (void *) &calib_M2lPriorRoutine_Data);
#else // Or you might want to use the default prior function provided by the MCMC tool
  uqDefault_M2lPriorRoutine_DataType<V,M> calib_M2lPriorRoutine_Data; // use default prior() routine
  V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uq_M2lProbDensity_Class<V,M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                             (void *) &calib_M2lPriorRoutine_Data);
#endif

  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************
  calib_MisfitLikelihoodRoutine_DataType<V,M> calib_MisfitLikelihoodRoutine_Data;
  V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_MisfitLikelihoodRoutine_Data.aVectorForInstance = &calib_ParamInitials;
  calib_MisfitLikelihoodRoutine_Data.aMatrixForInstance = NULL;

  uq_MisfitLikelihoodFunction_Class<V,M> calib_MisfitLikelihoodFunction_Obj(calib_MisfitLikelihoodRoutine<V,M>,
                                                                            (void *) &calib_MisfitLikelihoodRoutine_Data);

  //******************************************************
  // Step 4 of 6: Define the Markov chain generator.
  //******************************************************
  uqDRAM_MarkovChainGeneratorClass<V,M> mcg(env,
                                            "calib",
                                            calib_ParamSpace,
                                            calib_ObservableSpace,
                                            calib_M2lPriorProbDensity_Obj,
                                            calib_MisfitLikelihoodFunction_Obj);

  //******************************************************
  // Step 5 of 6: Compute the proposal covariance matrix.
  //******************************************************
  // Proposal covariance matrix will be computed internally
  // by the Markov chain generator.
  M* proposalCovMatrix = NULL;

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
