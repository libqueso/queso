#ifndef __UQ_NORMAL_EX_H__
#define __UQ_NORMAL_EX_H__

#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>
#include <uqCovCond.h>

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_MisfitLikelihoodRoutine_DataType
{
  const V* paramInitials;
  const M* matrix;
  bool     applyMatrixInvert;
};

template<class V, class M>
void calib_MisfitLikelihoodRoutine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const V& paramInitials     = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->paramInitials;
  const M& matrix            = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->matrix;
  bool     applyMatrixInvert =  ((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->applyMatrixInvert;

  V diffVec(paramValues - paramInitials);
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
template<class V, class M>
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
  uqParamSpaceClass<V,M>      calib_ParamSpace     (env,"calib");
  uqObservableSpaceClass<V,M> calib_ObservableSpace(env,"calib");

  //******************************************************
  // Step 2 of 6: Define the prior prob. density function object: -2*ln[prior]
  //******************************************************
  uqDefault_M2lPriorRoutine_DataType<V,M> calib_M2lPriorRoutine_Data; // use default prior() routine
  V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uq_M2lProbDensity_Class<V,M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                             (void *) &calib_M2lPriorRoutine_Data); 
  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************
  double condNumber = 100.0;
  V* direction = calib_ParamSpace.newVector();
  direction->cwSet(1.);
  M* precMatrix = calib_ParamSpace.newMatrix();
  M* covMatrix  = calib_ParamSpace.newMatrix();
  uqCovCond(condNumber,*direction,*covMatrix,*precMatrix);
  delete direction;

  calib_MisfitLikelihoodRoutine_DataType<V,M> calib_MisfitLikelihoodRoutine_Data;
  V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_MisfitLikelihoodRoutine_Data.paramInitials     = &calib_ParamInitials;
  calib_MisfitLikelihoodRoutine_Data.matrix            = precMatrix;
  calib_MisfitLikelihoodRoutine_Data.applyMatrixInvert = false;

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
  M proposalCovMatrix(*covMatrix);
  proposalCovMatrix *= (2.4*2.4/(double) calib_ParamSpace.dim());

  //******************************************************
  // Step 6 of 6: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //              In MATLAB, one should then run uqNormalEx.m,
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
