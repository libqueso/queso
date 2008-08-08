#ifndef __UQ_NORMAL_EX_H__
#define __UQ_NORMAL_EX_H__

#include <uqDRAM_MarkovChainGenerator.h>
#include <uqCovCond.h>

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> This application is using the default prior() routine provided by the MCMC tool.
//********************************************************
template<class V, class M>
double calib_PriorFunction_Routine(const V& paramValues, const void* functionDataPtr)
{
  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_PriorFunction_Routine()",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");
  return 0.;
}

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_LikelihoodFunction_DataType
{
  const V* paramInitials;
  const M* matrix;
  bool     applyMatrixInvert;
};

template<class V, class M>
void calib_LikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const V& paramInitials     = *((calib_LikelihoodFunction_DataType<V,M> *) functionDataPtr)->paramInitials;
  const M& matrix            = *((calib_LikelihoodFunction_DataType<V,M> *) functionDataPtr)->matrix;
  bool     applyMatrixInvert =  ((calib_LikelihoodFunction_DataType<V,M> *) functionDataPtr)->applyMatrixInvert;

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
  // Step 1 of 4: Define the finite dimensional linear spaces.
  //              Define the -2*ln[prior] and -2*ln[likelihood] function objects
  //              Define the Markov chain generator.
  //******************************************************
  uqParamSpaceClass <V,M> calib_ParamSpace (env,"calib");
  uqOutputSpaceClass<V,M> calib_OutputSpace(env,"calib");

  tmp_M2lPriorFunction_Class<V,M>      calib_PriorFunction_Obj     (uqDefault_M2lPriorFunction_Routine<V,M>); // use default prior() routine
  tmp_M2lLikelihoodFunction_Class<V,M> calib_LikelihoodFunction_Obj(calib_LikelihoodFunction_Routine<V,M>);
  uqDRAM_MarkovChainGeneratorClass<V,M> mcg(env,
                                            "calib",
                                            calib_ParamSpace,
                                            calib_OutputSpace,
                                            calib_PriorFunction_Obj,
                                            calib_LikelihoodFunction_Obj);

  //******************************************************
  // Step 2 of 4: Compute the proposal covariance matrix.
  //******************************************************
  double condNumber = 100.0;
  V* direction = calib_ParamSpace.newVector();
  direction->cwSet(1.);
  M* precMatrix = calib_ParamSpace.newMatrix();
  M* covMatrix  = calib_ParamSpace.newMatrix();
  uqCovCond(condNumber,*direction,*covMatrix,*precMatrix);
  delete direction;

  M proposalCovMatrix(*covMatrix);
  proposalCovMatrix *= (2.4*2.4/(double) calib_ParamSpace.dim());

  //******************************************************
  // Step 3 of 4: Prepare the data to be passed to
  //              calib_PriorFunction_Routine() and
  //              calib_LikelihoodFunction_Routine().
  //******************************************************
  uqDefault_M2lPriorFunction_DataType<V,M> calib_PriorFunction_Data; // use default prior() routine
  V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_PriorFunction_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_PriorFunction_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  calib_LikelihoodFunction_DataType<V,M> calib_LikelihoodFunction_Data;
  V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_LikelihoodFunction_Data.paramInitials     = &calib_ParamInitials;
  calib_LikelihoodFunction_Data.matrix            = precMatrix;
  calib_LikelihoodFunction_Data.applyMatrixInvert = false;

  //******************************************************
  // Step 4 of 4: Generate chains.
  //              In MATLAB, one should then run uqNormalEx.m,
  //              which calls 'mcmcplot.m'.
  //******************************************************
  mcg.generateChains(&proposalCovMatrix,
                     (void *) &calib_PriorFunction_Data, // use default prior() routine
                     (void *) &calib_LikelihoodFunction_Data,
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
