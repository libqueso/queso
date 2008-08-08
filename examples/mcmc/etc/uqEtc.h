#ifndef __UQ_TEMPLATE_EX_H__
#define __UQ_TEMPLATE_EX_H__

#include <uqDRAM_MarkovChainGenerator.h>

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> The application can use the default prior() routine provided by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_PriorFunction_DataType
{
  const V* maybeYouNeedAVector;
  const M* maybeYouNeedAMatrix;
};

template<class V, class M>
double calib_PriorFunction_Routine(const V& paramValues, const void* functionDataPtr)
{
  const V& v1 = *((calib_PriorFunction_DataType<V,M> *) functionDataPtr)->maybeYouNeedAVector;
  const M& m1 = *((calib_PriorFunction_DataType<V,M> *) functionDataPtr)->maybeYouNeedAMatrix;

  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_PriorFunction_Routine(), in uqEtcEx.h",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");

  return (m1 * v1).sumOfComponents();
}

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_LikelihoodFunction_DataType
{
  const V* aVectorForInstance;
  const M* aMatrixForInstance;
};

template<class V, class M>
void calib_LikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const V& v1 = *((calib_LikelihoodFunction_DataType<V,M> *) functionDataPtr)->aVectorForInstance;
  const M& m1 = *((calib_LikelihoodFunction_DataType<V,M> *) functionDataPtr)->aMatrixForInstance;

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

  //******************************************************
  // Step 1 of 4: Define the finite dimensional linear spaces.
  //              Define the -2*ln[prior] and -2*ln[likelihood] function objects
  //              Define the Markov chain generator.
  //******************************************************
  uqParamSpaceClass <V,M> calib_ParamSpace (env,"calib");
  uqOutputSpaceClass<V,M> calib_OutputSpace(env,"calib");

#if 1 // You might need to use your own prior function
  tmp_M2lPriorFunction_Class<V,M>      calib_PriorFunction_Obj(calib_PriorFunction_Routine<V,M>);
#else // Or you might want to use the default prior function provided by the MCMC tool
  tmp_M2lPriorFunction_Class<V,M>      calib_PriorFunction_Obj(uqDefault_M2lPriorFunction_Routine<V,M>); // use default prior() routine
#endif
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
  M* proposalCovMatrix = NULL;

  //******************************************************
  // Step 3 of 4: Prepare the data to be passed to
  //              calib_PriorFunction_Routine() and
  //              calib_LikelihoodFunction_Routine().
  //******************************************************
#if 1 // You might need to use your own prior function
  calib_PriorFunction_DataType<V,M> calib_PriorFunction_Data;
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_PriorFunction_Data.maybeYouNeedAVector = &calib_ParamPriorSigmas;
  calib_PriorFunction_Data.maybeYouNeedAMatrix = NULL;
#else // Or you might want to use the default prior function provided by the MCMC tool
  uqDefault_M2lPriorFunction_DataType<V,M> calib_PriorFunction_Data; // use default prior() routine
  V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_PriorFunction_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_PriorFunction_Data.paramPriorSigmas = &calib_ParamPriorSigmas;
#endif

  calib_LikelihoodFunction_DataType<V,M> calib_LikelihoodFunction_Data;
  V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_LikelihoodFunction_Data.aVectorForInstance = &calib_ParamInitials;
  calib_LikelihoodFunction_Data.aMatrixForInstance = NULL;

  //******************************************************
  // Step 4 of 4: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //******************************************************
  mcg.generateChains(proposalCovMatrix,
                     (void *) &calib_PriorFunction_Data,
                     (void *) &calib_LikelihoodFunction_Data);

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqEtcEx' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_TEMPLATE_EX_H__
