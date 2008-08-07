#ifndef __UQ_NORMAL_EX_H__
#define __UQ_NORMAL_EX_H__

#include <uqDRAM_MarkovChainGenerator.h>
#include <uqCovCond.h>

template<class V, class M>
struct
uqAppl_M2lLikelihoodFunction_DataType
{
  const V* paramInitials;
  const M* matrix;
  bool     applyMatrixInvert;
};

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
  //              Define the Markov chain generator.
  //******************************************************
  uqParamSpaceClass <V,M> calParamSpace (env,"cal");
  uqOutputSpaceClass<V,M> calOutputSpace(env,"cal");

  uq_M2lLikelihoodFunction_Class<V,M> uq_M2lLikelihoodFunction_Obj;
  uqDRAM_MarkovChainGeneratorClass<V,M> mcg(env,
                                            "cal",
                                            calParamSpace,
                                            calOutputSpace,
                                            NULL, // use default prior() routine
                                            uq_M2lLikelihoodFunction_Obj);

  //******************************************************
  // Step 2 of 4: Compute the proposal covariance matrix.
  //******************************************************
  double condNumber = 100.0;
  V* direction = calParamSpace.newVector();
  direction->cwSet(1.);
  M* precMatrix = calParamSpace.newMatrix();
  M* covMatrix  = calParamSpace.newMatrix();
  uqCovCond(condNumber,*direction,*covMatrix,*precMatrix);
  delete direction;

  M proposalCovMatrix(*covMatrix);
  proposalCovMatrix *= (2.4*2.4/(double) calParamSpace.dim());

  //******************************************************
  // Step 3 of 4: Prepare the data to be passed to
  //              uqAppl_M2lLikelihoodFunction_Routine().
  //******************************************************
  V paramInitials(calParamSpace.initialValues());

  uqAppl_M2lLikelihoodFunction_DataType<V,M> uqAppl_M2lLikelihoodFunction_Data;
  uqAppl_M2lLikelihoodFunction_Data.paramInitials     = &paramInitials;
  uqAppl_M2lLikelihoodFunction_Data.matrix            = precMatrix;
  uqAppl_M2lLikelihoodFunction_Data.applyMatrixInvert = false;

  //******************************************************
  // Step 4 of 4: Generate chains.
  //              In MATLAB, one should then run uqNormalEx.m,
  //              which calls 'mcmcplot.m'.
  //******************************************************
  mcg.generateChains(&proposalCovMatrix,
                     NULL, // use default prior() routine
                     (void *) &uqAppl_M2lLikelihoodFunction_Data,
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

template<class V, class M>
double uqAppl_M2lPriorFunction_Routine(const V& paramValues, const void* functionDataPtr)
{
  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "uqAppl_M2lPriorFunction_Routine()",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");
  return 0.;
}

template<class V, class M>
void uqAppl_M2lLikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const V& paramInitials     = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->paramInitials;
  const M& matrix            = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->matrix;
  bool     applyMatrixInvert =  ((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->applyMatrixInvert;

  V diffVec(paramValues - paramInitials);
  if (applyMatrixInvert) {
    resultValues[0] = scalarProduct(diffVec, matrix.invertMultiply(diffVec));
  }
  else {
    resultValues[0] = scalarProduct(diffVec, matrix * diffVec);
  }

  return;
}
#endif // __UQ_NORMAL_EX_H__
