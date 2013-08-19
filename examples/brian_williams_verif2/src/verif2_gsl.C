#include <verif2_gsl.h>
#include <uqStatisticalInverseProblem.h>

static unsigned int likelihoodCounter = 0;

int main(int argc, char* argv[])
{
  //***********************************************************************
  // Initialize MPI
  //***********************************************************************
  MPI_Init(&argc,&argv);

  //***********************************************************************
  // Initialize QUESO environment
  //***********************************************************************
  UQ_FATAL_TEST_MACRO((argc < 2),
                      UQ_UNAVAILABLE_RANK,
                      "main()",
                      "run as <executable> 'inputFileName'");
  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  //***********************************************************************
  // Run program
  //***********************************************************************
  solveSip(*env);

  //***********************************************************************
  // Finalize QUESO environment
  //***********************************************************************
  delete env;

  //***********************************************************************
  // Finalize MPI
  //***********************************************************************
  MPI_Finalize();

  return 0;
}

void solveSip(const uqFullEnvironmentClass& env)
{
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Entering solveSip()..."
                          << std::endl;
  }

  ////////////////////////////////////////////////////////
  // Step 1 of 5: Instantiate the parameter space
  ////////////////////////////////////////////////////////
  unsigned int p = 2;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paramSpace(env, "param_", p, NULL);
  uqGslVectorClass aVec(paramSpace.zeroVector());
  aVec[0] = 2.;
  aVec[1] = 5.;
  uqGslVectorClass xGiven(paramSpace.zeroVector());
  xGiven[0] = -1.;
  xGiven[1] =  7.;

  ////////////////////////////////////////////////////////
  // Step 2 of 5: Instantiate the parameter domain
  ////////////////////////////////////////////////////////
  //uqGslVectorClass paramMins    (paramSpace.zeroVector());
  //uqGslVectorClass paramMaxs    (paramSpace.zeroVector());
  //paramMins    [0] = -1.e+16;
  //paramMaxs    [0] =  1.e+16;
  //paramMins    [1] = -1.e+16;
  //paramMaxs    [1] =  1.e+16;
  //uqBoxSubsetClass<uqGslVectorClass,uqGslMatrixClass> paramDomain("param_",paramSpace,paramMins,paramMaxs);
  uqVectorSetClass<uqGslVectorClass,uqGslMatrixClass>* paramDomain = &paramSpace;

  ////////////////////////////////////////////////////////
  // Step 3 of 5: Instantiate the likelihood function object
  ////////////////////////////////////////////////////////
  unsigned int n = 5;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> dataSpace(env, "data_", n, NULL);
  
  uqGslVectorClass yMeanVec(dataSpace.zeroVector());
  double tmp = scalarProduct(aVec,xGiven);
  for (unsigned int i = 0; i < n; ++i) {
    yMeanVec[i] = tmp;
  }

  double sigmaEps = 2.1;
  uqGslMatrixClass yCovMat(dataSpace.zeroVector());
  tmp = sigmaEps*sigmaEps;
  for (unsigned int i = 0; i < n; ++i) {
    yCovMat(i,i) = tmp;
  }

  uqGslVectorClass ySamples(dataSpace.zeroVector());
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> yRv("y_", dataSpace, yMeanVec, yCovMat);
  yRv.realizer().realization(ySamples);

  double ySampleMean = 0.;
  for (unsigned int i = 0; i < n; ++i) {
    ySampleMean += ySamples[i];
  }
  ySampleMean /= ((double) n);

  struct likelihoodDataStruct likelihoodData;
  likelihoodData.aVec     = &aVec;
  likelihoodData.sigmaEps = sigmaEps;
  likelihoodData.ySamples = &ySamples;

  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>
    likelihoodFunctionObj("like_",
                          *paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodData,
                          true); // routine computes [ln(function)]

  ////////////////////////////////////////////////////////
  // Step 4 of 5: Instantiate the inverse problem
  ////////////////////////////////////////////////////////
  uqGslVectorClass xPriorMeanVec(paramSpace.zeroVector());
  xPriorMeanVec[0] = 0.;
  xPriorMeanVec[1] = 0.;
  uqGslMatrixClass sigma0Mat(paramSpace.zeroVector());
  sigma0Mat(0,0) = 1.e-3;
  sigma0Mat(0,1) = 0.;
  sigma0Mat(1,0) = 0.;
  sigma0Mat(1,1) = 1.e-3;
  uqGslMatrixClass sigma0MatInverse(paramSpace.zeroVector());
  sigma0MatInverse = sigma0Mat.inverse();
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> priorRv("prior_", *paramDomain, xPriorMeanVec, sigma0MatInverse);

  uqGenericVectorRVClass <uqGslVectorClass,uqGslMatrixClass> postRv ("post_", paramSpace);

  uqStatisticalInverseProblemClass<uqGslVectorClass,uqGslMatrixClass> sip("sip_", NULL, priorRv, likelihoodFunctionObj, postRv);

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSip():"
                          << "\n  p                = " << p
                          << "\n  xGiven           = " << xGiven
                          << "\n  sigma0Mat        = " << sigma0Mat
                          << "\n  sigma0MatInverse = " << sigma0MatInverse
                          << "\n  aVec             = " << aVec
                          << "\n  n                = " << n
                          << "\n  sigmaEps         = " << sigmaEps
                          << "\n  yMeanVec         = " << yMeanVec
                          << "\n  yCovMat          = " << yCovMat
                          << "\n  ySamples         = " << ySamples
                          << "\n  ySampleMean      = " << ySampleMean
                          << std::endl;
  }

  uqGslMatrixClass sigmaMatInverse(paramSpace.zeroVector());
  sigmaMatInverse = matrixProduct(aVec,aVec);
  sigmaMatInverse *= (((double) n)/sigmaEps/sigmaEps);
  sigmaMatInverse += sigma0Mat;
  uqGslMatrixClass sigmaMat(paramSpace.zeroVector());
  sigmaMat = sigmaMatInverse.inverse();

  uqGslVectorClass muVec(paramSpace.zeroVector());
  muVec = sigmaMat * aVec;
  muVec *= (((double) n) * ySampleMean)/sigmaEps/sigmaEps;

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSip():"
                          << "\n  muVec            = " << muVec
                          << "\n  sigmaMat         = " << sigmaMat
                          << "\n  sigmaMatInverse  = " << sigmaMatInverse
                          << std::endl;
  }

  ////////////////////////////////////////////////////////
  // Step 5 of 5: Solve the inverse problem
  ////////////////////////////////////////////////////////
  uqGslVectorClass initialValues(paramSpace.zeroVector());
  initialValues[0] = 25.;
  initialValues[1] = 25.;

  uqGslMatrixClass proposalCovMat(paramSpace.zeroVector());
  proposalCovMat(0,0) = 10.;
  proposalCovMat(0,1) = 0.;
  proposalCovMat(1,0) = 0.;
  proposalCovMat(1,1) = 10.;

  sip.solveWithBayesMetropolisHastings(NULL,initialValues,&proposalCovMat);

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Leaving solveSip()"
                          << std::endl;
  }

  return;
}

double likelihoodRoutine(
  const uqGslVectorClass& paramValues,
  const uqGslVectorClass* paramDirection,
  const void*             functionDataPtr,
  uqGslVectorClass*       gradVector,
  uqGslMatrixClass*       hessianMatrix,
  uqGslVectorClass*       hessianEffect)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  likelihoodCounter++;
  const uqBaseEnvironmentClass& env = paramValues.env();
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Entering likelihoodRoutine()..."
                          << ": likelihoodCounter = "       << likelihoodCounter
                          << ", params = "                  << paramValues
                          << ", env.subComm().NumProc() = " << env.subComm().NumProc()
                          << ", my subRank = "              << env.subRank()
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Entering likelihoodRoutine()"
              << ", likelihoodCounter = " << likelihoodCounter
              << std::endl;
#endif
  }

  //////////////////////////////////////////////////
  // Begin actual likelihood routine
  //////////////////////////////////////////////////
  double totalLnLikelihood = 0.;

  if (paramDirection  &&
      functionDataPtr &&
      gradVector      &&
      hessianMatrix   &&
      hessianEffect) {
    // Just to eliminate INTEL compiler warnings
  }

  struct likelihoodDataStruct* likelihoodData = (likelihoodDataStruct *) functionDataPtr; 
  uqGslVectorClass aVec(*(likelihoodData->aVec));
  unsigned int p = aVec.sizeLocal();
  double sigmaEps = likelihoodData->sigmaEps;
  uqGslVectorClass ySamples(*(likelihoodData->ySamples));
  unsigned int n = ySamples.sizeLocal();

  UQ_FATAL_TEST_MACRO(paramValues.sizeLocal() != p,
                      env.fullRank(),
                      "likelihoodRoutine()",
                      "invalid parameter vector size");

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
    *env.subDisplayFile() << "In likelihoodRoutine()"
                          << ": likelihoodCounter = " << likelihoodCounter
                          << ", params = "            << paramValues
                          << ", p = "                 << p
                          << ", aVec = "              << aVec
                          << ", sigmaEps = "          << sigmaEps
                          << ", n = "                 << n
                          << ", ySamples = "          << ySamples
                          << std::endl;
  }

  //******************************************************************************
  // Compute likelihood
  //******************************************************************************
  for (unsigned int i = 0; i < n; ++i) {
    double diff = (ySamples[i] - scalarProduct(aVec,paramValues))/sigmaEps;
    totalLnLikelihood -= 0.5 * diff * diff;
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
      *env.subDisplayFile() << "In likelihoodRoutine()"
                            << ": likelihoodCounter = " << likelihoodCounter
                            << ", params = "            << paramValues
                            << ", diff = "              << diff
                            << std::endl;
    }
  }

  //******************************************************************************
  // Prepare to return
  //******************************************************************************
  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Leaving likelihoodRoutine()"
                          << ": likelihoodCounter = " << likelihoodCounter
                          << ", params = "            << paramValues
                          << ", totalLnLikelihood = " << totalLnLikelihood
                          << " after "                << totalTime
                          << " seconds"
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Leaving likelihoodRoutine()"
              << ": likelihoodCounter = " << likelihoodCounter
              << ", params = "            << paramValues
              << ", totalLnLikelihood = " << totalLnLikelihood
              << " after "                << totalTime
              << " seconds"
              << std::endl;
#endif
  }

  env.subComm().Barrier();
  //exit(1);

  return totalLnLikelihood;
}
