#include <verif3_gsl.h>
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
  solveSipForX0(*env);

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

void solveSipForX0(const uqFullEnvironmentClass& env)
{
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Entering solveSipForX0()..."
                          << std::endl;
  }

  ////////////////////////////////////////////////////////
  // Step 1 of 5: Instantiate the parameter space
  ////////////////////////////////////////////////////////
  unsigned int p = 2;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paramSpace(env, "param_", p, NULL);
  uqGslVectorClass xGiven(paramSpace.zeroVector());
  xGiven[0] = -1.;
  xGiven[1] =  7.;

  uqGslMatrixClass sigmaMat(paramSpace.zeroVector());
  sigmaMat(0,0) = 36.;
  sigmaMat(0,1) = 0.;
  sigmaMat(1,0) = 0.;
  sigmaMat(1,1) = 4.;

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
  uqGslVectorClass zMeanVec(paramSpace.zeroVector());
  uqGslMatrixClass zCovMat (paramSpace.zeroVector());
  zCovMat(0,0) = 25.;
  zCovMat(0,1) = 0.;
  zCovMat(1,0) = 0.;
  zCovMat(1,1) = 16.;
  uqGslVectorClass zSample (paramSpace.zeroVector());
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> zRv("z_", paramSpace, zMeanVec, zCovMat);

  double sigmaEps = 2.1;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> oneDSpace(env, "oneD_", 1, NULL);
  uqGslVectorClass epsMeanVec(paramSpace.zeroVector());
  uqGslMatrixClass epsCovMat (paramSpace.zeroVector());
  epsCovMat(0,0) = sigmaEps*sigmaEps;
  uqGslVectorClass epsSample (paramSpace.zeroVector());
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> epsRv("eps_", paramSpace, epsMeanVec, epsCovMat);
 
  unsigned int n = 5;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> dataSpace(env, "data_", n, NULL);
  
  uqGslMatrixClass zMat(env,paramSpace.map(),n);

  uqGslVectorClass yMeanVec(dataSpace.zeroVector());
  for (unsigned int i = 0; i < n; ++i) {
    zRv.realizer().realization(zSample);
    epsRv.realizer().realization(epsSample);
    yMeanVec[i] = scalarProduct(zSample,xGiven) + epsSample[0];
  }

  uqGslMatrixClass yCovMat(dataSpace.zeroVector());
  double tmp = sigmaEps*sigmaEps;
  for (unsigned int i = 0; i < n; ++i) {
    yCovMat(i,i) = tmp;
  }

  uqGslVectorClass ySamples(dataSpace.zeroVector());
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> yRv("y_", dataSpace, yMeanVec, yCovMat);
  yRv.realizer().realization(ySamples);

  struct likelihoodDataStructForX0 likelihoodDataForX0;
  likelihoodDataForX0.zMat     = &zMat;
  likelihoodDataForX0.sigmaEps = sigmaEps;
  likelihoodDataForX0.ySamples = &ySamples;
  likelihoodDataForX0.sigmaMat = &sigmaMat;

  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>
    likelihoodFunctionObj("like_",
                          *paramDomain,
                          likelihoodRoutineForX0,
                          (void *) &likelihoodDataForX0,
                          true); // routine computes [ln(function)]

  ////////////////////////////////////////////////////////
  // Step 4 of 5: Instantiate the inverse problem
  ////////////////////////////////////////////////////////
  uqGslVectorClass mu0Vec(paramSpace.zeroVector());
  mu0Vec[0] = -1.;
  mu0Vec[1] =  7.;
  uqGslMatrixClass sigma0Mat(paramSpace.zeroVector());
  sigma0Mat(0,0) = 1.;
  sigma0Mat(0,1) = 0.;
  sigma0Mat(1,0) = 0.;
  sigma0Mat(1,1) = 4;
  uqGslMatrixClass sigma0MatInverse(paramSpace.zeroVector());
  sigma0MatInverse = sigma0Mat.inverse();
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> priorRv("prior_", *paramDomain, mu0Vec, sigma0MatInverse);

  uqGenericVectorRVClass <uqGslVectorClass,uqGslMatrixClass> postRv ("post_", paramSpace);

  uqStatisticalInverseProblemClass<uqGslVectorClass,uqGslMatrixClass> sip("sip_", NULL, priorRv, likelihoodFunctionObj, postRv);

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSipForX0():"
                          << "\n  p                = " << p
                          << "\n  xGiven           = " << xGiven
                          << "\n  sigma0Mat        = " << sigma0Mat
                          << "\n  sigma0MatInverse = " << sigma0MatInverse
                          << "\n  zMat             = " << zMat
                          << "\n  n                = " << n
                          << "\n  sigmaEps         = " << sigmaEps
                          << "\n  yMeanVec         = " << yMeanVec
                          << "\n  yCovMat          = " << yCovMat
                          << "\n  ySamples         = " << ySamples
                          << std::endl;
  }

#if 0
  uqGslMatrixClass sigmaMatInverse(paramSpace.zeroVector());
  sigmaMatInverse = matrixProduct(aVec,aVec);
  sigmaMatInverse *= (((double) n)/sigmaEps/sigmaEps);
  sigmaMatInverse += sigma0Mat;
  uqGslMatrixClass sigmaMat(paramSpace.zeroVector());
  sigmaMat = sigmaMatInverse.inverse();

  uqGslVectorClass muVec(paramSpace.zeroVector());
  muVec = sigmaMat * aVec;
  muVec *= 1.; //(((double) n) * ySampleMean)/sigmaEps/sigmaEps;

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSipForX0():"
                          << "\n  muVec            = " << muVec
                          << "\n  sigmaMat         = " << sigmaMat
                          << "\n  sigmaMatInverse  = " << sigmaMatInverse
                          << std::endl;
  }
#endif
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
    *env.subDisplayFile() << "Leaving solveSipForX0()"
                          << std::endl;
  }

  return;
}

double likelihoodRoutineForX0(
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
    *env.subDisplayFile() << "Entering likelihoodRoutineForX0()..."
                          << ": likelihoodCounter = "       << likelihoodCounter
                          << ", params = "                  << paramValues
                          << ", env.subComm().NumProc() = " << env.subComm().NumProc()
                          << ", my subRank = "              << env.subRank()
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Entering likelihoodRoutineForX0()"
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

  struct likelihoodDataStructForX0* likelihoodDataForX0 = (likelihoodDataStructForX0 *) functionDataPtr; 
  uqGslMatrixClass zMat(*(likelihoodDataForX0->zMat));
  unsigned int p = zMat.numRowsLocal();
  double sigmaEps = likelihoodDataForX0->sigmaEps;
  uqGslVectorClass ySamples(*(likelihoodDataForX0->ySamples));
  unsigned int n = ySamples.sizeLocal();

  UQ_FATAL_TEST_MACRO(paramValues.sizeLocal() != p,
                      env.fullRank(),
                      "likelihoodRoutineForX0()",
                      "invalid parameter vector size");

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
    *env.subDisplayFile() << "In likelihoodRoutineForX0()"
                          << ": likelihoodCounter = " << likelihoodCounter
                          << ", params = "            << paramValues
                          << ", p = "                 << p
                          << ", zMat = "              << zMat
                          << ", sigmaEps = "          << sigmaEps
                          << ", n = "                 << n
                          << ", ySamples = "          << ySamples
                          << std::endl;
  }

  //******************************************************************************
  // Compute likelihood
  //******************************************************************************
  for (unsigned int i = 0; i < n; ++i) {
    double diff = 0.;//(ySamples[i] - scalarProduct(aVec,paramValues))/sigmaEps;
    totalLnLikelihood -= 0.5 * diff * diff;
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
      *env.subDisplayFile() << "In likelihoodRoutineForX0()"
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
    *env.subDisplayFile() << "Leaving likelihoodRoutineForX0()"
                          << ": likelihoodCounter = " << likelihoodCounter
                          << ", params = "            << paramValues
                          << ", totalLnLikelihood = " << totalLnLikelihood
                          << " after "                << totalTime
                          << " seconds"
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Leaving likelihoodRoutineForX0()"
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
