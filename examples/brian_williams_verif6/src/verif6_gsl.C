#include <verif6_gsl.h>
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
                      "run as <executable> 'inputFileName' <'0 or 1'(useML; default 0)>");
  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);
  bool useML = false;
  if (argc >= 3) {
    useML = true;
  }

  //***********************************************************************
  // Run program
  //***********************************************************************
  solveSip(*env,useML);

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

void solveSip(const uqFullEnvironmentClass& env, bool useML)
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
  std::vector<double> as(3,0.);
  as[0] = 121370.92;
  as[1] =  17333.65;
  as[2] =  -9378.29;
  uqGslVectorClass bVec(paramSpace.zeroVector());
  bVec[0] = 50.45;
  bVec[1] = -2.99;

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
  unsigned int n = 60;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> dataSpace(env, "data_", n, NULL);

  std::vector<double> sigmas(4,0.);
  sigmas[0] = 4.148751;
  sigmas[1] = 4.170606;
  sigmas[2] = 4.190899;
  sigmas[3] = 4.209629;

  std::set<unsigned int> tmpSet;
  tmpSet.insert(env.subId());

  uqGslMatrixClass dataMat(env,dataSpace.map(),(unsigned int) 2);
  uqGslVectorClass tVec(dataSpace.zeroVector());
  uqGslVectorClass yVec(dataSpace.zeroVector());
  dataMat.subReadContents("input/dataPoints",
                          "m",
                          tmpSet);
  dataMat.getColumn(0,tVec);
  dataMat.getColumn(1,yVec);

  struct likelihoodDataStruct likelihoodData;
  likelihoodData.as       = &as;
  likelihoodData.bVec     = &bVec;
  likelihoodData.sigmas   = &sigmas;
  likelihoodData.tVec     = &tVec;
  likelihoodData.yVec     = &yVec;

  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>
    likelihoodFunctionObj("like_",
                          *paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodData,
                          true); // routine computes [ln(function)]

  ////////////////////////////////////////////////////////
  // Step 4 of 5: Instantiate the inverse problem
  ////////////////////////////////////////////////////////
  uqUniformVectorRVClass<uqGslVectorClass,uqGslMatrixClass> priorRv("prior_", *paramDomain);

  uqGenericVectorRVClass<uqGslVectorClass,uqGslMatrixClass> postRv ("post_", paramSpace);

  uqStatisticalInverseProblemClass<uqGslVectorClass,uqGslMatrixClass> sip("sip_", NULL, priorRv, likelihoodFunctionObj, postRv);

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSip():"
                          << "\n  p         = " << p
                          << "\n  as[0]     = " << as[0]
                          << "\n  as[1]     = " << as[1]
                          << "\n  as[2]     = " << as[2]
                          << "\n  bVec      = " << bVec
                          << "\n  n         = " << n
                          << "\n  sigmas[0] = " << sigmas[0]
                          << "\n  sigmas[1] = " << sigmas[1]
                          << "\n  sigmas[2] = " << sigmas[2]
                          << "\n  sigmas[3] = " << sigmas[3]
                          << "\n  tVec      = " << tVec
                          << "\n  yVec      = " << yVec
                          << "\n  useML     = " << useML
                          << std::endl;
  }

  ////////////////////////////////////////////////////////
  // Step 5 of 5: Solve the inverse problem
  ////////////////////////////////////////////////////////
  uqGslVectorClass initialValues(paramSpace.zeroVector());
  initialValues[0] = 0.;

  uqGslMatrixClass proposalCovMat(paramSpace.zeroVector());
  proposalCovMat(0,0) = 1.;

  if (useML) {
    sip.solveWithBayesMLSampling();
  }
  else {
    sip.solveWithBayesMetropolisHastings(NULL,initialValues,&proposalCovMat);
  }

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
  std::vector<double>& as = *(likelihoodData->as);
  uqGslVectorClass bVec(*(likelihoodData->bVec));
  unsigned int p = bVec.sizeLocal();
  std::vector<double>& sigmas = *(likelihoodData->sigmas);
  uqGslVectorClass tVec(*(likelihoodData->tVec));
  uqGslVectorClass yVec(*(likelihoodData->yVec));
  unsigned int n = yVec.sizeLocal();

  UQ_FATAL_TEST_MACRO(paramValues.sizeLocal() != p,
                      env.fullRank(),
                      "likelihoodRoutine()",
                      "invalid parameter vector size");

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
    *env.subDisplayFile() << "In likelihoodRoutine()"
                          << ": likelihoodCounter = " << likelihoodCounter
                          << ", params = "            << paramValues
                          << "\n  p         = " << p
                          << "\n  as[0]     = " << as[0]
                          << "\n  as[1]     = " << as[1]
                          << "\n  as[2]     = " << as[2]
                          << "\n  bVec      = " << bVec
                          << "\n  n         = " << n
                          << "\n  sigmas[0] = " << sigmas[0]
                          << "\n  sigmas[1] = " << sigmas[1]
                          << "\n  sigmas[2] = " << sigmas[2]
                          << "\n  sigmas[3] = " << sigmas[3]
                          << "\n  tVec      = " << tVec
                          << "\n  yVec      = " << yVec
                          << std::endl;
  }

  //******************************************************************************
  // Compute likelihood
  //******************************************************************************
  for (unsigned int i = 0; i < n; ++i) {
    double t = tVec[i];
    double diff = (yVec[i] - as[0] - as[1]*t - as[2]*t*t - scalarProduct(bVec,paramValues))/sigmas[i%4];
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
