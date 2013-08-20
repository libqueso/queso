#include <verif3_gsl.h>
#include <uqStatisticalInverseProblem.h>

static unsigned int likelihoodCounterForX0 = 0;
static unsigned int likelihoodCounterForX  = 0;
static unsigned int priorCounterForX       = 0;

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
  solveSips(*env);

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

void solveSips(const uqFullEnvironmentClass& env)
{
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Entering solveSips()..."
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
  uqGslMatrixClass sigmaMatInverse(paramSpace.zeroVector());
  sigmaMatInverse = sigmaMat.inverse();

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
  // Step 3 of 5: Instantiate the likelihood function objects
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
  
  uqGslMatrixClass zMat         (env,paramSpace.map(),n);
  uqGslMatrixClass zMatTranspose(env,dataSpace.map (),p);

  //uqGslVectorClass yMeanVec(dataSpace.zeroVector());
  uqGslVectorClass ySamples(dataSpace.zeroVector());
  for (unsigned int i = 0; i < n; ++i) {
    zRv.realizer().realization(zSample);
    zMat.setColumn(i,zSample);
    epsRv.realizer().realization(epsSample);
    ySamples[i]  = scalarProduct(zSample,xGiven) + epsSample[0];
    //yMeanVec[i] = 0.;
  }
  zMatTranspose.fillWithTranspose(0,0,zMat,true,true);

  //uqGslMatrixClass yCovMat(dataSpace.zeroVector());
  //double tmp = sigmaEps*sigmaEps;
  //for (unsigned int i = 0; i < n; ++i) {
  //  yCovMat(i,i) = tmp;
  //}

  //uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> yRv("y_", dataSpace, yMeanVec, yCovMat);
  //yRv.realizer().realization(ySamples);

  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> xRv("x_", paramSpace, paramSpace.zeroVector(), sigmaMatInverse); // In order to save computational time

  struct likelihoodDataStructForX0 likelihoodDataForX0;
  likelihoodDataForX0.zMat            = &zMat;
  likelihoodDataForX0.sigmaEps        = sigmaEps;
  likelihoodDataForX0.ySamples        = &ySamples;
  likelihoodDataForX0.sigmaMatInverse = &sigmaMatInverse;
  likelihoodDataForX0.xRv             = &xRv; // In order to save computational time

  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>
    likelihoodFunctionObjForX0("like_",
                               *paramDomain,
                               likelihoodRoutineForX0,
                               (void *) &likelihoodDataForX0,
                               true); // routine computes [ln(function)]

  struct likelihoodDataStructForX likelihoodDataForX;
  likelihoodDataForX.zMat            = &zMat;
  likelihoodDataForX.sigmaEps        = sigmaEps;
  likelihoodDataForX.ySamples        = &ySamples;
  likelihoodDataForX.sigmaMatInverse = &sigmaMatInverse;
  likelihoodDataForX.xRv             = &xRv; // In order to save computational time

  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>
    likelihoodFunctionObjForX("like_",
                               *paramDomain,
                               likelihoodRoutineForX,
                               (void *) &likelihoodDataForX,
                               true); // routine computes [ln(function)]

  ////////////////////////////////////////////////////////
  // Step 4 of 5: Instantiate the inverse problems
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

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSips():"
                          << "\n  p                = " << p
                          << "\n  xGiven           = " << xGiven
                          << "\n  mu0Vec           = " << mu0Vec
                          << "\n  sigma0Mat        = " << sigma0Mat
                          << "\n  sigma0MatInverse = " << sigma0MatInverse
                          << "\n  sigmaMat         = " << sigmaMat
                          << "\n  zMeanVec         = " << zMeanVec
                          << "\n  zCovMat          = " << zCovMat
                          << "\n  n                = " << n
                          << "\n  sigmaEps         = " << sigmaEps
                          << "\n  zMat             = " << zMat
      //<< "\n  yMeanVec         = " << yMeanVec
      //<< "\n  yCovMat          = " << yCovMat
                          << "\n  ySamples         = " << ySamples
                          << std::endl;
  }

  uqGslMatrixClass tmpZZt(paramSpace.zeroVector());
  tmpZZt = zMat * zMatTranspose;
  uqGslMatrixClass tmpZZtScaled(paramSpace.zeroVector());
  tmpZZtScaled = tmpZZt;
  tmpZZtScaled /= (sigmaEps * sigmaEps);
  uqGslMatrixClass tmpZZtInverse(paramSpace.zeroVector());
  tmpZZtInverse = tmpZZt.inverse();
  uqGslMatrixClass tmpZZtInverseScaled(paramSpace.zeroVector());
  tmpZZtInverseScaled = tmpZZtInverse;
  tmpZZtInverseScaled *= sigmaEps * sigmaEps;

  uqGslMatrixClass tmpMat1(paramSpace.zeroVector());
  tmpMat1 = tmpZZtInverseScaled;
  tmpMat1 += sigmaMatInverse;
  uqGslMatrixClass commonMat(paramSpace.zeroVector());
  commonMat = tmpMat1.inverse();
  uqGslMatrixClass sigmaPiMatInverse(paramSpace.zeroVector());
  sigmaPiMatInverse = sigma0Mat + commonMat;
  uqGslMatrixClass sigmaPiMat(paramSpace.zeroVector());
  sigmaPiMat = sigmaPiMatInverse.inverse();

  uqGslVectorClass xHatVec(paramSpace.zeroVector());
  xHatVec = tmpZZtInverse * (zMat * ySamples);
  uqGslVectorClass muPiVec(paramSpace.zeroVector());
  muPiVec = sigmaPiMat * ((commonMat * xHatVec) + (sigma0Mat * mu0Vec));

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSips():"
                          << "\n  xHatVec           = " << xHatVec
                          << "\n  muPiVec           = " << muPiVec
                          << "\n  sigmaPiMat        = " << sigmaPiMat
                          << "\n  sigmaPiMatInverse = " << sigmaPiMatInverse
                          << std::endl;
  }

  uqGslMatrixClass sigma0PlusSigmaMat(paramSpace.zeroVector());
  sigma0PlusSigmaMat = sigma0Mat + sigmaMat;
  uqGslMatrixClass sigma0PlusSigmaMatInverse(paramSpace.zeroVector());
  sigma0PlusSigmaMatInverse = sigma0PlusSigmaMat.inverse();

  uqGslVectorClass x0LimitVec(paramSpace.zeroVector());
  x0LimitVec = sigma0PlusSigmaMatInverse * ((sigmaMat * xGiven) + (sigma0Mat * mu0Vec));
  uqGslMatrixClass x0LimitMat(paramSpace.zeroVector());
  x0LimitMat = sigma0PlusSigmaMatInverse;

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSips():"
                          << "\n  x0LimitVec = " << x0LimitVec
                          << "\n  x0LimitMat = " << x0LimitMat
                          << std::endl;
  }

  uqGslMatrixClass sigma0InvPlusSigmaInvMat(paramSpace.zeroVector());
  sigma0InvPlusSigmaInvMat = sigma0MatInverse + sigmaMatInverse;
  uqGslMatrixClass sigma0InvPlusSigmaInvMatInverse(paramSpace.zeroVector());
  sigma0InvPlusSigmaInvMatInverse = sigma0InvPlusSigmaInvMat.inverse();

  uqGslMatrixClass sigmatMatInverse(paramSpace.zeroVector());
  sigmatMatInverse = sigma0InvPlusSigmaInvMatInverse + tmpZZtScaled;
  uqGslMatrixClass sigmatMat(paramSpace.zeroVector());
  sigmatMat = sigmatMatInverse.inverse();
  uqGslVectorClass tmpVec(paramSpace.zeroVector());
  tmpVec = zMat * ySamples;
  tmpVec /= (sigmaEps * sigmaEps);
  uqGslVectorClass mutVec(paramSpace.zeroVector());
  mutVec = sigmaMat * (tmpVec + (sigma0InvPlusSigmaInvMatInverse * mu0Vec));

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In solveSips():"
                          << "\n  mutVec           = " << mutVec
                          << "\n  sigmatMat        = " << sigmatMat
                          << "\n  sigmatMatInverse = " << sigmatMatInverse
                          << std::endl;
  }

  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass> priorRvForX0("prior_", *paramDomain, mu0Vec, sigma0MatInverse);
  uqGenericVectorRVClass <uqGslVectorClass,uqGslMatrixClass> postRvForX0 ("post_", paramSpace);
  uqStatisticalInverseProblemClass<uqGslVectorClass,uqGslMatrixClass> sipForX0("sipForX0_", NULL, priorRvForX0, likelihoodFunctionObjForX0, postRvForX0);

  struct dataStructForPriorRoutineForX dataForPriorRoutineForX;
  dataForPriorRoutineForX.priorRvForX0 = &priorRvForX0;
  dataForPriorRoutineForX.xRv          = &xRv;
  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass> genericScalarFunction("",*paramDomain,routinePriorPdfForX,&dataForPriorRoutineForX,true);
  uqGenericJointPdfClass      <uqGslVectorClass,uqGslMatrixClass> priorPdfForX         ("",genericScalarFunction);
//uqGenericVectorRealizerClass<uqGslVectorClass,uqGslMatrixClass> priorRealizerForX    ("",*paramDomain,1000,routinePriorRealizerForX,NULL);
  uqGenericVectorRVClass      <uqGslVectorClass,uqGslMatrixClass> priorRvForX          ("prior_", *paramDomain);
  priorRvForX.setPdf     (priorPdfForX);
//priorRvForX.setRealizer(priorRealizerForX);
  uqGenericVectorRVClass<uqGslVectorClass,uqGslMatrixClass> postRvForX ("post_", paramSpace);
  uqStatisticalInverseProblemClass<uqGslVectorClass,uqGslMatrixClass> sipForX("sipForX_", NULL, priorRvForX, likelihoodFunctionObjForX, postRvForX);

  ////////////////////////////////////////////////////////
  // Step 5 of 5: Solve the inverse problems
  ////////////////////////////////////////////////////////
  uqGslVectorClass initialValues(paramSpace.zeroVector());
  initialValues[0] = 25.;
  initialValues[1] = 25.;

  uqGslMatrixClass proposalCovMat(paramSpace.zeroVector());
  proposalCovMat(0,0) = 10.;
  proposalCovMat(0,1) = 0.;
  proposalCovMat(1,0) = 0.;
  proposalCovMat(1,1) = 10.;

  if (env.subRank() == 0) {
    std::cout << "Beginning to solve sipForX0" << std::enl;
  }
  sipForX0.solveWithBayesMetropolisHastings(NULL,initialValues,&proposalCovMat);
  if (env.subRank() == 0) {
    std::cout << "Finished solving sipForX0" << std::enl;
  }

  if (env.subRank() == 0) {
    std::cout << "Beginning to solve sipForX" << std::enl;
  }
  sipForX.solveWithBayesMetropolisHastings(NULL,initialValues,&proposalCovMat);
  if (env.subRank() == 0) {
    std::cout << "Finished solving sipForX" << std::enl;
  }

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Leaving solveSips()"
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

  likelihoodCounterForX0++;
  const uqBaseEnvironmentClass& env = paramValues.env();
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Entering likelihoodRoutineForX0()..."
                          << ": likelihoodCounterForX0 = "  << likelihoodCounterForX0
                          << ", params = "                  << paramValues
                          << ", env.subComm().NumProc() = " << env.subComm().NumProc()
                          << ", my subRank = "              << env.subRank()
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Entering likelihoodRoutineForX0()"
              << ", likelihoodCounterForX0 = " << likelihoodCounterForX0
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
  uqGslMatrixClass sigmaMatInverse(*(likelihoodDataForX0->sigmaMatInverse));
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>* xRv = likelihoodDataForX0->xRv;

  UQ_FATAL_TEST_MACRO(paramValues.sizeLocal() != p,
                      env.fullRank(),
                      "likelihoodRoutineForX0()",
                      "invalid parameter vector size");

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
    *env.subDisplayFile() << "In likelihoodRoutineForX0()"
                          << ": likelihoodCounterForX0 = " << likelihoodCounterForX0
                          << ", params = "                 << paramValues
                          << ", p = "                      << p
                          << ", zMat = "                   << zMat
                          << ", sigmaEps = "               << sigmaEps
                          << ", n = "                      << n
                          << ", ySamples = "               << ySamples
                          << std::endl;
  }

  //******************************************************************************
  // Compute likelihood
  //******************************************************************************
  xRv->updateLawExpVector(paramValues);
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paramSpace(env, "param_", p, NULL);
  uqGslVectorClass xMonteCarlo(paramSpace.zeroVector());
  uqGslVectorClass zCol(paramSpace.zeroVector());
  unsigned int numMonteCarloSamples = 256;
  double mcSum = 0.;
  for (unsigned int mcId = 0; mcId < numMonteCarloSamples; ++mcId) {
    xRv->realizer().realization(xMonteCarlo);
    double tmpSum = 0.;
    for (unsigned int i = 0; i < n; ++i) {
      zMat.getColumn(i,zCol);
      double diff = (ySamples[i] - scalarProduct(zCol,xMonteCarlo))/sigmaEps;
      tmpSum += std::exp(-0.5 * diff * diff);
      if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
        *env.subDisplayFile() << "In likelihoodRoutineForX0()"
                              << ": likelihoodCounterForX0 = " << likelihoodCounterForX0
                              << ", params = "                 << paramValues
                              << ", mcId = "                   << mcId
                              << ", i = "                      << i
                              << ", tmpSum = "                 << tmpSum
                              << std::endl;
      }
    }
    mcSum += tmpSum;
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
      *env.subDisplayFile() << "In likelihoodRoutineForX0()"
                            << ": likelihoodCounterForX0 = " << likelihoodCounterForX0
                            << ", params = "                 << paramValues
                            << ", mcId = "                   << mcId
                            << ", mcSum = "                  << mcSum
                            << std::endl;
    }
  }
  totalLnLikelihood = std::log(mcSum) - std::log((double) numMonteCarloSamples);

  //******************************************************************************
  // Prepare to return
  //******************************************************************************
  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Leaving likelihoodRoutineForX0()"
                          << ": likelihoodCounterForX0 = " << likelihoodCounterForX0
                          << ", params = "                 << paramValues
                          << ", totalLnLikelihood = "      << totalLnLikelihood
                          << " after "                     << totalTime
                          << " seconds"
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Leaving likelihoodRoutineForX0()"
              << ": likelihoodCounterForX0 = " << likelihoodCounterForX0
              << ", params = "                 << paramValues
              << ", totalLnLikelihood = "      << totalLnLikelihood
              << " after "                     << totalTime
              << " seconds"
              << std::endl;
#endif
  }

  env.subComm().Barrier();
  //exit(1);

  return totalLnLikelihood;
}

double likelihoodRoutineForX(
  const uqGslVectorClass& paramValues,
  const uqGslVectorClass* paramDirection,
  const void*             functionDataPtr,
  uqGslVectorClass*       gradVector,
  uqGslMatrixClass*       hessianMatrix,
  uqGslVectorClass*       hessianEffect)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  likelihoodCounterForX++;
  const uqBaseEnvironmentClass& env = paramValues.env();
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Entering likelihoodRoutineForX()..."
                          << ": likelihoodCounterForX = "   << likelihoodCounterForX
                          << ", params = "                  << paramValues
                          << ", env.subComm().NumProc() = " << env.subComm().NumProc()
                          << ", my subRank = "              << env.subRank()
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Entering likelihoodRoutineForX()"
              << ", likelihoodCounterForX = " << likelihoodCounterForX
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

  struct likelihoodDataStructForX* likelihoodDataForX = (likelihoodDataStructForX *) functionDataPtr; 
  uqGslMatrixClass zMat(*(likelihoodDataForX->zMat));
  unsigned int p = zMat.numRowsLocal();
  double sigmaEps = likelihoodDataForX->sigmaEps;
  uqGslVectorClass ySamples(*(likelihoodDataForX->ySamples));
  unsigned int n = ySamples.sizeLocal();
  uqGslMatrixClass sigmaMatInverse(*(likelihoodDataForX->sigmaMatInverse));

  UQ_FATAL_TEST_MACRO(paramValues.sizeLocal() != p,
                      env.fullRank(),
                      "likelihoodRoutineForX()",
                      "invalid parameter vector size");

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
    *env.subDisplayFile() << "In likelihoodRoutineForX()"
                          << ": likelihoodCounterForX = " << likelihoodCounterForX
                          << ", params = "                << paramValues
                          << ", p = "                     << p
                          << ", zMat = "                  << zMat
                          << ", sigmaEps = "              << sigmaEps
                          << ", n = "                     << n
                          << ", ySamples = "              << ySamples
                          << std::endl;
  }

  //******************************************************************************
  // Compute likelihood
  //******************************************************************************
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paramSpace(env, "param_", p, NULL);
  uqGslVectorClass zCol(paramSpace.zeroVector());
  for (unsigned int i = 0; i < n; ++i) {
    zMat.getColumn(i,zCol);
    double diff = (ySamples[i] - scalarProduct(zCol,paramValues))/sigmaEps;
    totalLnLikelihood -= 0.5 * diff * diff;
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
      *env.subDisplayFile() << "In likelihoodRoutineForX()"
                            << ": likelihoodCounterForX = " << likelihoodCounterForX
                            << ", params = "                << paramValues
                            << ", i = "                     << i
                            << ", diff = "                  << diff
                            << std::endl;
    }
  }

  //******************************************************************************
  // Prepare to return
  //******************************************************************************
  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Leaving likelihoodRoutineForX()"
                          << ": likelihoodCounterForX = " << likelihoodCounterForX
                          << ", params = "                << paramValues
                          << ", totalLnLikelihood = "     << totalLnLikelihood
                          << " after "                    << totalTime
                          << " seconds"
                          << std::endl;
  }

  if (env.subRank() == 0) {
#if 0
    std::cout << "Leaving likelihoodRoutineForX()"
              << ": likelihoodCounterForX = " << likelihoodCounterForX
              << ", params = "                << paramValues
              << ", totalLnLikelihood = "     << totalLnLikelihood
              << " after "                    << totalTime
              << " seconds"
              << std::endl;
#endif
  }

  env.subComm().Barrier();
  //exit(1);

  return totalLnLikelihood;
}

double routinePriorPdfForX(const uqGslVectorClass& domainVector,
                           const uqGslVectorClass* domainDirection,
                           const void*             routineDataPtr,
                           uqGslVectorClass*       gradVector,
                           uqGslMatrixClass*       hessianMatrix,
                           uqGslVectorClass*       hessianEffect)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  priorCounterForX++;

  const uqBaseEnvironmentClass& env = domainVector.env();
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Entering routinePriorPdfForX()..."
                          << ": domainVector = " << domainVector
                          << std::endl;
  }

  UQ_FATAL_TEST_MACRO((domainVector.sizeLocal() != 2),
                      env.fullRank(),
                      "routinePriorPdfForX()",
                      "domainVector has wrong size");

  double value = 0.;

  if (domainDirection &&
      routineDataPtr  &&
      gradVector      &&
      hessianMatrix   &&
      hessianEffect) {
    // Just to eliminate INTEL compiler warnings
  }

  struct dataStructForPriorRoutineForX* dataForPriorRoutineForX = (dataStructForPriorRoutineForX *) routineDataPtr;
  const uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>* priorRvForX0 = dataForPriorRoutineForX->priorRvForX0;
        uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>* xRv          = dataForPriorRoutineForX->xRv;
  uqGslVectorClass x0MonteCarlo(priorRvForX0->imageSet().vectorSpace().zeroVector());
  unsigned int numMonteCarloSamples = 256;
  double mcSum = 0.;
  for (unsigned int mcId = 0; mcId < numMonteCarloSamples; ++mcId) {
    priorRvForX0->realizer().realization(x0MonteCarlo);
    xRv->updateLawExpVector(x0MonteCarlo);
    double pdfValueForX  = xRv->pdf().actualValue(domainVector,NULL,NULL,NULL,NULL);
    double pdfValueForX0 = priorRvForX0->pdf().actualValue(x0MonteCarlo,NULL,NULL,NULL,NULL);
    mcSum += (pdfValueForX * pdfValueForX0);
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 4)) {
      *env.subDisplayFile() << "In routinePriorPdfForX()"
                            << ": priorCounterForX = " << priorCounterForX
                            << ", domainVector = "     << domainVector
                            << ", mcId = "             << mcId
                            << ", mcSum = "            << mcSum
                            << std::endl;
    }
  }
  value = std::log(mcSum) - std::log((double) numMonteCarloSamples);

  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
    *env.subDisplayFile() << "Leaving routinePriorPdfForX()"
                          << ": domainVector = " << domainVector
                          << ", value = "        << value
                          << " after "           << totalTime
                          << " seconds"
                          << std::endl;
  }

  return value;
}
