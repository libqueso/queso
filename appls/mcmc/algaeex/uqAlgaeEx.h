#ifndef __UQ_ALGAE_EX_H__
#define __UQ_ALGAE_EX_H__

double observationsOfAZP[] = {
   1.0000,   0.7674,   1.2794,   9.6170,
   5.0000,   0.6364,   1.9551,   7.7312,
   9.0000,   2.3456,   2.2325,   9.9886,
  13.0000,   1.4072,   1.6362,   8.8322,
  17.0000,   2.1013,   1.4278,   3.2399,
  21.0000,   3.3236,   1.6267,  12.9613,
  25.0000,  11.9567,   0.6698,  12.0581,
  29.0000,  13.1154,   3.4193,   6.3863,
  33.0000,  13.6240,   6.7866,   4.5111,
  37.0000,   9.3415,   8.3476,   5.1159,
  41.0000,   3.4367,  10.2362,   6.1919,
  45.0000,   0.6112,   8.2290,   9.6645,
  49.0000,   1.6024,   9.1845,  14.1957,
  53.0000,   2.6171,   5.5917,  11.4486,
  57.0000,   1.9495,   5.1942,  13.4226,
  61.0000,   2.5596,   3.8916,  13.3308,
  65.0000,   3.6307,   5.8413,   8.1629,
  69.0000,   4.1477,   4.7919,   8.9506,
  73.0000,   7.2978,   2.8895,   5.7932,
  77.0000,   5.2055,   2.6108,   6.4463,
  81.0000,   5.2294,   4.8424,   5.8734,
  85.0000,   3.0761,   3.8744,   7.3142,
  89.0000,   2.3220,   4.6127,   5.2613,
  93.0000,   3.5118,   3.7384,   7.3627,
  97.0000,   4.3230,   6.2917,   6.1683,
 101.0000,   3.4297,   6.5543,   5.5262,
 105.0000,   3.0706,   4.4377,   9.4023,
 109.0000,   3.2402,   4.2605,   9.6588,
 113.0000,   1.9054,   2.7570,  10.0952,
 117.0000,   5.3993,   6.3093,   6.5072
};

double observationsOfSomeParameters[] = {
    1.0000,   0.0380,   8.6756,  10.0000,
    2.0000,   0.0380,   9.0193,  10.0000,
    3.0000,   0.0380,   9.3630,  10.0000,
    4.0000,   0.0380,   9.7067,  10.0000,
    5.0000,   0.0380,  10.0504,  10.0000,
    6.0000,   0.0380,  10.3941,  10.0000,
    7.0000,   0.0380,  10.7378,  10.0000,
    8.0000,   0.0380,  11.0815,  10.0000,
    9.0000,   0.0380,  11.4252,  10.0000,
   10.0000,   0.0375,  11.7689,  10.0000,
   11.0000,   0.0375,  12.1126,  10.0000,
   12.0000,   0.0375,  12.4563,  10.0000,
   13.0000,   0.0355,  12.8000,  10.0000,
   14.0000,   0.0340,  12.7599,  10.0000,
   15.0000,   0.0340,  12.7198,  10.0000,
   16.0000,   0.0335,  12.6797,  10.0000,
   17.0000,   0.0320,  12.6396,  10.0000,
   18.0000,   0.0310,  12.5995,  10.0000,
   19.0000,   0.0285,  12.5594,  10.0000,
   20.0000,   0.0260,  12.5193, 100.0000,
   21.0000,   0.0260,  12.4793, 100.0000,
   22.0000,   0.0255,  12.4392, 100.0000,
   23.0000,   0.0255,  12.5255, 100.0000,
   24.0000,   0.0255,  12.6118, 100.0000,
   25.0000,   0.0250,  12.6982, 100.0000,
   26.0000,   0.0250,  12.7845,  10.0000,
   27.0000,   0.0240,  12.8708,  10.0000,
   28.0000,   0.0240,  12.9572,  10.0000,
   29.0000,   0.0250,  13.0435,  10.0000,
   30.0000,   0.0240,  13.1299,  10.0000,
   31.0000,   0.0240,  13.2162,  10.0000,
   32.0000,   0.0240,  13.3025,  10.0000,
   33.0000,   0.0240,  13.3889,  10.0000,
   34.0000,   0.0265,  13.6592,  10.0000,
   35.0000,   0.0265,  13.9296,  10.0000,
   36.0000,   0.0265,  14.2000,  10.0000,
   37.0000,   0.0265,  14.2364,  10.0000,
   38.0000,   0.0265,  14.2727,  10.0000,
   39.0000,   0.0265,  14.3091,  10.0000,
   40.0000,   0.0265,  14.3455,  10.0000,
   41.0000,   0.0270,  14.3818,  10.0000,
   42.0000,   0.0255,  14.4182,  10.0000,
   43.0000,   0.0230,  14.4545,  10.0000,
   44.0000,   0.0230,  14.4909,  10.0000,
   45.0000,   0.0230,  14.5273,  10.0000,
   46.0000,   0.0235,  14.5636,  10.0000,
   47.0000,   0.0250,  14.6000,  10.0000,
   48.0000,   0.0290,  14.8104,  10.0000,
   49.0000,   0.0335,  15.0207,  10.0000,
   50.0000,   0.0335,  15.2311,  10.0000,
   51.0000,   0.0335,  15.4415,  10.0000,
   52.0000,   0.0335,  15.6519,  10.0000,
   53.0000,   0.0335,  15.8622,  10.0000,
   54.0000,   0.0335,  16.0726,  10.0000,
   55.0000,   0.0340,  16.2830,  10.0000,
   56.0000,   0.0340,  16.4934,  10.0000,
   57.0000,   0.0350,  16.7037,  10.0000,
   58.0000,   0.0350,  16.9141,  10.0000,
   59.0000,   0.0350,  17.1245,  10.0000,
   60.0000,   0.0355,  17.3349,  10.0000,
   61.0000,   0.0355,  17.5452,  10.0000,
   62.0000,   0.0355,  17.7556,  10.0000,
   63.0000,   0.0340,  18.0529,  10.0000,
   64.0000,   0.0350,  18.0799,  10.0000,
   65.0000,   0.0350,  18.1068,  10.0000,
   66.0000,   0.0350,  18.1338,  10.0000,
   67.0000,   0.0350,  18.1607,  10.0000,
   68.0000,   0.0350,  18.1877,  10.0000,
   69.0000,   0.0355,  18.2146,  10.0000,
   70.0000,   0.0355,  18.2416,  10.0000,
   71.0000,   0.0355,  18.2685,  10.0000,
   72.0000,   0.0325,  18.2955,  10.0000,
   73.0000,   0.0290,  18.3224,  10.0000,
   74.0000,   0.0285,  18.3494,  10.0000,
   75.0000,   0.0285,  18.3763,  10.0000,
   76.0000,   0.0285,  18.4033,  10.0000,
   77.0000,   0.0285,  18.4303,  10.0000,
   78.0000,   0.0285,  18.3140,  10.0000,
   79.0000,   0.0285,  18.1977,  10.0000,
   80.0000,   0.0280,  18.0814,  10.0000,
   81.0000,   0.0280,  17.9651,  10.0000,
   82.0000,   0.0280,  17.8488,  10.0000,
   83.0000,   0.0280,  17.7326,  10.0000,
   84.0000,   0.0280,  17.6163,  10.0000,
   85.0000,   0.0280,  17.5000,  10.0000,
   86.0000,   0.0280,  17.6284,  10.0000,
   87.0000,   0.0280,  17.7568,  10.0000,
   88.0000,   0.0270,  17.8852,  10.0000,
   89.0000,   0.0270,  18.0136,  10.0000,
   90.0000,   0.0260,  18.1420,  10.0000,
   91.0000,   0.0255,  17.9926,  10.0000,
   92.0000,   0.0250,  17.8432,  10.0000,
   93.0000,   0.0250,  17.6939,  10.0000,
   94.0000,   0.0250,  17.5445,  10.0000,
   95.0000,   0.0250,  17.3951,  10.0000,
   96.0000,   0.0250,  17.2457,  10.0000,
   97.0000,   0.0250,  17.0963,  10.0000,
   98.0000,   0.0255,  16.9469,  10.0000,
   99.0000,   0.0255,  16.7975,  10.0000,
  100.0000,   0.0255,  16.6482,  10.0000,
  101.0000,   0.0250,  16.4988,  10.0000,
  102.0000,   0.0255,  16.3494,  10.0000,
  103.0000,   0.0255,  16.2000,  10.0000,
  104.0000,   0.0255,  16.1972,  10.0000,
  105.0000,   0.0255,  16.1945,  10.0000,
  106.0000,   0.0255,  16.1917,  10.0000,
  107.0000,   0.0260,  16.0581,  10.0000,
  108.0000,   0.0260,  15.9245,  10.0000,
  109.0000,   0.0260,  15.7909,  10.0000,
  110.0000,   0.0260,  15.6572,  10.0000,
  111.0000,   0.0260,  15.5236,  10.0000,
  112.0000,   0.0260,  15.3900,  10.0000,
  113.0000,   0.0265,  15.2564,  10.0000,
  114.0000,   0.0265,  15.1227,  10.0000,
  115.0000,   0.0265,  14.9891,  10.0000,
  116.0000,   0.0265,  14.8555,  10.0000,
  117.0000,   0.0265,  14.7219,  10.0000,
  118.0000,   0.0265,  14.5883,  10.0000,
  119.0000,   0.0265,  14.6385,  10.0000,
  120.0000,   0.0265,  14.6887,  10.0000
};

#include <uqAnyAppl.h>
#include <uqDRAM_MarkovChainGenerator.h>
#include <uq_defines.h>
#ifdef __APPL_USES_GSL__
#include <uqGslOdeSolver.h>
int uqGslAlgaeExStateDot(double t, const double currentState[], double stateDot[], void* infoForComputingStateDot);
#endif

template<class V, class M>
struct
uqAppl_M2lLikelihoodFunction_DataType
{
  const std::vector<double>* instantsOfAZPObservations;
  const std::vector<V*>*     observedEvolutionOfAZPConcentrations;
  const std::vector<double>* evolutionOfQpV;
  const std::vector<double>* evolutionOfT;
  const std::vector<double>* evolutionOfPin;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass*       gslOdeSolver;
#endif
};

template<class V, class M>
struct
uqAppl_StateDotFunction_DataType
{
  double                     muMax;
  double                     rhoA;
  double                     rhoZ;
  double                     k;
  double                     alpha;
  double                     th;
  const std::vector<double>* evolutionOfQpV;
  const std::vector<double>* evolutionOfT;
  const std::vector<double>* evolutionOfPin;
};

template<class V, class M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) std::cout << "Beginning run of 'uqAlgaeEx' example\n"
                                 << std::endl;

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //*****************************************************
  // Step 1 of 4: Define the finite dimensional linear spaces.
  //              Define the Markov chain generator.
  //*****************************************************
  uqParamSpaceClass <V,M> paramSpace (env);
  uqStateSpaceClass <V,M> stateSpace (env);
  uqOutputSpaceClass<V,M> outputSpace(env);

  uq_M2lLikelihoodFunction_Class<V,M> uq_M2lLikelihoodFunction_Obj;
  uqDRAM_MarkovChainGeneratorClass<V,M> mcg(env,
                                            paramSpace,
                                            outputSpace,
                                            NULL,  // use default prior() routine
                                            uq_M2lLikelihoodFunction_Obj);

  //******************************************************
  // Step 2 of 4: Compute the proposal covariance matrix.
  //******************************************************
  // Proposal covariance matrix will be computed internally
  // by the Markov chain generator.

  //******************************************************
  // Step 3 of 4: Prepare the data to be passed to
  //              uqAppl_M2lLikelihoodFunction_Routine().
  //******************************************************
  // The initial concentrations are also treated as parameters
  std::vector<double> instantsOfAZPObservations           (30,0.);
  std::vector<V*>     observedEvolutionOfAZPConcentrations(30,NULL);
  for (unsigned int i = 0; i < 30; ++i) {
    instantsOfAZPObservations[i]                  = observationsOfAZP[4*i+0];
    observedEvolutionOfAZPConcentrations[i]       = stateSpace.newVector();
    (*observedEvolutionOfAZPConcentrations[i])[0] = observationsOfAZP[4*i+1];
    (*observedEvolutionOfAZPConcentrations[i])[1] = observationsOfAZP[4*i+2];
    (*observedEvolutionOfAZPConcentrations[i])[2] = observationsOfAZP[4*i+3];
  }

  std::vector<double> evolutionOfQpV(120,0.);
  std::vector<double> evolutionOfT  (120,0.);
  std::vector<double> evolutionOfPin(120,0.);
  for (unsigned int i = 0; i < 120; ++i) {
    evolutionOfQpV[i] = observationsOfSomeParameters[4*i+1];
    evolutionOfT  [i] = observationsOfSomeParameters[4*i+2];
    evolutionOfPin[i] = observationsOfSomeParameters[4*i+3];
  }

#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass gslOdeSolver("gear2",stateSpace.dim());
#endif

  uqAppl_M2lLikelihoodFunction_DataType<V,M> uqAppl_M2lLikelihoodFunction_Data;
  uqAppl_M2lLikelihoodFunction_Data.instantsOfAZPObservations            = &instantsOfAZPObservations;
  uqAppl_M2lLikelihoodFunction_Data.observedEvolutionOfAZPConcentrations = &observedEvolutionOfAZPConcentrations;
  uqAppl_M2lLikelihoodFunction_Data.evolutionOfQpV                       = &evolutionOfQpV;
  uqAppl_M2lLikelihoodFunction_Data.evolutionOfT                         = &evolutionOfT;
  uqAppl_M2lLikelihoodFunction_Data.evolutionOfPin                       = &evolutionOfPin;
#ifdef __APPL_USES_GSL__
  uqAppl_M2lLikelihoodFunction_Data.gslOdeSolver                         = &gslOdeSolver;
#endif

  //******************************************************
  // Step 4 of 4: Generate chains.
  //              In MATLAB, one should then run uqAlgaeEx.m,
  //              which calls 'mcmcplot.m'.
  //******************************************************
  mcg.generateChains(NULL, // compute proposal cov matrix internally
                     NULL, // use default prior() routine
                     (void *) &uqAppl_M2lLikelihoodFunction_Data);

  //*****************************************************
  // Release memory before leaving routine.
  //*****************************************************
  for (unsigned int i = 0; i < observedEvolutionOfAZPConcentrations.size(); ++i) {
    if (observedEvolutionOfAZPConcentrations[i] != NULL) delete observedEvolutionOfAZPConcentrations[i];
  }

  if (env.rank() == 0) std::cout << "Finishing run of 'uqAlgaeEx' example"
                                 << std::endl;

  return;
}

template<class V, class M>
double
uqAppl_M2lPriorFunction_Routine(const V& paramValues, const void* functionDataPtr)
{
  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "uqAppl_M2lPriorFunction_Routine()",
                      "should not be here, since application is using the default prior() routine provided by ICES UQ library");
  return 0.;
}

template<class V, class M>
void
uqAppl_M2lLikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const std::vector<double>& instantsOfAZPObservations            = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->instantsOfAZPObservations;
  const std::vector<V*>&     observedEvolutionOfAZPConcentrations = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->observedEvolutionOfAZPConcentrations;
  const std::vector<double>& evolutionOfQpV                       = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->evolutionOfQpV;
  const std::vector<double>& evolutionOfT                         = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->evolutionOfT;
  const std::vector<double>& evolutionOfPin                       = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->evolutionOfPin;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass&       gslOdeSolver                         = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->gslOdeSolver;
#endif

  // The initial concentrations are also unknown and treated as extra parameters
  V initialConcentrations(*observedEvolutionOfAZPConcentrations[0]);
  initialConcentrations[0] = paramValues[6];
  initialConcentrations[1] = paramValues[7];
  initialConcentrations[2] = paramValues[8];

  uqAppl_StateDotFunction_DataType<V,M> uqAppl_StateDotFunction_Data;
  uqAppl_StateDotFunction_Data.muMax          = paramValues[0];
  uqAppl_StateDotFunction_Data.rhoA           = paramValues[1];
  uqAppl_StateDotFunction_Data.rhoZ           = paramValues[2];
  uqAppl_StateDotFunction_Data.k              = paramValues[3];
  uqAppl_StateDotFunction_Data.alpha          = paramValues[4];
  uqAppl_StateDotFunction_Data.th             = paramValues[5];
  uqAppl_StateDotFunction_Data.evolutionOfQpV = &evolutionOfQpV;
  uqAppl_StateDotFunction_Data.evolutionOfT   = &evolutionOfT;
  uqAppl_StateDotFunction_Data.evolutionOfPin = &evolutionOfPin;

  //if (initialConcentrations.env().rank() == 0) {
  //  std::cout << "In uqAppl_M2lLikelihoodFunction_Routine():"
  //            << "\n A0    = " << initialConcentrations[0]
  //            << "\n Z0    = " << initialConcentrations[1]
  //            << "\n P0    = " << initialConcentrations[2]
  //            << "\n muMax = " << uqAppl_StateDotFunction_Data.muMax
  //            << "\n rhoA  = " << uqAppl_StateDotFunction_Data.rhoA
  //            << "\n rhoZ  = " << uqAppl_StateDotFunction_Data.rhoZ
  //            << "\n k     = " << uqAppl_StateDotFunction_Data.k
  //           << "\n alpha = " << uqAppl_StateDotFunction_Data.alpha
  //            << "\n th    = " << uqAppl_StateDotFunction_Data.th
  //            << "\n"
  //            << std::endl;

    //std::cout << "In uqAppl_M2lLikelihoodFunction_Routine():"
    //          << " information on QpV, Pin, T"
    //          << std::endl;
    //for (unsigned int i = 0; i < evolutionOfQpV.size(); ++i) {
    //  std::cout << i                 << " "
    //            << evolutionOfQpV[i] << " "
    //            << evolutionOfT  [i] << " "
    //            << evolutionOfPin[i]
    //            << std::endl;
    //}
    //std::cout << std::endl;
  //}

  std::vector<V*> computedEvolutionOfConcentrations(0,NULL);
#ifdef __APPL_USES_GSL__
  gslOdeSolver.solveODE(uqGslAlgaeExStateDot,
                        NULL,
                        instantsOfAZPObservations,
                        initialConcentrations,
                        (void *)&uqAppl_StateDotFunction_Data, // data for stateDot computations
                        1.e-1,                                 // suggested time step
                        computedEvolutionOfConcentrations);
#endif
#ifdef __APPL_USES_TRILINOS__
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    paramValue.env().rank(),
                    "uqAppl_M2lLikelihoodFunction_Routine()",
                    "not yet implemented with Trilinos");
#endif
  //if (initialConcentrations.env().rank() == 0) {
  //  std::cout << "In uqAppl_M2lLikelihoodFunction_Routine()" << std::endl;
  //  for (unsigned int i = 0; i < instantsOfAZPObservations.size(); ++i) {
  //    std::cout << "  for t = " << instantsOfAZPObservations[i]
  //              << ", state = " << *(computedEvolutionOfConcentrations[i])
  //              << std::endl;
  //  }
  //}

  resultValues[0] = 0.;
  resultValues[1] = 0.;
  resultValues[2] = 0.;
  for (unsigned int i = 0; i < instantsOfAZPObservations.size(); ++i) {
    resultValues[0] += pow( (*computedEvolutionOfConcentrations[i])[0] - (*observedEvolutionOfAZPConcentrations[i])[0],2. );
    resultValues[1] += pow( (*computedEvolutionOfConcentrations[i])[1] - (*observedEvolutionOfAZPConcentrations[i])[1],2. );
    resultValues[2] += pow( (*computedEvolutionOfConcentrations[i])[2] - (*observedEvolutionOfAZPConcentrations[i])[2],2. );
  }

  for (unsigned int i = 0; i < computedEvolutionOfConcentrations.size(); ++i) {
    if (computedEvolutionOfConcentrations[i] != NULL) delete computedEvolutionOfConcentrations[i];
  }

  return;
}
#endif // __UQ_ALGAE_EX_H__
