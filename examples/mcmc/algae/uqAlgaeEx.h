/* uq/examples/mcmc/algae/uqAlgaeEx.h
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

#include <uqStateSpace.h>
#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>
#include <uqComputeQoIDistribution.h>

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> This application is using the default prior() routine provided by the MCMC tool.
//********************************************************
template<class V, class M>
double
calib_M2lPriorRoutine(const V& paramValues, const void* functionDataPtr)
{
  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_M2lPriorRoutine(), in uqAlgaeEx.h",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");
  return 0.;
}

//********************************************************
// The stateDot routine: provided by user and called by the likelihood routine.
//********************************************************
#ifdef __APPL_USES_GSL__
#include <uqGslOdeSolver.h>
int calib_StateDotRoutine_gsl(double t, const double currentState[], double stateDot[], void* infoForComputingStateDot);
#endif

template<class V, class M>
struct
calib_StateDotRoutine_DataType
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

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_MisfitLikelihoodRoutine_DataType
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
void
calib_StateEvolutionRoutine(const V& paramValues, const void* functionDataPtr, std::vector<V*>& states)
{
  const std::vector<double>& instantsOfAZPObservations            = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->instantsOfAZPObservations;
  const std::vector<V*>&     observedEvolutionOfAZPConcentrations = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->observedEvolutionOfAZPConcentrations;
  const std::vector<double>& evolutionOfQpV                       = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->evolutionOfQpV;
  const std::vector<double>& evolutionOfT                         = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->evolutionOfT;
  const std::vector<double>& evolutionOfPin                       = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->evolutionOfPin;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass&       gslOdeSolver                         = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->gslOdeSolver;
#endif

  // The initial concentrations are also unknown and treated as extra parameters
  V initialConcentrations(*observedEvolutionOfAZPConcentrations[0]);
  initialConcentrations[0] = paramValues[6];
  initialConcentrations[1] = paramValues[7];
  initialConcentrations[2] = paramValues[8];

  calib_StateDotRoutine_DataType<V,M> calib_StateDotRoutine_Data;
  calib_StateDotRoutine_Data.muMax          = paramValues[0];
  calib_StateDotRoutine_Data.rhoA           = paramValues[1];
  calib_StateDotRoutine_Data.rhoZ           = paramValues[2];
  calib_StateDotRoutine_Data.k              = paramValues[3];
  calib_StateDotRoutine_Data.alpha          = paramValues[4];
  calib_StateDotRoutine_Data.th             = paramValues[5];
  calib_StateDotRoutine_Data.evolutionOfQpV = &evolutionOfQpV;
  calib_StateDotRoutine_Data.evolutionOfT   = &evolutionOfT;
  calib_StateDotRoutine_Data.evolutionOfPin = &evolutionOfPin;

  //if (initialConcentrations.env().rank() == 0) {
  //  std::cout << "In calib_StateEvolutionRoutine():"
  //            << "\n A0    = " << initialConcentrations[0]
  //            << "\n Z0    = " << initialConcentrations[1]
  //            << "\n P0    = " << initialConcentrations[2]
  //            << "\n muMax = " << calib_StateDotRoutine_Data.muMax
  //            << "\n rhoA  = " << calib_StateDotRoutine_Data.rhoA
  //            << "\n rhoZ  = " << calib_StateDotRoutine_Data.rhoZ
  //            << "\n k     = " << calib_StateDotRoutine_Data.k
  //            << "\n alpha = " << calib_StateDotRoutine_Data.alpha
  //            << "\n th    = " << calib_StateDotRoutine_Data.th
  //            << "\n"
  //            << std::endl;
    //std::cout << "In calib_StateEvolutionRoutine():"
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

#ifdef __APPL_USES_GSL__
  gslOdeSolver.solveODE(calib_StateDotRoutine_gsl,
                        NULL,
                        instantsOfAZPObservations,
                        initialConcentrations,
                        (void *)&calib_StateDotRoutine_Data, // data for stateDot computations
                        1.e-1,                                 // suggested time step
                        states);
#endif
#ifdef __APPL_USES_TRILINOS__
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    paramValues.env().rank(),
                    "calib_StateEvolutionRoutine()",
                    "not yet implemented with Trilinos");
#endif
  //if (initialConcentrations.env().rank() == 0) {
  //  std::cout << "In calib_StateEvolutionRoutine()" << std::endl;
  //  for (unsigned int i = 0; i < instantsOfAZPObservations.size(); ++i) {
  //    std::cout << "  for t = " << instantsOfAZPObservations[i]
  //              << ", state = " << *(states[i])
  //              << std::endl;
  //  }
  //}

  return;
}

template<class V, class M>
void
calib_MisfitLikelihoodRoutine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  std::vector<V*> computedEvolutionOfConcentrations(0);//,NULL);
  calib_StateEvolutionRoutine<V,M>(paramValues, functionDataPtr, computedEvolutionOfConcentrations);

  const std::vector<double>& instantsOfAZPObservations = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->instantsOfAZPObservations;
  const std::vector<V*>&     observedEvolutionOfAZPConcentrations = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->observedEvolutionOfAZPConcentrations;
  resultValues[0] = 0.;
  resultValues[1] = 0.;
  resultValues[2] = 0.;
  for (unsigned int i = 0; i < instantsOfAZPObservations.size(); ++i) {
    resultValues[0] += pow( (*computedEvolutionOfConcentrations[i])[0] - (*observedEvolutionOfAZPConcentrations[i])[0],2. );
    resultValues[1] += pow( (*computedEvolutionOfConcentrations[i])[1] - (*observedEvolutionOfAZPConcentrations[i])[1],2. );
    resultValues[2] += pow( (*computedEvolutionOfConcentrations[i])[2] - (*observedEvolutionOfAZPConcentrations[i])[2],2. );
  }

  for (unsigned int i = 0; i < computedEvolutionOfConcentrations.size(); ++i) {
    if (computedEvolutionOfConcentrations[i]) delete computedEvolutionOfConcentrations[i];
  }
  computedEvolutionOfConcentrations.clear();

  return;
}

//********************************************************
// The QoI prediction routine: provided by user and called by the prediction tool.
//********************************************************
template<class V, class M>
struct
calib_QoIPredictionRoutine_DataType
{
  const std::vector<void*>* nothingYet;
};

template<class V, class M>
void
calib_QoIPredictionRoutine(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictions)
{
  // Nothing yet

  return;
}

//********************************************************
// The MCMC driving routine: called by main()
//********************************************************
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
  // Step 1 of 6: Define the finite dimensional linear spaces.
  //*****************************************************
  uqParamSpaceClass<V,M>      calib_ParamSpace     (env,"calib");
  uqStateSpaceClass<V,M>      calib_StateSpace     (env,"calib");
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
  // The initial concentrations are also treated as parameters
  std::vector<double> instantsOfAZPObservations           (30,0.);
  std::vector<V*>     observedEvolutionOfAZPConcentrations(30);//,NULL);
  for (unsigned int i = 0; i < 30; ++i) {
    instantsOfAZPObservations[i]                  = observationsOfAZP[4*i+0];
    observedEvolutionOfAZPConcentrations[i]       = calib_StateSpace.newVector();
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
  uqGslOdeSolverClass gslOdeSolver("gear2",calib_StateSpace.dim());
#endif

  calib_MisfitLikelihoodRoutine_DataType<V,M> calib_MisfitLikelihoodRoutine_Data;
  calib_MisfitLikelihoodRoutine_Data.instantsOfAZPObservations            = &instantsOfAZPObservations;
  calib_MisfitLikelihoodRoutine_Data.observedEvolutionOfAZPConcentrations = &observedEvolutionOfAZPConcentrations;
  calib_MisfitLikelihoodRoutine_Data.evolutionOfQpV                       = &evolutionOfQpV;
  calib_MisfitLikelihoodRoutine_Data.evolutionOfT                         = &evolutionOfT;
  calib_MisfitLikelihoodRoutine_Data.evolutionOfPin                       = &evolutionOfPin;
#ifdef __APPL_USES_GSL__
  calib_MisfitLikelihoodRoutine_Data.gslOdeSolver                         = &gslOdeSolver;
#endif

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

  //******************************************************
  // Step 6 of 6: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //******************************************************
  mcg.generateChains(NULL); // compute proposal cov matrix internally

  //******************************************************
  // Step 7 of 6: Compute distribution of quantity of interest
  //******************************************************
  if (0) {
  calib_QoIPredictionRoutine_DataType<V,M> calib_QoIPredictionRoutine_Data;
  uq_QoIPredictionFunction_Class<V,M>       calib_QoIPredictionFunction_Obj(calib_QoIPredictionRoutine<V,M>,
                                                                            (void *) &calib_QoIPredictionRoutine_Data);
  //uqComputeQoIDistribution(mcg.chain(),
  //                         mcg.misfitVarianceChain(),
  //                         500,
  //                         calib_QoIPredictionFunction_Obj,
  //                         mcg.outputFileName());
  }

  //*****************************************************
  // Release memory before leaving routine.
  //*****************************************************
  for (unsigned int i = 0; i < observedEvolutionOfAZPConcentrations.size(); ++i) {
    if (observedEvolutionOfAZPConcentrations[i]) delete observedEvolutionOfAZPConcentrations[i];
  }
  observedEvolutionOfAZPConcentrations.clear();

  if (env.rank() == 0) std::cout << "Finishing run of 'uqAlgaeEx' example"
                                 << std::endl;

  return;
}
#endif // __UQ_ALGAE_EX_H__
