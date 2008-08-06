#ifndef __UQ_HIMMEL_EX_H__
#define __UQ_HIMMEL_EX_H__

double observationsOfA[] = {
    0.,    0.02090,
    4.50,  0.01540,
    8.67,  0.01422,
   12.67,  0.01335,
   17.75,  0.01232,
   22.67,  0.01181,
   27.08,  0.01139,
   32.00,  0.01092,
   36.00,  0.01054,
   46.33,  0.00978,
   57.00,  0.009157,
   69.00,  0.008594,
   76.75,  0.008395,
   90.00,  0.007891,
  102.00,  0.007510,
  108.00,  0.007370,
  147.92,  0.006646,
  198.00,  0.005883,
  241.75,  0.005322,
  270.25,  0.004960,
  326.25,  0.004518,
  418.00,  0.004075,
  501.00,  0.003715
};

#include <uqDRAM_MarkovChainGenerator.h>
#ifdef __APPL_USES_GSL__
#include <uqGslOdeSolver.h>
int uqGslHimmelExStateDot(double t, const double currentState[], double stateDot[], void* infoForComputingStateDot);
#endif

template<class V, class M>
struct
uqAppl_M2lLikelihoodFunction_DataType
{
  const std::vector<double>* instantsOfAObservations;
  const std::vector<V*>*     observedEvolutionOfAConcentration;
  const V*                   initialConcentrations;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass*       gslOdeSolver;
#endif
};

template<class V, class M>
struct
uqAppl_StateDotFunction_DataType
{
  const V* concentrationRates;
};

template<class V, class M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) std::cout << "Beginning run of 'uqHimmelEx' example\n"
                                 << std::endl;

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //*****************************************************
  // Step 1 of 4: Define the finite dimensional linear spaces.
  //              Define the Markov chain generator.
  //*****************************************************
  uqParamSpaceClass<V,M>  paramSpace (env);
  uqStateSpaceClass<V,M>  stateSpace (env);
  uqOutputSpaceClass<V,M> outputSpace(env);

  uq_M2lLikelihoodFunction_Class<V,M> uq_M2lLikelihoodFunction_Obj;
  uqDRAM_MarkovChainGeneratorClass<V,M> mcg(env,
                                            paramSpace,
                                            outputSpace,
                                            NULL, // use default prior() routine
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
  std::vector<double> instantsOfAObservations          (23,0.);
  std::vector<V*>     observedEvolutionOfAConcentration(23,NULL);
  for (unsigned int i = 0; i < 23; ++i) {
    instantsOfAObservations[i]                 = observationsOfA[2*i+0];
    observedEvolutionOfAConcentration[i]       = stateSpace.newVector();
    (*observedEvolutionOfAConcentration[i])[0] = observationsOfA[2*i+1];
  }

  V* initialConcentrations = stateSpace.newVector();
  (*initialConcentrations)[0] = 0.02090;    // A0
  (*initialConcentrations)[1] = 0.02090/3.; // B0
  (*initialConcentrations)[2] = 0.;         // C0
  (*initialConcentrations)[3] = 0.;         // D0
  (*initialConcentrations)[4] = 0.;         // E0

#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass gslOdeSolver("45",stateSpace.dim());
#endif

  uqAppl_M2lLikelihoodFunction_DataType<V,M> uqAppl_M2lLikelihoodFunction_Data;
  uqAppl_M2lLikelihoodFunction_Data.instantsOfAObservations           = &instantsOfAObservations;
  uqAppl_M2lLikelihoodFunction_Data.observedEvolutionOfAConcentration = &observedEvolutionOfAConcentration;
  uqAppl_M2lLikelihoodFunction_Data.initialConcentrations             = initialConcentrations;
#ifdef __APPL_USES_GSL__
  uqAppl_M2lLikelihoodFunction_Data.gslOdeSolver                      = &gslOdeSolver;
#endif

  //******************************************************
  // Step 4 of 4: Generate chains.
  //              In MATLAB, one should then run uqHimmelEx.m,
  //******************************************************
  mcg.generateChains(NULL, // compute proposal cov matrix internally
                     NULL, // use default prior() routine
                     (void *) &uqAppl_M2lLikelihoodFunction_Data);

  //******************************************************
  // Step 5 of 5: Compute distribution of quantity of interest
  //******************************************************
  //mcg.computeQoiDistribution(NULL, // compute proposal cov matrix internally

  //*****************************************************
  // Release memory before leaving routine.
  //*****************************************************
  for (unsigned int i = 0; i < observedEvolutionOfAConcentration.size(); ++i) {
    if (observedEvolutionOfAConcentration[i]) delete observedEvolutionOfAConcentration[i];
  }
  observedEvolutionOfAConcentration.clear();
  delete initialConcentrations;

  if (env.rank() == 0) std::cout << "Finishing run of 'uqHimmelEx' example"
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
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");
  return 0.;
}

template<class V, class M>
void
uqAppl_M2lLikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  const std::vector<double>& instantsOfAObservations           = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->instantsOfAObservations;
  const std::vector<V*>&     observedEvolutionOfAConcentration = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->observedEvolutionOfAConcentration;
  const V&                   initialConcentrations             = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->initialConcentrations;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass&       gslOdeSolver                      = *((uqAppl_M2lLikelihoodFunction_DataType<V,M> *) functionDataPtr)->gslOdeSolver;
#endif

  uqAppl_StateDotFunction_DataType<V,M> uqAppl_StateDotFunction_Data;
  uqAppl_StateDotFunction_Data.concentrationRates = &paramValues;

  std::vector<V*> computedEvolutionOfConcentrations(0,NULL);
#ifdef __APPL_USES_GSL__
  gslOdeSolver.solveODE(uqGslHimmelExStateDot,
                        NULL,
                        instantsOfAObservations,
                        initialConcentrations,
                        (void *)&uqAppl_StateDotFunction_Data, // concentration rates
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
  //  for (unsigned int i = 0; i < instantsOfAObservations.size(); ++i) {
  //    std::cout << "  for t = " << instantsOfAObservations[i]
  //              << ", state = " << *(computedEvolutionOfConcentrations[i])
  //              << std::endl;
  //  }
  //}

  resultValues[0] = 0.;
  for (unsigned int i = 0; i < instantsOfAObservations.size(); ++i) {
    resultValues[0] += pow( (*computedEvolutionOfConcentrations[i])[0] - (*observedEvolutionOfAConcentration[i])[0],2. );
  }

  for (unsigned int i = 0; i < computedEvolutionOfConcentrations.size(); ++i) {
    if (computedEvolutionOfConcentrations[i]) delete computedEvolutionOfConcentrations[i];
  }
  computedEvolutionOfConcentrations.clear();

  return;
}
#endif // __UQ_HIMMEL_EX_H__
