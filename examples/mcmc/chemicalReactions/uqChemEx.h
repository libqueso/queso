/* uq/examples/mcmc/chemicalReactions/uqChemEx.h
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

#ifndef __UQ_CHEMICAL_REACTIONS_EX_H__
#define __UQ_CHEMICAL_REACTIONS_EX_H__

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
                      "calib_M2lPriorRoutine()",
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
  const V* concentrationRates;
};

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_MisfitLikelihoodRoutine_DataType
{
  const std::vector<double>* instantsOfAObservations;
  const std::vector<V*>*     observedEvolutionOfAConcentration;
  const V*                   initialConcentrations;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass*       gslOdeSolver;
#endif
};

template<class V, class M>
void
calib_MisfitLikelihoodRoutine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering calib_MisfitLikelihoodRoutine()" << std::endl;
  }

  const std::vector<double>& instantsOfAObservations           = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->instantsOfAObservations;
  const std::vector<V*>&     observedEvolutionOfAConcentration = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->observedEvolutionOfAConcentration;
  const V&                   initialConcentrations             = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->initialConcentrations;
#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass&       gslOdeSolver                      = *((calib_MisfitLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->gslOdeSolver;
#endif

  calib_StateDotRoutine_DataType<V,M> calib_StateDotRoutine_Data;
  calib_StateDotRoutine_Data.concentrationRates = &paramValues;

  std::vector<V*> computedEvolutionOfConcentrations(0);//,NULL);
#ifdef __APPL_USES_GSL__
  gslOdeSolver.solveODE(calib_StateDotRoutine_gsl,
                        NULL,
                        instantsOfAObservations,
                        initialConcentrations,
                        (void *)&calib_StateDotRoutine_Data, // concentration rates
                        1.e-1,                                // suggested time step
                        computedEvolutionOfConcentrations);
#endif
#ifdef __APPL_USES_TRILINOS__
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    paramValues.env().rank(),
                    "calib_MisfitLikelihoodRoutine()",
                    "not yet implemented with Trilinos");
#endif
  //if (initialConcentrations.env().rank() == 0) {
  //  std::cout << "In calib_MisfitLikelihoodRoutine()" << std::endl;
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

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving calib_MisfitLikelihoodRoutine()" << std::endl;
  }

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
  if (env.rank() == 0) std::cout << "Beginning run of 'uqChemEx' example\n"
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
  std::vector<double> instantsOfAObservations          (23,0.);
  std::vector<V*>     observedEvolutionOfAConcentration(23);//,NULL);
  for (unsigned int i = 0; i < 23; ++i) {
    instantsOfAObservations[i]                 = observationsOfA[2*i+0];
    observedEvolutionOfAConcentration[i]       = calib_StateSpace.newVector();
    (*observedEvolutionOfAConcentration[i])[0] = observationsOfA[2*i+1];
  }

  V* initialConcentrations = calib_StateSpace.newVector();
  (*initialConcentrations)[0] = 0.02090;    // A0
  (*initialConcentrations)[1] = 0.02090/3.; // B0
  (*initialConcentrations)[2] = 0.;         // C0
  (*initialConcentrations)[3] = 0.;         // D0
  (*initialConcentrations)[4] = 0.;         // E0

#ifdef __APPL_USES_GSL__
  uqGslOdeSolverClass gslOdeSolver("45",calib_StateSpace.dim());
#endif

  calib_MisfitLikelihoodRoutine_DataType<V,M> calib_MisfitLikelihoodRoutine_Data;
  calib_MisfitLikelihoodRoutine_Data.instantsOfAObservations           = &instantsOfAObservations;
  calib_MisfitLikelihoodRoutine_Data.observedEvolutionOfAConcentration = &observedEvolutionOfAConcentration;
  calib_MisfitLikelihoodRoutine_Data.initialConcentrations             = initialConcentrations;
#ifdef __APPL_USES_GSL__
  calib_MisfitLikelihoodRoutine_Data.gslOdeSolver                      = &gslOdeSolver;
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
  uq_QoIPredictionFunction_Class<V,M> calib_QoIPredictionFunction_Obj(calib_QoIPredictionRoutine<V,M>,
                                                                      (void *) &calib_QoIPredictionRoutine_Data);
  uqComputeQoIDistribution(mcg.chain(),
                           mcg.lrVarianceChain(),
                           500,
                           calib_QoIPredictionFunction_Obj,
                           mcg.outputFileName());
  }

  //*****************************************************
  // Release memory before leaving routine.
  //*****************************************************
  for (unsigned int i = 0; i < observedEvolutionOfAConcentration.size(); ++i) {
    if (observedEvolutionOfAConcentration[i]) delete observedEvolutionOfAConcentration[i];
  }
  observedEvolutionOfAConcentration.clear();
  delete initialConcentrations;

  if (env.rank() == 0) std::cout << "Finishing run of 'uqChemEx' example"
                                 << std::endl;

  return;
}
#endif // __UQ_CHEMICAL_REACTIONS_EX_H__
