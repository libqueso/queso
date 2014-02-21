#include "queso/StatisticalInverseProblem.h"
#include "queso/StatisticalInverseProblemOptions.h"
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"
#include "queso/Environment.h"
#include "queso/EnvironmentOptions.h"
#include "queso/Defines.h"
#include "queso/VectorSpace.h" 
#include "queso/JointPdf.h" 
#include "queso/GpmsaComputerModel.h"
#include "queso/SimulationModelOptions.h"

int main(int argc, char ** argv)
{
  ////////////////////////////////////////////////////////
  // Step 0: Instantiate the QUESO environment 
  ////////////////////////////////////////////////////////
  QUESO::EnvOptionsValues envOptionsValues;
  envOptionsValues.m_subDisplayFileName = "outputData/display";
  envOptionsValues.m_subDisplayAllowedSet.insert(0);
  envOptionsValues.m_subDisplayAllowedSet.insert(1);
  envOptionsValues.m_displayVerbosity = 3;
  envOptionsValues.m_seed = 1;
      
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);
 
  unsigned int numExperiments = 5;
  unsigned int numUncertainVars = 1;

  //***********************************************************************
  //  Step 01 of 09: Instantiate parameter space, parameter domain, and prior Rv
  // *********************************************************************
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    paramSpace(env, "param_", numUncertainVars, NULL);

  std::cout << "Got here 1" << std::endl;;

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  for (unsigned int i = 0; i < numUncertainVars; i++) {
    paramMins[i]=0.0;
    paramMaxs[i]=1.0;
  }
 
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  std::cout << "got here 2" << std::endl;;

#if 0
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>* paramPriorRvPtr = NULL;
  paramPriorRvPtr = new QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>("prior_",paramDomain);
#else
  QUESO::GslVector meanVec(paramSpace.zeroVector()); // prudenci_new_2013_09_06
  meanVec[0] = 0.5; // prudenci_new_2013_09_06
  meanVec[1] = 0.5; // prudenci_new_2013_09_06
  meanVec[2] = 0.5; // prudenci_new_2013_09_06
  QUESO::GslMatrix covMat(paramSpace.zeroVector()); // prudenci_new_2013_09_06
  covMat(0,0) = 10.*10.; // prudenci_new_2013_09_06
  covMat(1,1) = 10.*10.; // prudenci_new_2013_09_06
  covMat(2,2) = 10.*10.; // prudenci_new_2013_09_06

  QUESO::GaussianVectorRV<QUESO::GslVector, QUESO::GslMatrix>
    paramPriorRv("prior_", paramDomain, meanVec, covMat);
#endif
 

  std::cout << "Got here 3" << std::endl;

  //***********************************************************************
  // Step 02 of 09: Instantiate the 'scenario' and 'output' spaces
  //***********************************************************************

  // The only configuration variable is the height
  int num_config_vars = 1;
  std::cout << "num_config_vars " << num_config_vars << '\n';

  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> config_x_space(env, "scenario_", num_config_vars, NULL);
               
  // Number of responses the model is returning
  unsigned int num_eta = 1;

  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    n_eta_space(env, "output_", num_eta, NULL);

  std::cout << "num_eta " << num_eta << '\n';
  
  //***********************************************************************
  //  Step 03 of 09: Instantiate the simulation storage
  //  Regarding simulation scenario input values, QUESO will standardize them so that
  //  they exist inside a hypercube: this will be done in the 'SimulationModel'
  //  constructor, step 04 of 09.
  //  Regarding simulation output data, QUESO will transform it so that the mean is
  //       zero and the variance is one: this will be done in the 'SimulationModel'
  //        constructor, step 04 of 09.
  // ***********************************************************************
 
  // The number of expensive simulator runs we have: 25
  unsigned int num_simulations = 25;
  std::cout << "num_simulations " << num_simulations << '\n';

  // Simulation storage stores all the information about our simulation model
  // and simulations.  'Simulation model' is the expensive model
  QUESO::SimulationStorage<QUESO::GslVector, QUESO::GslMatrix,
    QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>
      simulationStorage(config_x_space,
                        paramSpace,
                        n_eta_space,
                        num_simulations);

 
  // Add simulations: what if none?
  // We have to dimension these properly. 
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    n_grid(env, "output_", num_eta, NULL);

  // std::vector containing all the points in scenario space where we have
  // simulations
  std::vector<QUESO::GslVector *> simulationScenarios(num_simulations, (QUESO::GslVector *) NULL);

  // std::vector containing all the points in parameter space where we have
  // simulations
  std::vector<QUESO::GslVector *> paramVecs(num_simulations, (QUESO::GslVector *) NULL);

  // std::vector containing all the points in scenario space where we have
  // simulations
  std::vector<QUESO::GslVector *> outputVecs(num_simulations, (QUESO::GslVector *) NULL);
            
  for (unsigned int i = 0; i < num_simulations; i++) {
    simulationScenarios[i] = new QUESO::GslVector(config_x_space.zeroVector());   // 'x_{i+1}^*' in paper
    paramVecs          [i] = new QUESO::GslVector(paramSpace.zeroVector());   // 't_{i+1}^*' in paper
    outputVecs         [i] = new QUESO::GslVector(n_eta_space.zeroVector()); // 'eta_{i+1}' in paper
  }
                           
  // All the positions in scenario space where simulations were run
  (*(simulationScenarios[0]))[0] = 0.0416667;
  (*(simulationScenarios[1]))[0] = 0.125;
  (*(simulationScenarios[2]))[0] = 0;
  (*(simulationScenarios[3]))[0] = 0.166667;
  (*(simulationScenarios[4]))[0] = 0.0833333;
  (*(simulationScenarios[5]))[0] = 0.333333;
  (*(simulationScenarios[6]))[0] = 0.208333;
  (*(simulationScenarios[7]))[0] = 0.291667;
  (*(simulationScenarios[8]))[0] = 0.375;
  (*(simulationScenarios[9]))[0] = 0.25;
  (*(simulationScenarios[10]))[0] = 0.583333;
  (*(simulationScenarios[11]))[0] = 0.416667;
  (*(simulationScenarios[12]))[0] = 0.541667;
  (*(simulationScenarios[13]))[0] = 0.5;
  (*(simulationScenarios[14]))[0] = 0.458333;
  (*(simulationScenarios[15]))[0] = 0.708333;
  (*(simulationScenarios[16]))[0] = 0.75;
  (*(simulationScenarios[17]))[0] = 0.666667;
  (*(simulationScenarios[18]))[0] = 0.625;
  (*(simulationScenarios[19]))[0] = 0.791667;
  (*(simulationScenarios[20]))[0] = 0.875;
  (*(simulationScenarios[21]))[0] = 1;
  (*(simulationScenarios[22]))[0] = 0.833333;
  (*(simulationScenarios[23]))[0] = 0.958333;
  (*(simulationScenarios[24]))[0] = 0.916667;

  // All the positions in parameter space where simulations were run
  (*(paramVecs[0]))[0] = 0.125;
  (*(paramVecs[1]))[0] = 0.375;
  (*(paramVecs[2]))[0] = 0.583333;
  (*(paramVecs[3]))[0] = 0.708333;
  (*(paramVecs[4]))[0] = 0.875;
  (*(paramVecs[5]))[0] = 0.0;
  (*(paramVecs[6]))[0] = 0.208333;
  (*(paramVecs[7]))[0] = 0.5;
  (*(paramVecs[8]))[0] = 0.791667;
  (*(paramVecs[9]))[0] = 0.916667;
  (*(paramVecs[10]))[0] = 0.0416667;
  (*(paramVecs[11]))[0] = 0.291667;
  (*(paramVecs[12]))[0] = 0.416667;
  (*(paramVecs[13]))[0] = 0.625;
  (*(paramVecs[14]))[0] = 1.0;
  (*(paramVecs[15]))[0] = 0.166667;
  (*(paramVecs[16]))[0] = 0.333333;
  (*(paramVecs[17]))[0] = 0.541667;
  (*(paramVecs[18]))[0] = 0.75;
  (*(paramVecs[19]))[0] = 0.833333;
  (*(paramVecs[20]))[0] = 0.0833333;
  (*(paramVecs[21]))[0] = 0.25;
  (*(paramVecs[22]))[0] = 0.458333;
  (*(paramVecs[23]))[0] = 0.666667;
  (*(paramVecs[24]))[0] = 0.958333;

  // Simulation output from sim_outputs file from matlab gpmsa implementation
  (*(outputVecs[0]))[0] = 1.60656419103786;
  (*(outputVecs[1]))[0] = 2.03091732612410;
  (*(outputVecs[2]))[0] = 1.63843806135107;
  (*(outputVecs[3]))[0] = 2.4638773088846;
  (*(outputVecs[4]))[0] = 2.21838509163166;
  (*(outputVecs[5]))[0] = 2.30939315250214;
  (*(outputVecs[6]))[0] = 2.17300595474166;
  (*(outputVecs[7]))[0] = 2.71440886744846;
  (*(outputVecs[8]))[0] = 3.45738385598637;
  (*(outputVecs[9]))[0] = 3.08913132945656;
  (*(outputVecs[10]))[0] = 2.8542802908597;
  (*(outputVecs[11]))[0] = 2.82284112252142;
  (*(outputVecs[12]))[0] = 3.35111130604648;
  (*(outputVecs[13]))[0] = 3.62106706299426;
  (*(outputVecs[14]))[0] = 4.27778530159255;
  (*(outputVecs[15]))[0] = 3.28874614306084;
  (*(outputVecs[16]))[0] = 3.71321395867495;
  (*(outputVecs[17]))[0] = 3.97952570579624;
  (*(outputVecs[18]))[0] = 4.37705002709939;
  (*(outputVecs[19]))[0] = 5.31766210668879;
  (*(outputVecs[20]))[0] = 3.43505963946617;
  (*(outputVecs[21]))[0] = 4.04494232394982;
  (*(outputVecs[22]))[0] = 4.23694028613074;
  (*(outputVecs[23]))[0] = 5.31781675076092;
  (*(outputVecs[24]))[0] = 6.39189142027432;

  // Finally, add information to the simulation storage
  for (unsigned int i = 0; i < num_simulations; i++) {
    simulationStorage.addSimulation(*(simulationScenarios[i]),
        *(paramVecs[i]), *(outputVecs[i]));
  }

  std::cout << "got here 4 " << std::endl;

  //***********************************************************************
  // Step 04 of 09: Instantiate the simulation model
  //***********************************************************************
  QUESO::SmOptionsValues *smVarOptions = NULL;
  smVarOptions = new QUESO::SmOptionsValues();
  smVarOptions->m_dataOutputFileName = ".";
  smVarOptions->m_dataOutputAllowAll = 0;
  smVarOptions->m_dataOutputAllowedSet.insert(0);

  // 1 GP (scalar use-case)
  smVarOptions->m_p_eta = 1;

  smVarOptions->m_zeroRelativeSingularValue = 0.;
  smVarOptions->m_cdfThresholdForPEta = 0.;
  //simOptionsValues->m_ov.cdfThresholdForPEta = 0.95;

  // DM: These should be set properly for the tower problem.
  smVarOptions->m_a_w = 5.;
  smVarOptions->m_b_w = 5.;
  smVarOptions->m_a_rho_w = 1.;
  smVarOptions->m_b_rho_w = 0.1;
  smVarOptions->m_a_eta = 5.;
  smVarOptions->m_b_eta = 0.005;
  smVarOptions->m_a_s = 1.;
  smVarOptions->m_b_s = 0.0001;

  QUESO::SimulationModel<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector,
    QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>
      simulationModel("", smVarOptions, simulationStorage);

  // number of simulation basis; 'p_eta' in paper; 1 in tower example
  unsigned int num_bases_eta = simulationModel.numBasis();

  std::cout << "num_bases_eta " << num_bases_eta << std::endl;
  std::cout << "got here 5" << std::endl;

  //***********************************************************************
  // Step 05 of 09: Instantiate the experiment storage
  // Regarding experimental scenario input values, QUESO will standardize them so that
  //    they exist inside a hypercube: this will be done in the 'ExperimentModel'
  //    constructor, step 06 of 09.
  // Regarding experimental data, the user has to provide it already in transformed
  //    format, that is, with mean zero and standard deviation one.
  //***********************************************************************

  // number of experiments; 'n' in paper; 6 in tower example, numExperiments

  QUESO::ExperimentStorage<QUESO::GslVector, QUESO::GslMatrix,
    QUESO::GslVector, QUESO::GslMatrix> experimentStorage(config_x_space,
        numExperiments);
 
  std::cout << "Got here 6" << std::endl;

  // Add experiments
  std::vector<QUESO::GslVector *> experimentScenarios_original(numExperiments,
      (QUESO::GslVector *) NULL);

  std::vector<QUESO::GslVector *> experimentScenarios_standard(numExperiments,
      (QUESO::GslVector *) NULL);

  std::vector<unsigned int> experimentDims(numExperiments, 0);
  std::vector<unsigned int> experimentDimsEach(numExperiments, 0);

  std::vector<QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> * >
    experimentSpaces(numExperiments, (QUESO::VectorSpace<QUESO::GslVector,
          QUESO::GslMatrix> *) NULL);

  std::vector<QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> * >
    eachResponseSpace(numExperiments, (QUESO::VectorSpace<QUESO::GslVector,
          QUESO::GslMatrix> *) NULL);

  std::vector<QUESO::GslVector *> experimentVecs_original(numExperiments,
      (QUESO::GslVector *) NULL);

  std::vector<QUESO::GslVector *> experimentVecs_auxMean(numExperiments,
      (QUESO::GslVector *) NULL);

  std::vector<QUESO::GslVector *> experimentVecs_transformed(numExperiments,
      (QUESO::GslVector *) NULL);

  std::vector<QUESO::GslMatrix *> experimentMats_original(numExperiments,
      (QUESO::GslMatrix *) NULL);

  std::vector<QUESO::GslMatrix *> experimentMats_transformed(numExperiments,
      (QUESO::GslMatrix *) NULL);

  std::vector<QUESO::GslMatrix *>
    experimentMats_transformed_inv(numExperiments, (QUESO::GslMatrix *) NULL);

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentDims[i] = 1;
    experimentDimsEach[i] = 1;
  }

  std::cout << "Got here 7" << std::endl;
 
  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios_original[i] =
      new QUESO::GslVector(config_x_space.zeroVector()); // 'x_{i+1}' in paper

    experimentScenarios_standard[i] =
      new QUESO::GslVector(config_x_space.zeroVector());

    experimentSpaces[i] =
      new QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>(env,
          "expSpace", experimentDims[i], NULL);

    eachResponseSpace[i] =
      new QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>(env,
          "expSpace", experimentDimsEach[i], NULL);

    experimentVecs_original[i] =
      new QUESO::GslVector(experimentSpaces[i]->zeroVector());

    experimentVecs_auxMean[i] =
      new QUESO::GslVector(experimentSpaces[i]->zeroVector());

    experimentVecs_transformed[i] =
      new QUESO::GslVector(experimentSpaces[i]->zeroVector()); // 'y_{i+1}' in paper

    experimentMats_original   [i] =
      new QUESO::GslMatrix(experimentSpaces[i]->zeroVector());

    experimentMats_transformed[i] =
      new QUESO::GslMatrix(experimentSpaces[i]->zeroVector()); // 'W_{i+1}' in paper

    experimentMats_transformed_inv[i] =
      new QUESO::GslMatrix(experimentSpaces[i]->zeroVector()); // 'W_{i+1}' in paper
  }

  (*(experimentScenarios_original[0]))[0] = 0.0;
  (*(experimentScenarios_original[1]))[0] = 0.2;
  (*(experimentScenarios_original[2]))[0] = 0.4;
  (*(experimentScenarios_original[3]))[0] = 0.6;
  (*(experimentScenarios_original[4]))[0] = 0.8;
  (*(experimentScenarios_original[5]))[0] = 1.0;

  // These are the expeimental observations
  (*(experimentVecs_original[0]))[0] = 1.59294826065229;
  (*(experimentVecs_original[1]))[0] = 2.17696977275016;
  (*(experimentVecs_original[2]))[0] = 2.87061286591332;
  (*(experimentVecs_original[3]))[0] = 3.8330395105599;
  (*(experimentVecs_original[4]))[0] = 4.59654198432239;
  (*(experimentVecs_original[5]))[0] = 4.75087857533489;

  for (unsigned int i = 0; i < numExperiments; i++) {
    std::cout << *(experimentScenarios_original[i]) << std::endl;
    *(experimentScenarios_standard[i])  = *(experimentScenarios_original[i]);
    std::cout << *(experimentVecs_original[i]) << std::endl;
  }

  std::cout << "Got here 8" << std::endl;

  // We do this because GslVector doesn't support operator-(GslVector, double)
  double all_sim_means = simulationModel.etaSeq_allMean();
  QUESO::GslVector all_sim_means_vec(experimentSpaces[0]->zeroVector());
  all_sim_means_vec.cwSet(all_sim_means);

  // We do this because GslVector doesn't support operator-(GslVector, double)
  double all_sim_stds = simulationModel.etaSeq_allStd();
  QUESO::GslVector all_sim_stds_vec(experimentSpaces[0]->zeroVector());
  all_sim_stds_vec.cwSet(all_sim_stds);

  for (unsigned int i = 0; i < numExperiments; i++) {
    // DM: Why obs standardised with statistics from the simulation output??
    // 'y_1' in paper
    *(experimentVecs_transformed[i]) = (*(experimentVecs_original[i]) -
        all_sim_means_vec) / all_sim_stds_vec;
  }

  std::cout << "Got here 9" << std::endl;

  for (unsigned int i = 0; i < numExperiments; i++) {
    for (unsigned int j = 0; j < experimentDims[i] ; j++) {
      (*(experimentMats_original[i]))(j,j) = 0.0001;
    }
  }

  for (unsigned int i = 0; i < numExperiments; i++) {
    for (unsigned int j = 0; j < experimentDims[i] ; j++) {
      (*(experimentMats_transformed[i]))(j,j) =
        (*(experimentMats_original[i]))(j,j) / simulationModel.etaSeq_allStd();
    }
  }

  for (unsigned int i = 0; i < numExperiments; i++) {
    for (unsigned int j = 0; j < experimentDims[i] ; j++) {
      (*(experimentMats_transformed_inv[i]))(j,j) =
        1.0 / (*(experimentMats_transformed[i]))(j,j);                                                                            
    }
  }

  //***********************************************************************
  // Finally, add information to the experiment storage
  //***********************************************************************
  for (unsigned int i = 0; i < numExperiments; i++) {
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
      *env.subDisplayFile() << "In compute(), step 05"
                            << ": calling experimentStorage.addExperiment() for experiment of id '" << i << "'..."
                            << std::endl;
    }

    experimentStorage.addExperiment(*(experimentScenarios_standard[i]),
        *(experimentVecs_transformed[i]),
        *(experimentMats_transformed_inv[i]));

    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
      *env.subDisplayFile() << "In compute(), step 05"
                            << ": returned from experimentStorage.addExperiment() for experiment of id '" << i << "'"
                            << std::endl;
    }
  }
  std::cout << "Got here 10" << std::endl;

  //***********************************************************************
  // Step 06 of 09: Instantiate the experiment model
  // User has to provide 'D' matrices for discrepencies
  // User has to interpolate 'K_eta' matrix in order to form 'K_i' matrices
  // 'K_eta' is 'Ksim' in the GPMSA tower example document (page 9)
  //  'K_i' is 'Kobs' in the GPMSA tower example document (page 9)
  //***********************************************************************

  // number of experiment basis; 'p_delta' in paper; 1 in tower example
  unsigned int num_delta_bases = 1;

  double kernelSigma = 0.1;
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> oneDSpace(env,
      "oneDSpace", 1, NULL);

  QUESO::GslVector oneDVec(oneDSpace.zeroVector());
  QUESO::GslMatrix oneDMat(oneDSpace.zeroVector());

  oneDMat(0,0) = kernelSigma * kernelSigma;

  //***********************************************************************
  // Form and compute 'DsimMat'
  // Not mentioned in the paper
  // 'Dsim' in the GPMSA tower example document (page 11)
  //***********************************************************************
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> delta_space(env,
      "delta_space", num_delta_bases, NULL);

  QUESO::GslVector kernelSigmas(delta_space.zeroVector());
  kernelSigmas[0] = 1.0;

  QUESO::GslVector kernelCenters(delta_space.zeroVector());
  kernelCenters[0] = 0.0;

  std::cout << "Before we do 2nd subgroup of DsimMat" << '\n';

  std::cout << "Got here 11" << std::endl;

  std::cout << "Got here 12" << std::endl;
  
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "In compute()"
                           << ": finished computing 'K_i' matrices"
                           << std::endl;
  }

  QUESO::EmOptionsValues emVarOptions;
  emVarOptions.m_Gvalues.resize(1,0.);
  emVarOptions.m_Gvalues[0] = 14;
  emVarOptions.m_a_v = 1.;
  emVarOptions.m_b_v = 0.0001;
  emVarOptions.m_a_rho_v = 1.;
  emVarOptions.m_b_rho_v = 0.1;
  emVarOptions.m_a_y = 1000.;
  emVarOptions.m_b_y = 1000.;

  QUESO::ExperimentModel<QUESO::GslVector, QUESO::GslMatrix,
                         QUESO::GslVector, QUESO::GslMatrix>
                           experimentModel("",  // prefix
                           &emVarOptions,
                           experimentStorage);

  // experimentModel = new QUESO::ExperimentModel<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>("",   // prefix
  //                   emVarOptions, // options values
  //                   experimentStorage,
  //                   DobsMats,
  //                   Kmats_interp);

  //***********************************************************************
  // Step 07 of 09: Instantiate the GPMSA computer model
  //***********************************************************************
 
  QUESO::GcmOptionsValues gcmVarOptions;
  gcmVarOptions.m_checkAgainstPreviousSample = 0;
  gcmVarOptions.m_dataOutputFileName = ".";
  gcmVarOptions.m_dataOutputAllowAll = 0;
  gcmVarOptions.m_dataOutputAllowedSet.insert(0);
  gcmVarOptions.m_dataOutputAllowedSet.insert(1);
  gcmVarOptions.m_priorSeqNumSamples = 0;
  gcmVarOptions.m_priorSeqDataOutputFileName = "outputData/priorSeq";
  gcmVarOptions.m_priorSeqDataOutputFileType = "m";
  gcmVarOptions.m_priorSeqDataOutputAllowAll = 0;
  gcmVarOptions.m_priorSeqDataOutputAllowedSet.insert(0);
  gcmVarOptions.m_nuggetValueForBtWyB = 0.0001;       // IMPORTANT
  gcmVarOptions.m_nuggetValueForBtWyBInv = 1e-06;     // IMPORTANT
  gcmVarOptions.m_formCMatrix = 0;
  gcmVarOptions.m_useTildeLogicForRankDefficientC = 0;
  gcmVarOptions.m_predVUsBySamplingRVs = 0;
  gcmVarOptions.m_predVUsBySummingRVs = 0;
  gcmVarOptions.m_predVUsAtKeyPoints = 1;
  gcmVarOptions.m_predWsBySamplingRVs = 0;
  gcmVarOptions.m_predWsBySummingRVs = 1;
  gcmVarOptions.m_predWsAtKeyPoints = 0;
  gcmVarOptions.m_predLag = 15;

  QUESO::GpmsaComputerModel<QUESO::GslVector, QUESO::GslMatrix,
                            QUESO::GslVector, QUESO::GslMatrix,
                            QUESO::GslVector, QUESO::GslMatrix,
                            QUESO::GslVector, QUESO::GslMatrix>
                            gcm("",
                                &gcmVarOptions,
                                simulationStorage,
                                simulationModel,
                                &experimentStorage, // pass "NULL" if there is no experimental data available
                                &experimentModel,   // pass "NULL" if there is no experimental data available
                                &paramPriorRv);

   //***********************************************************************
   // Step 08 of 09: Calibrate the computer model
   //***********************************************************************
  QUESO::GslVector totalInitialVec(gcm.totalSpace().zeroVector());
  gcm.totalPriorRv().realizer().realization(totalInitialVec);
  // totalInitialVec[ 0] = 126.11626143185;  // lambda_eta = lamWOs
  // totalInitialVec[ 1] = 0.701343955862486;  // lambda_w_1 = lamUz
  // totalInitialVec[ 2] = 0.708829052004483;  // lambda_w_2 =
  // totalInitialVec[ 3] = 0.709682686862828;  // lambda_w_3 =
  // totalInitialVec[ 4] = 0.997799643317481;  // rho_w_1_1 = exp(-model.betaU.*(0.5^2));
  // totalInitialVec[ 5] = 0.751419041988001;  // rho_w_1_2 =
  // totalInitialVec[ 6] = 0.810538180186257;  // rho_w_1_3 =
  // totalInitialVec[ 7] = 0.998404593547895;  // rho_w_1_4 =
  // totalInitialVec[ 8] = 0.999507021484067;  // rho_w_2_1 =
  // totalInitialVec[ 9] = 0.00715444906006562;  // rho_w_2_2 =
  // totalInitialVec[10] = 0.000231553394124439;  // rho_w_2_3 =
  // totalInitialVec[11] = 0.388519721938543;  // rho_w_2_4 =
  // totalInitialVec[12] = 0.999482367093646;  // rho_w_3_1 =
  // totalInitialVec[13] = 0.000666947922032456;  // rho_w_3_2 =
  // totalInitialVec[14] = 9.18617635486179e-07;  // rho_w_3_3 =
  // totalInitialVec[15] = 0.156278430012145;  // rho_w_3_4 =
  // totalInitialVec[16] =  221.39241555032;  // lambda_s_1 = lamWs
  // totalInitialVec[17] = 1956.58311200054;  // lambda_s_2 =
  // totalInitialVec[18] = 8843.41472347505;  // lambda_s_2 =
  // totalInitialVec[19] = 0.71728115298871;  // lambda_y = lamOs
  // totalInitialVec[20] = 0.14404321630243;  // lambda_v_1 = lamVz
  // totalInitialVec[21] = 0.999908091519677;  // rho_v_1_1 = betaV prudenci 2013-08-24
  // totalInitialVec[22] = 0.000429516943004888;  // theta_1 = theta
  // totalInitialVec[23] = 0.0101666957753953;  // theta_2 = theta
  // totalInitialVec[24] = 0.0601876182547638;  // theta_3 = theta
 
  totalInitialVec[0] = 6.11;  // \theta
  totalInitialVec[1] = 0.99;  // \mu
  totalInitialVec[2] = 1.39;  // \lambda_\eta
  totalInitialVec[3] = 0.71;  // \rho_\eta_1
  totalInitialVec[4] = 0.14;  // \rho_\eta_2
  totalInitialVec[5] = 0.99;  // \lambda_\delta
  totalInitialVec[6] = 0.00;  // \rho_\delta

  QUESO::GslVector diagVec(gcm.totalSpace().zeroVector());
  diagVec.cwSet(0.01);

  QUESO::GslMatrix totalInitialProposalCovMatrix(diagVec);
  std::cout << "Got here 13" << std::endl; 

  unsigned int numSamples = 1000;
  
  QUESO::MhOptionsValues mhVarOptions;
  mhVarOptions.m_dataOutputFileName = ".";
  mhVarOptions.m_dataOutputAllowAll = 0;
  mhVarOptions.m_dataOutputAllowedSet.insert(0);
  mhVarOptions.m_dataOutputAllowedSet.insert(1);
  mhVarOptions.m_totallyMute = 0;
  mhVarOptions.m_initialPositionDataInputFileName = ".";
  mhVarOptions.m_initialPositionDataInputFileType = "m";
  mhVarOptions.m_initialProposalCovMatrixDataInputFileName = ".";
  mhVarOptions.m_initialProposalCovMatrixDataInputFileType = "m";
  mhVarOptions.m_rawChainDataInputFileName = ".";
  mhVarOptions.m_rawChainDataInputFileType = "m";
  mhVarOptions.m_rawChainSize = numSamples;                                  
   // IMPORTANT
  mhVarOptions.m_rawChainGenerateExtra = 0;
  mhVarOptions.m_rawChainDisplayPeriod = 1000;
  mhVarOptions.m_rawChainMeasureRunTimes = 1;
  mhVarOptions.m_rawChainDataOutputPeriod = 1000;
  mhVarOptions.m_rawChainDataOutputFileName = "outputData/rawChain_mh";
   // IMPORTANT
  mhVarOptions.m_rawChainDataOutputFileType = "m";                     
   // IMPORTANT
  mhVarOptions.m_rawChainDataOutputAllowAll = 0;                       
   // IMPORTANT
  mhVarOptions.m_rawChainDataOutputAllowedSet.insert(0);               
   // IMPORTANT
  mhVarOptions.m_filteredChainGenerate = 1;                            
   // IMPORTANT
  mhVarOptions.m_filteredChainDiscardedPortion = 0.1;                 
   // IMPORTANT
  mhVarOptions.m_filteredChainLag = 6;
   // IMPORTANT
   mhVarOptions.m_filteredChainDataOutputFileName = "outputData/filtChain_mh";
   mhVarOptions.m_filteredChainDataOutputFileType = "m";                
  mhVarOptions.m_filteredChainDataOutputAllowAll = 0;                  
  mhVarOptions.m_filteredChainDataOutputAllowedSet.insert(0);          
  mhVarOptions.m_displayCandidates = 0;
  mhVarOptions.m_putOutOfBoundsInChain = 0;
  mhVarOptions.m_tkUseLocalHessian = 0;
  mhVarOptions.m_tkUseNewtonComponent = 1;
  mhVarOptions.m_drMaxNumExtraStages = 0;          // IMPORTANT, typically 1
  mhVarOptions.m_drScalesForExtraStages.resize(1);
  //mhVarOptions.m_drScalesForExtraStages[0] = 10; // IMPORTANT
 // mhVarOptions.m_drScalesForExtraStages[1] = 7; // IMPORTANT
 // mhVarOptions.m_drScalesForExtraStages[2] = 10; // IMPORTANT
 // mhVarOptions.m_drScalesForExtraStages[3] = 20; // IMPORTANT
  mhVarOptions.m_drDuringAmNonAdaptiveInt = 1;     // IMPORTANT
  mhVarOptions.m_amKeepInitialMatrix = 0;          // IMPORTANT, typically 0
  mhVarOptions.m_amInitialNonAdaptInterval = 5000;  // IMPORTANT, typically 1
  mhVarOptions.m_amAdaptInterval = 5000;            // IMPORTANT
  mhVarOptions.m_amAdaptedMatricesDataOutputPeriod = 5000;
  mhVarOptions.m_amAdaptedMatricesDataOutputFileName = ".";
  mhVarOptions.m_amAdaptedMatricesDataOutputFileType = "m";
  mhVarOptions.m_amAdaptedMatricesDataOutputAllowAll = 0;
  //mhVarOptions.m_am_adaptedMatrices_dataOutputAllowedSet.insert(0);
  mhVarOptions.m_amEta = 2.4*2.4/((double) totalInitialVec.sizeLocal()); // IMPORTANT // prudenci 2013-08-25
  mhVarOptions.m_amEpsilon = 1e-05;                // IMPORTANT
  mhVarOptions.m_enableBrooksGelmanConvMonitor = 0;
  mhVarOptions.m_BrooksGelmanLag = 100;

  //gcm.calibrateWithBayesMetropolisHastings(mhVarOptions,totalInitialVec,&totalInitialProposalCovMatrix);
  gcm.calibrateWithBayesMLSampling();
 
  std::cout << "Got here 14" << std::endl; 
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "In compute()"
                          << ": finished calibrating 'gcm'"
                          << std::endl;
  }
 
  //***********************************************************************
  // Step 09 of 09: Make predictions with the calibrated computer model
  //***********************************************************************
  //QUESO::GslVector newExperimentScenarioVec(config_x_space.zeroVector());              // todo_rr
  //QUESO::GslMatrix newKmat_interp          (env,n_eta_space.map(),num_bases_eta);   // todo_rr
  //QUESO::GslMatrix newDmat                 (env,n_eta_space.map(),num_delta_bases); // todo_rr
  //QUESO::GslVector simulationOutputVec     (n_eta_space.zeroVector()); // Yes, size of simulation, since it is a prediction using the emulator
  //QUESO::GslVector discrepancyVec          (n_eta_space.zeroVector());
  //gcm->predictExperimentResults(newExperimentScenarioVec,newKmat_interp,newDmat,simulationOutputVec,discrepancyVec);
 
  //QUESO::GslVector newSimulationScenarioVec (config_x_space.zeroVector()); // todo_rr
  //QUESO::GslVector newSimulationParameterVec(paramSpace.zeroVector()); // todo_rr
  //QUESO::GslVector newSimulationOutputVec   (n_eta_space.zeroVector());
  //gcm->predictSimulationOutputs(newSimulationScenarioVec,newSimulationParameterVec,newSimulationOutputVec);

  //MPI_Finalize();

  return 0;
}
