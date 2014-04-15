#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/VectorSet.h>
#include <queso/GaussianProcessEmulator.h>

int main(int argc, char ** argv) {
  // Step 0: Set up some variables
  unsigned int numExperiments = 6;  // Number of experiments
  unsigned int numUncertainVars = 1;  // Number of things to calibrate
  unsigned int numSimulations = 25;  // Number of simulations
  unsigned int numConfigVars = 1;  // Dimension of configuration space
  unsigned int numEta = 1;  // Number of responses the model is returning
  unsigned int experimentSize = 1;  // Size of each experiment

  MPI_Init(&argc, &argv);

  // Step 1: Set up QUESO environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  // Step 2: Set up prior for calibration parameters
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", numUncertainVars, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(0.0);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(1.0);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorRv("prior_",
      paramDomain);

  // Step 3: Instantiate the 'scenario' and 'output' spaces for simulation
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> configSpace(env,
      "scenario_", numConfigVars, NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> nEtaSpace(env,
      "output_", numEta, NULL);

  // Step 4: Instantiate the 'output' space for the experiments
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> experimentSpace(env,
      "experimentspace_", experimentSize, NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> totalExperimentSpace(env,
      "experimentspace_", experimentSize * numExperiments, NULL);

  // Step 5: Instantiate the Gaussian process emulator object
  //
  // Regarding simulation scenario input values, QUESO should standardise them
  // so that they exist inside a hypercube.
  //
  // Regarding simulation output data, QUESO should transform it so that the
  // mean is zero and the variance is one.
  //
  // Regarding experimental scenario input values, QUESO should standardize
  // them so that they exist inside a hypercube.
  //
  // Regarding experimental data, QUESO should transformed it so that it has
  // zero mean and variance one.

  // GaussianProcessEmulator stores all the information about our simulation
  // data and experimental data.  It also stores default information about the
  // hyperparameter distributions.
  QUESO::GaussianProcessFactory<QUESO::GslVector, QUESO::GslMatrix>
    gpFactory("gp_",
              priorRv,
              configSpace,
              paramSpace,
              nEtaSpace,
              experimentSpace,
              numSimulations,
              numExperiments);

  // std::vector containing all the points in scenario space where we have
  // simulations
  std::vector<QUESO::GslVector *> simulationScenarios(numSimulations,
      (QUESO::GslVector *) NULL);

  // std::vector containing all the points in parameter space where we have
  // simulations
  std::vector<QUESO::GslVector *> paramVecs(numSimulations,
      (QUESO::GslVector *) NULL);

  // std::vector containing all the simulation output data
  std::vector<QUESO::GslVector *> outputVecs(numSimulations,
      (QUESO::GslVector *) NULL);

  // std::vector containing all the points in scenario space where we have
  // experiments
  std::vector<QUESO::GslVector *> experimentScenarios(numExperiments,
      (QUESO::GslVector *) NULL);

  // std::vector containing all the experimental output data
  std::vector<QUESO::GslVector *> experimentVecs(numExperiments,
      (QUESO::GslVector *) NULL);

  // The experimental output data observation error covariance matrix
  QUESO::GslMatrix experimentMat(totalExperimentSpace.zeroVector());

  // Instantiate each of the simulation points/outputs
  for (unsigned int i = 0; i < numSimulations; i++) {
    simulationScenarios[i] = new QUESO::GslVector(configSpace.zeroVector());  // 'x_{i+1}^*' in paper
    paramVecs          [i] = new QUESO::GslVector(paramSpace.zeroVector());  // 't_{i+1}^*' in paper
    outputVecs         [i] = new QUESO::GslVector(nEtaSpace.zeroVector());  // 'eta_{i+1}' in paper
  }

  // All the positions in scenario space where simulations were run
  // This should probably be read from a file
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
  // This should probably be read from a file
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
  // This should probably be read from a file
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

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i] = new QUESO::GslVector(configSpace.zeroVector()); // 'x_{i+1}' in paper
    experimentVecs[i] = new QUESO::GslVector(experimentSpace.zeroVector());
  }

  // All the positions in scenario space where we have experimental data
  // This should probably be read from a file
  (*(experimentScenarios[0]))[0] = 0.0;
  (*(experimentScenarios[1]))[0] = 0.2;
  (*(experimentScenarios[2]))[0] = 0.4;
  (*(experimentScenarios[3]))[0] = 0.6;
  (*(experimentScenarios[4]))[0] = 0.8;
  (*(experimentScenarios[5]))[0] = 1.0;

  // These are the experimental data
  // This should probably be read from a file
  (*(experimentVecs[0]))[0] = 1.59294826065229;
  (*(experimentVecs[1]))[0] = 2.17696977275016;
  (*(experimentVecs[2]))[0] = 2.87061286591332;
  (*(experimentVecs[3]))[0] = 3.8330395105599;
  (*(experimentVecs[4]))[0] = 4.59654198432239;
  (*(experimentVecs[5]))[0] = 4.75087857533489;

  for (unsigned int i = 0; i < 5; i++) {
    experimentMat(i, i) = 0.075 * 0.075;
  }

  // Add simulation and experimental data
  gpFactory.addSimulations(simulationScenarios, paramVecs, outputVecs);
  gpFactory.addExperiments(experimentScenarios, experimentVecs, &experimentMat);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix> postRv(
      "post_",
      gpFactory.prior().imageSet().vectorSpace());
  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> ip("",
      NULL, gpFactory, postRv);

  QUESO::GslVector paramInitials(
      gpFactory.prior().imageSet().vectorSpace().zeroVector());

  // Initial condition of the chain
  for (unsigned int i = 0; i < paramInitials.sizeLocal(); i++) {
    paramInitials[i] = 0.5;
  }

  QUESO::GslMatrix proposalCovMatrix(
      gpFactory.prior().imageSet().vectorSpace().zeroVector());
  for (unsigned int i = 0; i < proposalCovMatrix.numRowsLocal(); i++) {
    proposalCovMatrix(i, i) = 0.01;
  }
  std::cout << "got here 25" << std::endl;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
  std::cout << "got here 26" << std::endl;

  MPI_Finalize();
  std::cout << "got here 27" << std::endl;

  return 0;
}
