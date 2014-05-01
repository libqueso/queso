#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/VectorSet.h>
#include <queso/GaussianProcessEmulator.h>

#include <cstdio>

// Read in data files
void readExpData(const std::vector<QUESO::GslVector *> & experimentScenarios,
    const std::vector<QUESO::GslVector *> & experimentOutputs) {
  FILE * fp_in = fopen("gp/scalar/ctf_dat.txt", "r");
  unsigned int i = 0;
  double pressure;
  
  while (fscanf(fp_in, "%lf\n", &pressure) != EOF) {
    (*(experimentOutputs[i]))[0] = pressure;
    i++;
  }

  fclose(fp_in);
}

void readSimData(const std::vector<QUESO::GslVector *> & simulationScenarios,
    const std::vector<QUESO::GslVector *> & simulationParameters,
    const std::vector<QUESO::GslVector *> & simulationOutputs) {
  FILE * fp_in = fopen("gp/scalar/dakota_pstudy.dat", "r");
  unsigned int i, id, size = 512;
  double k_tmasl, k_tmoml, k_tnrgl, k_xkwlx, k_cd, pressure;
  char line[size];

  // First line is a header, so we ignore it
  fgets(line, size, fp_in);

  i = 0;
  while (fscanf(fp_in, "%d %lf %lf %lf %lf %lf %lf\n", &id, &k_tmasl, &k_tmoml,
        &k_tnrgl, &k_xkwlx, &k_cd, &pressure) != EOF) {
    (*(simulationScenarios[i]))[0] = 0.5;

    (*(simulationParameters[i]))[0] = k_tmasl;
    (*(simulationParameters[i]))[1] = k_tmoml;
    (*(simulationParameters[i]))[2] = k_tnrgl;
    (*(simulationParameters[i]))[3] = k_xkwlx;
    (*(simulationParameters[i]))[4] = k_cd;

    (*(simulationOutputs[i]))[0] = pressure;
    i++;
  }

    fclose(fp_in);
}

int main(int argc, char ** argv) {
  // Step 0: Set up some variables
  unsigned int numExperiments = 10;  // Number of experiments
  unsigned int numUncertainVars = 5;  // Number of things to calibrate
  unsigned int numSimulations = 50;  // Number of simulations
  unsigned int numConfigVars = 1;  // Dimension of configuration space
  unsigned int numEta = 1;  // Number of responses the model is returning
  unsigned int experimentSize = 1;  // Size of each experiment

  MPI_Init(&argc, &argv);

  // Step 1: Set up QUESO environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  // Step 2: Set up prior for calibration parameters
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", numUncertainVars, NULL);

  // Parameter (theta) bounds:
  //   descriptors   'k_tmasl' 'k_xkle' 'k_xkwew' 'k_xkwlx' 'k_cd'
  //   upper_bounds   1.05      1.1      1.1       1.1       1.1
  //   lower_bounds   0.95      0.9      0.9       0.9       0.9
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  paramMins.cwSet(0.9);
  paramMaxs.cwSet(1.1);

  paramMins[0] = 0.95;
  paramMaxs[0] = 1.05;

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

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i] = new QUESO::GslVector(configSpace.zeroVector()); // 'x_{i+1}' in paper
    experimentVecs[i] = new QUESO::GslVector(experimentSpace.zeroVector());
  }

  readExpData(experimentScenarios, experimentVecs);
  readSimData(simulationScenarios, paramVecs, outputVecs);

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentMat(i, i) = 0.025 * 0.025;
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
  // for (unsigned int i = 0; i < paramInitials.sizeLocal(); i++) {
  //   paramInitials[i] = 0.95;
  // }

  // 'Best' initial conditions (from cobra report) are:
  paramInitials[0] = 1.05;  // k_tmasl
  paramInitials[1] = 0.9;  // k_tmoml
  paramInitials[2] = 0.9888;  // k_tnrgl
  paramInitials[3] = 0.995;  // k_xkwlx
  paramInitials[4] = 1.0178;  // k_cd
  paramInitials[5] = 0.0;  // Emulator mean

  for (unsigned int i = 6; i < 15; i++) {
    paramInitials[i] = 0.01;
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
