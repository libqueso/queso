#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/VectorSet.h>
#include <queso/GaussianProcessEmulator.h>

#include <cstdio>

// Read in data files
void readData(const std::vector<QUESO::GslVector *> & simulationScenarios,
    const std::vector<QUESO::GslVector *> & simulationParameters,
    const std::vector<QUESO::GslVector *> & simulationOutputs,
    const std::vector<QUESO::GslVector *> & experimentScenarios,
    const std::vector<QUESO::GslVector *> & experimentOutputs) {
  FILE * fp_in = fopen("gp/scalar/dakota_pstudy.dat", "r");
  unsigned int i, id, size = 512;
  double k_tmasl, k_tmoml, k_tnrgl, k_xkwlx, k_cd, pressure;
  char line[size];

  double mins[] = {0.95, 0.9, 0.9, 0.9, 0.9};
  double maxs[] = {1.05, 1.1, 1.1, 1.1, 1.1};

  double meansim = 0;
  double m2sim = 0;

  // First line is a header, so we ignore it
  fgets(line, size, fp_in);

  i = 0;
  while (fscanf(fp_in, "%d %lf %lf %lf %lf %lf %lf\n", &id, &k_tmasl, &k_tmoml,
        &k_tnrgl, &k_xkwlx, &k_cd, &pressure) != EOF) {
    (*(simulationScenarios[i]))[0] = 0.5;

    (*(simulationParameters[i]))[0] = (k_tmasl - mins[0]) / (maxs[0] - mins[0]);
    (*(simulationParameters[i]))[1] = (k_tmoml - mins[1]) / (maxs[1] - mins[1]);
    (*(simulationParameters[i]))[2] = (k_tnrgl - mins[2]) / (maxs[2] - mins[2]);
    (*(simulationParameters[i]))[3] = (k_xkwlx - mins[3]) / (maxs[3] - mins[3]);
    (*(simulationParameters[i]))[4] = (k_cd    - mins[4]) / (maxs[4] - mins[4]);

    (*(simulationOutputs[i]))[0] = pressure;
    i++;

    double delta = pressure - meansim;
    meansim += delta / i;
    m2sim += delta * (pressure - meansim);
  }

  double stdsim = m2sim / (double)(i - 1);

  for (unsigned int j = 0; j < i; j++) {
    (*(simulationOutputs[j]))[0] -= meansim;
    (*(simulationOutputs[j]))[0] /= stdsim;
  }

  fclose(fp_in);  // Done with simulation data

  // Read in experimental data
  fp_in = fopen("gp/scalar/ctf_dat.txt", "r");
  i = 0;
  while (fscanf(fp_in, "%lf\n", &pressure) != EOF) {
    (*(experimentOutputs[i]))[0] = pressure;
    i++;
  }

  for (unsigned int j = 0; j < i; j++) {
    (*(experimentOutputs[j]))[0] -= meansim;
    (*(experimentOutputs[j]))[0] /= stdsim;
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

  paramMins.cwSet(0.0);
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

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i] = new QUESO::GslVector(configSpace.zeroVector()); // 'x_{i+1}' in paper
    experimentVecs[i] = new QUESO::GslVector(experimentSpace.zeroVector());
  }

  readData(simulationScenarios,
      paramVecs,
      outputVecs,
      experimentScenarios,
      experimentVecs);

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentMat(i, i) = 1.0;
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

  paramInitials[5]  = 0.00;  // Emulator mean

  for (unsigned int i = 6; i < 15; i++) {
    paramInitials[i] = 0.01;
  }

  QUESO::GslMatrix proposalCovMatrix(
      gpFactory.prior().imageSet().vectorSpace().zeroVector());

  double scale = 1000.0;
  proposalCovMatrix(0, 0)   = 0.2;
  proposalCovMatrix(1, 1)   = 0.2;
  proposalCovMatrix(2, 2)   = 0.2;
  proposalCovMatrix(3, 3)   = 0.2;
  proposalCovMatrix(4, 4)   = 0.2;
  proposalCovMatrix(5, 5)   = 0.1 / scale;  // not used
  proposalCovMatrix(6, 6)   = 5.0 / scale;
  proposalCovMatrix(7, 7)   = 0.1 / scale;
  proposalCovMatrix(8, 8)   = 0.1 / scale;
  proposalCovMatrix(9, 9)   = 0.1 / scale;
  proposalCovMatrix(10, 10) = 0.1 / scale;
  proposalCovMatrix(11, 11) = 0.1 / scale;
  proposalCovMatrix(12, 12) = 0.1 / scale;
  proposalCovMatrix(13, 13) = 10.0 / scale;
  proposalCovMatrix(14, 14) = 0.1 / scale;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
  std::cout << "got here 26" << std::endl;

  MPI_Finalize();
  std::cout << "got here 27" << std::endl;

  return 0;
}
