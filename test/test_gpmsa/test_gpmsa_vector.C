#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/VectorSet.h>
#include <queso/GPMSA.h>

#include <cstdio>
#include <cstdlib>

// Read in data files.  Return std deviations for each output.
std::vector<double>
readData(const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationScenarios,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationParameters,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationOutputs,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & experimentScenarios,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & experimentOutputs) {

  std::string simulationsFileName = "test_Regression/dakota_pstudy.dat";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    simulationsFileName = test_srcdir + ('/' + simulationsFileName);

  FILE * fp_in = fopen(simulationsFileName.c_str(), "r");
  if (!fp_in)
    queso_error_msg("Cannot find " << simulationsFileName);

  unsigned int i, id, size = 512;
  double k_tmasl, k_tmoml, k_tnrgl, k_xkwlx, k_cd, pressure;
  char line[size];

  // The user should know what the bounds on the data are
  double mins[] = {0.95, 0.9, 0.9, 0.9, 0.9};
  double maxs[] = {1.05, 1.1, 1.1, 1.1, 1.1};

  double meansim = 0;
  double m2sim = 0;

  double meanexpsim = 0;
  double m2expsim = 0;

  // First line is a header, so we ignore it
  char * gotline = fgets(line, size, fp_in);
  if (!gotline)
    queso_error_msg(simulationsFileName << " was unreadable");

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

    double exp_pressure = std::exp(pressure);
    (*(simulationOutputs[i]))[1] = exp_pressure;

    i++;

    double delta = pressure - meansim;
    meansim += delta / i;
    m2sim += delta * (pressure - meansim);

    double delta2 = exp_pressure - meanexpsim;
    meanexpsim += delta2 / i;
    m2expsim += delta2 * (exp_pressure - meanexpsim);
  }

  double stdsim = m2sim / (double)(i - 1);
  double stdexpsim = m2expsim / (double)(i - 1);

  // The user was required to standardise the experimental and
  // simulation data in the original GPMSA code.
  //
  // This is no longer required, but we'll do a rescaling *without*
  // recentering now to make sure that we can replicate our older
  // results.
  for (unsigned int j = 0; j < i; j++) {
    (*(simulationOutputs[j]))[0] -= meansim;
    (*(simulationOutputs[j]))[0] /= stdsim;
    (*(simulationOutputs[j]))[0] += meansim;

    (*(simulationOutputs[j]))[1] -= meanexpsim;
    (*(simulationOutputs[j]))[1] /= stdexpsim;
    (*(simulationOutputs[j]))[1] += meanexpsim;
  }

  fclose(fp_in);  // Done with simulation data

  // Read in experimental data
  std::string experimentsFileName = "test_Regression/ctf_dat.txt";
  if (test_srcdir)
    experimentsFileName = test_srcdir + ('/' + experimentsFileName);

  fp_in = fopen(experimentsFileName.c_str(), "r");
  if (!fp_in)
    queso_error_msg("Cannot find " << experimentsFileName);

  i = 0;
  while (fscanf(fp_in, "%lf\n", &pressure) != EOF) {
    (*(experimentOutputs[i]))[0] = pressure;
    (*(experimentOutputs[i]))[1] = std::exp(pressure);
    i++;
  }

  // We used to be required to standardise experimental data too
  for (unsigned int j = 0; j < i; j++) {
    (*(experimentOutputs[j]))[0] -= meansim;
    (*(experimentOutputs[j]))[0] /= stdsim;
    (*(experimentOutputs[j]))[0] += meansim;

    (*(experimentOutputs[j]))[1] -= meanexpsim;
    (*(experimentOutputs[j]))[1] /= stdexpsim;
    (*(experimentOutputs[j]))[1] += meanexpsim;
  }

  fclose(fp_in);

  // Returning standard deviation of data (simulation data is treated as
  // plentiful and the standard deviation of this data is treated as the 'true'
  // standard deviation of the data.  This can be argued.)
  std::vector<double> allstds;
  allstds.push_back(stdsim);
  allstds.push_back(stdexpsim);
  return allstds;
}

int main(int argc, char ** argv) {
  // Step 0: Set up some variables
  unsigned int numExperiments = 10;  // Number of experiments
  unsigned int numUncertainVars = 5;  // Number of things to calibrate
  unsigned int numSimulations = 50;  // Number of simulations
  unsigned int numConfigVars = 1;  // Dimension of configuration space

  // Number of outputs the experiment is returning.
  // For simplicity, this test will be isomorphic to the GPMSA scalar
  // test: Eta0 will be the same, and Eta1 will be exp(Eta0)
  unsigned int experimentSize = 2;      // Size of each experiment
  unsigned int numEta = experimentSize;

  if (argc < 2)
    queso_error_msg("Usage: " << argv[0] << " input_filename");

  std::string inputFileName = argv[1];

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);

  // Step 1: Set up QUESO environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);
#else
  QUESO::FullEnvironment env(inputFileName, "", NULL);
#endif

  // Step 2: Set up prior for calibration parameters
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", numUncertainVars, NULL);

  // Parameter (theta) bounds:
  //   descriptors   'k_tmasl' 'k_xkle' 'k_xkwew' 'k_xkwlx' 'k_cd'
  //   upper_bounds   1.05      1.1      1.1       1.1       1.1
  //   lower_bounds   0.95      0.9      0.9       0.9       0.9
  //
  // These bounds are dealt with when reading in the data
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
  // Regarding simulation scenario input values, the user should standardise
  // them so that they exist inside a hypercube.
  //
  // Regarding simulation output data, the user should transform it so that the
  // mean is zero and the variance is one.
  //
  // Regarding experimental scenario input values, the user should standardize
  // them so that they exist inside a hypercube.
  //
  // Regarding experimental data, the user should transformed it so that it has
  // zero mean and variance one.

  // GPMSA stores all the information about our simulation
  // data and experimental data.  It also stores default information about the
  // hyperparameter distributions.
  QUESO::GPMSAFactory<QUESO::GslVector, QUESO::GslMatrix>
    gpmsaFactory(env,
                 NULL,
                 priorRv,
                 configSpace,
                 paramSpace,
                 nEtaSpace,
                 experimentSpace,
                 numSimulations,
                 numExperiments);

  // We want to calibrate for the observation precision
  gpmsaFactory.options().m_calibrateObservationalPrecision = true;

  // Print our final options
  gpmsaFactory.options().print(std::cout);

  // std::vector containing all the points in scenario space where we have
  // simulations
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type>
    simulationScenarios(numSimulations);

  // std::vector containing all the points in parameter space where we have
  // simulations
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type>
    paramVecs(numSimulations);

  // std::vector containing all the simulation output data
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type>
    outputVecs(numSimulations);

  // std::vector containing all the points in scenario space where we have
  // experiments
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type>
    experimentScenarios(numExperiments);

  // std::vector containing all the experimental output data
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> experimentVecs(numExperiments);

  // The experimental output data observation error covariance matrix
  QUESO::SharedPtr<QUESO::GslMatrix>::Type experimentMat
    (new QUESO::GslMatrix(totalExperimentSpace.zeroVector()));

  // Instantiate each of the simulation points/outputs
  for (unsigned int i = 0; i < numSimulations; i++) {
    simulationScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector()));  // 'x_{i+1}^*' in paper
    paramVecs          [i].reset(new QUESO::GslVector(paramSpace.zeroVector()));  // 't_{i+1}^*' in paper
    outputVecs         [i].reset(new QUESO::GslVector(nEtaSpace.zeroVector()));  // 'eta_{i+1}' in paper
  }

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector())); // 'x_{i+1}' in paper
    experimentVecs[i].reset(new QUESO::GslVector(experimentSpace.zeroVector()));
  }

  // Read in data and store the standard deviations of the simulation
  // data.  We will need the standard deviation when we pass the
  // experiment error covariance matrix to QUESO.
  std::vector<double> stdsim =
    readData(simulationScenarios,
             paramVecs,
             outputVecs,
             experimentScenarios,
             experimentVecs);

  for (unsigned int i = 0; i < numExperiments; i++) {
      for (unsigned int j = 0; j < experimentSize; j++) {
          // Passing in error of experiments (standardised).
          // The "magic number" here will in practice be an estimate
          // of experimental error.  In this test example it is pulled
          // from our posterior.
          (*experimentMat)(experimentSize*i+j, experimentSize*i+j) = 1;
//          (*experimentMat)(experimentSize*i+j, experimentSize*i+j) =
//            (0.025 / stdsim[j]) * (0.025 / stdsim[j]);
        }
    }

  // Add simulation and experimental data
  gpmsaFactory.addSimulations(simulationScenarios, paramVecs, outputVecs);
  gpmsaFactory.addExperiments(experimentScenarios, experimentVecs, experimentMat);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix> postRv(
      "post_",
      gpmsaFactory.prior().imageSet().vectorSpace());
  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> ip("",
      NULL, gpmsaFactory, postRv);

  QUESO::GslVector paramInitials(
      gpmsaFactory.prior().imageSet().vectorSpace().zeroVector());

  // Initial condition of the chain

  // Start with the mean of the prior
  gpmsaFactory.prior().pdf().distributionMean(paramInitials);

  // These we'll override, not because we have to, but because the
  // regression gold standard predates distributionMean()
  paramInitials[5]  = 0.4; // emulator precision
  paramInitials[6]  = 0.4; // weights0 precision
  paramInitials[7]  = 0.4; // weights1 precision
  paramInitials[8]  = 0.97; // emulator corr str
  paramInitials[9]  = 0.97; // emulator corr str
  paramInitials[10] = 0.97; // emulator corr str
  paramInitials[11] = 0.97; // emulator corr str
  paramInitials[12] = 0.20; // emulator corr str
  paramInitials[13] = 0.80; // emulator corr str
  paramInitials[14] = 10.0; // discrepancy precision
  paramInitials[15] = 0.97; // discrepancy corr str
  paramInitials[16] = 8000.0; // emulator data precision
  // paramInitials[17]  = 1.0;  // observation error precision

  QUESO::GslMatrix proposalCovMatrix(
      gpmsaFactory.prior().imageSet().vectorSpace().zeroVector());

  // Start with the covariance matrix for the prior.  We won't
  // actually use these values, since we want to match a gold standard
  // that was generated from hand-twiddled values, but this will at
  // least get a touch more test coverage on distributionVariance()
  gpmsaFactory.prior().pdf().distributionVariance(proposalCovMatrix);

  // Setting the proposal covariance matrix by hand.  This requires great
  // forethough, and can generally be referred to as a massive hack.  These
  // values were taken from the gpmsa matlab code and fiddled with.
  double scale = 3000.0;
  proposalCovMatrix(0, 0)   = 3.1646 / 10.0;  // param 1
  proposalCovMatrix(1, 1)   = 3.1341 / 10.0;  // param 2
  proposalCovMatrix(2, 2)   = 3.1508 / 10.0;  // param 3
  proposalCovMatrix(3, 3)   = 0.3757 / 10.0;  // param 4
  proposalCovMatrix(4, 4)   = 0.6719 / 10.0;  // param 5
  proposalCovMatrix(5, 5)   = 0.4953 / scale;  // emulator precision
  proposalCovMatrix(6, 6)   = 0.4953 / scale;  // weights0 precision
  proposalCovMatrix(7, 7)   = 0.4953 / scale;  // weights1 precision
  proposalCovMatrix(8, 8)   = 0.6058 / scale;  // emulator corr str
  proposalCovMatrix(9, 9)   = 7.6032e-04 / scale;  // emulator corr str
  proposalCovMatrix(10, 10) = 8.3815e-04 / scale;  // emulator corr str
  proposalCovMatrix(11, 11) = 7.5412e-04 / scale;  // emulator corr str
  proposalCovMatrix(12, 12) = 0.2682 / scale;  // emulator corr str
  proposalCovMatrix(13, 13) = 0.0572 / scale;  // emulator corr str
  proposalCovMatrix(14, 14) = 1.3417 / scale;  // discrepancy precision
  proposalCovMatrix(15, 15) = 0.3461 / scale;  // discrepancy corr str
  proposalCovMatrix(16, 16) = 495.3 / scale;  // emulator data precision
  proposalCovMatrix(17, 17) = 0.4953 / scale;  // observation error precision

  // Square to get variances
  for (unsigned int i = 0; i < 15; i++) {
    proposalCovMatrix(i, i) = proposalCovMatrix(i, i) * proposalCovMatrix(i, i);
  }

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
