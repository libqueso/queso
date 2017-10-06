#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/VectorSet.h>
#include <queso/GPMSA.h>
#include <queso/TensorProductMesh.h>

#include <cstdio>
#include <cstdlib>

#define LINE_SIZE 512

const double t_length = 2;   // We'll have a grid in time t
const int n_grid_points = 5; // With 5 points
const int max_experiment_points = 3; // But each experiment will have 3 evaluations or fewer

// Read in data files.  Return std deviations for each output.
void
readData(const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationScenarios,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationParameters,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationOutputs,
         const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & /* experimentScenarios */,
         std::vector<std::vector<QUESO::SimulationOutputPoint> > & experimentPoints,
         std::vector<QUESO::SharedPtr<QUESO::Map>::Type> & experimentMaps,
         std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & experimentOutputs) {

  std::string simulationsFileName = "test_Regression/dakota_pstudy.dat";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    simulationsFileName = test_srcdir + ('/' + simulationsFileName);

  FILE * fp_in = fopen(simulationsFileName.c_str(), "r");
  if (!fp_in)
    queso_error_msg("Cannot find " << simulationsFileName);

  unsigned int id;
  double pressure;
  char line[LINE_SIZE];

  // First line is a header, so we ignore it
  char * gotline = fgets(line, LINE_SIZE, fp_in);
  if (!gotline)
    queso_error_msg(simulationsFileName << " was unreadable");

  // Read in simulation data.  Fake functional data by just adding a
  // spatially-dependent term to the real data.
  for (unsigned int i=0; i != simulationParameters.size(); ++i) {
    fscanf(fp_in, "%u %lf %lf %lf %lf %lf %lf\n", &id,
           &((*(simulationParameters[i]))[0]),
           &((*(simulationParameters[i]))[1]),
           &((*(simulationParameters[i]))[2]),
           &((*(simulationParameters[i]))[3]),
           &((*(simulationParameters[i]))[4]),
           &pressure);
    (*(simulationScenarios[i]))[0] = 0.5;

    for (unsigned int p = 0; p != n_grid_points; ++p) {
      (*(simulationOutputs[i]))[p] = pressure + ((double)p / n_grid_points);
    }
  }

  queso_assert_equal_to(fgetc(fp_in), EOF); // Done with simulation data?
  fclose(fp_in);                            // Done with simulation data.

  // Read in experimental data.  Fake one or two point locations for
  // each experiment and add the (consistently fake) output offset for
  // each point.
  std::string experimentsFileName = "test_Regression/ctf_dat.txt";
  if (test_srcdir)
    experimentsFileName = test_srcdir + ('/' + experimentsFileName);

  fp_in = fopen(experimentsFileName.c_str(), "r");
  if (!fp_in)
    queso_error_msg("Cannot find " << experimentsFileName);

  unsigned int i = 0;
  double t_loc = 0;
  while (fscanf(fp_in, "%lf\n", &pressure) != EOF) {
    unsigned int n_points = (i % (max_experiment_points - 1)) + 1;

    experimentMaps.push_back
      (std::shared_ptr<QUESO::Map>
        (new QUESO::Map 
          (n_points, 0, simulationOutputs[0]->map().Comm())));
    experimentOutputs.push_back
      (std::shared_ptr<QUESO::GslVector>
        (new QUESO::GslVector
          (simulationOutputs[0]->env(), *experimentMaps.back())));
    experimentPoints.push_back(std::vector<QUESO::SimulationOutputPoint>(n_points));

    for (unsigned int p = 0; p != n_points; ++p) {
      t_loc += 1.6;
      t_loc = std::fmod(t_loc, t_length);
      (*(experimentOutputs[i]))[p] = pressure + (t_loc/t_length);
      experimentPoints[i][p].t() = t_loc;
    }
    
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
  unsigned int numEta = n_grid_points; // Number of simulation outputs

  // Number of outputs the experiment is returning.

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

  // Step 4: Instantiate the Gaussian process emulator object
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
                 numSimulations,
                 numExperiments);

  // We want to calibrate for the observation precision
  gpmsaFactory.options().m_calibrateObservationalPrecision = true;

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

  // std::vector containers of all the experimental output structure
  // and data.  We'll just let readData handle filling empty
  // containers for these.
  //
  // We no longer have a single 'output' space for the experiments;
  // each experiment may have a different number of outputs, which
  // with the QUESO::Vector interface means each experiment output
  // vector needs a different Map.

  std::vector<QUESO::SharedPtr<QUESO::Map>::Type>       experimentMaps;
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> experimentVecs;
  std::vector<std::vector<QUESO::SimulationOutputPoint> > experimentPoints;

  // Instantiate each of the simulation points/outputs
  for (unsigned int i = 0; i < numSimulations; i++) {
    simulationScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector()));  // 'x_{i+1}^*' in paper
    paramVecs          [i].reset(new QUESO::GslVector(paramSpace.zeroVector()));  // 't_{i+1}^*' in paper
    outputVecs         [i].reset(new QUESO::GslVector(nEtaSpace.zeroVector()));  // 'eta_{i+1}' in paper
  }

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector())); // 'x_{i+1}' in paper
  }

  // Read in data and store the standard deviations of the simulation
  // data.  We will need the standard deviation when we pass the
  // experiment error covariance matrix to QUESO.
  readData(simulationScenarios,
           paramVecs,
           outputVecs,
           experimentScenarios,
           experimentPoints,
           experimentMaps,
           experimentVecs);

  // Let QUESO know that every experiment output is for variable 0.
  // We should probably make this an automatic default...
  std::vector<std::vector<unsigned int> > experimentVars(numExperiments);
  for (unsigned int i=0; i != numExperiments; ++i)
    experimentVars[i].resize(experimentPoints[i].size(), 0);

  // std::vector containing error covariance matrices
  std::vector<QUESO::SharedPtr<QUESO::GslMatrix>::Type> experimentErrorMats;

  for (unsigned int i = 0; i < numExperiments; i++) {
      experimentErrorMats.push_back
        (QUESO::SharedPtr<QUESO::GslMatrix>::Type
          (new QUESO::GslMatrix(*experimentVecs[i], 0)));
      for (unsigned int j = 0; j < experimentVecs[i]->sizeGlobal(); j++) {
          // Passing in error of experiments (standardised).
          // The "magic number" here will in practice be an estimate
          // of experimental error.  In this test example it is pulled
          // from our posterior.
          (*experimentErrorMats[i])(j, j) = 1;
        }
    }

  // Add simulation mesh.  We'll place our cooked-up outputs on an
  // equispaced grid from t=0 to t=t_length.
  std::vector<double> gridPoints(n_grid_points,0);
  const double n_elem = n_grid_points-1; // double to force upconversion later
  for (unsigned int i=0; i != n_grid_points; ++i)
    gridPoints[i] = i/n_elem*t_length;

  // Our current APIs force us to use shared pointers for this too,
  // which is a bit awkward...
  typedef QUESO::TensorProductMesh<QUESO::GslVector> Mesh;
  typedef QUESO::SharedPtr<Mesh>::Type MeshPtr;
  MeshPtr simulationMesh (MeshPtr(new Mesh));
  simulationMesh->set_t_coordinates(gridPoints);

  // Add simulation and experimental data
  gpmsaFactory.addSimulationMesh(simulationMesh);
  gpmsaFactory.addSimulations(simulationScenarios, paramVecs, outputVecs);

  gpmsaFactory.addExperiments
    (experimentScenarios, experimentVecs, experimentErrorMats,
     &experimentPoints, &experimentVars);

  // Print our final options.  Note that this had to come *after* the
  // data was added, because only then is the GPMSAFactory able to
  // finally determine gaussian discrepancy autogeneration parameters.
  gpmsaFactory.options().print(std::cout);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix> postRv(
      "post_",
      gpmsaFactory.prior().imageSet().vectorSpace());
  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> ip("",
      NULL, gpmsaFactory, postRv);

  QUESO::GslVector paramInitials(
      gpmsaFactory.prior().imageSet().vectorSpace().zeroVector());

  // Start with the mean of the prior as our initial sample
  gpmsaFactory.prior().pdf().distributionMean(paramInitials);

  // Take the covariance matrix for the prior, scaled down to try to
  // keep samples within sane bounds, as our proposal covariance.
  QUESO::GslMatrix proposalCovMatrix(
      gpmsaFactory.prior().imageSet().vectorSpace().zeroVector());
  gpmsaFactory.prior().pdf().distributionVariance(proposalCovMatrix); 
  proposalCovMatrix *= 0.001;

  // GPMSA *really* dislikes too-low discrepancy precision values...
  proposalCovMatrix(13, 13) *= 0.001;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
