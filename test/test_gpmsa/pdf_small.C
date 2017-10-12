//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// Bayes linear verification problems due to Brian Williams
//
// Described in: [1] Adams, B.M., Coleman, Kayla, Hooper, Russell W.,
//   Khuwaileh, Bassam A., Lewis, Allison, Smith, Ralph C., Swiler,
//   Laura P., Turinsky, Paul J., Williams, Brian J., "User Guidelines
//   and Best Practices for CASL VUQ Analysis Using Dakota, Sandia
//   Technical Report 2016-11614 (also CASL-U-2016-1233-000),
//   September 2016.
//
// Original QUESO/GPMSA example created by Brian Adams, 20170505

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/VectorSet.h>
#include <queso/GPMSA.h>

#include <cstdio>

#define TOL 1e-5

void open_data_file(const std::string& data_filename, std::ifstream& data_stream) {
  std::string example_dir = "gp/linear_verify";
  data_stream.open(example_dir + "/" + data_filename);
  if (!data_stream.good()) {
    data_stream.open(data_filename);
    if (!data_stream.good())
      queso_error_msg(std::string("Cannot open data file: ") + data_filename);
  }
  data_stream.exceptions(std::ifstream::badbit | std::ifstream::eofbit |
                         std::ifstream::failbit);
}

// Read in data files
void readData
(const std::string& sim_data_filename, const std::string& exp_data_filename,
 const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationScenarios,
 const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationParameters,
 const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & simulationOutputs,
 const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & experimentScenarios,
 const std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type> & experimentOutputs,
 QUESO::GslMatrix& experimentMat)
{
  unsigned int num_config = simulationScenarios[0]->sizeGlobal();
  unsigned int num_params = simulationParameters[0]->sizeGlobal();
  unsigned int num_responses = simulationOutputs[0]->sizeGlobal();
  unsigned int num_simulations = simulationOutputs.size();

  // File containing x, theta, y (vertical concatenation of lhs.txt, y_mod.txt)
  std::ifstream sim_data;
  open_data_file(sim_data_filename, sim_data);
  try {
    for (unsigned int i = 0; i < num_simulations; ++i) {
      for (unsigned int j = 0; j < num_config; ++j)
        sim_data >> (*(simulationScenarios[i]))[j];
      for (unsigned int j = 0; j < num_params; ++j)
        sim_data >> (*(simulationParameters[i]))[j];
      for (unsigned int j = 0; j < num_responses; ++j) {
        sim_data >> (*(simulationOutputs[i]))[j];
      }
    }
  }
  catch (const std::ifstream::failure& e) {
    queso_error_msg(std::string("Error reading ") + sim_data_filename + ": " +
                    e.what());
  }

  // Read in experimental data
  unsigned int num_experiments = experimentOutputs.size();
  std::ifstream exp_data;
  open_data_file(exp_data_filename, exp_data);
  try {
    for (unsigned int i = 0; i < num_experiments; ++i) {
      for (unsigned int j = 0; j < num_config; ++j)
        exp_data >> (*(experimentScenarios[i]))[j];
      for (unsigned int j = 0; j < num_responses; ++j)
        exp_data >> (*(experimentOutputs[i]))[j];
    }
  }
  catch (const std::ifstream::failure& e) {
    queso_error_msg(std::string("Error reading ") + exp_data_filename + ": " +
                    e.what());
  }
}

void run_scalar(const QUESO::FullEnvironment& env)
{
  // Step 0: Set up some variables
  unsigned int numExperiments = 1;  // Number of experiments
  unsigned int numUncertainVars = 1;  // Number of things to calibrate (3 beta)
  unsigned int numSimulations = 2;  // Number of simulations
  unsigned int numConfigVars = 1;  // Dimension of configuration space (3 x)
  unsigned int numEta = 1;  // Number of responses the model is returning
  unsigned int experimentSize = 1;  // Size of each experiment

  // Step 2: Set up prior for calibration parameters
  QUESO::VectorSpace<> paramSpace(env, "param_", numUncertainVars, NULL);

  // Parameter (theta) bounds:
  //   descriptors   'beta0'
  //   upper_bounds    0.45
  //   lower_bounds   -0.1
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = -0.1;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] =  0.45;

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);
  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  // Step 3: Instantiate the 'scenario' and 'output' spaces for simulation

  // Config space:
  //         Min  Max
  // x1      0.9   1.1
  QUESO::VectorSpace<> configSpace(env, "scenario_", numConfigVars, NULL);
  QUESO::VectorSpace<> nEtaSpace(env, "output_", numEta, NULL);

  // Step 4: Instantiate the 'output' space for the experiments
  QUESO::VectorSpace<> experimentSpace(env, "experimentspace_", experimentSize,
      NULL);

  QUESO::VectorSpace<> totalExperimentSpace(env, "experimentspace_",
      experimentSize * numExperiments, NULL);

  // Step 5: Instantiate the Gaussian process emulator object

  // GPMSA stores all the information about our simulation data and experimental
  // data, and if the users opts to, will normalise this to have zero mean and
  // unit variance.  It also stores default information about the hyperparameter
  // distributions.
  QUESO::GPMSAFactory<> gpmsaFactory(env,
                                     NULL,
                                     priorRv,
                                     configSpace,
                                     paramSpace,
                                     nEtaSpace,
                                     numSimulations,
                                     numExperiments);

  QUESO::GPMSAOptions& gp_opts = gpmsaFactory.options();

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
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type>
    experimentVecs(numExperiments);

  // The experimental output data observation error covariance matrix
  QUESO::SharedPtr<QUESO::GslMatrix>::Type experimentMat
    (new QUESO::GslMatrix(totalExperimentSpace.zeroVector()));

  // Instantiate each of the simulation points/outputs
  for (unsigned int i = 0; i < numSimulations; i++) {
    simulationScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector()));  // 'x_{i+1}^*' in paper
    paramVecs          [i].reset(new QUESO::GslVector(paramSpace.zeroVector()));  // 't_{i+1}^*' in paper
    outputVecs         [i].reset(new QUESO::GslVector(nEtaSpace.zeroVector()));  // 'eta_{i+1}' in paper
  }

  // Must set scaling before adding experiments due to construct order
  for (unsigned int i = 0; i < numUncertainVars; i++)
    gp_opts.set_autoscale_minmax_uncertain_parameter(i);
  for (unsigned int i = 0; i < numConfigVars; i++)
    gp_opts.set_autoscale_minmax_scenario_parameter(i);
  for (unsigned int i = 0; i < numEta; i++)
    gp_opts.set_autoscale_meanvar_output(i);

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector())); // 'x_{i+1}' in paper
    experimentVecs     [i].reset(new QUESO::GslVector(experimentSpace.zeroVector()));
  }

  // Experiment observation error (also needs to be scaled by
  // simulation variance)
  for (unsigned int i = 0; i < numExperiments; i++)
    (*experimentMat)(i, i) = 1.0;

  // Read in data and store the standard deviation of the simulation
  // data (ignored for now).
  const char * test_srcdir = std::getenv("srcdir");

  std::string simInputFileName = "test_gpmsa/sim_scalar_small.dat";
  std::string expInputFileName = "test_gpmsa/y_exp_scalar_small.txt";
  std::string solInputFileName = "test_gpmsa/regression_solution_pdf_small.txt";

  if (test_srcdir) {
    simInputFileName = test_srcdir + ('/' + simInputFileName);
    expInputFileName = test_srcdir + ('/' + expInputFileName);
    solInputFileName = test_srcdir + ('/' + solInputFileName);
  }

  readData(simInputFileName,
           expInputFileName,
           simulationScenarios,
           paramVecs,
           outputVecs,
           experimentScenarios,
           experimentVecs,
           *experimentMat);

  // Add simulation and experimental data
  gpmsaFactory.addSimulations(simulationScenarios, paramVecs, outputVecs);
  gpmsaFactory.addExperiments(experimentScenarios, experimentVecs, experimentMat);

  QUESO::GenericVectorRV<> postRv( "post_",
      gpmsaFactory.prior().imageSet().vectorSpace());

  QUESO::GslVector point(
      gpmsaFactory.prior().imageSet().vectorSpace().zeroVector());

  unsigned int num_lines = 10000;
  std::vector<double> expected_log_likelihoods(num_lines);
  std::vector<double> expected_log_priors(num_lines);
  std::vector<double> expected_log_posteriors(num_lines);
  std::vector<double> computed_log_likelihoods(num_lines);
  std::vector<double> computed_log_priors(num_lines);
  std::vector<double> computed_log_posteriors(num_lines);

  // File containing x, theta, y (vertical concatenation of lhs.txt, y_mod.txt)
  std::ifstream solution_data;
  open_data_file(solInputFileName, solution_data);
  for (unsigned int i = 0; i < num_lines; i++) {
    for (unsigned int j = 0; j < point.sizeLocal(); ++j)
      solution_data >> point[j];

    double log_pdf;
    solution_data >> log_pdf;
    expected_log_likelihoods[i] = log_pdf;

    solution_data >> log_pdf;
    expected_log_priors[i] = log_pdf;

    solution_data >> log_pdf;
    expected_log_posteriors[i] = log_pdf;

    log_pdf = (gpmsaFactory.prior().pdf().lnValue(point));
    computed_log_priors[i] = log_pdf;

    log_pdf = (gpmsaFactory.getGPMSAEmulator().lnValue(point, NULL, NULL, NULL, NULL));
    computed_log_likelihoods[i] = log_pdf;
  }
  solution_data.close();

  // We subtract off the first difference in log priors because the
  // normalisation constants between our implementation and Matlab's
  // implementation are not the same.  We expect them to only be the same up to
  // an additive constant.
  double initial_diff_prior = expected_log_priors[0] - computed_log_priors[0];
  for (unsigned int i = 0; i < expected_log_priors.size(); i++) {
    double diff_prior = expected_log_priors[i] - computed_log_priors[i] - initial_diff_prior;
    queso_require_less_equal_msg(std::abs(diff_prior), TOL, "computed log prior differs too much from expected");
  }

  // We don't subtract off the initial because the normalisation constants
  // between our implementation and Matlab's implementation appear to be the
  // same
  for (unsigned int i = 0; i < expected_log_likelihoods.size(); i++) {
    double diff_likelihood = expected_log_likelihoods[i] - computed_log_likelihoods[i];
    queso_require_less_equal_msg(std::abs(diff_likelihood), TOL, "computed log likelihood differs too much from expected");
  }
}

void run_multivariate(const QUESO::FullEnvironment& env)
{
  // Step 0: Set up some variables
  unsigned int numExperiments = 1;  // Number of experiments
  unsigned int numUncertainVars = 1;  // Number of things to calibrate (1 beta)
  unsigned int numSimulations = 2;  // Number of simulations
  unsigned int numConfigVars = 1;  // Dimension of configuration space (1 x)
  unsigned int numEta = 2;  // Number of responses the model is returning
  unsigned int experimentSize = 2;  // Size of each experiment

  // Step 2: Set up prior for calibration parameters
  QUESO::VectorSpace<> paramSpace(env, "param_", numUncertainVars, NULL);

  // Parameter (theta) bounds:
  //   descriptors   'beta'
  //   upper_bounds  -0.1
  //   lower_bounds  -0.5
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = -0.5;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = -0.1;

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  // Step 3: Instantiate the 'scenario' and 'output' spaces for simulation

  // Config space:
  //         Min  Max
  // x      -1.5   1.0

  QUESO::VectorSpace<> configSpace(env, "scenario_", numConfigVars, NULL);

  QUESO::VectorSpace<> nEtaSpace(env, "output_", numEta, NULL);

  // Step 4: Instantiate the 'output' space for the experiments
  QUESO::VectorSpace<> experimentSpace(env, "experimentspace_", experimentSize,
      NULL);

  QUESO::VectorSpace<> totalExperimentSpace(env, "experimentspace_",
      experimentSize * numExperiments, NULL);

  // Step 5: Instantiate the Gaussian process emulator object

  // GPMSA stores all the information about our simulation data and experimental
  // data, and if the users opts to, will normalise this to have zero mean and
  // unit variance.  It also stores default information about the hyperparameter
  // distributions.
  QUESO::GPMSAFactory<> gpmsaFactory(env,
                                     NULL,
                                     priorRv,
                                     configSpace,
                                     paramSpace,
                                     nEtaSpace,
                                     numSimulations,
                                     numExperiments);

  QUESO::GPMSAOptions& gp_opts = gpmsaFactory.options();

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
  std::vector<QUESO::SharedPtr<QUESO::GslVector>::Type>
    experimentVecs(numExperiments);

  // The experimental output data observation error covariance matrix
  QUESO::SharedPtr<QUESO::GslMatrix>::Type experimentMat
    (new QUESO::GslMatrix(totalExperimentSpace.zeroVector()));

  // Instantiate each of the simulation points/outputs
  for (unsigned int i = 0; i < numSimulations; i++) {
    simulationScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector()));  // 'x_{i+1}^*' in paper
    paramVecs          [i].reset(new QUESO::GslVector(paramSpace.zeroVector()));  // 't_{i+1}^*' in paper
    outputVecs         [i].reset(new QUESO::GslVector(nEtaSpace.zeroVector()));  // 'eta_{i+1}' in paper

    // Must set scaling before adding experiments due to construct order
    // As of this implementation autoscale only affects params/scenarios
    gp_opts.set_autoscale_minmax_uncertain_parameter(i);
    gp_opts.set_autoscale_minmax_scenario_parameter(i);
    gp_opts.set_autoscale_meanvar_output(i);
  }

  for (unsigned int i = 0; i < numExperiments; i++) {
    experimentScenarios[i].reset(new QUESO::GslVector(configSpace.zeroVector())); // 'x_{i+1}' in paper
    experimentVecs     [i].reset(new QUESO::GslVector(experimentSpace.zeroVector()));
  }


  // Experiment observation error

  // True observation error for each experiment uses this R among responses
  //   0.0025   0.002
  //   0.002    0.0025
  QUESO::GslMatrix covarianceR(experimentSpace.zeroVector());
  for (unsigned int i = 0; i < 2; ++i)
    covarianceR(i, i) = 0.0025;
  for (unsigned int i = 0; i < 1; ++i)
    covarianceR(i, i+1) = covarianceR(i+1, i) = 0.002;

  // Populate the totalExperimentSpace covariance matrix
  std::vector<const QUESO::GslMatrix* > vec_covmat_ptrs(numExperiments,
                                                        &covarianceR);
  experimentMat->fillWithBlocksDiagonally(0, 0, vec_covmat_ptrs, true, true);

  // Read in data and store the standard deviation of the simulation
  // data (ignored for now).
  const char * test_srcdir = std::getenv("srcdir");

  std::string simInputFileName = "test_gpmsa/sim_mv_small.dat";
  std::string expInputFileName = "test_gpmsa/y_exp_mv_small.txt";
  std::string solInputFileName = "test_gpmsa/regression_solution_pdf_mv_small.txt";

  if (test_srcdir) {
    simInputFileName = test_srcdir + ('/' + simInputFileName);
    expInputFileName = test_srcdir + ('/' + expInputFileName);
    solInputFileName = test_srcdir + ('/' + solInputFileName);
  }

  readData(simInputFileName,
           expInputFileName,
           simulationScenarios,
           paramVecs,
           outputVecs,
           experimentScenarios,
           experimentVecs,
           *experimentMat);

  // Add simulation and experimental data
  gpmsaFactory.addSimulations(simulationScenarios, paramVecs, outputVecs);
  gpmsaFactory.addExperiments(experimentScenarios, experimentVecs, experimentMat);

  QUESO::GenericVectorRV<> postRv( "post_",
      gpmsaFactory.prior().imageSet().vectorSpace());

  QUESO::GslVector point(
      gpmsaFactory.prior().imageSet().vectorSpace().zeroVector());

  std::cout << "Point has size: " << point.sizeLocal() << std::endl;
  queso_require_equal_to(point.sizeLocal(), 10);

  unsigned int num_lines = 10000;
  std::vector<double> expected_log_likelihoods(num_lines);
  std::vector<double> expected_log_priors(num_lines);
  std::vector<double> expected_log_posteriors(num_lines);
  std::vector<double> computed_log_likelihoods(num_lines);
  std::vector<double> computed_log_priors(num_lines);
  std::vector<double> computed_log_posteriors(num_lines);

  // File containing x, theta, y (vertical concatenation of lhs.txt, y_mod.txt)
  std::ifstream solution_data;
  open_data_file(solInputFileName, solution_data);
  for (unsigned int i = 0; i < num_lines; i++) {
    for (unsigned int j = 0; j < point.sizeLocal(); ++j)
      solution_data >> point[j];

    double log_pdf;
    solution_data >> log_pdf;
    expected_log_likelihoods[i] = log_pdf;

    solution_data >> log_pdf;
    expected_log_priors[i] = log_pdf;

    solution_data >> log_pdf;
    expected_log_posteriors[i] = log_pdf;

    log_pdf = (gpmsaFactory.prior().pdf().lnValue(point));
    computed_log_priors[i] = log_pdf;

    log_pdf = (gpmsaFactory.getGPMSAEmulator().lnValue(point, NULL, NULL, NULL, NULL));
    computed_log_likelihoods[i] = log_pdf;
  }
  solution_data.close();

  // We subtract off the first difference in log priors because the
  // normalisation constants between our implementation and Matlab's
  // implementation are not the same.  We expect them to only be the same up to
  // an additive constant.
  double initial_diff_prior = expected_log_priors[0] - computed_log_priors[0];
  for (unsigned int i = 0; i < expected_log_priors.size(); i++) {
    double diff_prior = expected_log_priors[i] - computed_log_priors[i] - initial_diff_prior;
    queso_require_less_equal_msg(std::abs(diff_prior), TOL, "computed log prior differs too much from expected");
  }

  // We don't subtract off the initial because the normalisation constants
  // between our implementation and Matlab's implementation appear to be the
  // same
  for (unsigned int i = 0; i < expected_log_likelihoods.size(); i++) {
    double diff_likelihood = expected_log_likelihoods[i] - computed_log_likelihoods[i];
    queso_require_less_equal_msg(std::abs(diff_likelihood), TOL, "computed log likelihood differs too much from expected");
  }
}

int main(int argc, char ** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: argv[0] gpmsa_<case>.txt\n";
    return 1;
  }

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);

  // Step 1: Set up QUESO environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);
#else
  QUESO::FullEnvironment env(argv[1], "", NULL);
#endif

  std::string input_file(argv[1]);
  if (input_file.find("scalar") != std::string::npos)
    run_scalar(env);
  else if  (input_file.find("mv") != std::string::npos)
    run_multivariate(env);
  else {
    std::cerr << "Unknown case with input: " << input_file <<"\n";
    return 1;
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
