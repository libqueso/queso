#include <uqEnvironment.h>
#include <uqStatisticalInverseProblem.h>

int main(int argc, char **argv)
{
  int num_times = 40;
  double obs_var = 1e-2;
  double burger_speed = 1.0;
  double dt = 0.01;

  MPI_Init(&argc, &argv);

  // Set up the environment
  uqFullEnvironmentClass *env = new uqFullEnvironmentClass(MPI_COMM_WORLD,
      "my_input_file.in", "", NULL);

  // Set up the prior -- only Gaussians are supported for now
  // A gaussian prior is specified by a mean function and a covariance operator

  // Mean function
  FUNCTION *mu = new FUNCTION(something here);
  // Maybe the default constructor should default to the zero function?

  // Inverse covariance operator: e.g., the inverse laplacian
  OPERATOR *C = new OPERATOR("Laplacian");

  // Set up the Gaussian measure. s is the power that determines regularity > 0
  uqGaussianRF *prior = new GAUSSIAN_MEASURE(mu, C, s);

  // Set up the likelihood. It takes a field
  Likelihood *likelihood = new Likelihood(*env, true_initial_condition,
      burger_speed, num_times, dt, obs_var);
  likelihood->synthesize_data();

  // Set up the inverse problem
  uqSipOptionsValuesClass opts;
  uqStatisticalInverseProblemClass *ip = new
    uqStatisticalInverseProblemClass("", &opts, *prior, *likelihood);

  // MCMC seed
  FUNCTION initialVals();  // Zero as the MCMC seed
  ip->solveWithBayesMetropolisHastings (NULL, initialVals);

  delete env;
  delete mu;
  delete C;
  delete prior;
  delete likelihood;
  MPI_Finalize();

  return 0;
}
