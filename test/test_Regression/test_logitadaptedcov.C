#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/UniformVectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/StatisticalInverseProblem.h>

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD,
      "test_Regression/adaptedcov_input.txt", "", NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 2, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(0.0);

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(1.0);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> prior("prior_",
      paramDomain);

  QUESO::GslVector mean(paramSpace.zeroVector());
  mean.cwSet(0.5);

  QUESO::GslMatrix cov(paramSpace.zeroVector());
  cov(0, 0) = 0.05 * 0.05;
  cov(0, 1) = 0.001;
  cov(1, 0) = 0.001;
  cov(1, 1) = 0.05 * 0.05;

  QUESO::GaussianVectorRV<QUESO::GslVector, QUESO::GslMatrix> likelihood(
      "likelihood_", paramDomain, mean, cov);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix> posterior(
      "posterior_", paramDomain);

  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> ip("",
      NULL, prior, likelihood.pdf(), posterior);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials.cwSet(0.5);

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0, 0) = 0.5;
  proposalCovMatrix(0, 1) = 0.0;
  proposalCovMatrix(1, 0) = 0.0;
  proposalCovMatrix(1, 1) = 0.5;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  // ip.sequenceGenerator().transitionKernel().rv(0);

  MPI_Finalize();

  return 0;
}
