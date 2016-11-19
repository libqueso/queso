#include <queso/GenericScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>

template<class V, class M>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain)
    : QUESO::BaseScalarFunction<V, M>(prefix, domain)
  {
    // Setup here
  }

  virtual ~Likelihood()
  {
    // Deconstruct here
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    if (gradVector != NULL) {
      (*gradVector)[0] = -domainVector[0];
    }

    return -0.5 * domainVector[0] * domainVector[0];
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }

private:
  // Maybe store the observed data, y, here.
};

int main(int argc, char ** argv) {
  std::string inputFileName = "test_algorithms/input_test_mala.txt";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir) {
    inputFileName = test_srcdir + ('/' + inputFileName);
  }

  MPI_Init(&argc, &argv);

  // Step 0 of 5: Set up environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);

  unsigned int dim = 1;

  // Step 1 of 5: Instantiate the parameter space
  QUESO::VectorSpace<> paramSpace(env, "param_", dim, NULL);

  double min_val = -10000.0;
  double max_val =  10000.0;

  // Step 2 of 5: Set up the prior
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(min_val);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  // Uniform prior here.  Could be a different prior.
  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  // Step 3 of 5: Set up the likelihood using the class above
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood("llhd_", paramDomain);

  // Step 4 of 5: Instantiate the inverse problem
  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

  // Step 5 of 5: Solve the inverse problem
  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  // Initial condition of the chain
  for (unsigned int i = 0; i < dim; i++) {
    paramInitials[i] = 0.0;
  }

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  for (unsigned int i = 0; i < dim; i++) {
    // Might need to tweak this
    proposalCovMatrix(i, i) = 1.0;
  }

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  QUESO::GslVector draw(paramSpace.zeroVector());
  QUESO::GslVector mean(paramSpace.zeroVector());
  QUESO::GslVector sumsq(paramSpace.zeroVector());
  QUESO::GslVector delta(paramSpace.zeroVector());
  QUESO::GslVector var(paramSpace.zeroVector());

  unsigned int num_samples = 10000;
  for (unsigned int i = 1; i < num_samples + 1; i++) {
    postRv.realizer().realization(draw);
    for (unsigned int j = 0; j < dim; j++) {
      delta[j] = draw[j] - mean[j];
      mean[j] += (double) delta[j] / i;
      sumsq[j] += delta[j] * (draw[j] - mean[j]);
    }
  }

  for (unsigned int j = 0; j < dim; j++) {
    // This is the sample variance
    var[j] = sumsq[j] / (num_samples - 1);
  }

  double sigmasq = 1.0;
  double mean_min = -3.0 * sqrt(sigmasq) / sqrt(num_samples);
  double mean_max =  3.0 * sqrt(sigmasq) / sqrt(num_samples);

  int return_val = 0;
  for (unsigned int j = 0; j < dim; j++) {
    if (mean[j] < mean_min || mean[j] > mean_max) {
      return_val = 1;
      break;
    }
  }

  double var_min = sigmasq - 3.0 * sqrt(2.0 * sigmasq * sigmasq / (num_samples - 1));
  double var_max = sigmasq + 3.0 * sqrt(2.0 * sigmasq * sigmasq / (num_samples - 1));

  // var[j] should be approximately ~ N(sigma^2, 2 sigma^4 / (num_samples - 1))
  for (unsigned int j = 0; j < dim; j++) {
    if (var[j] < var_min || var[j] > var_max) {
      std::cout << "var" << std::endl;
      std::cout << var[j] << std::endl;
      return_val = 1;
      break;
    }
  }

  MPI_Finalize();

  return return_val;
}
