#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>

template<class V = QUESO::GslVector, class M = QUESO::GslMatrix>
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

  virtual double lnValue(const V & domainVector) const
  {
    // 1) Run the forward code at the point domainVector
    //    domainVector[0] is the first element of the parameter vector
    //    damainVector[1] is the second element of the parameter vector
    //    and so on
    //
    // 2) Compare to data, y
    //    Usually we compute something like
    //    || MyModel(domainVector) - y ||^2 / (sigma * sigma)
    //
    // 3) Return below

    double misfit = 0.0;

    return -0.5 * misfit;
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector));
  }

private:
  // Maybe store the observed data, y, here.
};

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  // Step 0 of 5: Set up environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  // Step 1 of 5: Instantiate the parameter space
  QUESO::VectorSpace<> paramSpace(env, "param_", 1, NULL);

  double min_val = 0.0;
  double max_val = 1.0;

  // Step 2 of 5: Set up the prior
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(min_val);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  // Uniform prior here.  Could be a different prior.
  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  // Step 3 of 5: Set up the likelihood using the class above
  Likelihood<> lhood("llhd_", paramDomain);

  // Step 4 of 5: Instantiate the inverse problem
  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

  // Step 5 of 5: Solve the inverse problem
  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  // Initial condition of the chain
  paramInitials[0] = 0.0;
  paramInitials[1] = 0.0;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  for (unsigned int i = 0; i < 1; i++) {
    // Might need to tweak this
    proposalCovMatrix(i, i) = 0.1;
  }

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  MPI_Finalize();

  return 0;
}
