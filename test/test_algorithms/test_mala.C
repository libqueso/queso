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
  MPI_Init(&argc, &argv);

  // Step 0 of 5: Set up environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  unsigned int dim = 1;

  // Step 1 of 5: Instantiate the parameter space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", dim, NULL);

  double min_val = -10000.0;
  double max_val =  10000.0;

  // Step 2 of 5: Set up the prior
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(min_val);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  // Uniform prior here.  Could be a different prior.
  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorRv("prior_",
      paramDomain);

  // Step 3 of 5: Set up the likelihood using the class above
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood("llhd_", paramDomain);

  // Step 4 of 5: Instantiate the inverse problem
  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix>
    postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix>
    ip("", NULL, priorRv, lhood, postRv);

  // Step 5 of 5: Solve the inverse problem
  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  // Initial condition of the chain
  for (unsigned int i = 0; i < dim; i++) {
    paramInitials[i] = 0.0;
  }

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  for (unsigned int i = 0; i < dim; i++) {
    // Might need to tweak this
    proposalCovMatrix(i, i) = atof(argv[2]);
  }

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  MPI_Finalize();

  return 0;
}
