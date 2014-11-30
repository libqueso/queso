#include <queso/Environment.h>
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
  }

  virtual ~Likelihood()
  {
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    double misfit = 1.0;

    return -0.5 * misfit;
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }
};

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 1, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  double min_val = 0.0;
  double max_val = 1.0;
  paramMins.cwSet(min_val);
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorRv("prior_",
      paramDomain);

  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood("llhd_", paramDomain);
  
  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix>
    postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix>
    ip("", NULL, priorRv, lhood, postRv);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.0;
  
  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0, 0) = 0.1;
  
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  MPI_Finalize();
 
  return 0;
}
