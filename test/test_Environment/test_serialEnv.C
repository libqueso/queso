#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>

#include <cstdlib>

template<class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain)
    : QUESO::BaseScalarFunction<V, M>(prefix, domain) {}

  virtual ~Likelihood() {}

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
  std::string inputFileName = "test_Environment/input_test_serialEnv";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

  QUESO::FullEnvironment env(inputFileName, "", NULL);

  QUESO::VectorSpace<> paramSpace(env, "param_", 1, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  double min_val = 0.0;
  double max_val = 1.0;
  paramMins.cwSet(min_val);
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  Likelihood<> lhood("llhd_", paramDomain);

  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.0;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0, 0) = 0.1;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  return 0;
}
