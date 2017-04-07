#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalInverseProblemOptions.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
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

  virtual double lnValue(const V & domainVector) const
  {
    double x1 = domainVector[0];
    double x2 = domainVector[1];

    return -0.5 * (x1 * x1 + x2 * x2);
  }

  virtual double actualValue(const V & domainVector, const V *, V *, M *,
      V *) const
  {
    return std::exp(this->lnValue(domainVector));
  }
};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues envOptions;
  envOptions.m_numSubEnvironments = 1;
  envOptions.m_subDisplayFileName = ".";
  envOptions.m_subDisplayAllowAll = 0;
  envOptions.m_displayVerbosity = 0;
  envOptions.m_seed = 0;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &envOptions);
#else
  QUESO::FullEnvironment env("", "", &envOptions);
#endif

  unsigned int dim = 2;
  QUESO::VectorSpace<> paramSpace(env, "param_", dim, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  double min_val = -10.0;
  double max_val = 10.0;
  paramMins.cwSet(min_val);
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::GslVector mean(paramSpace.zeroVector());

  QUESO::GslMatrix cov(paramSpace.zeroVector());
  cov(0, 0) = 1.0;
  cov(0, 1) = 0.0;
  cov(1, 0) = 0.0;
  cov(1, 1) = 1.0;

  QUESO::GaussianVectorRV<> prior("prior_", paramDomain, mean, cov);

  Likelihood<> lhood("llhd_", paramDomain);

  QUESO::GenericVectorRV<> post("post_", paramSpace);

  QUESO::SipOptionsValues sipOptions;
  sipOptions.m_computeSolution = 1;
  sipOptions.m_dataOutputFileName = ".";
  sipOptions.m_dataOutputAllowedSet.clear();
  sipOptions.m_dataOutputAllowedSet.insert(0);

  QUESO::StatisticalInverseProblem<> ip1("", &sipOptions, prior, lhood, post);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.0;
  paramInitials[1] = 0.0;

  QUESO::GslMatrix proposalCovMatrix(paramInitials);

  proposalCovMatrix(0, 0) = 1.0;
  proposalCovMatrix(0, 1) = 0.0;
  proposalCovMatrix(1, 0) = 0.0;
  proposalCovMatrix(1, 1) = 1.0;

  ip1.solveWithBayesMetropolisHastings(NULL, paramInitials,
      &proposalCovMatrix);


  // Now do the second IP
  env.resetSeed(0);

  QUESO::GenericVectorRV<> post2("post_", paramSpace);
  QUESO::StatisticalInverseProblem<> ip2("", &sipOptions, prior, lhood, post2);

  ip2.solveWithBayesMetropolisHastings();

  QUESO::GslVector sample1(paramSpace.zeroVector());
  QUESO::GslVector sample2(paramSpace.zeroVector());

  for (unsigned int i = 0; i < 100; i++) {
    ip1.postRv().realizer().realization(sample1);
    ip2.postRv().realizer().realization(sample2);

    queso_require_equal_to(sample1, sample2);
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
