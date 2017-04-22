#include <iostream>
#include <cmath>
#include <queso/asserts.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/ScalarFunction.h>
#include <queso/GenericVectorRV.h>
#include <queso/GslOptimizer.h>
#include <queso/StatisticalInverseProblem.h>

template <class V, class M>
class Likelihood : public QUESO::BaseScalarFunction<V, M> {
public:
  Likelihood(const char * prefix,
      const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet) {
      // Do nothing
    }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const {
    return std::exp(this->lnValue(domainVector));
  }

  virtual double lnValue(const V & domainVector) const {
    return -domainVector[0] * domainVector[0];
  }

  virtual double lnValue(const V & domainVector, V & gradVector) const {
    // Need to check if NULL because QUESO will somtimes call this with a
    // NULL pointer
    gradVector[0] = -2.0 * domainVector[0];

    return this->lnValue(domainVector);
  }
};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);
#else
  QUESO::FullEnvironment env("", "", NULL);
#endif

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "space_", 1, NULL);

  QUESO::GslVector minBound(paramSpace.zeroVector());
  minBound[0] = -10.0;

  QUESO::GslVector maxBound(paramSpace.zeroVector());
  maxBound[0] = 10.0;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> domain("", paramSpace,
      minBound, maxBound);

  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> prior("", domain);

  Likelihood<QUESO::GslVector, QUESO::GslMatrix> likelihood("", domain);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix> posterior("",
      domain);

  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> ip("",
      NULL, prior, likelihood, posterior);

  QUESO::GslVector initialValues(paramSpace.zeroVector());
  initialValues[0] = 9.0;

  QUESO::GslMatrix proposalCovarianceMatrix(paramSpace.zeroVector());
  proposalCovarianceMatrix(0, 0) = 1.0;

  ip.seedWithMAPEstimator();
  ip.solveWithBayesMetropolisHastings(NULL, initialValues,
      &proposalCovarianceMatrix);

  // The first sample should be the seed
  QUESO::GslVector first_sample(paramSpace.zeroVector());
  posterior.realizer().realization(first_sample);

  // Tight tolerance for analytical derivative
  if (std::abs(first_sample[0]) > 1e-10) {
    std::cerr << "seedWithMAPEstimator failed.  Seed was: " << first_sample[0]
              << std::endl;
    std::cerr << "Actual seed should be 0.0" << std::endl;
    queso_error();
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
