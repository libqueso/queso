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
    return std::exp(this->actualValue(domainVector, domainDirection,
          gradVector, hessianMatrix, hessianEffect));
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const {

    // Need to check if NULL because QUESO will somtimes call this with a
    // NULL pointer
    if (gradVector != NULL) {
      (*gradVector)[0] = -2.0 * domainVector[0];
    }

    return -(domainVector[0] * domainVector[0]);
  }
};

int main(int argc, char ** argv) {
  std::string inputFileName = "test_optimizer/input_test_optimizer_input_parameters";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir) {
    inputFileName = test_srcdir + ('/' + inputFileName);
  }

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);
#else
  QUESO::FullEnvironment env(inputFileName, "", NULL);
#endif

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "space_", 1, NULL);

  QUESO::GslVector minBound(paramSpace.zeroVector());
  minBound[0] = -10.0;

  QUESO::GslVector maxBound(paramSpace.zeroVector());
  maxBound[0] = 10.0;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> domain("", paramSpace,
      minBound, maxBound);

  Likelihood<QUESO::GslVector, QUESO::GslMatrix> likelihood("", domain);

  QUESO::OptimizerOptions options(&env);

  QUESO::GslOptimizer optimizer(options, likelihood);

  queso_require_equal_to_msg(optimizer.getMaxIterations(), 1, "m_maxIterations read incorrectly");
  queso_require_equal_to_msg(optimizer.getFiniteDifferenceStepSize(), 2.0, "m_finiteDifferenceStepSize read incorrectly");
  queso_require_equal_to_msg(optimizer.getTolerance(), 3.0, "m_tolerance read incorrectly");

  // // Tight tolerance for analytical derivative
  // if (std::abs(first_sample[0]) > 1e-10) {
  //   std::cerr << "seedWithMAPEstimator failed.  Seed was: " << first_sample[0]
  //             << std::endl;
  //   std::cerr << "Actual seed should be 0.0" << std::endl;
  //   queso_error();
  // }

  return 0;
}
