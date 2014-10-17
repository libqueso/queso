#include <iostream>
#include <cmath>
#include <queso/asserts.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/ScalarFunction.h>
#include <queso/GslOptimizer.h>
#include <queso/OptimizerMonitor.h>

template <class V, class M>
class ObjectiveFunction : public QUESO::BaseScalarFunction<V, M> {
public:
  ObjectiveFunction(const char * prefix,
      const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet) {
      // Do nothing
    }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const {
    return std::exp(this->lnValue(domainVector, domainDirection,
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
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "space_", 1, NULL);

  QUESO::GslVector minBound(paramSpace.zeroVector());
  minBound[0] = -10.0;

  QUESO::GslVector maxBound(paramSpace.zeroVector());
  maxBound[0] = 10.0;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> domain("", paramSpace,
      minBound, maxBound);

  ObjectiveFunction<QUESO::GslVector, QUESO::GslMatrix> objectiveFunction(
      "", domain);

  QUESO::GslVector initialPoint(paramSpace.zeroVector());
  initialPoint[0] = 9.0;

  QUESO::GslOptimizer optimizer(objectiveFunction);
  optimizer.setInitialPoint(initialPoint);
  optimizer.minimize();

  if (std::abs((optimizer.minimizer())[0]) > 1e-10) {
    std::cerr << "GslOptimize failed.  Found minimizer at: "
              << (optimizer.minimizer())[0]
              << std::endl;
    std::cerr << "Actual minimizer is 0.0" << std::endl;
    queso_error();
  }

  optimizer.setInitialPoint(initialPoint);
  optimizer.set_solver_type(QUESO::GslOptimizer::NELDER_MEAD2);

  QUESO::OptimizerMonitor monitor;
  monitor.set_display_output(true,true);

  optimizer.minimize(&monitor);

  if (std::abs((optimizer.minimizer())[0]) > 1e-10) {
    std::cerr << "GslOptimize failed.  Found minimizer at: "
              << (optimizer.minimizer())[0]
              << std::endl;
    std::cerr << "Actual minimizer is 0.0" << std::endl;
    queso_error();
  }

  return 0;
}
