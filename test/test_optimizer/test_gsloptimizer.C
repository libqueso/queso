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
    return std::exp(this->lnValue(domainVector));
  }

  virtual double lnValue(const V & domainVector, V & gradVector) const {
    gradVector[0] = -2.0 * (domainVector[0]-1);
    gradVector[1] = -2.0 * (domainVector[1]-2);
    gradVector[2] = -2.0 * (domainVector[2]-3);

    // Mean = (1,2)
    return this->lnValue(domainVector);
  }

  virtual double lnValue(const V & domainVector) const {
    // Mean = (1,2)
    return -( (domainVector[0]-1)*(domainVector[0]-1) +
              (domainVector[1]-2)*(domainVector[1]-2) +
              (domainVector[2]-3)*(domainVector[2]-3) );
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
      "space_", 3, NULL);

  QUESO::GslVector minBound(paramSpace.zeroVector());
  minBound[0] = -10.0;
  minBound[1] = -10.0;
  minBound[2] = -10.0;

  QUESO::GslVector maxBound(paramSpace.zeroVector());
  maxBound[0] = 10.0;
  maxBound[1] = 10.0;
  maxBound[2] = 10.0;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> domain("", paramSpace,
      minBound, maxBound);

  ObjectiveFunction<QUESO::GslVector, QUESO::GslMatrix> objectiveFunction(
      "", domain);

  QUESO::GslVector initialPoint(paramSpace.zeroVector());
  initialPoint[0] = 9.0;
  initialPoint[1] = -9.0;
  initialPoint[1] = -1.0;

  QUESO::GslOptimizer optimizer(objectiveFunction);

  double tol = 1.0e-10;
  optimizer.setTolerance(tol);
  optimizer.set_solver_type(QUESO::GslOptimizer::STEEPEST_DESCENT);

  QUESO::OptimizerMonitor monitor(env);
  monitor.set_display_output(true,true);

  std::cout << "Solving with Steepest Decent" << std::endl;
  optimizer.minimize(&monitor);

  if (std::abs( optimizer.minimizer()[0] - 1.0) > tol) {
    std::cerr << "GslOptimize failed.  Found minimizer at: " << optimizer.minimizer()[0]
              << std::endl;
    std::cerr << "Actual minimizer is 1.0" << std::endl;
    queso_error();
  }

  std::string nm = "nelder_mead2";
  optimizer.set_solver_type(nm);
  monitor.reset();
  monitor.set_display_output(true,true);

  std::cout << std::endl << "Solving with Nelder Mead" << std::endl;
  optimizer.minimize(&monitor);

  monitor.print(std::cout,false);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
