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

// Common function for both objectives
template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
double
f(const V & domainVector)
{
  // Mean = (1, 2, 3)
  return -((domainVector[0] - 1)*(domainVector[0] - 1) +
           (domainVector[1] - 2)*(domainVector[1] - 2) +
           (domainVector[2] - 3)*(domainVector[2] - 3));
}

// Function without a gradient (falls back to finite difference approximation
template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class ObjectiveFunctionWithoutGradient : public QUESO::BaseScalarFunction<V, M> {
public:
  ObjectiveFunctionWithoutGradient(const char * prefix,
      const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet) {
      // Do nothing
    }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const {
    return std::exp(this->lnValue(domainVector, domainDirection,
          gradVector, hessianMatrix, hessianEffect));
  }

  virtual double lnValue(const V & domainVector) const {
    return f(domainVector);
  }

  using QUESO::BaseScalarFunction<V, M>::lnValue;
};

// Objective function with analytical gradient (same function basically)
template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class ObjectiveFunctionWithGradient : public QUESO::BaseScalarFunction<V, M> {
public:
  ObjectiveFunctionWithGradient(const char * prefix,
      const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet) {
      // Do nothing
    }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const {
    return std::exp(this->lnValue(domainVector, domainDirection,
          gradVector, hessianMatrix, hessianEffect));
  }

  virtual double lnValue(const V & domainVector) const {
    return f(domainVector);
  }

  virtual double lnValue(const V & domainVector, V & gradVector) const {
    gradVector[0] = -2.0 * (domainVector[0]-1);
    gradVector[1] = -2.0 * (domainVector[1]-2);
    gradVector[2] = -2.0 * (domainVector[2]-3);

    return f(domainVector);
  }

  using QUESO::BaseScalarFunction<V, M>::lnValue;
};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);
#else
  QUESO::FullEnvironment env("", "", NULL);
#endif

  QUESO::VectorSpace<> paramSpace(env, "space_", 3, NULL);

  QUESO::GslVector minBound(paramSpace.zeroVector());
  minBound[0] = -10.0;
  minBound[1] = -10.0;
  minBound[2] = -10.0;

  QUESO::GslVector maxBound(paramSpace.zeroVector());
  maxBound[0] = 10.0;
  maxBound[1] = 10.0;
  maxBound[2] = 10.0;

  QUESO::BoxSubset<> domain("", paramSpace, minBound, maxBound);

  ObjectiveFunctionWithoutGradient<> objectiveFunctionWithoutGradient("", domain);
  ObjectiveFunctionWithGradient<> objectiveFunctionWithGradient("", domain);

  QUESO::GslVector point(paramSpace.zeroVector());
  point[0] = 9.0;
  point[1] = -9.0;
  point[1] = -1.0;

  double value1;
  value1 = objectiveFunctionWithoutGradient.lnValue(point);

  double value2;
  value2 = objectiveFunctionWithGradient.lnValue(point);

  double tol = 1e-13;
  if (std::abs(value1 - value2) > tol) {
    std::string msg;

    msg += "Objective function values with and without gradients do not match";
    queso_error_msg(msg);
  }

  QUESO::GslVector grad1(paramSpace.zeroVector());
  QUESO::GslVector grad2(paramSpace.zeroVector());

  objectiveFunctionWithoutGradient.lnValue(point, grad1);
  objectiveFunctionWithGradient.lnValue(point, grad2);

  std::cout << "grad1: " << std::endl;
  std::cout << grad1 << std::endl;

  std::cout << "grad2: " << std::endl;
  std::cout << grad2 << std::endl;

  grad1 -= grad2;

  std::cout << "diff: " << std::endl;
  std::cout << grad1 << std::endl;

  double norm2 = grad1.norm2();  // Should be close to zero?

  std::cout << "norm2: " << std::endl;
  std::cout << norm2 << std::endl;

  tol = 1e-5;  // Because this is a finite difference approximation, we use a
               // much smaller tolerance for the following comparison.

  if (std::abs(norm2) > tol) {
    std::string msg;

    msg += "Objective function gradients (finite diff and analytical) do not";
    msg += " match.";

    queso_error_msg(msg);
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
