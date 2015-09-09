#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/VectorSpace.h>
#include <queso/GaussianVectorRV.h>
#include <queso/GslMatrix.h>

#define PI 3.14159265358979323846

#define QUESO_REQUIRE_CLOSE(a, b, c) do { if (!require_close(a, b, c)) { \
                                            std::cerr << "FAILED: " << a \
                                                      << " and " << b \
                                                      << " differ by " << c \
                                                      << " in the relative " \
                                                      << "sense." \
                                                      << std::endl; \
                                            queso_error(); \
                                          } \
                                        } while (0)

using namespace QUESO;

int require_close(double a, double b, double tol) {
  return (std::abs(a - b) / std::abs(b) > tol) ? 0 : 1;
}

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  EnvOptionsValues envOptionsValues;
#ifdef QUESO_HAS_MPI
  FullEnvironment env(MPI_COMM_WORLD, "", "", &envOptionsValues);
#else
  FullEnvironment env("", "", &envOptionsValues);
#endif

  VectorSpace<GslVector, GslMatrix> imageSpace(env, "test_space", 2, NULL);
  Map eMap(2, 0, env.fullComm());

  GslVector imageMinVal(env, eMap, -INFINITY);
  GslVector imageMaxVal(env, eMap,  INFINITY);

  BoxSubset<GslVector, GslMatrix> domain("domain", imageSpace, imageMinVal,
      imageMaxVal);

  GslVector initExpectedValues(env, eMap, 0.0);
  GslMatrix initCovMatrix(env, eMap, 1.0);

  GslVector finalExpectedValues(env, eMap, 1.0);
  GslMatrix finalCovMatrix(env, eMap, 3.0);

  GslVector testValues(env, eMap, 0.0);

  GaussianVectorRV<GslVector, GslMatrix> gaussianRV("test_rv", domain,
      initExpectedValues, initCovMatrix);
  double normalisingConstant = 1.0 / (2.0 * PI);

  double tolClose = 1e-13;

  //***********************************************************************
  // Test pdf
  //***********************************************************************

  // mean = [0; 0], var = [1; 1], testValues = [0; 0]
  QUESO_REQUIRE_CLOSE(gaussianRV.pdf().actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant, tolClose);

  // change mean and check that new pdf is correct
  gaussianRV.updateLawExpVector(finalExpectedValues);
  QUESO_REQUIRE_CLOSE(gaussianRV.pdf().actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);

  //***********************************************************************
  // Test realizer
  // NOTE: just calls it... doesn't check values
  //***********************************************************************
  GslVector myRealization(testValues);
  gaussianRV.realizer().realization(myRealization);

  std::cout << myRealization;

  // finalize
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
