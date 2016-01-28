#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/VectorSpace.h>
#include <queso/GaussianJointPdf.h>
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
  // Initialize
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  EnvOptionsValues envOptionsValues;
#ifdef QUESO_HAS_MPI
  FullEnvironment env(MPI_COMM_WORLD, "", "", &envOptionsValues);
#else
  FullEnvironment env("", "", &envOptionsValues);
#endif

  VectorSpace<GslVector, GslMatrix> domainSpace(env, "test_space", 2, NULL);
  Map eMap(2, 0, env.fullComm());

  GslVector domainMinVal(env, eMap, -1e30);
  GslVector domainMaxVal(env, eMap,  1e30);

  BoxSubset<GslVector, GslMatrix> domain("domain", domainSpace, domainMinVal, domainMaxVal);

  GaussianJointPdf<GslVector, GslMatrix>* gaussianPdf;
  double tolClose = 1e-13;

  //***********************************************************************
  // Tests for diagonal covariance matrix
  //***********************************************************************

  // mean = [0; 0], var = [1; 1]
  GslVector expectedVal(env, eMap, 0.0);
  GslVector varianceVal(env, eMap, 1.0);
  GslVector testValues(env, eMap, 0.0);

  gaussianPdf = new GaussianJointPdf<GslVector, GslMatrix>("test_pdf", domain, expectedVal, varianceVal);
  double normalisingConstant = 1.0 / (2.0 * PI);
  double logNormalisingConstant = std::log(normalisingConstant);

  testValues[0] = testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant, tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  delete gaussianPdf;

  // mean = [0; 0], var = [0.25; 0.5]
  varianceVal[0] = 0.25; varianceVal[1] = 0.5;

  gaussianPdf = new GaussianJointPdf<GslVector, GslMatrix>("test_pdf", domain, expectedVal, varianceVal);
  normalisingConstant = 1.0 / (2.0 * PI * std::sqrt(0.125));
  logNormalisingConstant = std::log(normalisingConstant);

  testValues[0] = testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant, tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-3.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 3.0, tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-3.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 3.0, tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-2.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 2.0, tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-3.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 3.0, tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-3.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 3.0, tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-2.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 2.0, tolClose);

  // mean = [1.0; -0.5], var = [0.25; 0.5]
  expectedVal[0] = 1.0; expectedVal[1] = -0.5;

  gaussianPdf->updateLawExpVector(expectedVal); // just update expected value (don't reallocate everything)

  testValues[0] = 0.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-2.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 2.25, tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-2.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 2.25, tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-4.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 4.25, tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-10.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 10.25, tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-8.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 8.25, tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-8.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 8.25, tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-2.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 2.25, tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.25, tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.25), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.25, tolClose);

  delete gaussianPdf;

  //***********************************************************************
  // Tests for general covariance matrix
  //***********************************************************************

  // mean = [0; 0], covar = [1, 0; 0, 1], i.e. same as first case for diagonal matrices
  expectedVal[0] = expectedVal[1] = 0.0;

  GslMatrix covMatrix(env, eMap, 1.0);  // actually diagonal

  gaussianPdf = new GaussianJointPdf<GslVector, GslMatrix>("test_pdf", domain, expectedVal, covMatrix);
  normalisingConstant = 1.0 / (2.0 * PI);
  logNormalisingConstant = std::log(normalisingConstant);

  testValues[0] = testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant, tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.5), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.5, tolClose);

  //delete gaussianPdf;

  // mean = [0, 0], covar = [2, 1; 1, 2];
  covMatrix(0,0) = 2.0; covMatrix(0,1) = 1.0;
  covMatrix(1,0) = 1.0; covMatrix(1,1) = 2.0;

  gaussianPdf->updateLawCovMatrix(covMatrix);
  normalisingConstant = 1.0 / (2.0 * PI * std::sqrt(3.0));
  logNormalisingConstant = std::log(normalisingConstant);

  testValues[0] = testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant, tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0/3.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - (1.0/3.0), tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0/3.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - (1.0/3.0), tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.0), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.0, tolClose);

  // mean = [1.0, -0.5], covar = [2, 1; 1, 2];
  expectedVal[0] = 1.0; expectedVal[1] = -0.5;

  gaussianPdf->updateLawExpVector(expectedVal); // just update expected value (don't reallocate everything)

  testValues[0] = testValues[1] = 0.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-5.833333333333333e-01), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 5.833333333333333e-01, tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-0.75), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 0.75, tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-1.583333333333333e+00), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 1.583333333333333e+00, tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  QUESO_REQUIRE_CLOSE(gaussianPdf->actualValue(testValues, NULL, NULL, NULL, NULL), normalisingConstant * std::exp(-3.083333333333333e+00), tolClose);
  QUESO_REQUIRE_CLOSE(gaussianPdf->lnValue(testValues, NULL, NULL, NULL, NULL), logNormalisingConstant - 3.083333333333333e+00, tolClose);

  delete gaussianPdf;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
