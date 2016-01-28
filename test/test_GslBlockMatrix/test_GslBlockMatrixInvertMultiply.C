#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/VectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GslBlockMatrix.h>

#define TOL 1e-10

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment env("", "", &options);
#endif

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 3, NULL);

  // Example RHS
  QUESO::GslVector b(paramSpace.zeroVector());
  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = 3.0;

  // Set up block sizes for observation covariance matrix
  std::vector<unsigned int> blockSizes(2);
  blockSizes[0] = 1;  // First block is 1x1 (scalar)
  blockSizes[1] = 2;  // Second block is 2x2

  // Set up block (identity) matrix with specified block sizes
  QUESO::GslBlockMatrix covariance(blockSizes, b, 1.0);

  // The matrix [[1, 0, 0], [0, 1, 2], [0, 2, 8]]
  // has inverse 0.25 * [[1, 0, 0], [0, 2, -0.5], [0, -0.5, 0.25]]
  covariance.getBlock(0)(0, 0) = 1.0;
  covariance.getBlock(1)(0, 0) = 1.0;
  covariance.getBlock(1)(0, 1) = 2.0;
  covariance.getBlock(1)(1, 0) = 2.0;
  covariance.getBlock(1)(1, 1) = 8.0;

  // Compute solution
  QUESO::GslVector x(paramSpace.zeroVector());
  covariance.invertMultiply(b, x);

  // This is the analytical solution
  QUESO::GslVector sol(paramSpace.zeroVector());
  sol[0] = 1.0;
  sol[1] = 2.5;
  sol[2] = -0.25;

  // So if the solve worked, this sucker should be the zero vector
  sol -= x;

  // So its norm should be zero
  if (sol.norm2() > TOL) {
    std::cerr << "Computed solution:" << std::endl;
    std::cerr << b << std::endl;
    std::cerr << "Actual solution:" << std::endl;
    std::cerr << sol << std::endl;
    queso_error_msg("TEST: GslBlockMatrix::invertMultiply failed.");
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
