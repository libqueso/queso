#include <uqEnvironment.h>
#include <uqVectorSpace.h>
#include <uqGslVector.h>
#include <uqGslMatrix.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

#define TOL 1e-10

using namespace std;

int matrixIsDiag(const QUESO::uqGslMatrixClass &M, double diagValue) {
  // M is assumed to be square
  unsigned int i, j;
  unsigned int n;

  n = M.numRowsLocal();

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) {
        if (std::abs(M(i, i) - diagValue) > TOL) {
          std::cerr << "Matrix not diagonal" << std::endl;
          return 0;
        }
      }
      else {
        if (std::abs(M(i, j)) > TOL) {
          std::cerr << "Matrix not diagonal" << std::endl;
          return 0;
        }
      }
    }
  }
  return 1;
}

void fill2By2Matrix(QUESO::uqGslMatrixClass &M) {
  M(0, 0) = 2.0; M(0, 1) = 3.0;
  M(1, 0) = 2.0; M(1, 1) = 2.0;
}

int main(int argc, char **argv) {
  unsigned int i, j;
  double diagValue = 1.5;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::uqEnvOptionsValuesClass options;
  options.m_numSubEnvironments = 1;

  QUESO::uqFullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new QUESO::uqFullEnvironmentClass(MPI_COMM_WORLD, "", "", &options);
#else
    new QUESO::uqFullEnvironmentClass(0, "", "", &options);
#endif

  QUESO::uqVectorSpaceClass<QUESO::uqGslVectorClass, QUESO::uqGslMatrixClass> *param_space =
    new QUESO::uqVectorSpaceClass<QUESO::uqGslVectorClass, QUESO::uqGslMatrixClass>(*env, "param_", 3, NULL);

  QUESO::uqGslVectorClass v(param_space->zeroVector());
  QUESO::uqGslMatrixClass M1(*env, v.map(), diagValue);
  QUESO::uqGslMatrixClass M2(v, diagValue);

  if (!matrixIsDiag(M1, diagValue)) {
    std::cerr << "matrix not diagonal" << std::endl;
    return 1;
  }

  if (!matrixIsDiag(M2, diagValue)) {
    std::cerr << "matrix not diagonal" << std::endl;
    return 1;
  }

  M1 /= diagValue;
  if (!matrixIsDiag(M1, 1.0)) {
    std::cerr << "operator /= failed" << std::endl;
    return 1;
  }

  M2 -= M1;
  if (!matrixIsDiag(M2, 0.5)) {
    std::cerr << "operator -= failed" << std::endl;
    return 1;
  }

  if (std::abs(M2.normFrob() - std::sqrt(3.0 * 0.5 * 0.5)) > TOL) {
    std::cerr << "frobenius norm failed" << "  " << M2.normFrob() << std::endl;
    return 1;
  }

  if (std::abs(M2.normMax() - 0.5) > TOL) {
    std::cerr << "max norm failed " << M2.normMax() << std::endl;
    return 1;
  }

  if (std::abs(M2.max() - 0.5) > TOL) {
    std::cerr << "max failed" << std::endl;
    return 1;
  }

  M2.cwSet(1.0);
  for (i = 0; i < M2.numRowsLocal(); i++) {
    for (j = 0; j < M2.numCols(); j++) {
      if (std::abs(M2(i, j) - 1.0) > TOL) {
        std::cerr << "cwSet failed" << std::endl;
        return 1;
      }
    }
  }

  M2.zeroLower(false);
  for (i = 0; i < M2.numRowsLocal(); i++) {
    for (j = 0; j < M2.numCols(); j++) {
      if (i > j) {
        if (std::abs(M2(i, j)) > TOL) {
          std::cerr << "zero lower failed 1" << std::endl;
          return 1;
        }
      }
      else {
        if (std::abs(M2(i, j) - 1.0) > TOL) {
          std::cerr << "zero lower failed 2" << std::endl;
        }
      }
    }
  }

  M2.zeroLower(true);
  for (i = 0; i < M2.numRowsLocal(); i++) {
    for (j = 0; j < M2.numCols(); j++) {
      if (i >= j) {
        if (std::abs(M2(i, j)) > TOL) {
          std::cerr << "zero lower failed 1" << std::endl;
          return 1;
        }
      }
      else {
        if (std::abs(M2(i, j) - 1.0) > TOL) {
          std::cerr << "zero lower failed 2" << std::endl;
        }
      }
    }
  }

  M2.cwSet(1.0);
  M2.zeroUpper(true);
  for (i = 0; i < M2.numRowsLocal(); i++) {
    for (j = 0; j < M2.numCols(); j++) {
      if (i <= j) {
        if (std::abs(M2(i, j)) > TOL) {
          std::cerr << "zero upper failed 1" << std::endl;
          return 1;
        }
      }
      else {
        if (std::abs(M2(i, j) - 1.0) > TOL) {
          std::cerr << "zero lower failed 2" << std::endl;
        }
      }
    }
  }

  if (std::abs(M2.determinant()) > TOL) {
    std::cerr << "determinant failed" << std::endl;
    return 1;
  }

  M2(0, 0) = 1.0;
  M2(1, 1) = 1.0;
  M2(2, 2) = 1.0;
  if (std::abs(M2.lnDeterminant()) > TOL) {
    std::cerr << "ln determinant failed" << std::endl;
    return 1;
  }

  if (M2.rank(TOL, TOL) != 3) {
    std::cerr << "rank failed" << std::endl;
    return 1;
  }

  QUESO::uqVectorSpaceClass<QUESO::uqGslVectorClass, QUESO::uqGslMatrixClass> space(*env, "", 2, NULL);

  QUESO::uqGslVectorClass v2(space.zeroVector());
  QUESO::uqGslMatrixClass M3(v2, 0.0);

  fill2By2Matrix(M3);
  M3.filterSmallValues(2.5);
  if (std::abs(M3(0, 1) - 3.0) > TOL ||
      std::abs(M3(0, 0)) > TOL ||
      std::abs(M3(1, 0)) > TOL ||
      std::abs(M3(1, 1)) > TOL) {
    std::cerr << "filter small values failed" << std::endl;
    return 1;
  }

  fill2By2Matrix(M3);
  M3.filterLargeValues(2.5);
  if (std::abs(M3(0, 1)) > TOL ||
      std::abs(M3(0, 0) - 2.0) > TOL ||
      std::abs(M3(1, 0) - 2.0) > TOL ||
      std::abs(M3(1, 1) - 2.0) > TOL) {
    std::cerr << "filter large values failed" << std::endl;
    return 1;
  }

  QUESO::uqGslMatrixClass M4(v2, 0.0);
  fill2By2Matrix(M3);
  M4 = M3.inverse();
  if (std::abs(M4(0, 0) + 1.0) > TOL ||
      std::abs(M4(0, 1) - 1.5) > TOL ||
      std::abs(M4(1, 0) - 1.0) > TOL ||
      std::abs(M4(1, 1) + 1.0) > TOL) {
    std::cerr << "inverse failed" << std::endl;
    return 1;
  }

  M4 = M3;
  QUESO::uqGslMatrixClass I(M3.invertMultiply(M4));
  if (!matrixIsDiag(I, 1.0)) {
    std::cerr << "invert multiply failed" << std::endl;
    return 1;
  }

  fill2By2Matrix(M3);
  M3(1,0) = 3.0;
  M3.eigen(v2, &M4);
  double v00 = -1.0 / std::sqrt(2.0);
  double v01 = 1.0 / std::sqrt(2.0);
  double v10 = 1.0 / std::sqrt(2.0);
  double v11 = 1.0 / std::sqrt(2.0);
  double l1 = -1.0;
  double l2 = 5.0;
  if (std::abs(v2[0] - l1) > TOL ||
      std::abs(v2[1] - l2) > TOL ||
      std::abs(M4(0, 0) - v00) > TOL ||
      std::abs(M4(0, 1) - v01) > TOL ||
      std::abs(M4(1, 0) - v10) > TOL ||
      std::abs(M4(1, 1) - v11) > TOL) {
    std::cerr << "eigen failed" << std::endl;
    return 1;
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
