#include <vector>
#include <set>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>

#define TOL 1e-10

int checkLinearSpacing(const QUESO::GslVector &v, double d1, double d2) {
  unsigned int i;
  double alpha, elt;

  for (i = 0; i < v.sizeLocal(); i++) {
    alpha = (double) i / (v.sizeLocal() - 1);
    elt = (1.0 - alpha) * d1 + alpha * d2;
    if (std::abs(elt - v[i]) >= TOL) {
      return 1;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  unsigned int i;
  double d1 = 0.0;
  double d2 = 1.0;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment("", "", &options);
#endif

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> *param_space =
    new QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>(*env, "param_", 3, NULL);

  QUESO::GslVector v1(param_space->zeroVector());
  v1.cwSet(2.0);

  QUESO::GslVector *v2 = new QUESO::GslVector(*env, d1, d2, v1.map());
  if (checkLinearSpacing(*v2, d1, d2) == 1) {
    std::cerr << "Linear spacing test 1 failed" << std::endl;
  }
  delete v2;

  v2 = new QUESO::GslVector(v1, d1, d2);
  if (checkLinearSpacing(*v2, d1, d2) == 1) {
    std::cerr << "Linear spacing test 2 failed" << std::endl;
  }

  QUESO::GslVector v3(*v2);
  *v2 /= v1;

  for (i = 0; i < v2->sizeLocal(); i++) {
    if (std::abs((*v2)[i] - ((double) v3[i] / v1[i])) > TOL) {
      std::cerr << "Self-divide test failed" << std::endl;
      return 1;
    }
  }

  // At this point, v3 is the vector (0.0, 0.5, 1.0)
  if (std::abs(v3.norm1() - 1.5) > TOL) {
    std::cerr << "1-norm test failed. You should probably call it a day." << std::endl;
    return 1;
  }

  if (std::abs(v3.normInf() - 1.0) > TOL) {
    std::cerr << "inf-norm test failed" << std::endl;
    return 1;
  }

  // Testing concatenate so we need a bigger param space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> *big_param_space =
    new QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>(*env, "", 6, NULL);
  QUESO::GslVector v4(big_param_space->zeroVector());
  v4.cwSetConcatenated(v1, v3);

  for (i = 0; i < v1.sizeLocal(); i++) {
    if (std::abs(v4[i] - v1[i]) > TOL) {
      std::cerr << "Concatenate (first half) test failed" << std::endl;
    }
  }
  for (i = 0; i < v3.sizeLocal(); i++) {
    if (std::abs(v4[i+3] - v3[i]) > TOL) {
      std::cerr << "Concatenate (second half) test failed" << std::endl;
      return 1;
    }
  }

  std::vector<const QUESO::GslVector*> vecs;
  vecs.push_back(&v1);
  vecs.push_back(&v3);
  QUESO::GslVector v5(big_param_space->zeroVector());
  v5.cwSetConcatenated(vecs);
  for (i = 0; i < v5.sizeLocal(); i++) {
    if (std::abs(v5[i] - v4[i]) > TOL) {
      std::cerr << "Concatenate test 2 failed" << std::endl;
      return 1;
    }
  }

  v5.cwSet(0.0);
  v5.cwSet(3, v3);
  for (i = 0; i < v5.sizeLocal(); i++) {
    if (i < 3) {
      if (std::abs(v5[i] - 0.0) > TOL) {
        std::cerr << "cwSet with initial position failed (first half)" << std::endl;
        return 1;
      }
    }
    else if (i >= 3) {
      if (std::abs(v5[i] - v3[i-3]) > TOL) {
        std::cerr << "cwSet with intial position failed (second half)" << std::endl;
        return 1;
      }
    }
  }


  v5.cwExtract(0, v3);
  for (i = 0; i < v3.sizeLocal(); i++) {
    if (std::abs(v3[i]) > TOL) {
      std::cerr << "cwExtract test failed" << std::endl;
      return 1;
    }
  }

  v3[0] = 0.25;
  v3[1] = 0.50;
  v3[2] = 1.00;
  v3.cwInvert();
  if (std::abs(v3[0] - 4.0) > TOL ||
      std::abs(v3[1] - 2.0) > TOL ||
      std::abs(v3[2] - 1.0 > TOL)) {
    std::cerr << "cwInvert test failed" << std::endl;
    return 1;
  }

  v3.matlabDiff(0, 0.0, v1);
  if (std::abs(v1[0] + 2.0) > TOL ||
      std::abs(v1[1] + 1.0) > TOL ||
      std::abs(v1[2]) > TOL) {
    std::cerr << "matlabDiff test failed" << std::endl;
    return 1;
  }

  v3.matlabDiff(1, 0.0, v1);
  if (std::abs(v1[0]) > TOL ||
      std::abs(v1[1] + 2.0) > TOL ||
      std::abs(v1[2] + 1.0) > TOL) {
    std::cerr << "matlabDiff 2 test failed" << std::endl;
    return 1;
  }

  v1.sort();
  if (std::abs(v1[0] + 2.0) > TOL ||
      std::abs(v1[1] + 1.0) > TOL ||
      std::abs(v1[2]) > TOL) {
    std::cerr << "sort test failed" << std::endl;
    return 1;
  }

  // There is a shell
  std::set<unsigned int> ids;
  ids.insert(0);
  v1.subWriteContents("gslvector", "gslvector_out", UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, ids);
  v3.subReadContents("gslvector_out_sub0", UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, ids);
  std::cerr << v3[0] << "  " << v3[1] << "  " << v3[2] << std::endl;
  for (i = 0; i < v3.sizeLocal(); i++) {
    if (std::abs(v1[i] - v3[i]) > TOL) {
      std::cerr << "read/write test failed" << std::endl;
      return 1;
    }
  }

  v3[2] = -2.0;
  if (!v3.atLeastOneComponentSmallerOrEqualThan(v1)) {
    std::cerr << "smaller test failed" << std::endl;
    return 1;
  }

  v3[2] = 0.0;
  if (!v3.atLeastOneComponentBiggerOrEqualThan(v1)) {
    std::cerr << "bigger test failed" << std::endl;
    return 1;
  }

  if (std::abs(v3.getMaxValue()) > TOL) {
    std::cerr << "Max value test failed" << std::endl;
    return 1;
  }

  if (v3.getMinValueIndex() != 0) {
    std::cerr << "Min value index test failed" << std::endl;
    return 1;
  }

  v3[0] = 1.0;
  v3[1] = 2.0;
  v3[2] = 3.0;
  v1 = 1.0 / v3;
  v3.cwInvert();
  for (i = 0; i < v3.sizeLocal(); i++) {
    if (std::abs(v1[i] - v3[i]) > TOL) {
      std::cerr << "Operator / test failed" << std::endl;
      return 1;
    }
  }

  if (!(v1 == v3)) {
    std::cerr << "operator == test failed" << std::endl;
    return 1;
  }

  QUESO::GslVector ones(v1, 1.0, 1.0);
  if (!(ones == (v1 / v3))) {
    std::cerr << "division test failed" << std::endl;
    return 1;
  }
  delete big_param_space;
  delete param_space;
  delete env;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
