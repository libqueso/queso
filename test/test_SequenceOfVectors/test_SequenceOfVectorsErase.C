#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/SequenceOfVectors.h>

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = "outputData/test_SequenceOfVectorsErase";
  options.m_subDisplayAllowAll = 0;
  options.m_subDisplayAllowedSet.insert(0);
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 55;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment env("", "", &options);
#endif

  // Create a vector space
  std::vector<std::string> names(2);
  names[0] = "my_name_0";
  names[1] = "my_name_1";
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> vec_space(env,
      "vec_prefix", 2, &names);

  // Create some things to put in the sequence
  QUESO::GslVector v(vec_space.zeroVector());

  // Create a sequence of vectors
  QUESO::SequenceOfVectors<QUESO::GslVector, QUESO::GslMatrix> vec_seq(
      vec_space, 13, "vec_seq");

  for (unsigned int i = 0; i < 13; i++) {
    v[0] = i;
    v[1] = i + 1;
    vec_seq.setPositionValues(i, v);
  }

  // Now erase
  vec_seq.erasePositions(2, 5);

  QUESO::GslVector expected(vec_space.zeroVector());
  QUESO::GslVector computed(vec_space.zeroVector());
  for (unsigned int i = 0; i < 8; i++) {
    if (i < 2) {
      expected[0] = i;
      expected[1] = i + 1;
    }
    else {
      expected[0] = i + 5;
      expected[1] = i + 6;
    }

    vec_seq.getPositionValues(i, computed);
    queso_require_equal_to(expected[0], computed[0]);
    queso_require_equal_to(expected[1], computed[1]);
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
