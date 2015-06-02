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
  MPI_Init(&argc, &argv);

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = "outputData/test_SequenceOfVectorsErase";
  options.m_subDisplayAllowAll = 0;
  options.m_subDisplayAllowedSet.insert(0);
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 55;

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &options);

  // Create a vector space
  std::vector<std::string> names(1);
  names[0] = "my_name";
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> vec_space(env,
      "vec_prefix", 1, &names);

  // Create some things to put in the sequence
  QUESO::GslVector v1(vec_space.zeroVector());
  QUESO::GslVector v2(vec_space.zeroVector());
  v1[0] = 0.0;
  v2[0] = 1.0;

  // Create a sequence of vectors
  QUESO::SequenceOfVectors<QUESO::GslVector, QUESO::GslMatrix> vec_seq(
      vec_space, 2, "vec_seq");

  vec_seq.setPositionValues(0, v1);
  vec_seq.setPositionValues(1, v2);

  // Now erase
  vec_seq.erasePositions(0, 1);

  // Should generate an error
  try {
    vec_seq.getPositionValues(0, v1);
  }
  catch (...) {
    return 0;
  }

  MPI_Finalize();
  return 1;
}
