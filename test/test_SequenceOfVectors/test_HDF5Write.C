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
#ifndef QUESO_HAS_HDF5
  // If we don't have HDF5, skip this test
  return 77;
#else
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = ".";
  options.m_subDisplayAllowAll = 0;
  options.m_subDisplayAllowedSet.insert(0);
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 0;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment env("", "", &options);
#endif

  // Create a vector space
  QUESO::VectorSpace<> vec_space(env, "", 2, NULL);

  // Create some things to put in the sequence
  QUESO::GslVector v1(vec_space.zeroVector());
  v1[0] = 1.0;
  v1[1] = 2.0;

  QUESO::GslVector v2(vec_space.zeroVector());
  v2[0] = 3.0;
  v2[1] = 4.0;

  // Create a sequence of vectors
  QUESO::SequenceOfVectors<> vec_seq(vec_space, 2, "");

  vec_seq.setPositionValues(0, v1);
  vec_seq.setPositionValues(1, v2);

  // Now write
  if (env.fullRank() == 0) {
    std::set<unsigned int> allowedIds;
    allowedIds.insert(0);
    vec_seq.subWriteContents(0, 2, "output_test_hdf5", "h5", allowedIds);
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
#endif  // QUESO_HAS_HDF
}
