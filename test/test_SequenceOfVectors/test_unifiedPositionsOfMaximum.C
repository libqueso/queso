#include <string>

#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/ScalarSequence.h>
#include <queso/SequenceOfVectors.h>

int main(int argc, char **argv) {
#ifndef QUESO_HAS_MPI
  return 77;
#else
  MPI_Init(&argc, &argv);

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 2;
  options.m_subDisplayFileName = "outputData/test_SequenceOfVectorsMax";
  options.m_subDisplayAllowAll = 1;
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 55;

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &options);

  // Create a vector space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env, "", 1,
      NULL);

  QUESO::SequenceOfVectors<QUESO::GslVector, QUESO::GslMatrix> pretendChain(
      paramSpace, 3, "pretendChain");

  QUESO::GslVector v(paramSpace.zeroVector());
  v[0] = 0.0;

  // Create a scalar sequence on each processor
  std::string name = "name";
  QUESO::ScalarSequence<double> scalarSequence(env, 3, name);
  pretendChain.setPositionValues(0, v);
  pretendChain.setPositionValues(1, v);
  pretendChain.setPositionValues(2, v);

  if (env.inter0Rank() == 0) {
    scalarSequence[0] = 10.0;
    scalarSequence[1] = 11.0;
    scalarSequence[2] = 12.0;
  }
  else if (env.inter0Rank() == 1) {
    scalarSequence[0] = 0.0;
    scalarSequence[1] = 1.0;
    scalarSequence[2] = 2.0;
  }
  else {
    queso_error_msg("Test should not get here!");
  }

  // Create a sequence of vectors
  QUESO::SequenceOfVectors<QUESO::GslVector, QUESO::GslMatrix> maxs(paramSpace,
      0, "name2");

  pretendChain.unifiedPositionsOfMaximum(scalarSequence, maxs);

  // Should not fail
  QUESO::GslVector tmpVec(paramSpace.zeroVector());

  // The loop should only actually do anything on process zero
  for (unsigned int i = 0; i < maxs.subSequenceSize(); i++) {
    maxs.getPositionValues(0, tmpVec);
  }

  MPI_Finalize();

  return 0;
#endif
}
