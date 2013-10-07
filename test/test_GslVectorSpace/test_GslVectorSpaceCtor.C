#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;

  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "",
      "", &options);

  try {
    QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> vec_space;
  }
  catch (...) {
    std::cerr << "Caught QUESO exception" << std::endl;
    return 0;
  }
  std::cerr << "Did not catch QUESO exception" << std::endl;

  delete env;
  MPI_Finalize();
  return 1;
}
