#include <queso/Environment.h>
#include <mpi.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "test_Environment/copy_env", "",
      NULL);

  try {
    QUESO::FullEnvironment another_env(env);
  }
  catch (...) {
    std::cerr << "Caught QUESO exception" << std::endl;
    return 0;
  }
  std::cerr << "Did not catch QUESO exception" << std::endl;

  MPI_Finalize();
  return 1;
}
