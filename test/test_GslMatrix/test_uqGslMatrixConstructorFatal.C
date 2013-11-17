#include <queso/GslMatrix.h>
#include <mpi.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  try {
    QUESO::GslMatrix matrix;
  }
  catch (...) {
    std::cerr << "Caught QUESO exception!" << std::endl;
    return 0;
  }
  std::cerr << "Did not catch QUESO exception!" << std::endl;

  MPI_Finalize();
  return 1;
}
