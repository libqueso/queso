#include <stdio.h>
#include <string.h>
#include <cmath>
#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>
#include <queso/VectorRV.h>

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  QUESO::FullEnvironment *env =
#ifdef QUESO_HAS_MPI
    new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", NULL);
#else
    new QUESO::FullEnvironment(0, "", "", NULL);
#endif

  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  /*
   * This code should never get here. If it does, the bash script that wraps
   * around it negates the return value, making this a failure
   */
  return 0;
}
