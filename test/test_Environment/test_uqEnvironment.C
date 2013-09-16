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
  // This should fail because we didn't call MPI_Init
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "",
      "", NULL);

  delete env;
  MPI_Finalize();

  /*
   * This code should never get here. If it does, the bash script that wraps
   * around it negates the return value, making this a failure
   */
  return 0;
}
