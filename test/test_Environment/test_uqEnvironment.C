#include <stdio.h>
#include <string.h>
#include <cmath>
#include <uqEnvironment.h>
#include <uqGslVector.h>
#include <uqGslMatrix.h>
#include <uqVectorSpace.h>
#include <uqVectorSubset.h>
#include <uqVectorRV.h>

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  QUESO::uqFullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new QUESO::uqFullEnvironmentClass(MPI_COMM_WORLD, "", "", NULL);
#else
    new QUESO::uqFullEnvironmentClass(0, "", "", NULL);
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
