//
// regression that tests exception handling 
//
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <uqEnvironment.h>
#include <uqGslVector.h>
#include <uqGslMatrix.h>
#include <uqVectorSpace.h>
#include <uqVectorSubset.h>
#include <uqVectorRV.h>

//
#include <uq_asserts.h>

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  uqFullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new uqFullEnvironmentClass(MPI_COMM_WORLD, "", "", NULL);
#else
    new uqFullEnvironmentClass(0, "", "", NULL);
#endif

  delete env;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  //QUESO_THROW('error');
  queso_error();

  return 0;
}
