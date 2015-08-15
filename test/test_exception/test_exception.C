//
// regression that tests exception handling
//
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>
#include <queso/VectorRV.h>

// queso error handling
#include <queso/asserts.h>

using namespace std;

int main(int argc, char **argv)
{
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  // hit with exception
  try
    {
      queso_error_msg("Executing error");
    }
  catch(...)
    {
      printf("Caught QUESO exception!\n");
#ifdef QUESO_HAS_MPI
      MPI_Finalize();
#endif
      return 0;
    }

  return 1;
}
