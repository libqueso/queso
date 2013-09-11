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
  MPI_Init(&argc, &argv);

  // hit with exception
  try
    {
      printf("Executing error\n");
      queso_error();
    }
  catch(...)
    {
      printf("Caught QUESO exception!\n");
      MPI_Finalize();
      return 0;
    }

  return 1;
}
