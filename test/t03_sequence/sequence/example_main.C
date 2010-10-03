//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
// Brief description of this file: 
// 
// This is an example of how to define and use some basic classes and
// algorithms of QUESO, mainly 'uqGaussianVectorRVClass',
// 'uqVectorSequenceClass', realizations and statistical computations.
// The code itself is in the routine 'compute(*env)'. This routine is
// called right after the initialization of the MPI environment and of
// the QUESO environment and is available in file 'example_compute.C'.
//
//--------------------------------------------------------------------------


#include <example_compute.h>

int main(int argc, char* argv[])
{
  // Initialize environment
  MPI_Init(&argc,&argv);
  uqFullEnvironmentClass* env =
    new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // Compute
  compute(*env);

  // Finalize environment
  delete env;
  MPI_Finalize();

  return 0;
}
