//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id:$
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
