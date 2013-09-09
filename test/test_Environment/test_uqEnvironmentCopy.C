#include <uqEnvironment.h>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  QUESO::FullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new QUESO::FullEnvironmentClass(MPI_COMM_WORLD, "copy_env.inp", "", NULL);
#else
    new QUESO::FullEnvironmentClass(0, "copy_env.inp", "", NULL);
#endif

  QUESO::FullEnvironmentClass another_env(*env);

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
