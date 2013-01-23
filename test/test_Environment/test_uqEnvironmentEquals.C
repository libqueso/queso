#include <uqEnvironment.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  uqFullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new uqFullEnvironmentClass(MPI_COMM_WORLD, "copy_env.inp", "", NULL);
#else
    new uqFullEnvironmentClass(0, "copy_env.inp", "", NULL);
#endif

#ifdef QUESO_HAS_MPI
  uqFullEnvironmentClass another_env(MPI_COMM_WORLD, "copy_env.inp", "", NULL);
#else
  uqFullEnvironmentClass another_env(0, "copy_env.inp", "", NULL);
#endif

  another_env = *env;

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
