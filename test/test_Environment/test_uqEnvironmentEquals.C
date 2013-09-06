#include <uqEnvironment.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  QUESO::uqFullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new QUESO::uqFullEnvironmentClass(MPI_COMM_WORLD, "copy_env.inp", "", NULL);
#else
    new QUESO::uqFullEnvironmentClass(0, "copy_env.inp", "", NULL);
#endif

#ifdef QUESO_HAS_MPI
  QUESO::uqFullEnvironmentClass another_env(MPI_COMM_WORLD, "copy_env.inp", "", NULL);
#else
  QUESO::uqFullEnvironmentClass another_env(0, "copy_env.inp", "", NULL);
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
