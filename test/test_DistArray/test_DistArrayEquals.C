#include <mpi.h>
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/MpiComm.h>
#include <queso/Map.h>
#include <queso/DistArray.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;

  QUESO::FullEnvironment *env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", &options);

  const QUESO::MpiComm & comm = env->fullComm();

  QUESO::Map map(1, 1, comm);

  QUESO::DistArray<std::string> d(map, 1);
  QUESO::DistArray<std::string> e(map, 1);

  try {
    e = d;
  }
  catch (...) {
    delete env;
    return 0;
  }

  MPI_Finalize();
  return 1;
}
