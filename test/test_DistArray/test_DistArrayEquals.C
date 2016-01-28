#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/MpiComm.h>
#include <queso/Map.h>
#include <queso/DistArray.h>

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment *env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment *env =
    new QUESO::FullEnvironment("", "", &options);
#endif

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

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 1;
}
