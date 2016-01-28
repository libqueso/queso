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

  d(0, 0) = "1";
  const std::string & one = d(0, 0);
  if (!(one == "1")) {
    std::cerr << "operator= failed" << std::endl;
  }

  if (d.RowSize() != 1) {
    std::cerr << "row size test failed" << std::endl;
    return 1;
  }

  std::cerr << "d is: " << std::endl;
  std::cerr << d << std::endl;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
