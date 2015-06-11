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

  MPI_Finalize();
  return 0;
}
