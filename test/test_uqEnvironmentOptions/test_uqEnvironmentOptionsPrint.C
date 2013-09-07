/* #include <set> */
#include <uqEnvironment.h>
#include <uqDefines.h>
#include <uqEnvironmentOptions.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = "";
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_platformName = "my_platform";
  options.m_displayVerbosity = 0;
  options.m_identifyingString = "my_identifying_string";

  std::set<unsigned int> subDisplayAllowed;
  subDisplayAllowed.insert(1);

  options.m_subDisplayAllowedSet = subDisplayAllowed;
  options.m_subDisplayAllowAll = 0;

  QUESO::FullEnvironment *env =
#ifdef QUESO_HAS_MPI
    new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", &options);
#else
    new QUESO::FullEnvironment(0, "", "", &options);
#endif

  QUESO::EnvironmentOptions env_options(*env, "", options);
  std::cout << env_options;

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
