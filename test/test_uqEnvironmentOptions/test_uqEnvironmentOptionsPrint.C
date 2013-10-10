/* #include <set> */
#include <queso/Environment.h>
#include <queso/Defines.h>
#include <queso/EnvironmentOptions.h>

#include <mpi.h>

using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

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

  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", &options);

  QUESO::EnvironmentOptions env_options(*env, "", options);
  std::cout << env_options;

  delete env;
  MPI_Finalize();
  /*
   * This code should never get here. If it does, the bash script that wraps
   * around it negates the return value, making this a failure
   */
  return 0;
}
