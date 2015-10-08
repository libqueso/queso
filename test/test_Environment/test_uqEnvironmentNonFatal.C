/* #include <set> */
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/Defines.h>

using namespace std;

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = "debug_output";
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_platformName = "my_platform";
  options.m_displayVerbosity = 20;

  /* std::set<unsigned int> subDisplayAllowed; */
  /* subDisplayAllowed.insert(0); */

  /* options.m_subDisplayAllowedSet = subDisplayAllowed; */
  /* options.m_subDisplayAllowAll = 0; */

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment("", "", &options);
#endif

  if (!env->fullEnvIsReady()) {
    std::cerr << "Full env ready test failed" << std::endl;
    return 1;
  }

  if (env->displayVerbosity() > 10) {
    std::cout << "have high enough display verbosity" << std::endl;
  }
  else {
    std::cout << "not high enough" << std::endl;
  }

  if (env->worldRank() != 0) {
    std::cerr << "World rank test failed" << std::endl;
    return 1;
  }

  if (env->subDisplayFileName() != "debug_output") {
    std::cerr << "subDisplayFileName test failed" << std::endl;
    return 1;
  }

  env->setOptionsInputFileAccessState(false);
  if (env->optionsInputFileName() != "") {
    std::cerr << "Input file acces state test failed" << std::endl;
    return 1;
  }
  env->setOptionsInputFileAccessState(true);

  if (env->checkingLevel() != options.m_checkingLevel) {
    std::cerr << "Checking level test failed" << std::endl;
    return 1;
  }

  if (env->seed() != options.m_seed) {
    std::cerr << "Seed test failed" << std::endl;
    return 1;
  }

#ifdef QUESO_USES_NEW_RNG_CLASS
  env->resetSeed(2);
  if (env->seed() != 2) {
    std::cerr << "Second seed test failed" << std::endl;
    return 1;
  }

  env->resetSeed(-2);
  if (env->seed() != (2 + env->worldRank())) {
    std::cerr << "Third seed test failed" << std::endl;
    return 1;
  }
#else
  env->resetGslSeed(2);
  if (gsl_rng_default_seed != 2) {
    std::cerr << "Second seed test failed" << std::endl;
    return 1;
  }

  env->resetGslSeed(-2);
  if (gsl_rng_default_seed != (2 + env->worldRank())) {
    std::cerr << "Third seed test failed" << std::endl;
    return 1;
  }
#endif
  if (env->platformName() != options.m_platformName) {
    std::cerr << "Platform name test failed" << std::endl;
    return 1;
  }

  env->resetIdentifyingString("my_identifying_string");
  if (env->identifyingString() != "my_identifying_string") {
    std::cerr << "Identifying string test failed" << std::endl;
    return 1;
  }

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
