#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = "outputData/debug_output";
  options.m_subDisplayAllowAll = 0;
  options.m_subDisplayAllowedSet.insert(0);
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 20;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "",
      "", &options);
#else
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment("",
      "", &options);
#endif

  std::vector<std::string> names(1);
  names[0] = "my_name";
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> vec_space(*env,
      "vec_prefix", 1, &names);

  // 1x1 diagonal matrix containing the value 3.0 on the diagonal
  QUESO::GslMatrix * diag_matrix = vec_space.newDiagMatrix(3.0);

  if ((*diag_matrix)(0,0) != 3.0) {
    std::cerr << "newDiagMatrix test failed" << std::endl;
    return 1;
  }

  if (vec_space.globalIdOfFirstComponent() != 0) {
    std::cerr << "failed globalIdOfFirstComponent test" << std::endl;
    return 1;
  }

  // Create a vector
  QUESO::GslVector v1(vec_space.zeroVector());
  v1.cwSet(2.0);

  // Now create a new vector
  QUESO::GslVector * v2 = vec_space.newVector(v1);

  if ((*v2)[0] != v1[0]) {
    std::cerr << "newVector test failed" << std::endl;
    return 1;
  }

  const QUESO::DistArray<std::string> *names_array = vec_space.componentsNamesArray();

  std::cout << "Test print: componentsNamesArray: "
            << *names_array
            << std::endl;

  std::cout << "Test print: localComponentName: "
            << vec_space.localComponentName(0)
            << std::endl;

  std::cout << "Test print: printComponentsNames horizontally: ";
  vec_space.printComponentsNames(std::cout, true);
  std::cout << std::endl;

  std::cout << "Test print: printComponentsNames vertically: ";
  vec_space.printComponentsNames(std::cout, false);
  std::cout << std::endl;

  std::cout << "Test print: print: ";
  vec_space.print(std::cout);
  std::cout << std::endl;

  delete diag_matrix;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
