#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/IntersectionSubset.h>

int main(int argc, char **argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = "outputData/testIntersectionSubsetContains";
  options.m_subDisplayAllowAll = 0;
  options.m_subDisplayAllowedSet.insert(0);
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 55;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment(MPI_COMM_WORLD, "",
            "", &options);
#else
  QUESO::FullEnvironment *env = new QUESO::FullEnvironment("",
            "", &options);
#endif

  std::vector<std::string> names(1);
  names[0] = "my_name";

  // Create a vector space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> vec_space(*env,
      "vec_prefix", 1, &names);

  // Create two vector sets
  QUESO::GslVector min1(vec_space.zeroVector());
  min1[0] = 0.0;
  QUESO::GslVector max1(vec_space.zeroVector());
  max1[0] = 1.0;
  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> set1("set1", vec_space,
      min1, max1);

  // Now for the second one
  QUESO::GslVector min2(vec_space.zeroVector());
  min2[0] = 0.5;
  QUESO::GslVector max2(vec_space.zeroVector());
  max2[0] = 1.5;
  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> set2("set1", vec_space,
      min2, max2);

  // Create their intersection
  QUESO::IntersectionSubset<QUESO::GslVector, QUESO::GslMatrix> intersection(
      "intersection", vec_space, 1.0, set1, set2);

  // Test the containment
  bool does_contain;
  QUESO::GslVector test_vec(vec_space.zeroVector());

  // Should be true
  test_vec[0] = 0.75;
  does_contain = intersection.contains(test_vec);
  if (!does_contain) {
    std::cerr << "First IntersectionSubset contains test failed" << std::endl;
    return 1;
  }

  // Should be false
  test_vec[0] = 2.0;
  does_contain = intersection.contains(test_vec);
  if (does_contain) {
    std::cerr << "Second contains test failed" << std::endl;
    return 1;
  }

  // Print out some info
  intersection.print(std::cout);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
