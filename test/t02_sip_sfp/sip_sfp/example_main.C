//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------
 * routine 'compute(*env)'. This routine is called right after the initialization
 * of the MPI environment and of the QUESO environment and is available in
 * file 'example_compute.C'.
 *
 * It is easier to understand this example after undertanding two other simpler
 * examples:
 * (1) examples/statisticalInverseProblem1/src/exStatisticalInverseProblem1_gsl.C
 * (2) examples/statisticalForwardProblem/src/exStatisticalForwardProblem_gsl.C
 *
 * Example (1) focuses on the templated class uqStatisticalInverseProblemClass<P_V,P_M>.
 * Example (2) focuses on the templated class uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>.
 * 
 * This example uses both statistical problem templated classes, once each,
 * since it uses the solution of the inverse problem as input to the forward problem.
 * 
 * This example is explained in detail in the user's manual (pdf file)
 *
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <example_compute.h>

int main(int argc, char* argv[])
{
  // Initialize environment
  MPI_Init(&argc,&argv);

  UQ_FATAL_TEST_MACRO(argc != 2,
                      UQ_UNAVAILABLE_RANK,
                      "main()",
                      "input file must be specified in command line as argv[1], just after executable argv[0]");
  uqFullEnvironmentClass* env =
    new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // Compute
  compute(*env);

  // Finalize environment
  delete env;
  MPI_Finalize();

  return 0;
}
