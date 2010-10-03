//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
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
