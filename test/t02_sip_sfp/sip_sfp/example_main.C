/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 *
 * This is an example of how to define and run both types of statistical problem, inverse and forward,
 * using QUESO classes and algorithms. The code itself is in the
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
  uqFullEnvironmentClass* env =
    new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"");

  // Compute
  compute(*env);

  // Finalize environment
  delete env;
  MPI_Finalize();

  return 0;
}
