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
 * Usage: ./example_gsl example.inp  <(int) numModes>
 *-------------------------------------------------------------------------- */

#include <example_compute.h>

int main(int argc, char* argv[])
{
  // Initialize environment
  MPI_Init(&argc,&argv);
 
  UQ_FATAL_TEST_MACRO(argc != 3,
                      QUESO::UQ_UNAVAILABLE_RANK,
                      "main()",
                      "after executable argv[0], input file must be specified in command line as argv[1], then numModes (1 or 2) must be specified as argv[2]");
                      
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
    
  // Compute
  unsigned int numModes = (unsigned int) atoi(argv[2]);
  
  compute(*env,numModes);

  // Finalize environment
  delete env;
  MPI_Finalize();

  std::cout << std::endl << "FIM!" << std::endl << std::endl; 
  
  return 0;
}
