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
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <exPhysics1_validation.h>
#include <uqGslMatrix.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  MPI_Init(&argc,&argv);

  UQ_FATAL_TEST_MACRO(argc != 2,
                      UQ_UNAVAILABLE_RANK,
                      "main()",
                      "input file must be specified in command line as argv[1], just after executable argv[0]");

  uqFullEnvironment* env = new uqFullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);

  //************************************************
  // Check for the invalidation of model(s): either 'invalid' or 'not invalid'
  //************************************************
  exPhysics1Validation<uqGslVector, // type for parameter vectors
                            uqGslMatrix, // type for parameter matrices
                            uqGslVector, // type for qoi vectors
                            uqGslMatrix  // type for qoi matrices
                           > physics1Validation(*env,"physics1_");
  physics1Validation.run();

#if 0
  uqPhysics2Validation<uqGslVector, // type for parameter vectors
                            uqGslMatrix, // type for parameter matrices
                            uqGslVector, // type for qoi vectors
                            uqGslMatrix  // type for qoi matrices
                           > physics2Validation(*env);
  physics2Validation.run();

  uqPhysics1And2Validation<uqGslVector, // type for parameter vectors
                                uqGslMatrix, // type for parameter matrices
                                uqGslVector, // type for qoi vectors
                                uqGslMatrix  // type for qoi matrices
                               > physics1And2Validation(physics1Validation,physics2Validation);
  physics1And2Validation.run();
#endif

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
  MPI_Finalize();
  return 0;
}
