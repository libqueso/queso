/* uq/examples/queso/uqPyramid_gsl.C
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqTgaValidation.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  MPI_Init(&argc,&argv);
  uqEnvironmentClass* env = new uqEnvironmentClass(argc,argv);

  UQ_FATAL_TEST_MACRO(env->isThereInputFile() == false,
                      env->rank(),
                      "main()",
                      "input file must be specified in command line, after the '-i' option");

  //************************************************
  // Check for the invalidation of model(s): either 'invalid' or 'not invalid'
  //************************************************
  uqTgaValidationClass<uqGslVectorClass, // type for parameter vectors
                       uqGslMatrixClass, // type for parameter matrices
                       uqGslVectorClass, // type for qoi vectors
                       uqGslMatrixClass  // type for qoi matrices
                      > tgaValidation(*env,"tga_");
  tgaValidation.run();

#if 0
  uqTurValidationClass<uqGslVectorClass, // type for parameter vectors
                       uqGslMatrixClass, // type for parameter matrices
                       uqGslVectorClass, // type for qoi vectors
                       uqGslMatrixClass  // type for qoi matrices
                      > turbValidation(*env);
  turbValidation.run();

  uqTgaTurValidationClass<uqGslVectorClass, // type for parameter vectors
                          uqGslMatrixClass, // type for parameter matrices
                          uqGslVectorClass, // type for qoi vectors
                          uqGslMatrixClass  // type for qoi matrices
                         > tgaTurValidation(tgaValidation,turValidation);
  tgaTurbValidation.run();
#endif

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
  MPI_Finalize();

  return 0;
}
