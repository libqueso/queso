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
 * $Id: exInfoTheory_gsl.C 14596 2011-02-03 15:10:42Z gabriel $
 *
 * Brief description of this file: 
 * 
 * This is an example of how to use the information theoretic measures such as
 * entropy, mutual information and Kullback-Leibler divergenece available in 
 * QUESO. The code itself is in the templated routine 'uqAppl(*env)'. This 
 * routine is called right after the initialization of the MPI environment and 
 * of the QUESO environment and is available in file 'exInfoTheory_appl.h'
 *
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <exInfoTheory_appl.h>
#include <uqGslMatrix.h>
#include <config_queso.h>

int main(int argc, char* argv[])
{

#ifdef HAVE_ANN

  //************************************************
  // Initialize environment
  //************************************************
  MPI_Init(&argc,&argv);

  UQ_FATAL_TEST_MACRO(argc != 2,
                      UQ_UNAVAILABLE_RANK,
                      "main()",
                      "input file must be specified in command line as argv[1], just after executable argv[0]");
  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  //************************************************
  // Call application
  //************************************************
  uqAppl<uqGslVectorClass, // type for parameter vectors
         uqGslMatrixClass, // type for parameter matrices
         uqGslVectorClass, // type for qoi vectors
         uqGslMatrixClass  // type for qoi matrices
        >(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
  MPI_Finalize();

  return 0;


#else // HAVE_ANN

  std::cout << "Error: the ANN has not been enabled. Please recompile with --enable-infotheory" << std::endl;
  exit(1);

#endif // HAVE_ANN

}
