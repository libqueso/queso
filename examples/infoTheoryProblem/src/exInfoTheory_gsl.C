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
#include <queso/asserts.h>
#include <queso/GslMatrix.h>

int main(int argc, char* argv[])
{
#ifdef QUESO_HAS_ANN
  queso_require_equal_to_msg(argc, 2, "Please specify an input file on the command line");

  // Initialize environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(MPI_COMM_WORLD, argv[1], "", NULL);
#else
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(argv[1], "", NULL);
#endif

  //************************************************
  // Call application
  //************************************************
  uqAppl<QUESO::GslVector, // type for parameter vectors
         QUESO::GslMatrix, // type for parameter matrices
         QUESO::GslVector, // type for qoi vectors
         QUESO::GslMatrix  // type for qoi matrices
        >(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;

#else // QUESO_HAS_ANN
  queso_error_msg("ANN has not been enabled. Please recompile QUESO with --enable-ann=yes");
#endif // QUESO_HAS_ANN
}
