//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <multiple_rhs_matrix_solve.h>
#include <queso/GslMatrix.h>
#include <queso/VectorRV.h>

int main(int argc, char* argv[])
{
  int return_flag = 0;

   // Initialize environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(MPI_COMM_WORLD,"input","",NULL);
#else
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment("input","",NULL);
#endif

  return_flag = actualChecking(env);

  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  // Fin.
  return return_flag;
}

/* Separated this out into a function because we want
   the destructor for paramSpace to be called before
   we delete env in main. */
int actualChecking(const QUESO::FullEnvironment* env)
{
  int return_flag = 0;

  const double fp_tol = 1.0e-14;

  // Instantiate the parameter space
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    paramSpace( (*env), "param_", 2, NULL);

  // Populate Matrix.
  QUESO::GslMatrix* Matrix = paramSpace.newMatrix();
  (*Matrix)(0,0) = 4.; (*Matrix)(0,1) = 3.;
  (*Matrix)(1,0) = 5.; (*Matrix)(1,1) = 7.;

  // Conduct test.
  QUESO::GslMatrix Result( (*Matrix) );

  Matrix->invertMultiply( (*Matrix), Result );

  if( fabs(Result(0,0) - 1.0) > fp_tol ||
      fabs(Result(1,1) - 1.0) > fp_tol ||
      fabs(Result(1,0))       > fp_tol ||
      fabs(Result(0,1))       > fp_tol    )
    {
      return_flag = 1;
    }

  // Deallocate pointers we created.
  delete Matrix;

  return return_flag;
}
