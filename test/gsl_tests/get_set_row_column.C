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

#include <uqEnvironment.h>
#include <uqGslMatrix.h>
#include <uqVectorRV.h>

int actualChecking(const uqFullEnvironmentClass* env);

int main(int argc, char* argv[]) 
{
  int return_flag = 0;

   // Initialize environment
  MPI_Init(&argc,&argv);
  uqFullEnvironmentClass* env =
    new uqFullEnvironmentClass(MPI_COMM_WORLD,"input","");

  return_flag = actualChecking(env);

  delete env;
  MPI_Finalize();

  // Fin.
  return return_flag;
}

/* Separated this out into a function because we want
   the destructor for paramSpace to be called before
   we delete env in main. */
int actualChecking(const uqFullEnvironmentClass* env)
{
  int return_flag = 0;

  // Instantiate the parameter space
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace( (*env), "param_", 2, NULL);

  // Set matrix values
  uqGslMatrixClass* Matrix = paramSpace.newMatrix();
  (*Matrix)(0,0) = 4.; (*Matrix)(0,1) = 3.;
  (*Matrix)(1,0) = 5.; (*Matrix)(1,1) = 7.;

  // Conduct tests
  uqGslVectorClass row(Matrix->getRow(0));
  if( row[0] != (*Matrix)(0,0) || row[1] != (*Matrix)(0,1) )
    {
      return_flag = 1;
    }

  row = Matrix->getRow(1);
  if( row[0] != (*Matrix)(1,0) || row[1] != (*Matrix)(1,1) )
    {
      return_flag = 1;
    }

  uqGslVectorClass column(Matrix->getColumn(0));
  if( column[0] != (*Matrix)(0,0) || column[1] != (*Matrix)(1,0) )
    {
      return_flag = 1;
    }

  column = Matrix->getColumn(1);
  if( column[0] != (*Matrix)(0,1) || column[1] != (*Matrix)(1,1) )
    {
      return_flag = 1;
    }

  row[0] = 9.;
  row[1] = 15.;

  Matrix->setRow( 0, row );
  row = Matrix->getRow(0);
  if( row[0] != (*Matrix)(0,0) || row[1] != (*Matrix)(0,1) )
    {
      return_flag = 1;
    }
  
  Matrix->setColumn( 1, row );
  column = Matrix->getColumn(1);
  if( column[0] != (*Matrix)(0,1) || column[1] != (*Matrix)(1,1) )
    {
      return_flag = 1;
    }


  // Deallocate pointers we created
  delete Matrix;

  return return_flag;
}
