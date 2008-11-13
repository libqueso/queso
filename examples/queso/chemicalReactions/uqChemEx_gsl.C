/* uq/examples/queso/chemicalReactions/uqChemEx_gsl.C
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

#include <uqChemEx.h>
#include <uqApplRoutines.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  MPI_Init(&argc,&argv);
  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(argc,argv);

  //************************************************
  // Call application
  //************************************************
  uqAppl<uqGslVectorClass, // type for parameter vectors
         uqGslMatrixClass, // type for parameter matrices
         uqGslVectorClass, // type for state vectors
         uqGslMatrixClass, // type for state matrices
         uqGslVectorClass, // type for likelihood vectors
         uqGslMatrixClass  // type for likelihood matrices
        >(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
  MPI_Finalize();

  return 0;
}

int
calib_StateDotRoutine_gsl( // Compute state dot = dConcentrations/dt
  double       t,
  const double currentState[],           // current concentrations
  double       stateDot[],               // dConcentrations/dt
  void*        infoForComputingStateDot) // concentration rates
{
  const uqGslVectorClass* k = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->concentrationRates;

  double k1 = (*k)[0];
  double k2 = (*k)[1];
  double k3 = (*k)[2]; 
  double A = currentState[0];
  double B = currentState[1];
  double C = currentState[2];
  double D = currentState[3];
  double E = currentState[4];
  E = E; // E is not used. Add this line just to avoid compiler warning;

  stateDot[0] = -k1*A*B - k2*A*C - k3*A*D;
  stateDot[1] = -k1*A*B;
  stateDot[2] =  k1*A*B - k2*A*C;
  stateDot[3] =           k2*A*C - k3*A*D;
  stateDot[4] =                    k3*A*D;

  return GSL_SUCCESS;
}
