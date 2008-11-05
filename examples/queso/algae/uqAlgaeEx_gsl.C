/* uq/examples/queso/algae/uqAlgaeEx_gsl.C
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

#include <uqAlgaeEx.h>
#include <uqApplRoutines.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  MPI_Init(&argc,&argv);
  uqEnvironmentClass* env = new uqEnvironmentClass(argc,argv);

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
  double muMax = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->muMax;
  double rhoA  = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->rhoA;
  double rhoZ  = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->rhoZ;
  double k     = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->k;
  double alpha = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->alpha;
  double th    = ((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->th;

  const std::vector<double>& evolutionOfQpV = *((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->evolutionOfQpV;
  const std::vector<double>& evolutionOfT   = *((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->evolutionOfT;
  const std::vector<double>& evolutionOfPin = *((calib_StateDotRoutine_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->evolutionOfPin;

  bool bRC = ((evolutionOfQpV.size() == evolutionOfPin.size()) &&
              (evolutionOfQpV.size() == evolutionOfT.size()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      UQ_UNAVAILABLE_RANK,
                      "calib_StateDotRoutine_gsl()",
                      "evolution arrays are not all of the same size");

  // E.g.: if t == 0.  then instantIndex = 0
  //       if t == 0.1 then instantIndex = 0
  //       if t == 1.  then instantIndex = 0
  //       if t == 1.1 then instantIndex = 1
  unsigned int instantIndex = (unsigned int) ceil(t);
  if (instantIndex > 0) instantIndex--;
  //std::cout << "t = " << t << ", instantIndex = " << instantIndex << std::endl;
  UQ_FATAL_TEST_MACRO(instantIndex >= evolutionOfQpV.size(),
                      UQ_UNAVAILABLE_RANK,
                      "calib_StateDotRoutine_gsl()",
                      "instantIndex is too big, that is, requested instant 't' is incompatible with available data");

  double qpv  = evolutionOfQpV[instantIndex];
  double temp = evolutionOfT  [instantIndex];
  double pin  = evolutionOfPin[instantIndex];

  double A = currentState[0];
  double Z = currentState[1];
  double P = currentState[2];

  double mu = muMax*pow(th,temp-20)*P/(k+P);

  stateDot[0] = (mu - rhoA - qpv - alpha*Z)*A;
  stateDot[1] = alpha*Z*A - rhoZ*Z;
  stateDot[2] = -qpv*(P-pin) + (rhoA-mu)*A + rhoZ*Z;

  //std::cout << "In calib_StateDotRoutine_gsl():"
  //          << "\n t     = " << t
  //          << "\n A     = " << A
  //          << "\n Z     = " << Z
  //          << "\n P     = " << P
  //          << "\n qpv   = " << qpv
  //          << "\n pin   = " << pin
  //          << "\n temp  = " << temp
  //          << "\n muMax = " << muMax
  //          << "\n rhoA  = " << rhoA
  //          << "\n rhoZ  = " << rhoZ
  //          << "\n k     = " << k
  //          << "\n alpha = " << alpha
  //          << "\n th    = " << th
  //          << "\n mu    = " << mu
  //          << "\n y'[0] = " << stateDot[0]
  //          << "\n y'[1] = " << stateDot[1]
  //          << "\n y'[2] = " << stateDot[2]
  //          << "\n"
  //          << std::endl;

  return GSL_SUCCESS;
}
