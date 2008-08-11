#include <uqChemEx.h>
#include <uqApplRoutines.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  uqEnvironmentClass* env = new uqEnvironmentClass(argc,argv);

  //************************************************
  // Call application
  //************************************************
  uqAppl<uqGslVectorClass,uqGslMatrixClass>(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;

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
