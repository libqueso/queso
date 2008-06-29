#include <uqAlgaeEx.h>
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
uqGslAlgaeExStateDot( // Compute state dot = dConcentrations/dt
  double       t,
  const double currentState[],           // current concentrations
  double       stateDot[],               // dConcentrations/dt
  void*        infoForComputingStateDot) // concentration rates
{
  double muMax = ((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->muMax;
  double rhoA  = ((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->rhoA;
  double rhoZ  = ((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->rhoZ;
  double k     = ((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->k;
  double alpha = ((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->alpha;
  double th    = ((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->th;

  const std::vector<double>& evolutionOfQpV = *((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->evolutionOfQpV;
  const std::vector<double>& evolutionOfT   = *((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->evolutionOfT;
  const std::vector<double>& evolutionOfPin = *((uqAppl_StateDotFunction_DataType<uqGslVectorClass,uqGslMatrixClass> *)infoForComputingStateDot)->evolutionOfPin;

  bool bRC = ((evolutionOfQpV.size() == evolutionOfPin.size()) &&
              (evolutionOfQpV.size() == evolutionOfT.size()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      UQ_UNAVAILABLE_RANK,
                      "uqGslAlgaeExStateDot()",
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
                      "uqGslAlgaeExStateDot()",
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

  //std::cout << "In uqGslAlgaeExStateDot():"
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
