#include <uqEtc.h>
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
