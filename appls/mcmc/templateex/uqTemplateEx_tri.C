#include <uqTemplateEx.h>
#include <uqTrilinosMatrix.h>

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
  uqAppl<uqTrilinosVectorClass,uqTrilinosMatrixClass>(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
  MPI_Finalize();

  return 0;
}
