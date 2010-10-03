//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
// This is an example very similar to
// examples/validationCycle/src/exTgaValidationCycle_gsl.C.  It just
// generates smaller chains and so is faster and more suitable to
// serve as an automatic regression test. Indeed, it was set by Karl
// as a first example to the QUESO team on how to run regression tests
// with buildbot and how to setup them in the configure & make package
// of QUESO.
//
//--------------------------------------------------------------------------


#include <TgaValidationCycle_appl.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
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
  // Run application
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
}
