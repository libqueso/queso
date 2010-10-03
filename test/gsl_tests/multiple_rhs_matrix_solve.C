//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <multiple_rhs_matrix_solve.h>
#include <uqGslMatrix.h>
#include <uqVectorRV.h>

int main(int argc, char* argv[]) 
{
  int return_flag = 0;

   // Initialize environment
  MPI_Init(&argc,&argv);
  uqFullEnvironmentClass* env =
    new uqFullEnvironmentClass(MPI_COMM_WORLD,"input","",NULL);

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

  const double fp_tol = 1.0e-14;

  // Instantiate the parameter space
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace( (*env), "param_", 2, NULL);

  // Populate Matrix.
  uqGslMatrixClass* Matrix = paramSpace.newMatrix();
  (*Matrix)(0,0) = 4.; (*Matrix)(0,1) = 3.;
  (*Matrix)(1,0) = 5.; (*Matrix)(1,1) = 7.;

  // Conduct test.
  uqGslMatrixClass Result( (*Matrix) );

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
