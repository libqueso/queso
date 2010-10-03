//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <get_min_max_vec.h>
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

  // Deallocate pointers we created.
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

  // Instantiate the parameter domain
  uqGslVectorClass Vector( paramSpace.zeroVector() );

  Vector[0] = -4.;
  Vector[1] =  3.;

  // Conduct tests.
  double value;
  int index;
  Vector.getMinValueAndIndex( value, index );

  if( index != 0 || value != Vector[0] )
    {
      return_flag = 1;
    }

  Vector.getMaxValueAndIndex( value, index );
  if( index != 1 || value != Vector[1] )
    {
      return_flag = 1;
    }

  return return_flag;
}
