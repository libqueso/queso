//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <get_set_row_column.h>
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
