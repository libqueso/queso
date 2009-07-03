#include <stdlib.h>
#include <stdio.h>
#include <queso.h>

//----------------------------------------
// A simple Calibration example using the
// basic QUESO interface..
//---------------------------------------

double my_likelihood(double *params);

int main(int argc, char *argv[])
{

  MPI_Init(&argc,&argv);

  printf("--> Initializing QUESO Environment...\n");

  /* Initialize QUESO environment and define input file */

  QUESO_init("queso.inp");	                   

 /* Register an application likelihood function with QUESO and run a
  * statistical inversion problem using MCMC sampling */

  QUESO_statistical_inversion(my_likelihood);     

  /* Finalize the analysis and output statistics */

  QUESO_finalize();
  
  MPI_Finalize();
  return 0;

}

double my_likelihood(double *params)
{
  static double needle   = 42.;   
  static double variance = 1.*1.;  /* sigma^2 */

  return( (params[0] - needle)*(params[0] - needle)/(variance) );
}
