#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/PeriodicBoxSubset.h>
#include <queso/BoxSubset.h>


int main(int argc, char ** argv) {

 
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 3, NULL);

  //////////////////////////////////////////////////////////////////
  // Testing the normal functionality of the periodic domain here.
  //////////////////////////////////////////////////////////////////

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  // Setting it up so that we have domains [0,1], [0,2], [0,3].
  for ( int i=0; i < 3; i++ ){
   paramMaxs[i]= double(i) + 1.; 
  }

  QUESO::PeriodicBoxSubset<QUESO::GslVector, QUESO::GslMatrix> normalPeriodicDomain("param_",
     paramSpace, paramMins, paramMaxs);


  QUESO::GslVector truePeriods(paramMaxs);
  double periodsNorm = truePeriods.norm2();
  truePeriods -= normalPeriodicDomain.periods();
   double periodsRelErr= truePeriods.norm2()/periodsNorm;
   if ( periodsRelErr > 1E-15 ){
     queso_error_msg( "Failure to compute correct periods for periodic box subset." );
   }

   double trueVolume = paramMaxs[0] * paramMaxs[1] * paramMaxs[2];
   double volumeRelErr = abs( trueVolume - normalPeriodicDomain.volume() ) / trueVolume;
   if ( volumeRelErr > 1E-15 ){
     queso_error_msg( "Failure to compute correct volume for periodic box subset." );
   }


  //////////////////////////////////////////////////////////////////
  // Testing that the periodic domain catches infinite domains
  // with an error here.
  //////////////////////////////////////////////////////////////////
  
  paramMaxs[0]=INFINITY;

  try{
  QUESO::PeriodicBoxSubset<QUESO::GslVector, QUESO::GslMatrix> infiniteDomain("param_",
     paramSpace, paramMins, paramMaxs);
  printf("Error for infinite maximum value in domain was not caught by PeriodicBoxSubset.C\n");
  return 1;
  }
  catch(...){
    printf("\nSuccessfully caught error for infinite maximum value in domain.\n");
  }

  paramMaxs[0]=1.;
  paramMins[0]=-INFINITY;

  try{
  QUESO::PeriodicBoxSubset<QUESO::GslVector, QUESO::GslMatrix> infiniteDomain("param_",
     paramSpace, paramMins, paramMaxs);
  printf("Error for infinite minimum value in domain was not caught by PeriodicBoxSubset.C\n"); 
    return 1;
  }
  catch(...){
    printf("\nSuccessfully caught error for infinite minimum value in domain.\n");
  }

  MPI_Finalize();
  
  return 0;
}
