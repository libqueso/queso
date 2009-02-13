#define BOOST_TEST_MODULE $Id: uqTestVectorRV_gsl.C $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <uqVectorSpace.h>
#include <uqVectorRV.h>
#include <uqGslMatrix.h>

BOOST_AUTO_TEST_CASE( test_uqGaussianVectorRVClass )
{
  // Initialize
  MPI_Init(NULL, NULL);
  uqFullEnvironmentClass env;
  uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass> imageSpace(env, "test_space", 2, NULL);
  Epetra_Map eMap(2, 0, env.comm());

  uqGslVectorClass imageMinVal(env, eMap, -INFINITY);
  uqGslVectorClass imageMaxVal(env, eMap,  INFINITY);

  uqGslVectorClass initExpectedValues(env, eMap, 0.0);
  uqGslMatrixClass initCovMatrix(env, eMap, 1.0); 

  uqGslVectorClass finalExpectedValues(env, eMap, 1.0);
  uqGslMatrixClass finalCovMatrix(env, eMap, 3.0);

  uqGslVectorClass testValues(env, eMap, 0.0);

  uqGaussianVectorRVClass<uqGslVectorClass, uqGslMatrixClass> gaussianRV("test_rv", imageSpace, imageMinVal, imageMaxVal,
									 initExpectedValues, initCovMatrix);
  double tolClose = 1e-13, tolSmall = 1e-16;
  
  //***********************************************************************
  // Test pdf
  // NOTE: pdf is not normalized
  //***********************************************************************

  // mean = [0; 0], var = [1; 1], testValues = [0; 0]
  BOOST_REQUIRE_CLOSE(gaussianRV.pdf().actualDensity(testValues), 1.0, tolClose);

  // change mean and check that new pdf is correct
  gaussianRV.updateExpectedValues(finalExpectedValues);
  BOOST_REQUIRE_CLOSE(gaussianRV.pdf().actualDensity(testValues), exp(-1.0), tolClose);

  //***********************************************************************
  // Test realizer
  // NOTE: just calls it... doesn't check values
  //***********************************************************************
  uqGslVectorClass myRealization(testValues);
  gaussianRV.realizer().realization(myRealization);

  std::cout << myRealization;

  // finalize
  MPI_Finalize();
}
