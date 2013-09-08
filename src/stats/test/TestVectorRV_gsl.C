#define BOOST_TEST_MODULE $Id: uqTestVectorRV_gsl.C $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <uqVectorSpace.h>
#include <uqVectorRV.h>
#include <uqGslMatrix.h>

using namespace QUESO;

BOOST_AUTO_TEST_CASE( test_uqGaussianVectorRV )
{
  // Initialize
  MPI_Init(NULL, NULL);
  uqFullEnvironment env;
  uqVectorSpace<uqGslVector, uqGslMatrix> imageSpace(env, "test_space", 2, NULL);
  uqMap eMap(2, 0, env.comm());

  uqGslVector imageMinVal(env, eMap, -INFINITY);
  uqGslVector imageMaxVal(env, eMap,  INFINITY);

  uqGslVector initExpectedValues(env, eMap, 0.0);
  uqGslMatrix initCovMatrix(env, eMap, 1.0); 

  uqGslVector finalExpectedValues(env, eMap, 1.0);
  uqGslMatrix finalCovMatrix(env, eMap, 3.0);

  uqGslVector testValues(env, eMap, 0.0);

  uqGaussianVectorRV<uqGslVector, uqGslMatrix> gaussianRV("test_rv", imageSpace, imageMinVal, imageMaxVal,
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
  BOOST_REQUIRE_CLOSE(gaussianRV.pdf().actualDensity(testValues), std::exp(-1.0), tolClose);

  //***********************************************************************
  // Test realizer
  // NOTE: just calls it... doesn't check values
  //***********************************************************************
  uqGslVector myRealization(testValues);
  gaussianRV.realizer().realization(myRealization);

  std::cout << myRealization;

  // finalize
  MPI_Finalize();
}
