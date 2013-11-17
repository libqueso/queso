#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <queso/VectorSpace.h>
#include <queso/VectorPdf.h>
#include <queso/GslMatrix.h>

#define PI 3.14159265358979323846

using namespace QUESO;

BOOST_AUTO_TEST_CASE( test_uqGaussianVectorPdf )
{
  // Initialize
  MPI_Init(NULL, NULL);
  uqFullEnvironment env;
  uqVectorSpace<uqGslVector, uqGslMatrix> domainSpace(env, "test_space", 2, NULL);
  uqMap eMap(2, 0, env.comm());

  uqGslVector domainMinVal(env, eMap, -1e30);
  uqGslVector domainMaxVal(env, eMap,  1e30);

  uqGaussianVectorPdf<uqGslVector, uqGslMatrix>* gaussianPdf;
  double tolClose = 1e-13, tolSmall = 1e-16;
  
  //***********************************************************************
  // Tests for diagonal covariance matrix
  // NOTE: distribution is not normalized
  //***********************************************************************

  // mean = [0; 0], var = [1; 1]
  uqGslVector expectedVal(env, eMap, 0.0);
  uqGslVector varianceVal(env, eMap, 1.0);
  uqGslVector testValues(env, eMap, 0.0);

  gaussianPdf = new uqGaussianVectorPdf<uqGslVector, uqGslMatrix>("test_pdf", domainSpace, domainMinVal, 
										 domainMaxVal, expectedVal, varianceVal);

  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), 1.0, tolClose);
  BOOST_REQUIRE_SMALL(gaussianPdf->minus2LnDensity(testValues),    tolSmall); // can't check close b/c exact val = 0

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  delete gaussianPdf;

  // mean = [0; 0], var = [0.25; 0.5]
  varianceVal[0] = 0.25; varianceVal[1] = 0.5;

  gaussianPdf = new uqGaussianVectorPdf<uqGslVector, uqGslMatrix>("test_pdf", domainSpace, domainMinVal, 
										 domainMaxVal, expectedVal, varianceVal);
  
  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), 1.0, tolClose);
  BOOST_REQUIRE_SMALL(gaussianPdf->minus2LnDensity(testValues),    tolSmall); // can't check close b/c exact val = 0

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         6.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         6.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-2.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         4.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         6.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         6.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-2.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         4.0 , tolClose);

  // mean = [1.0; -0.5], var = [0.25; 0.5]
  expectedVal[0] = 1.0; expectedVal[1] = -0.5;

  gaussianPdf->updateExpectedValues(expectedVal); // just update expected value (don't reallocate everything)

  testValues[0] = 0.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-2.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         4.5  , tolClose); 

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-2.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         4.5 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-4.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         8.5  , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-10.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         20.5  , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-8.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),        16.5  , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-8.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),        16.5  , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-2.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         4.5  , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         0.5 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         0.5  , tolClose);

  delete gaussianPdf;

  //***********************************************************************
  // Tests for general covariance matrix
  // NOTE: distribution is not normalized
  //***********************************************************************
  
  // mean = [0; 0], covar = [1, 0; 0, 1], i.e. same as first case for diagonal matrices
  expectedVal[0] = expectedVal[1] = 0.0;

  uqGslMatrix covMatrix(env, eMap, 1.0);  // actually diagonal

  gaussianPdf = new uqGaussianVectorPdf<uqGslVector, uqGslMatrix>("test_pdf", domainSpace, domainMinVal, 
										 domainMaxVal, expectedVal, covMatrix);

  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), 1.0, tolClose);
  BOOST_REQUIRE_SMALL(gaussianPdf->minus2LnDensity(testValues),    tolSmall); // can't check close b/c exact val = 0

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.0 , tolClose);

  //delete gaussianPdf;

  // mean = [0, 0], covar = [2, 1; 1, 2];
  covMatrix(0,0) = 2.0; covMatrix(0,1) = 1.0;
  covMatrix(1,0) = 1.0; covMatrix(1,1) = 2.0;

//   gaussianPdf = new uqGaussianVectorPdf<uqGslVector, uqGslMatrix>("test_pdf", domainSpace, domainMinVal, 
// 										 domainMaxVal, expectedVal, covMatrix);

  gaussianPdf->updateCovMatrix(covMatrix);

  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), 1.0, tolClose);
  BOOST_REQUIRE_SMALL(gaussianPdf->minus2LnDensity(testValues),    tolSmall); // can't check close b/c exact val = 0

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0/3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0/3.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0/3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0/3.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         2.0 , tolClose);

  // mean = [1.0, -0.5], covar = [2, 1; 1, 2];
  expectedVal[0] = 1.0; expectedVal[1] = -0.5;

  gaussianPdf->updateExpectedValues(expectedVal); // just update expected value (don't reallocate everything)

  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-5.833333333333333e-01), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.166666666666667e+00 , tolClose);

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-0.75), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         1.5 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-1.583333333333333e+00), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         3.166666666666667e+00 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), std::exp(-3.083333333333333e+00), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),         6.166666666666666e+00 , tolClose);

  delete gaussianPdf;

  // Clean up
  MPI_Finalize();
}
