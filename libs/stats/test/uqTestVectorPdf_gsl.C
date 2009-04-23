#define BOOST_TEST_MODULE $Id: uqTestVectorPdf_gsl.C $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <uqVectorSpace.h>
#include <uqVectorPdf.h>
#include <uqGslMatrix.h>

#define PI 3.14159265358979323846

BOOST_AUTO_TEST_CASE( test_uqGaussianVectorPdfClass )
{
  // Initialize
  MPI_Init(NULL, NULL);
  uqFullEnvironmentClass env;
  uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass> domainSpace(env, "test_space", 2, NULL);
  Epetra_Map eMap(2, 0, env.comm());

  uqGslVectorClass domainMinVal(env, eMap, -1e30);
  uqGslVectorClass domainMaxVal(env, eMap,  1e30);

  uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>* gaussianPdf;
  double tolClose = 1e-13, tolSmall = 1e-16;
  
  //***********************************************************************
  // Tests for diagonal covariance matrix
  // NOTE: distribution is not normalized
  //***********************************************************************

  // mean = [0; 0], var = [1; 1]
  uqGslVectorClass expectedVal(env, eMap, 0.0);
  uqGslVectorClass varianceVal(env, eMap, 1.0);
  uqGslVectorClass testValues(env, eMap, 0.0);

  gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
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

  gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
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

  uqGslMatrixClass covMatrix(env, eMap, 1.0);  // actually diagonal

  gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
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

//   gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
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
