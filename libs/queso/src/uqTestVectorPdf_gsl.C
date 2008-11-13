#define BOOST_TEST_MODULE $Id: uqTestVectorPdf_gsl.C $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <uqVectorSpace.h>
#include <uqVectorPdf.h>
#include <uqGslMatrix.h>

#define PI 3.14159265358979323846

BOOST_AUTO_TEST_CASE( test_uqGaussianVectorPdfClass_diagonalCovariance )
{
  // Initialize
  MPI_Init(NULL, NULL);
  uqFullEnvironmentClass env;
  uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass> domainSpace(env, "test_space", 2, NULL);
  Epetra_Map eMap(2, 0, env.comm());

  uqGslVectorClass domainMinVal(env, eMap, -1e30);
  uqGslVectorClass domainMaxVal(env, eMap,  1e30);
  uqGslVectorClass expectedVal(env, eMap, 0.0);

  uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>* gaussianPdf;

  // Tests
  double tolClose = 1e-14, tolSmall = 1e-16;
  
  // mean = [0; 0], var = [1; 1]
  uqGslVectorClass varianceVal(env, eMap, 1.0);
  uqGslVectorClass testValues(env, eMap, 0.0);

  gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
										 domainMaxVal, expectedVal, varianceVal);


  // NOTE: distribution is not normalized

  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), 1.0, tolClose);
  BOOST_REQUIRE_SMALL(gaussianPdf->minus2LnDensity(testValues),    tolSmall); // can't check close b/c exact val = 0

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    2.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    1.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    2.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    1.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    2.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    1.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    2.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-0.5), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    1.0 , tolClose);

  delete gaussianPdf;

  // mean = [0; 0], var = [0.25; 0.5]
  varianceVal[0] = 0.25; varianceVal[1] = 0.5;

  gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
										 domainMaxVal, expectedVal, varianceVal);

  testValues[0] = testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), 1.0, tolClose);
  BOOST_REQUIRE_SMALL(gaussianPdf->minus2LnDensity(testValues),    tolSmall); // can't check close b/c exact val = 0

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    6.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    2.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    6.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-2.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    4.0 , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    6.0 , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-1.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    2.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-3.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    6.0 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-2.0), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    4.0 , tolClose);

  delete gaussianPdf;

  // mean = [1.0; -0.5], var = [0.25; 0.5]
  expectedVal[0] = 1.0; expectedVal[1] = -0.5;

  gaussianPdf = new uqGaussianVectorPdfClass<uqGslVectorClass, uqGslMatrixClass>("test_pdf", domainSpace, domainMinVal, 
										 domainMaxVal, expectedVal, varianceVal);

  testValues[0] = 0.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-2.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    4.5  , tolClose); 

  testValues[0] = 1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-2.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    4.5 , tolClose);

  testValues[0] = 0.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-4.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    8.5  , tolClose);

  testValues[0] = -1.0; testValues[1] = 1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-10.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    20.5  , tolClose);

  testValues[0] = -1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-8.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),   16.5  , tolClose);

  testValues[0] = -1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-8.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),   16.5  , tolClose);

  testValues[0] = 0.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-2.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    4.5  , tolClose);

  testValues[0] = 1.0; testValues[1] = -1.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-0.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    0.5 , tolClose);

  testValues[0] = 1.0; testValues[1] = 0.0;
  BOOST_REQUIRE_CLOSE(gaussianPdf->actualDensity(testValues), exp(-0.25), tolClose);
  BOOST_REQUIRE_CLOSE(gaussianPdf->minus2LnDensity(testValues),    0.5  , tolClose);

  delete gaussianPdf;

  // Clean up
  MPI_Finalize();
}
