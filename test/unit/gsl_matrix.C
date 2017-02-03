//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include "config_queso.h"

#ifdef QUESO_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <queso/EnvironmentOptions.h>
#include <queso/GslMatrix.h>
#include <queso/VectorRV.h>

namespace QUESOTesting
{
  class GslMatrixTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( GslMatrixTest );

    CPPUNIT_TEST( test_get_set_row_column );
    CPPUNIT_TEST( test_inverse_power_method );
    CPPUNIT_TEST( test_power_method );
    CPPUNIT_TEST( test_multiple_rhs_matrix_solve );
    CPPUNIT_TEST( test_cw_extract );
    CPPUNIT_TEST( test_svd );

    CPPUNIT_TEST_SUITE_END();

    // yes, this is necessary
  public:
    void setUp()
    {
      _env.reset( new QUESO::FullEnvironment("","",&_options) );
    }

    void test_get_set_row_column()
    {
      QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
        paramSpace( (*_env), "param_", 2, NULL);

      typename QUESO::ScopedPtr<QUESO::GslMatrix>::Type
        matrix( paramSpace.newMatrix() );

      (*matrix)(0,0) = 4.; (*matrix)(0,1) = 3.;
      (*matrix)(1,0) = 5.; (*matrix)(1,1) = 7.;

      QUESO::GslVector row((*matrix).getRow(0));
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,0), row[0]);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,1), row[1]);

      row = (*matrix).getRow(1);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(1,0), row[0]);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(1,1), row[1]);

      QUESO::GslVector column((*matrix).getColumn(0));
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,0), column[0]);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(1,0), column[1]);

      column = (*matrix).getColumn(1);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,1), column[0]);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(1,1), column[1]);

      row[0] = 9.;
      row[1] = 15.;
      (*matrix).setRow(0,row);

      row = (*matrix).getRow(0);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,0), row[0]);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,1), row[1]);

      row[0] = 12.;
      row[1] = 2.;
      (*matrix).setColumn(1,row);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(0,1), row[0]);
      CPPUNIT_ASSERT_EQUAL( (*matrix)(1,1), row[1]);
    }

    void test_inverse_power_method()
    {
      QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
        paramSpace( (*_env), "param_", 2, NULL);

      typename QUESO::ScopedPtr<QUESO::GslMatrix>::Type
        matrix( paramSpace.newMatrix() );
      (*matrix)(0,0) = 4.; (*matrix)(0,1) = 3.;
      (*matrix)(1,0) = 5.; (*matrix)(1,1) = 7.;

      double eValue = 0.0;
      QUESO::GslVector eVector( (*paramSpace.newVector()) );

      matrix->smallestEigen( eValue, eVector );

      // MATLAB reports the following for the eigenvalues and eigenvectors
      // of the matrix:
      // A = [ 4, 3; 5, 7];
      // [V,D] = eig(A)
      // V =
      //
      //  -0.749062754969087  -0.468750367387953
      //   0.662499048390352  -0.883330681610041
      //
      // D =
      //
      //   1.346688068540963                   0
      //                   0   9.653311931459037

      double eValueExact = 1.346688068540963;

      // Note the minus sign doesn't matter (as long as the sign change is
      // consistent within the vector). We just need to make sure
      // that we get the right values in the vector.
      QUESO::GslVector eVectorExact( (*paramSpace.newVector() ) );
      eVectorExact[0] =  0.749062754969087;
      eVectorExact[1] = -0.662499048390352;

      // MATLAB returns normalized vectors while the Power method
      // function does not, so we must normalize.
      double norm = eVector.norm2();
      eVector /= norm;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(eValueExact,eValue,1.0e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(eVectorExact[0],eVector[0],1.0e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(eVectorExact[1],eVector[1],1.0e-13);
    }

    void test_power_method()
    {
      QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
        paramSpace( (*_env), "param_", 2, NULL);

      typename QUESO::ScopedPtr<QUESO::GslMatrix>::Type
        matrix( paramSpace.newMatrix() );
      (*matrix)(0,0) = 4.; (*matrix)(0,1) = 3.;
      (*matrix)(1,0) = 5.; (*matrix)(1,1) = 7.;

      double eValue = 0.0;
      QUESO::GslVector eVector( (*paramSpace.newVector()) );

      matrix->largestEigen( eValue, eVector );

      // MATLAB reports the following for the eigenvalues and eigenvectors
      // of the matrix:
      // A = [ 4, 3; 5, 7];
      // [V,D] = eig(A)
      // V =
      //
      //  -0.749062754969087  -0.468750367387953
      //   0.662499048390352  -0.883330681610041
      //
      // D =
      //
      //   1.346688068540963                   0
      //                   0   9.653311931459037

      double eValueExact = 9.653311931459037;

      QUESO::GslVector eVectorExact( (*paramSpace.newVector() ) );
      eVectorExact[0] = 0.468750367387953;
      eVectorExact[1] = 0.883330681610041;

      // MATLAB returns normalized vectors while the Power method
      // function does not, so we must normalize.
      double norm = eVector.norm2();
      eVector /= norm;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(eValueExact,eValue,1.0e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(eVectorExact[0],eVector[0],1.0e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(eVectorExact[1],eVector[1],1.0e-13);
    }

    void test_multiple_rhs_matrix_solve()
    {
      QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
        paramSpace( (*_env), "param_", 2, NULL);

      typename QUESO::ScopedPtr<QUESO::GslMatrix>::Type
        matrix( paramSpace.newMatrix() );
      (*matrix)(0,0) = 4.; (*matrix)(0,1) = 3.;
      (*matrix)(1,0) = 5.; (*matrix)(1,1) = 7.;

      QUESO::GslMatrix result( (*matrix) );

      matrix->invertMultiply( (*matrix), result );

      // We should be getting back the identity matrix in result
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,result(0,0),1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,result(1,1),1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,result(0,1),1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,result(1,0),1.0e-14);
    }

    void test_cw_extract()
    {
      QUESO::VectorSpace<> space4(*_env, "", 4, NULL);
      QUESO::GslMatrix mat4(space4.zeroVector());

      for (unsigned int i = 0; i < 4; i++) {
        for (unsigned int j = 0; j < 4; j++) {
          mat4(i,j) = 4*i+j;
        }
      }

      QUESO::VectorSpace<> space2(*_env, "", 2, NULL);
      QUESO::GslMatrix mat2(space2.zeroVector());

      mat4.cwExtract(1, 1, mat2);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, mat2(0,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, mat2(0,1), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0, mat2(1,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, mat2(1,1), 1.0e-14);
    }

    void test_svd()
    {
      QUESO::VectorSpace<> space(*_env, "", 2, NULL);
      QUESO::GslMatrix M(space.zeroVector());
      QUESO::GslMatrix U(space.zeroVector());
      QUESO::GslMatrix Vs(space.zeroVector());
      QUESO::GslVector S(space.zeroVector());

      M(0,0) = 1.0;
      M(0,1) = 2.0;
      M(1,0) = -2.0;
      M(1,1) = 1.0;

      M.svd(U, S, Vs);

      // Sanity check the decomposition
      QUESO::GslMatrix M_computed(space.zeroVector());
      QUESO::GslMatrix temp(space.zeroVector());
      QUESO::GslMatrix S_matrix(S);

      U.multiply(S_matrix, temp);
      temp.multiply(Vs, M_computed);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(M(0,0), M_computed(0,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(M(0,1), M_computed(0,1), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(M(1,0), M_computed(1,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(M(1,1), M_computed(1,1), 1.0e-14);

      // Check solves
      QUESO::GslMatrix X(space.zeroVector());
      QUESO::GslMatrix X_computed(space.zeroVector());
      QUESO::GslMatrix B(space.zeroVector());

      B(0,0) = 11.0;
      B(0,1) = 17.0;
      B(1,0) = -2.0;
      B(1,1) = -4.0;

      X(0,0) = 3.0;
      X(0,1) = 5.0;
      X(1,0) = 4.0;
      X(1,1) = 6.0;

      M.svdSolve(B, X_computed);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(X(0,0), X_computed(0,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(X(0,1), X_computed(0,1), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(X(1,0), X_computed(1,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(X(1,1), X_computed(1,1), 1.0e-14);

      // Test getters
      CPPUNIT_ASSERT_DOUBLES_EQUAL(U(0,0), M.svdMatU()(0,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(U(0,1), M.svdMatU()(0,1), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(U(1,0), M.svdMatU()(1,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(U(1,1), M.svdMatU()(1,1), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(Vs(0,0), M.svdMatV()(0,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(Vs(0,1), M.svdMatV()(0,1), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(Vs(1,0), M.svdMatV()(1,0), 1.0e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(Vs(1,1), M.svdMatV()(1,1), 1.0e-14);
    }

  private:
    QUESO::EnvOptionsValues _options;
    typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslMatrixTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
