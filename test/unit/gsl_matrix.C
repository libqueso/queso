//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
    CPPUNIT_TEST( inverse_power_method );

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

    void inverse_power_method()
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

  private:

    QUESO::EnvOptionsValues _options;
    typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslMatrixTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
