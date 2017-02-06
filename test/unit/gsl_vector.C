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
#include <queso/GslVector.h>
#include <queso/VectorRV.h>

namespace QUESOTesting
{
  class GslVectorTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( GslVectorTest );
    CPPUNIT_TEST( test_min_max );
    CPPUNIT_TEST( test_beta );
    CPPUNIT_TEST( test_gamma );
    CPPUNIT_TEST( test_inverse_gamma );
    CPPUNIT_TEST_SUITE_END();

    // yes, this is necessary
  public:

    void setUp()
    {
      _env.reset( new QUESO::FullEnvironment("","",&_options) );
    }

    void test_min_max()
    {
      // Instantiate the parameter space
      QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
        paramSpace( (*_env), "param_", 2, NULL);

      // Instantiate the parameter domain
      QUESO::GslVector vec( paramSpace.zeroVector() );

      vec[0] = -4.;
      vec[1] =  3.;

      double value;
      int index;
      vec.getMinValueAndIndex( value, index );

      CPPUNIT_ASSERT_EQUAL(0,index);
      CPPUNIT_ASSERT_EQUAL(vec[0],value);

      vec.getMaxValueAndIndex( value, index );

      CPPUNIT_ASSERT_EQUAL(1,index);
      CPPUNIT_ASSERT_EQUAL(vec[1],value);
    }

    void test_beta()
    {
      unsigned int num_samples = 1000000;

      // Instantiate the parameter space
      QUESO::VectorSpace<> paramSpace(*_env, "", num_samples, NULL);

      // Instantiate the parameter domain
      QUESO::GslVector vec(paramSpace.zeroVector());
      QUESO::GslVector a(paramSpace.zeroVector());
      QUESO::GslVector b(paramSpace.zeroVector());
      a.cwSet(2.0);
      b.cwSet(2.0);
      vec.cwSetBeta(a, b);

      double mean = vec.sumOfComponents() / num_samples;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, mean, 1e-2);
    }

    void test_gamma()
    {
      unsigned int num_samples = 1000000;

      // Instantiate the parameter space
      QUESO::VectorSpace<> paramSpace(*_env, "", num_samples, NULL);

      // Instantiate the parameter domain
      QUESO::GslVector vec(paramSpace.zeroVector());
      QUESO::GslVector a(paramSpace.zeroVector());
      QUESO::GslVector b(paramSpace.zeroVector());
      a.cwSet(10.0);
      b.cwSet(0.5);
      vec.cwSetGamma(a, b);

      double mean = vec.sumOfComponents() / num_samples;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, mean, 1e-2);
    }

    void test_inverse_gamma()
    {
      unsigned int num_samples = 1000000;

      // Instantiate the parameter space
      QUESO::VectorSpace<> paramSpace(*_env, "", num_samples, NULL);

      // Instantiate the parameter domain
      QUESO::GslVector vec(paramSpace.zeroVector());
      QUESO::GslVector a(paramSpace.zeroVector());
      QUESO::GslVector b(paramSpace.zeroVector());
      a.cwSet(11.0);
      b.cwSet(1.0);
      vec.cwSetInverseGamma(a, b);

      double mean = vec.sumOfComponents() / num_samples;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, mean, 1e-3);
    }

  private:

    QUESO::EnvOptionsValues _options;
    typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslVectorTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
