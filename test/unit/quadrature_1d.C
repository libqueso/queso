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

#include <queso/1DQuadrature.h>

#include <quadrature_testing_helper.h>

#include <cmath>
#include <limits>

namespace QUESOTesting
{
  class Quadrature1DTestBase : public CppUnit::TestCase
  {
  public:

    void test_quadrature_rule( QUESO::Base1DQuadrature & quad_rule,
                               OneDQuadratureFunction & func,
                               double min_domain_value, double max_domain_value,
                               double tol )
    {
      const std::vector<double> & x = quad_rule.positions();
      const std::vector<double> & w = quad_rule.weights();

      CPPUNIT_ASSERT_EQUAL(x.size(),w.size());

      unsigned int n_points = x.size();

      double int_value = 0.0;
      for( unsigned int q = 0; q < n_points; q++ )
        int_value += func.f(x[q])*w[q];

      double exact_int_value = func.int_f(min_domain_value,max_domain_value);

      double rel_error = std::abs( (int_value-exact_int_value)/exact_int_value );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, rel_error, tol );
    }
  };

  class LegendreQuadrature1DTest : public Quadrature1DTestBase,
                                   public LegendreQuadratureTestingHelper
  {
  public:
    CPPUNIT_TEST_SUITE( LegendreQuadrature1DTest );

    CPPUNIT_TEST( test_1d_legendre_quadrature_unscaled_interval );
    CPPUNIT_TEST( test_1d_legendre_quadrature_scaled_interval );

    CPPUNIT_TEST_SUITE_END();

    // yes, this is necessary
  public:

    void test_1d_legendre_quadrature_unscaled_interval()
    {
      // This test only checks on the interval \pm 1
      double min_domain_value=-1.0;
      double max_domain_value=1.0;

      this->test_1d_legendre_quadrature(min_domain_value,max_domain_value);
    }

    void test_1d_legendre_quadrature_scaled_interval()
    {
      // This test only checks on the interval \pm 1
      double min_domain_value=-2.7183;
      double max_domain_value=3.1415;

      this->test_1d_legendre_quadrature(min_domain_value,max_domain_value);
    }

  private:

    void test_1d_legendre_quadrature(double min_domain_value, double max_domain_value)
    {
      std::vector<unsigned int> testing_orders;
      this->testing_orders( testing_orders );

      for( std::vector<unsigned int>::const_iterator it = testing_orders.begin();
           it < testing_orders.end(); ++it )
      {
        unsigned int quad_order = *it;
        QUESO::UniformLegendre1DQuadrature quad_rule(min_domain_value,max_domain_value,quad_order,false);

        // We should be able to integrate polynomials of order 2*n-1 exactly
        // order = quad_order + 1 by QUESO convention
        PolynomialFunction func( 2*(quad_rule.order()+1)-1 );

        double tol = std::numeric_limits<double>::epsilon()*20;

        // 7th order quadrature only has positions/weights out to 8 digits...
        if( quad_order == 7 )
          tol = 1.0e-7;

        this->test_quadrature_rule( quad_rule, func,
                                    min_domain_value, max_domain_value, tol );
      }
    }
  };

  class HermiteQuadrature1DTest : public Quadrature1DTestBase,
                                  public HermiteQuadratureTestingHelper
  {
  public:
    CPPUNIT_TEST_SUITE( HermiteQuadrature1DTest );

    CPPUNIT_TEST( test_1d_hermite_quadrature_erf );
    CPPUNIT_TEST( test_1d_hermite_quadrature_x2erf );
    CPPUNIT_TEST( test_1d_hermite_quadrature_x4erf );

    CPPUNIT_TEST_SUITE_END();

    // yes, this is necessary
  public:

    void test_1d_hermite_quadrature_erf()
    {
      Erf func;
      std::vector<unsigned int> testing_orders;
      this->testing_orders(testing_orders);

      double tol = std::numeric_limits<double>::epsilon()*50000;

      this->test_1d_hermite_quadrature(func,testing_orders,tol);
    }

    void test_1d_hermite_quadrature_x2erf()
    {
      X2Erf func;
      std::vector<unsigned int> testing_orders;
      this->testing_orders(testing_orders);

      double tol = 1.0e-10;

      this->test_1d_hermite_quadrature(func,testing_orders,tol);
    }

    void test_1d_hermite_quadrature_x4erf()
    {
      std::cout << std::endl << "Beginning Hermite X4Erf test!" << std::endl;
      X4Erf func;
      std::vector<unsigned int> testing_orders;
      this->testing_orders(testing_orders);

      // We need at least a 3-point rule, so get rid of the first entry
      testing_orders.erase(testing_orders.begin());

      double tol = 2.0e-10;

      this->test_1d_hermite_quadrature(func,testing_orders,tol);
    }

  private:

    void test_1d_hermite_quadrature(OneDQuadratureFunction & func,
                                    const std::vector<unsigned int> & testing_orders,
                                    double tol)
    {
      for( std::vector<unsigned int>::const_iterator it = testing_orders.begin();
           it < testing_orders.end(); ++it )
      {
        unsigned int quad_order = *it;

        // 0 mean, unit variance
        QUESO::GaussianHermite1DQuadrature quad_rule(0.0,1.0,quad_order);

        this->test_quadrature_rule( quad_rule, func,
                                    -INFINITY, INFINITY, tol );
      }
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( LegendreQuadrature1DTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( HermiteQuadrature1DTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
