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

    virtual double f( double x ) =0;

    virtual double int_f( double lower, double upper ) =0;

    void test_quadrature_rule( QUESO::Base1DQuadrature & quad_rule,
                               double min_domain_value, double max_domain_value,
                               double tol )
    {
      const std::vector<double> & x = quad_rule.positions();
      const std::vector<double> & w = quad_rule.weights();

      CPPUNIT_ASSERT_EQUAL(x.size(),w.size());

      unsigned int n_points = x.size();

      double int_value = 0.0;
      for( unsigned int q = 0; q < n_points; q++ )
        int_value += this->f(x[q])*w[q];

      double exact_int_value = this->int_f(min_domain_value,max_domain_value);

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

    // e.g. 4*x^3 + 3*x^2 + 2*x + 1
    virtual double f( double x )
    {
      CPPUNIT_ASSERT(_order>=0);
      double value = 0.0;

      for( int i = 0; i <= _order; i++ )
        value += (i+1)*std::pow( x, i );

      return value;
    }

    // Integral of the function in f()
    virtual double int_f( double lower, double upper )
    {
      CPPUNIT_ASSERT(_order>=0);
      double value = 0.0;

      for( int i = 0; i <= _order; i++ )
        value += (std::pow(upper,i+1.0) - std::pow(lower,i+1.0) );

      return value;
    }

    void setUp()
    {
      // Setup for error checking, in order to make sure this gets reset later.
      _order = -1;
    }

  private:

    int _order;

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
        _order = 2*(quad_rule.order()+1)-1;

        double tol = std::numeric_limits<double>::epsilon()*20;

        // 7th order quadrature only has positions/weights out to 8 digits...
        if( quad_order == 7 )
          tol = 1.0e-7;

        this->test_quadrature_rule( quad_rule, min_domain_value, max_domain_value, tol );
      }
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( LegendreQuadrature1DTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
