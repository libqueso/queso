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

#include <queso/1DQuadrature.h>

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

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
