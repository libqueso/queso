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
#include <queso/BoxSubset.h>
#include <queso/MultiDQuadratureBase.h>

namespace QUESOTesting
{
  class LegendreQuadratureTestingHelper
  {
  public:

    void testing_orders( std::vector<unsigned int> & orders )
    {
      // These are the valid listed orders for Legendre quadrature in QUESO;
      // TODO: With C++11, we can initialize this with array syntax
      orders.resize(11);
      orders[0] = 1;
      orders[1] = 2;
      orders[2] = 3;
      orders[3] = 4;
      orders[4] = 5;
      orders[5] = 6;
      orders[6] = 7;
      orders[7] = 10;
      orders[8] = 11;
      orders[9] = 12;
      orders[10] = 16;
    }

  };

  class HermiteQuadratureTestingHelper
  {
  public:

    void testing_orders( std::vector<unsigned int> & orders )
    {
      // These are the valid listed orders for Legendre quadrature in QUESO;
      // TODO: With C++11, we can initialize this with array syntax
      orders.resize(10);
      orders[0] = 1;
      orders[1] = 2;
      orders[2] = 3;
      orders[3] = 4;
      orders[4] = 5;
      orders[5] = 6;
      orders[6] = 7;
      orders[7] = 8;
      orders[8] = 9;
      orders[9] = 19;
    }

  };

  class OneDQuadratureFunction
  {
  public:

    virtual double f( double x ) =0;

    virtual double int_f( double lower, double upper ) =0;
  };

  class PolynomialFunction : public OneDQuadratureFunction
  {
  public:
    PolynomialFunction(int order )
      : _order(order)
    {}

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

  private:

    int _order;
  };


  // For integrating -\infty,\infty using Gauss-Hermite
  // The weighting function is not included in the function evaluation
  class Erf : public OneDQuadratureFunction
  {
  public:

    virtual double f( double /*x*/ )
    {
      return 1.0;
    }

    virtual double int_f( double /*lower*/, double /*upper*/ )
    {
      return std::sqrt(M_PI);
    }
  };

  // For integrating -\infty,\infty using Gauss-Hermite
  // The weighting function is not included in the function evaluation
  class X2Erf : public OneDQuadratureFunction
  {
  public:

    virtual double f( double x )
    {
      return x*x;
    }

    virtual double int_f( double /*lower*/, double /*upper*/ )
    {
      return std::sqrt(M_PI)/2.0;
    }
  };

  // For integrating -\infty,\infty using Gauss-Hermite
  // The weighting function is not included in the function evaluation
  class X4Erf : public OneDQuadratureFunction
  {
  public:

    virtual double f( double x )
    {
      return x*x*x*x;
    }

    virtual double int_f( double /*lower*/, double /*upper*/ )
    {
      return 3.0*std::sqrt(M_PI)/4.0;
    }
  };

  template <class V, class M>
  class MultiDQuadratureFunction
  {
  public:

    virtual double f( const V & x ) const =0;

    virtual double int_f( const QUESO::BoxSubset<V,M> & domain ) const =0;
  };

  template <class V, class M>
  class MultiDLinearFunction : public MultiDQuadratureFunction<V,M>
  {
  public:

    //! e.g. in 3D: f =  1 + x + y + z
    virtual double f( const V & x ) const
    {
      double value = 1.0;

      unsigned int dim = x.sizeGlobal();

      for( unsigned int i = 0; i < dim; i++ )
        value += x[i];

      return value;
    }

    //! e.g. in 3D: \int f = x*y*z + (x^2/2)*y*z + x*(y^2/2)*z + x*y*(z^2/2)
    virtual double int_f( const QUESO::BoxSubset<V,M> & domain ) const
    {
      unsigned int dim = domain.vectorSpace().dimGlobal();
      unsigned int n_terms = dim+1;

      const V & min_values = domain.minValues();
      const V & max_values = domain.maxValues();

      double value = 0.0;

      for( unsigned int t = 0; t < n_terms; t++ )
        {
          double term = 1.0;

          if( t == 0 )
            {
              for( unsigned int i = 0; i < dim; i++ )
                term *= (max_values[i] - min_values[i]);
            }
          else
            {
              for( unsigned int i = 0; i < dim; i++ )
                {
                  if( i == (t-1) )
                    term *= (max_values[i]*max_values[i] - min_values[i]*min_values[i])/2.0;
                  else
                    term *= (max_values[i] - min_values[i]);
                }
            }

          value += term;
        }

      return value;
    }

  };

  //! Base class for testing multi-dimensional quadrature rules
  template <class V, class M>
  class QuadratureMultiDTestBase : public CppUnit::TestCase
  {
  public:

    double integrate_func( const MultiDQuadratureFunction<V,M> & func,
                           QUESO::MultiDQuadratureBase<V,M> & quad_rule,
                           unsigned int n_points )
    {
      const std::vector<typename QUESO::SharedPtr<V>::Type> & x = quad_rule.positions();
      const std::vector<double> & w = quad_rule.weights();

      CPPUNIT_ASSERT_EQUAL(x.size(),w.size());
      CPPUNIT_ASSERT( n_points <= x.size());

      double int_value = 0.0;

      for( unsigned int q = 0; q < n_points; q++ )
        int_value += func.f( (*x[q]) )*w[q];

      return int_value;
    }

    void exact_quadrature_rule_test( const MultiDQuadratureFunction<V,M> & func,
                                     const QUESO::BoxSubset<V,M> & domain,
                                     QUESO::MultiDQuadratureBase<V,M> & quad_rule,
                                     double tol )
    {
      double int_value = this->integrate_func(func,quad_rule,quad_rule.positions().size());

      double exact_int_value = func.int_f( domain );

      double rel_error = std::abs( (int_value-exact_int_value)/exact_int_value );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, rel_error, tol );
    }

    protected:

      QUESO::EnvOptionsValues _options;

      typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;

      void init_env()
      {
        _env.reset( new QUESO::FullEnvironment("","",&_options) );
      }
  };

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
