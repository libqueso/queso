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

#include <quadrature_testing_helper.h>

#include <queso/MonteCarloQuadrature.h>
#include <queso/TensorProductQuadrature.h>
#include <queso/1DQuadrature.h>

#include <convergence_rate_helper.h>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <cmath>
#include <limits>

namespace QUESOTesting
{
  template <class V, class M>
  class TensorProductLegendreQuadratureRuleTestBase :
    public QuadratureMultiDTestBase<V,M>,
    public LegendreQuadratureTestingHelper
  {
  public:

    void setUp()
    {
      this->init_env();
    }

    void test_legendre_linear_func_2d()
    {
      std::vector<unsigned int> orders;
      this->testing_orders(orders);

      double tol = std::numeric_limits<double>::epsilon()*10;
      this->test_legendre_linear_func_dim_wrapper(2,tol,orders);
    }

    void test_legendre_linear_func_3d()
    {
      std::vector<unsigned int> orders;
      this->testing_orders(orders);

      double tol = std::numeric_limits<double>::epsilon()*50;
      this->test_legendre_linear_func_dim_wrapper(3,tol,orders);
    }

    void test_legendre_linear_func_4d()
    {
      std::vector<unsigned int> orders;
      this->testing_orders(orders);

      double tol = std::numeric_limits<double>::epsilon()*500;
      this->test_legendre_linear_func_dim_wrapper(4,tol,orders);
    }

    void test_legendre_linear_func_5d()
    {
      std::vector<unsigned int> orders;

      // To keep testing time reasonable, only test orders < 10
      // Currently, Legendre missing orders 8,9,10.
      for( unsigned int i = 1; i < 7; i++ )
        orders.push_back(i);

      double tol = std::numeric_limits<double>::epsilon()*5000;
      this->test_legendre_linear_func_dim_wrapper(5,tol,orders);
    }

  private:

    void test_legendre( unsigned int dim,
                        double min_domain_value,
                        double max_domain_value,
                        unsigned int order,
                        const MultiDQuadratureFunction<V,M> & func,
                        double tol )
    {
      // Instantiate the parameter space
      QUESO::VectorSpace<V,M> param_space( (*this->_env), "param_", dim, NULL);

      typename QUESO::ScopedPtr<V>::Type min_values( param_space.newVector(min_domain_value) );
      typename QUESO::ScopedPtr<V>::Type max_values( param_space.newVector(max_domain_value) );

      QUESO::BoxSubset<V,M> param_domain( "param_domain_", param_space, (*min_values), (*max_values) );

      QUESO::SharedPtr<QUESO::Base1DQuadrature>::Type qrule_1d( new QUESO::UniformLegendre1DQuadrature(min_domain_value, max_domain_value, order, false) );

      std::vector<QUESO::SharedPtr<QUESO::Base1DQuadrature>::Type> all_1d_qrules(dim,qrule_1d);

      QUESO::TensorProductQuadrature<V,M> tp_qrule(param_domain,all_1d_qrules);

      this->exact_quadrature_rule_test( func,
                                        param_domain,
                                        tp_qrule,
                                        tol );
    }

    void test_legendre_linear_func_dim_wrapper(unsigned int dim,
                                               double tol,
                                               const std::vector<unsigned int> & orders)
    {
      MultiDLinearFunction<V,M> func;

      for( std::vector<unsigned int>::const_iterator oit = orders.begin();
           oit < orders.end(); ++oit )
        {
          // "standard" bounds
          this->test_legendre(dim,-1,1,*oit,func,tol);

          // "non-standard" bounds
          this->test_legendre(dim,-3.14,2.71,*oit,func,tol);
        }
    }

  };

  class TensorProductLegendreQuadratureRuleGslTest :
    public TensorProductLegendreQuadratureRuleTestBase<QUESO::GslVector,QUESO::GslMatrix>
  {
  public:

    CPPUNIT_TEST_SUITE( TensorProductLegendreQuadratureRuleGslTest );

    CPPUNIT_TEST( test_legendre_linear_func_2d );
    CPPUNIT_TEST( test_legendre_linear_func_3d );
    CPPUNIT_TEST( test_legendre_linear_func_4d );
    CPPUNIT_TEST( test_legendre_linear_func_5d );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( TensorProductLegendreQuadratureRuleGslTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
