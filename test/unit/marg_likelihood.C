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
#include <queso/TensorProductQuadrature.h>
#include <queso/LikelihoodBase.h>
#include <queso/UniformVectorRV.h>
#include <queso/1DQuadrature.h>
#include <queso/GaussianVectorRV.h>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <cmath>
#include <limits>

namespace QUESOTesting
{
  template <class V, class M>
  class TestlingLikelihoodBase : public QUESO::LikelihoodBase<V,M>
  {
  public:

    TestlingLikelihoodBase(const char * prefix,
                          const QUESO::VectorSet<V, M> & domainSet,
                          const V & observations,
                          typename QUESO::SharedPtr<QUESO::BaseVectorRV<V,M> >::Type & marg_param_pdf,
                          typename QUESO::SharedPtr<QUESO::MultiDQuadratureBase<V,M> >::Type & marg_integration,
                          bool marg_pdf_is_weight_func)
      : QUESO::LikelihoodBase<V,M>(prefix,domainSet,observations,marg_param_pdf,marg_integration,marg_pdf_is_weight_func)
    {}

    //! Evaluate the exact ln value of the marginalized likelihood
    virtual double lnMargLikelihood( const V & domainVector ) const =0;

    void testLikelihoodValue( const V & domainVector,
                              const QUESO::LikelihoodBase<V,M> & likelihood,
                              double tol ) const
    {
      double computed_value = likelihood.lnValue(domainVector);

      double exact_value = this->lnMargLikelihood(domainVector);

      double rel_error = std::abs( (computed_value-exact_value)/exact_value );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, rel_error, tol );
    }

    protected:

      virtual double lnLikelihood(const V & /*domainVector*/,
                                  V & modelOutput) const
      {
        queso_require_equal_to(1,modelOutput.sizeGlobal());
        return std::log(modelOutput[0]);
      }

  };

  template <class V, class M>
  class LinearModelLikelihood : public TestlingLikelihoodBase<V,M>
  {
  public:

    LinearModelLikelihood(const char * prefix,
                          const QUESO::VectorSet<V, M> & domainSet,
                          const V & observations,
                          typename QUESO::SharedPtr<QUESO::BaseVectorRV<V,M> >::Type & marg_param_pdf,
                          typename QUESO::SharedPtr<QUESO::MultiDQuadratureBase<V,M> >::Type & marg_integration,
                          bool marg_pdf_is_weight_func)
      : TestlingLikelihoodBase<V,M>(prefix,domainSet,observations,marg_param_pdf,marg_integration,marg_pdf_is_weight_func)
    {}

    virtual void evaluateModel(const V & domainVector,
                               const V & marginalVector,
                               V & modelOutput) const
    {
      queso_require_equal_to(1,domainVector.sizeGlobal());
      queso_require_equal_to(3,marginalVector.sizeGlobal());
      queso_require_equal_to(1,modelOutput.sizeGlobal());

      modelOutput[0] = 3.14*domainVector[0] + 1.1*marginalVector[0]
        + 2.2*marginalVector[1] + 3.3*marginalVector[2];
    }

    virtual double lnMargLikelihood( const V & domainVector ) const
    {
      queso_require_equal_to(1,domainVector.sizeGlobal());

      return std::log(3.14*domainVector[0] + 1.1/2.0 + 2.2/2.0 + 3.3/2.0);
    }

  };

  template <class V, class M>
  class GaussModelLikelihood : public TestlingLikelihoodBase<V,M>
  {
  public:

    GaussModelLikelihood(const char * prefix,
                         const QUESO::VectorSet<V, M> & domainSet,
                         const V & observations,
                         typename QUESO::SharedPtr<QUESO::BaseVectorRV<V,M> >::Type & marg_param_pdf,
                         typename QUESO::SharedPtr<QUESO::MultiDQuadratureBase<V,M> >::Type & marg_integration,
                         bool marg_pdf_is_weight_func)
      : TestlingLikelihoodBase<V,M>(prefix,domainSet,observations,marg_param_pdf,marg_integration,marg_pdf_is_weight_func)
    {}

    virtual void evaluateModel(const V & domainVector,
                               const V & marginalVector,
                               V & modelOutput) const
    {
      queso_require_equal_to(1,domainVector.sizeGlobal());
      queso_require_equal_to(2,marginalVector.sizeGlobal());
      queso_require_equal_to(1,modelOutput.sizeGlobal());

      double q1 = marginalVector[0];
      double q2 = marginalVector[1];
      double m = domainVector[0];

      modelOutput[0] = 3.14*m + 1.1*q1 + 2.2*q2*q2;
    }

    virtual double lnMargLikelihood( const V & domainVector ) const
    {
      queso_require_equal_to(1,domainVector.sizeGlobal());

      double m = domainVector[0];

      return std::log(3.14*m*M_PI + 2.2*M_PI/2.0);
    }

  };

  template <class V, class M>
  class MarginalLikelihoodTestBase : public CppUnit::TestCase
  {
  public:

    void setUp()
    {
      this->init_env();
    }

    void test_linear_func_uniform_marg_space()
    {
      // Instantiate the parameter space
      unsigned int param_dim = 1;
      unsigned int marg_dim = 3;
      QUESO::VectorSpace<V,M> param_space( (*this->_env), "param_", param_dim, NULL);
      QUESO::VectorSpace<V,M> marg_space( (*this->_env), "marg_", marg_dim, NULL);

      double param_min_domain_value = 0.0;
      double param_max_domain_value = 1.0;

      double marg_min_domain_value = 0.0;
      double marg_max_domain_value = 1.0;

      typename QUESO::ScopedPtr<V>::Type param_min_values( param_space.newVector(param_min_domain_value) );
      typename QUESO::ScopedPtr<V>::Type param_max_values( param_space.newVector(param_max_domain_value) );

      typename QUESO::ScopedPtr<V>::Type marg_min_values( marg_space.newVector(marg_min_domain_value) );
      typename QUESO::ScopedPtr<V>::Type marg_max_values( marg_space.newVector(marg_max_domain_value) );

      QUESO::BoxSubset<V,M> param_domain( "param_domain_", param_space, (*param_min_values), (*param_max_values) );
      QUESO::BoxSubset<V,M> marg_domain( "marg_domain_", marg_space, (*marg_min_values), (*marg_max_values) );

      typename QUESO::SharedPtr<QUESO::BaseVectorRV<V,M> >::Type
        marg_param_rv( new QUESO::UniformVectorRV<V,M>("marg_param_rv_", marg_domain) );

      const V & data = param_space.zeroVector();

      unsigned int int_order = 1;

      QUESO::SharedPtr<QUESO::Base1DQuadrature>::Type
        qrule_1d( new QUESO::UniformLegendre1DQuadrature(marg_min_domain_value, marg_max_domain_value,
                                                         int_order, false) );

      std::vector<QUESO::SharedPtr<QUESO::Base1DQuadrature>::Type> all_1d_qrules(marg_dim,qrule_1d);

      typename QUESO::SharedPtr<QUESO::MultiDQuadratureBase<V,M> >::Type
        marg_integration( new QUESO::TensorProductQuadrature<V,M>(marg_domain,all_1d_qrules) );

      LinearModelLikelihood<V,M> likelihood( "likelihood_test_",
                                             param_domain,
                                             data,
                                             marg_param_rv,
                                             marg_integration,
                                             false );

      // Test marginalized likelihood over a range of points in the parameter domain
      unsigned int n_intervals = 20;
      double tol = std::numeric_limits<double>::epsilon()*10;

      this->test_likelihood_values_range( n_intervals, param_min_domain_value, param_max_domain_value,
                                          param_space, tol, likelihood );
    }

    void test_linear_func_gaussian_marg_space()
    {
      // Instantiate the parameter space
      unsigned int param_dim = 1;
      unsigned int marg_dim = 2;
      QUESO::VectorSpace<V,M> param_space( (*this->_env), "param_", param_dim, NULL);
      QUESO::VectorSpace<V,M> marg_space( (*this->_env), "marg_", marg_dim, NULL);

      double param_min_domain_value = 0.0;
      double param_max_domain_value = 1.0;

      double marg_min_domain_value = -INFINITY;
      double marg_max_domain_value = INFINITY;

      typename QUESO::ScopedPtr<V>::Type param_min_values( param_space.newVector(param_min_domain_value) );
      typename QUESO::ScopedPtr<V>::Type param_max_values( param_space.newVector(param_max_domain_value) );

      typename QUESO::ScopedPtr<V>::Type marg_min_values( marg_space.newVector(marg_min_domain_value) );
      typename QUESO::ScopedPtr<V>::Type marg_max_values( marg_space.newVector(marg_max_domain_value) );

      QUESO::BoxSubset<V,M> param_domain( "param_domain_", param_space, (*param_min_values), (*param_max_values) );
      QUESO::BoxSubset<V,M> marg_domain( "marg_domain_", marg_space, (*marg_min_values), (*marg_max_values) );

      typename QUESO::ScopedPtr<V>::Type var_vec( param_space.newVector(1.0) );

      // Zero mean, unit variance Gaussian
      typename QUESO::SharedPtr<QUESO::BaseVectorRV<V,M> >::Type
        marg_param_rv( new QUESO::GaussianVectorRV<V,M>("marg_param_rv_", marg_domain,
                                                        marg_space.zeroVector(),
                                                        *var_vec) );

      const V & data = param_space.zeroVector();

      unsigned int int_order = 1;

      QUESO::SharedPtr<QUESO::Base1DQuadrature>::Type
        qrule_1d( new QUESO::GaussianHermite1DQuadrature(0.0,1.0,int_order) );

      std::vector<QUESO::SharedPtr<QUESO::Base1DQuadrature>::Type> all_1d_qrules(marg_dim,qrule_1d);

      typename QUESO::SharedPtr<QUESO::MultiDQuadratureBase<V,M> >::Type
        marg_integration( new QUESO::TensorProductQuadrature<V,M>(marg_domain,all_1d_qrules) );

      GaussModelLikelihood<V,M> likelihood( "likelihood_test_",
                                             param_domain,
                                             data,
                                             marg_param_rv,
                                             marg_integration,
                                             true );

      // Test marginalized likelihood over a range of points in the parameter domain
      unsigned int n_intervals = 20;
      double tol = std::numeric_limits<double>::epsilon()*50;

      this->test_likelihood_values_range( n_intervals, param_min_domain_value, param_max_domain_value,
                                          param_space, tol, likelihood );
    }

  protected:

    QUESO::EnvOptionsValues _options;

    typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;

    void init_env()
    {
      _env.reset( new QUESO::FullEnvironment("","",&_options) );
    }

    void test_likelihood_values_range( unsigned int n_intervals, double param_min_domain_value,
                                       double param_max_domain_value, const QUESO::VectorSpace<V,M> & param_space,
                                       double tol, const TestlingLikelihoodBase<V,M> & likelihood )
    {
      for( unsigned int i = 0; i < n_intervals; i++ )
        {
          double x = param_min_domain_value + i*(param_max_domain_value-param_min_domain_value)/(double)n_intervals;

          typename QUESO::ScopedPtr<V>::Type param_vec( param_space.newVector(x) );

          likelihood.testLikelihoodValue(*param_vec,likelihood,tol);
        }
    }
  };

  class MarginalLikelihoodGslTest : public MarginalLikelihoodTestBase<QUESO::GslVector,QUESO::GslMatrix>
  {
  public:
    CPPUNIT_TEST_SUITE( MarginalLikelihoodGslTest );

    CPPUNIT_TEST( test_linear_func_uniform_marg_space );
    CPPUNIT_TEST( test_linear_func_gaussian_marg_space );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( MarginalLikelihoodGslTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
