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
#include <convergence_rate_helper.h>

#include <queso/MonteCarloQuadrature.h>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <cmath>
#include <limits>

namespace QUESOTesting
{
  template <class V, class M>
  class MonteCarloQuadratureHypersphereVolumeTestBase : public MultiDQuadratureFunction<V,M>,
                                                        public QuadratureMultiDTestBase<V,M>,
                                                        public ConvergenceRateHelper
  {
  public:

    void setUp()
    {
      this->init_env();
    }

    //! Heavyside function on hypersphere
    virtual double f( const V & x ) const
    {
      CPPUNIT_ASSERT_EQUAL(x.sizeGlobal(),this->_dim);

      double value = 0.0;

      double radius = 0.0;
      for( unsigned int i = 0; i < this->_dim; i++ )
        radius += x[i]*x[i];

      if( radius <= 1.0 )
        value = 1.0;

      return value;
    }

    virtual double int_f( const QUESO::BoxSubset<V,M> & /*domain*/ ) const
    {
      return this->unit_hypersphere_volume(this->_dim);
    }

    void test_mc_2d_hyperball()
    {
      this->mc_hyperball_volume_test_impl(2);
    }

    void test_mc_3d_hyperball()
    {
      this->mc_hyperball_volume_test_impl(3);
    }

    void test_mc_4d_hyperball()
    {
      this->mc_hyperball_volume_test_impl(4);
    }

    void test_mc_5d_hyperball()
    {
      this->mc_hyperball_volume_test_impl(5);
    }

  protected:

    unsigned int _dim;

    double unit_hypersphere_volume( unsigned int dim ) const
    {
      double pi = 3.14159265358979323846;
      double value = 0.0;
      switch(dim)
        {
        case(0):
          {
            value = 1.0;
          }
          break;

        case(1):
          {
            value = 2.0;
          }
          break;

        case(2):
          {
            value = pi;
          }
          break;

        default:
          {
            value = 2.0*pi/(double)(dim)*this->unit_hypersphere_volume(dim-2);
          }
          break;
        }

      return value;
    }

    void mc_hyperball_volume_test_impl( unsigned int dim )
    {
      this->_dim = dim;

      // Instantiate the parameter space
      QUESO::VectorSpace<V,M> param_space( (*this->_env), "param_", dim, NULL);

      typename QUESO::ScopedPtr<V>::Type min_values( param_space.newVector(-1.0) );
      typename QUESO::ScopedPtr<V>::Type max_values( param_space.newVector(1.0) );

      QUESO::BoxSubset<V,M> param_domain( "param_domain_", param_space, (*min_values), (*max_values) );

      std::vector<unsigned int> samples(5);
      samples[0] = 1e1;
      samples[1] = 1e2;
      samples[2] = 1e3;
      samples[3] = 1e4;
      samples[4] = 1e5;

      QUESO::MonteCarloQuadrature<V,M> quad_rule( param_domain, samples.back() );

      std::vector<double> log_samples(samples.size());

      std::vector<double> int_error(samples.size());

      double exact = this->int_f(param_domain);

      for( unsigned int i = 0; i < samples.size(); i++ )
        {
          // We pregenerate the samples/weights, but the weights will have V/N
          // where N is the number of requested samples at construction time. But
          // we're using fewer points than N to compute the integral, so we need
          // to rescale for this test.
          double factor = (double)(samples.back())/(double)(samples[i]);

          double value = factor*this->integrate_func( (*this), quad_rule, samples[i] );

          double rel_error = std::abs( (value - exact)/exact );

          int_error[i] = std::log10( rel_error );
          log_samples[i] = std::log10( samples[i] );
        }

      double rate = this->compute_convergence_rate( log_samples, int_error, 0, samples.size()-1 );

      // Monte Carlo convergence is very noisy for this problem, so we
      // can't really meaningfully check for a rate = -0.5.
      // But we can at least make sure the rate is negative as a general trend.
      CPPUNIT_ASSERT( rate < 0.0 );
    }

  };


  class MonteCarloQuadratureHypersphereVolumeGslTest :
    public MonteCarloQuadratureHypersphereVolumeTestBase<QUESO::GslVector,QUESO::GslMatrix>
  {
  public:

    CPPUNIT_TEST_SUITE( MonteCarloQuadratureHypersphereVolumeGslTest );

    CPPUNIT_TEST(test_mc_2d_hyperball);
    CPPUNIT_TEST(test_mc_3d_hyperball);
    CPPUNIT_TEST(test_mc_4d_hyperball);
    CPPUNIT_TEST(test_mc_5d_hyperball);

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( MonteCarloQuadratureHypersphereVolumeGslTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
