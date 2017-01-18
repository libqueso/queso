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

#include <cmath>
#include <limits>

namespace QUESOTesting
{
  //! Provides helper functions for fitting data to a linear curve.
  /*! Intended use case is for fitting data from convergence tests to ascertain
      whether or not the correct convergence rate was achieved. x is the indepdent
      variable and y is the dependent variable.
   */
  class ConvergenceRateHelper
  {
  public:

    //! Check the convergence rate of the supplied data to within the given tolerance
    /*! Relative error is checked */
    void check_convergence_rate( const std::vector<double> & x,
                                 const std::vector<double> & y,
                                 unsigned int start_index,
                                 unsigned int end_index,
                                 double expected_rate,
                                 double tol )
    {
      CPPUNIT_ASSERT(tol > 0.0);

      double computed_rate = this->compute_convergence_rate( x, y, start_index, end_index );

      double rel_error = std::abs( (computed_rate - expected_rate)/expected_rate);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,rel_error,tol);
    }

    double compute_convergence_rate( const std::vector<double> & x,
                                     const std::vector<double> & y,
                                     unsigned int start_index,
                                     unsigned int end_index )
    {
      CPPUNIT_ASSERT_EQUAL(x.size(), y.size());
      CPPUNIT_ASSERT(start_index < end_index);
      CPPUNIT_ASSERT(end_index < x.size());

      return this->slope(x,y,start_index,end_index);
    }

  private:

    //! Slope \f$ m \f$ for linear regression
    double slope( const std::vector<double> & x,
                  const std::vector<double> & y,
                  unsigned int start_index,
                  unsigned int end_index )
    {
      unsigned int n = x.size();

      double sum_y, sum_x2, sum_x, sum_xy;

      this->compute_fit_quantities( x, y, sum_y, sum_x, sum_x2, sum_xy, start_index, end_index );

      return ( n*sum_xy - sum_x*sum_y )/( n*sum_x2 - (sum_x*sum_x) );
    }

    //! y-intercept \f$ b \f$ for linear regression
    double intercept( const std::vector<double> & x,
                      const std::vector<double> & y,
                      unsigned int start_index,
                      unsigned int end_index )
    {
      unsigned int n = x.size();

      double sum_y, sum_x2, sum_x, sum_xy;

      this->compute_fit_quantities( x, y, sum_y, sum_x, sum_x2, sum_xy, start_index, end_index );

      return ( sum_y*sum_x2 - sum_x*sum_xy )/( n*sum_x2 - (sum_x*sum_x) );
    }

    void compute_fit_quantities( const std::vector<double> & x,
                                 const std::vector<double> & y,
                                 double & sum_y, double & sum_x, double & sum_x2, double & sum_xy,
                                 unsigned int start_index, unsigned int end_index )
    {
      CPPUNIT_ASSERT_EQUAL(x.size(),y.size());
      CPPUNIT_ASSERT(start_index < end_index);

      unsigned int n = x.size();
      CPPUNIT_ASSERT(end_index < n);

      sum_y = 0.0;
      sum_x2 = 0.0;
      sum_x = 0.0;
      sum_xy = 0.0;

      for (unsigned int i=start_index; i<=end_index; i++)
        {
          double xi = x[i];
          double yi = y[i];

          sum_y  += yi;
          sum_x2 += xi*xi;
          sum_x  += xi;
          sum_xy += xi*yi;
        }
    }

  };
}

#endif // QUESO_HAVE_CPPUNIT
