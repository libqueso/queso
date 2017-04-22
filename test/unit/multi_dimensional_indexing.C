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

#include <queso/MultiDimensionalIndexing.h>

namespace QUESOTesting
{
  class MultiDimensionalIndexingTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( MultiDimensionalIndexingTest );

    CPPUNIT_TEST( test_local_to_global );
    CPPUNIT_TEST( test_global_to_local );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      _n_points.resize(5);
      _n_points[0] = 11;
      _n_points[1] = 21;
      _n_points[2] = 31;
      _n_points[3] = 41;
      _n_points[4] = 51;

      _indices.resize(5);
      _indices[0] = 1;
      _indices[1] = 2;
      _indices[2] = 3;
      _indices[3] = 4;
      _indices[4] = 5;

      _global_exact =
        _indices[0] +
        _indices[1]*_n_points[0] +
        _indices[2]*_n_points[0]*_n_points[1] +
        _indices[3]*_n_points[0]*_n_points[1]*_n_points[2] +
        _indices[4]*_n_points[0]*_n_points[1]*_n_points[2]*_n_points[3];
    }

    void test_local_to_global()
    {
      unsigned int global =
        QUESO::MultiDimensionalIndexing::coordToGlobal( _indices, _n_points );

      CPPUNIT_ASSERT_EQUAL(_global_exact,global);
    }

    void test_global_to_local()
    {
      std::vector<unsigned int> indices_test;
      QUESO::MultiDimensionalIndexing::globalToCoord( _global_exact,
                                                      _n_points,
                                                      indices_test );

      for( unsigned int d = 0; d < 5; d++ )
        CPPUNIT_ASSERT_EQUAL(_indices[d],indices_test[d]);
    }

  private:

    std::vector<unsigned int> _n_points;
    std::vector<unsigned int> _indices;
    unsigned int _global_exact;
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( MultiDimensionalIndexingTest );

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
