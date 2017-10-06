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

#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/ScopedPtr.h>
#include <queso/SimulationOutputPoint.h>
#include <queso/TensorProductMesh.h>

namespace
{
  const int n_elem      = 5;
  const int n_tests     = 7;
  const double base_len = 1.5;

  const std::size_t n_points    = n_elem + 1;
}

namespace QUESOTesting
{
  class TensorProductMeshTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( TensorProductMeshTest );
    CPPUNIT_TEST( test_1D );
    CPPUNIT_TEST_SUITE_END();

  public:
    void setUp()
    {
      env.reset(new QUESO::FullEnvironment("","",NULL));
    }

    void test_1D()
    {
      QUESO::Map map(n_points, 0, env->fullComm());

      QUESO::GslVector coefs(*env, map);

      QUESO::TensorProductMesh<> mesh;

      // Test coords going out of scope
      {
        std::vector<double> coords;
        const double h = 1.0 / n_elem;
        for (unsigned int i=0; i != n_points; ++i)
          {
            const double frac = h * i;
            coords.push_back(base_len * frac * frac);
          }

        mesh.set_t_coordinates(coords);

        CPPUNIT_ASSERT_EQUAL(mesh.n_outputs(), n_points);

        // Linear solution, so we expect exact values
        for (unsigned int i=0; i != n_points; ++i)
          coefs[i] = 3.0 - 2.0 * coords[i];
      }

      for (unsigned int i=0; i != n_tests; ++i)
        {
          const double t = base_len / n_tests * i;
          QUESO::SimulationOutputPoint p;
          p.t() = t;

          const double mesh_u = mesh.interpolateOutput(coefs, p);
          const double true_u = 3.0 - 2.0 * t;
          CPPUNIT_ASSERT_DOUBLES_EQUAL(mesh_u, true_u, 1.e-12);
        }
    }

  private:
    typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( TensorProductMeshTest );
} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
