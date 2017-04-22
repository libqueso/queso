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
#include <queso/ScopedPtr.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>

namespace QUESOTesting
{

class VectorSpaceTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(VectorSpaceTest);
  CPPUNIT_TEST(test_centroid);
  CPPUNIT_TEST(test_moments);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));
    space.reset(new QUESO::VectorSpace<>(*env, "", 2, NULL));
  }

  void test_centroid()
  {
    QUESO::GslVector centroid(space->zeroVector());
    space->centroid(centroid);

    // Is this peligroso?
    CPPUNIT_ASSERT_EQUAL(INFINITY, (float)centroid[0]);
    CPPUNIT_ASSERT_EQUAL(INFINITY, (float)centroid[1]);
  }

  void test_moments()
  {
    QUESO::GslMatrix moments(space->zeroVector());
    space->moments(moments);

    // Is this peligroso?
    CPPUNIT_ASSERT_EQUAL(INFINITY, (float)moments(0,0));
    CPPUNIT_ASSERT_EQUAL(0.0f, (float)moments(0,1));
    CPPUNIT_ASSERT_EQUAL(0.0f, (float)moments(1,0));
    CPPUNIT_ASSERT_EQUAL(INFINITY, (float)moments(1,1));
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> >::Type space;
};

CPPUNIT_TEST_SUITE_REGISTRATION(VectorSpaceTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
