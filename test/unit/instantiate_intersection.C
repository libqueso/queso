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
#include <queso/VectorSpace.h>
#include <queso/InstantiateIntersection.h>
#include <queso/GslVector.h>
#include <queso/BoxSubset.h>
#include <queso/IntersectionSubset.h>

namespace QUESOTesting
{

class InstantiateIntersectionTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(InstantiateIntersectionTest);
  CPPUNIT_TEST(test_intersect1);
  CPPUNIT_TEST(test_intersect2);
  CPPUNIT_TEST(test_intersect3);
  CPPUNIT_TEST(test_intersect4);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));
  }

  void test_intersect1()
  {
    QUESO::VectorSpace<> space1(*env, "", 1, NULL);
    QUESO::VectorSpace<> space2(*env, "", 2, NULL);

    QUESO::ScopedPtr<QUESO::VectorSet<> >::Type vectorSetPtr;
    vectorSetPtr.reset(QUESO::InstantiateIntersection<>(space1, space2));

    unsigned int dim = 1;
    CPPUNIT_ASSERT_EQUAL(vectorSetPtr->vectorSpace().dimGlobal(), dim);
  }

  void test_intersect2()
  {
    QUESO::VectorSpace<> space1(*env, "", 1, NULL);
    QUESO::VectorSpace<> space2(*env, "", 1, NULL);

    QUESO::ScopedPtr<QUESO::VectorSet<> >::Type vectorSetPtr;
    vectorSetPtr.reset(QUESO::InstantiateIntersection<>(space1, space2));

    unsigned int dim = 1;
    CPPUNIT_ASSERT_EQUAL(vectorSetPtr->vectorSpace().dimGlobal(), dim);
  }

  void test_intersect3()
  {
    QUESO::VectorSpace<> space1(*env, "", 2, NULL);
    QUESO::VectorSpace<> space2(*env, "", 1, NULL);

    QUESO::ScopedPtr<QUESO::VectorSet<> >::Type vectorSetPtr;
    vectorSetPtr.reset(QUESO::InstantiateIntersection<>(space1, space2));

    unsigned int dim = 1;
    CPPUNIT_ASSERT_EQUAL(vectorSetPtr->vectorSpace().dimGlobal(), dim);
  }

  void test_intersect4()
  {
    QUESO::VectorSpace<> space1(*env, "", 1, NULL);
    QUESO::GslVector min1(space1.zeroVector());
    min1[0] = -2;
    QUESO::GslVector max1(space1.zeroVector());
    max1[0] = 2;
    QUESO::BoxSubset<> box1("", space1, min1, max1);

    QUESO::VectorSpace<> space2(*env, "", 1, NULL);

    QUESO::ScopedPtr<QUESO::VectorSet<> >::Type vectorSetPtr;
    vectorSetPtr.reset(QUESO::InstantiateIntersection<>(space2, box1));

    // It should be a box
    QUESO::BoxSubset<> & box2 = dynamic_cast<QUESO::BoxSubset<> &>(*vectorSetPtr);

    CPPUNIT_ASSERT_EQUAL(box1.volume(), box2.volume());

    vectorSetPtr.reset(QUESO::InstantiateIntersection<>(box1, space2));

    // It should be a box
    QUESO::BoxSubset<> & box3 = dynamic_cast<QUESO::BoxSubset<> &>(*vectorSetPtr);

    CPPUNIT_ASSERT_EQUAL(box1.volume(), box3.volume());
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
};

CPPUNIT_TEST_SUITE_REGISTRATION(InstantiateIntersectionTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
