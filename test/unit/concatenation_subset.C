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
#include <queso/ConcatenationSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorRV.h>

namespace QUESOTesting
{

class ConcatenationSubsetTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(ConcatenationSubsetTest);
  CPPUNIT_TEST(test_contains);
  CPPUNIT_TEST(test_moments);
  CPPUNIT_TEST(test_print);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));

    space.reset(new QUESO::VectorSpace<>(*env, "", 1, NULL));

    min1.reset(new QUESO::GslVector(space->zeroVector()));
    (*min1)[0] = 1;
    max1.reset(new QUESO::GslVector(space->zeroVector()));
    (*max1)[0] = 3;

    min2.reset(new QUESO::GslVector(space->zeroVector()));
    (*min2)[0] = 2;
    max2.reset(new QUESO::GslVector(space->zeroVector()));
    (*max2)[0] = 14;

    domain1.reset(new QUESO::BoxSubset<>("", *space, *min1, *max1));
    domain2.reset(new QUESO::BoxSubset<>("", *space, *min2, *max2));

    double volume = domain1->volume() * domain2->volume();

    sets.reset(new std::vector<const QUESO::VectorSet<> *>(2));
    (*sets)[0] = domain1.get();
    (*sets)[1] = domain2.get();

    concat_subset.reset(new QUESO::ConcatenationSubset<>("", *space, volume, *sets));
  }

  void test_contains()
  {
    QUESO::VectorSpace<> big_space(*env, "", 2, NULL);
    QUESO::GslVector vec(big_space.zeroVector());

    vec[0] = 1.5;
    vec[1] = 4.1;

    CPPUNIT_ASSERT(concat_subset->contains(vec));
  }

  void test_moments()
  {
    QUESO::VectorSpace<> big_space(*env, "", 2, NULL);
    QUESO::GslMatrix moments(big_space.zeroVector());

    concat_subset->moments(moments);

    CPPUNIT_ASSERT_EQUAL(moments(0, 0), 2.0/3.0);
    CPPUNIT_ASSERT_EQUAL(moments(1, 1), 144.0);
  }

  void test_print()
  {
    concat_subset->print(std::cout);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::VectorSpace<> >::Type space;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type min1;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type max1;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type min2;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type max2;
  typename QUESO::ScopedPtr<QUESO::BoxSubset<> >::Type domain1;
  typename QUESO::ScopedPtr<QUESO::BoxSubset<> >::Type domain2;
  typename QUESO::ScopedPtr<std::vector<const QUESO::VectorSet<> *> >::Type sets;
  typename QUESO::ScopedPtr<QUESO::ConcatenationSubset<> >::Type concat_subset;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ConcatenationSubsetTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
