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
#include <queso/GslVector.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/ConstantVectorFunction.h>

namespace QUESOTesting
{

class ConstantVectorFunctionTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(ConstantVectorFunctionTest);
  CPPUNIT_TEST(test_compute);
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

    domain.reset(new QUESO::BoxSubset<>("", *space, *min1, *max1));
    image.reset(new QUESO::BoxSubset<>("", *space, *min1, *max1));

    val.reset(new QUESO::GslVector(space->zeroVector()));
    (*val)[0] = 2.0;

    fn.reset(new QUESO::ConstantVectorFunction<>("", *domain, *image, *val));
  }

  void test_compute()
  {
    QUESO::GslVector point(space->zeroVector());
    point[0] = 1.5;

    QUESO::GslVector result(space->zeroVector());
    fn->compute(point, NULL, result, NULL, NULL, NULL);

    CPPUNIT_ASSERT_EQUAL(result[0], (*val)[0]);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::VectorSpace<> >::Type space;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type min1;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type max1;
  typename QUESO::ScopedPtr<QUESO::BoxSubset<> >::Type domain;
  typename QUESO::ScopedPtr<QUESO::BoxSubset<> >::Type image;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type val;
  typename QUESO::ScopedPtr<QUESO::ConstantVectorFunction<> >::Type fn;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ConstantVectorFunctionTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
