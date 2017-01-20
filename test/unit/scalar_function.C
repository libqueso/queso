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
#include <queso/BoxSubset.h>
#include <queso/VectorSet.h>
#include <queso/ScalarFunction.h>

namespace QUESOTesting
{

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet)
  {
  }

  virtual double lnValue(const V & domainVector) const
  {
    return domainVector[0] * domainVector[0];
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, NULL, NULL, NULL, NULL));
  }

  using QUESO::BaseScalarFunction<V, M>::lnValue;
};

class BaseScalarFunctionTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(BaseScalarFunctionTest);
  CPPUNIT_TEST(test_fd);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));
  }

  void test_fd()
  {
    QUESO::VectorSpace<> space(*env, "", 1, NULL);

    QUESO::GslVector min(space.zeroVector());
    min[0] = -INFINITY;
    QUESO::GslVector max(space.zeroVector());
    max[0] = INFINITY;

    QUESO::BoxSubset<> domain("", space, min, max);

    Likelihood<> lhood("", domain);

    QUESO::GslVector point(space.zeroVector());
    point[0] = 10.0;

    QUESO::GslVector grad1(space.zeroVector());
    QUESO::GslVector grad2(space.zeroVector());

    lhood.setFiniteDifferenceStepSize(1e-4);
    lhood.lnValue(point, grad1);

    lhood.setFiniteDifferenceStepSize(0, 1e-4);
    lhood.lnValue(point, grad2);

    // So the gradients should be the same
    CPPUNIT_ASSERT_DOUBLES_EQUAL(grad1[0], grad2[0], 1e-4);

    // And it should be 20.0
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20.0, grad1[0], 1e-4);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
};

CPPUNIT_TEST_SUITE_REGISTRATION(BaseScalarFunctionTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
