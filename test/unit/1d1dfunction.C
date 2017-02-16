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

#include <queso/1D1DFunction.h>

namespace QUESOTesting
{

double value(double domainValue, const void * routinesDataPtr)
{
  return 2.0;
}

double deriv(double domainValue, const void * routinesDataPtr)
{
  return 0.0;
}

class Test1D1DFunction : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(Test1D1DFunction);
  CPPUNIT_TEST(test_constant);
  CPPUNIT_TEST(test_generic);
  CPPUNIT_TEST(test_linear);
  CPPUNIT_TEST(test_piecewise_linear);
  CPPUNIT_TEST(test_quadratic);
  CPPUNIT_TEST(test_sampled);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
  }

  void test_constant()
  {
    QUESO::Constant1D1DFunction constant(-1.0, 1.0, 2.0);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, constant.minDomainValue(), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, constant.maxDomainValue(), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, constant.value(0.0), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, constant.deriv(0.0), 1e-14);
  }

  void test_generic()
  {
    QUESO::Generic1D1DFunction generic(-1.0, 1.0, value, deriv, NULL);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, generic.value(0.0), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, generic.deriv(0.0), 1e-14);
  }

  void test_linear()
  {
    QUESO::Linear1D1DFunction linear(-1.0, 1.0, 0.0, 1.0, 2.0);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, linear.value(0.0), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, linear.deriv(0.0), 1e-14);
  }

  void test_piecewise_linear()
  {
    unsigned int num_pieces = 3;

    std::vector<double> ref_domain_vals(num_pieces);
    ref_domain_vals[0] = 0.0;
    ref_domain_vals[1] = 1.0;
    ref_domain_vals[2] = 2.0;

    std::vector<double> rates(num_pieces);
    rates[0] = 1.0;
    rates[1] = 2.0;
    rates[2] = 3.0;

    QUESO::PiecewiseLinear1D1DFunction piecewise_linear(0.0, 3.0,
        ref_domain_vals, 3.0, rates);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5, piecewise_linear.value(0.5), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, piecewise_linear.value(1.5), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.5, piecewise_linear.value(2.5), 1e-14);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, piecewise_linear.deriv(0.5), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, piecewise_linear.deriv(1.5), 1e-14);
  }

  void test_quadratic()
  {
    // x(1 - x)
    QUESO::Quadratic1D1DFunction quadratic(0.0, 1.0, -1.0, 1.0, 0.0);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, quadratic.value(0.0), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, quadratic.value(0.5), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, quadratic.value(1.0), 1e-14);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, quadratic.deriv(0.0), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, quadratic.deriv(0.5), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, quadratic.deriv(1.0), 1e-14);
  }

  void test_sampled()
  {
    std::vector<double> domain_vals;
    domain_vals.push_back(0.0);
    domain_vals.push_back(1.0);

    std::vector<double> image_vals;
    image_vals.push_back(1.0);
    image_vals.push_back(2.0);

    QUESO::Sampled1D1DFunction sampled(domain_vals, image_vals);

    domain_vals.push_back(2.0);
    image_vals.push_back(3.0);

    sampled.set(domain_vals, image_vals);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, sampled.value(0.0), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, sampled.value(0.5), 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5, sampled.value(1.5), 1e-14);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, sampled.domainValues()[0], 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, sampled.domainValues()[1], 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, sampled.domainValues()[2], 1e-14);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, sampled.imageValues()[0], 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, sampled.imageValues()[1], 1e-14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, sampled.imageValues()[2], 1e-14);

    CPPUNIT_ASSERT_EQUAL(true, sampled.domainValueMatchesExactly(1.0));

    // I don't care about this.  I just print to hit lines in the gcov output.
    QUESO::FullEnvironment env("", "", NULL);
    std::ofstream out_stream("/dev/null");
    sampled.printForMatlab(env, out_stream, "");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(Test1D1DFunction);

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
