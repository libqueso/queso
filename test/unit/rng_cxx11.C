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

#ifdef QUESO_HAVE_CXX11
#include <queso/RngCXX11.h>

namespace QUESOTesting
{

class RngCXX11Test : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(RngCXX11Test);
  CPPUNIT_TEST(test_beta);
  CPPUNIT_TEST(test_gamma);
  CPPUNIT_TEST(test_uniform);
  CPPUNIT_TEST(test_gaussian);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    _env.reset(new QUESO::FullEnvironment("","",NULL));
  }

  void test_beta()
  {
    QUESO::RngCXX11 rng_cxx11(0, 0);
    rng_cxx11.resetSeed(0);

    double mean = 0.0;
    unsigned int num_samples = 1000000;

    for (unsigned int i = 0; i < num_samples; i++) {
      mean += rng_cxx11.betaSample(2.0, 2.0) / num_samples;
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, mean, 1e-2);
  }

  void test_gamma()
  {
    QUESO::RngCXX11 rng_cxx11(0, 0);

    double mean = 0.0;
    unsigned int num_samples = 1000000;

    for (unsigned int i = 0; i < num_samples; i++) {
      mean += rng_cxx11.gammaSample(10.0, 0.5) / num_samples;
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, mean, 1e-2);
  }

  void test_uniform()
  {
    QUESO::RngCXX11 rng_cxx11(0, 0);

    double mean = 0.0;
    unsigned int num_samples = 1000000;

    for (unsigned int i = 0; i < num_samples; i++) {
      mean += rng_cxx11.uniformSample() / num_samples;
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, mean, 1e-2);
  }

  void test_gaussian()
  {
    QUESO::RngCXX11 rng_cxx11(0, 0);

    double mean = 0.0;
    unsigned int num_samples = 1000000;

    for (unsigned int i = 0; i < num_samples; i++) {
      mean += rng_cxx11.gaussianSample(0.1) / num_samples;
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, mean, 1e-2);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;
};

CPPUNIT_TEST_SUITE_REGISTRATION(RngCXX11Test);

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CXX11
#endif // QUESO_HAVE_CPPUNIT
