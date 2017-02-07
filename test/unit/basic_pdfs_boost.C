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

#include <queso/BasicPdfsBoost.h>
#include <cmath>

namespace QUESOTesting
{

class BasicPdfsBoostTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(BasicPdfsBoostTest);
  CPPUNIT_TEST(test_beta);
  CPPUNIT_TEST(test_gamma);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    _env.reset(new QUESO::FullEnvironment("","",NULL));
  }

  void test_beta()
  {
    QUESO::BasicPdfsBoost pdf_boost(0);

    double value = pdf_boost.betaPdfActualValue(0.5, 2.0, 2.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, value, 1e-14);
  }

  void test_gamma()
  {
    QUESO::BasicPdfsBoost pdf_boost(0);

    double value = pdf_boost.gammaPdfActualValue(1.5, 4.0, 0.5);
    double actualValue = 1.5 * 1.5 * 1.5 * std::exp(-3.0) / (3.0 * 0.5 * 0.5 * 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualValue, value, 1e-14);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type _env;
};

CPPUNIT_TEST_SUITE_REGISTRATION(BasicPdfsBoostTest);

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
