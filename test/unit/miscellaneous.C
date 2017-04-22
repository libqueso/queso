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

#include <queso/Miscellaneous.h>
#include <queso/RngGsl.h>

namespace QUESOTesting
{

class MiscellaneousTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(MiscellaneousTest);
  CPPUNIT_TEST(test_words_from_string);
  CPPUNIT_TEST(test_string_and_double_from_file);
  CPPUNIT_TEST(test_chars_and_double_from_file);
  CPPUNIT_TEST(test_debug_msgs);
  CPPUNIT_TEST(test_hamming);
  CPPUNIT_TEST(test_gammar);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
  }

  void test_words_from_string()
  {
    std::string inputString("these are some words");
    std::vector<std::string> words;

    QUESO::MiscReadWordsFromString(inputString, words);

    CPPUNIT_ASSERT_EQUAL(std::string("these"), words[0]);
    CPPUNIT_ASSERT_EQUAL(std::string("are"), words[1]);
    CPPUNIT_ASSERT_EQUAL(std::string("some"), words[2]);
    CPPUNIT_ASSERT_EQUAL(std::string("words"), words[3]);
  }

  void test_string_and_double_from_file()
  {
    std::string fname("test_string_and_double_from_file");
    std::ofstream ofile(fname);

    // Write some data
    double actualValue = 5.32;
    if (ofile.good()) {
      ofile << actualValue << std::endl;
    }
    ofile.close();

    // Read that data
    std::ifstream ifile(fname);
    std::string output;
    double value;
    QUESO::MiscReadStringAndDoubleFromFile(ifile, output, &value);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualValue, value, 1e-14);

    int returnval = std::remove(fname.c_str());
    CPPUNIT_ASSERT_EQUAL(0, returnval);
  }

  void test_chars_and_double_from_file()
  {
    std::string fname("test_string_and_double_from_file");
    std::ofstream ofile(fname);

    // Write some data
    double actualValue = 5.32;
    if (ofile.good()) {
      ofile << actualValue << std::endl;
    }
    ofile.close();

    // Read that data
    std::ifstream ifile(fname);
    std::string output;
    double value;
    bool eol;
    QUESO::MiscReadCharsAndDoubleFromFile(ifile, output, &value, eol);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualValue, value, 1e-14);

    int returnval = std::remove(fname.c_str());
    CPPUNIT_ASSERT_EQUAL(0, returnval);
  }

  void test_debug_msgs()
  {
    unsigned int val1 = 1;
    int val2 = -1;
    double val3 = -1.5;

    QUESO::MiscUintDebugMessage(val1, "uint debug msg");
    QUESO::MiscIntDebugMessage(val2, "int debug msg");
    QUESO::MiscDoubleDebugMessage(val3, "debug msg");
  }

  void test_hamming()
  {
    double value = 0.53836 - 0.46164;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(value, QUESO::MiscHammingWindow(1, 0), 1e-14);
  }

  void test_gammar()
  {
    // I don't know what MiscGammar is supposed to do, so this is a regression
    // test
    QUESO::RngGsl rng_gsl(1, 0);

    double value = QUESO::MiscGammar(0.5, 1.0, &rng_gsl);

    double actualValue = 0.06164141389893609824;  // Regression solution
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualValue, value, 1e-14);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(MiscellaneousTest);

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
