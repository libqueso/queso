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
#include <queso/ScalarSequence.h>

namespace QUESOTesting
{

class ScalarSequenceTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(ScalarSequenceTest);
  CPPUNIT_TEST(test_setup);
  CPPUNIT_TEST(test_erase);
  CPPUNIT_TEST(test_sub_min_plain);
  CPPUNIT_TEST(test_unified_min_plain);
  CPPUNIT_TEST(test_sub_mean_plain);
  CPPUNIT_TEST(test_sub_median_plain);
  CPPUNIT_TEST(test_unified_median_plain);
  CPPUNIT_TEST(test_sub_sample_variance_plain);
  CPPUNIT_TEST(test_unified_sample_variance_plain);
  CPPUNIT_TEST(test_set_gaussian);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));
    sequence.reset(new QUESO::ScalarSequence<double>(*env, 13, "my_sequence"));

    // Set *sequence to (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    for (unsigned int i = 0; i < 13; i++) {
      (*sequence)[i] = i;
    }
  }

  void test_setup()
  {
    // Some dumb tests just to get coverage up
    const QUESO::BaseEnvironment * env_ptr = env.get();
    CPPUNIT_ASSERT_EQUAL(env_ptr, &(sequence->env()));

    CPPUNIT_ASSERT_EQUAL(std::string("my_sequence"), sequence->name());
  }

  void test_erase()
  {
    unsigned int expectedSize = 13;
    CPPUNIT_ASSERT_EQUAL(expectedSize, sequence->subSequenceSize());

    // Remove (0-based) elements 2 to 7 (not including 7)
    sequence->erasePositions(2, 5);

    // *sequence should now be (0, 1, 7, 8, 9, 10, 11, 12)
    for (unsigned int i = 0; i < 8; i++) {
      double expected;
      if (i < 2) {
        expected = i;
      }
      else {
        expected = i + 5;
      }
      CPPUNIT_ASSERT_EQUAL(expected, (*sequence)[i]);
    }
  }

  void test_sub_min_plain()
  {
    double subMinPlain = sequence->subMinPlain();

    CPPUNIT_ASSERT_EQUAL(0.0, subMinPlain);
  }

  void test_unified_min_plain()
  {
    double unifiedMinPlain = sequence->unifiedMinPlain(false);

    CPPUNIT_ASSERT_EQUAL(0.0, unifiedMinPlain);
  }

  void test_sub_mean_plain()
  {
    double actualMean = (1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12) / 13.0;

    double subMeanPlain = sequence->subMeanPlain();

    CPPUNIT_ASSERT_EQUAL(actualMean, subMeanPlain);
  }

  void test_sub_median_plain()
  {
    double actualMedian = 6;
    double subMedianPlain = sequence->subMedianPlain();

    CPPUNIT_ASSERT_EQUAL(actualMedian, subMedianPlain);
  }

  void test_unified_median_plain()
  {
    double actualMedian = 6;
    double unifiedMedianPlain = sequence->unifiedMedianPlain(false);

    CPPUNIT_ASSERT_EQUAL(actualMedian, unifiedMedianPlain);
  }

  void test_sub_sample_variance_plain()
  {
    double actualMean = 6.0;

    double runningSumSq = 0.0;
    for (unsigned int i = 0; i < 13; i++) {
      double diff = (*sequence)[i] - actualMean;
      runningSumSq += diff * diff;
    }

    double actualVar = runningSumSq / 12.0;  // Sample variance, not population

    double subVariancePlain = sequence->subSampleVariancePlain();

    CPPUNIT_ASSERT_EQUAL(actualVar, subVariancePlain);
  }

  void test_unified_sample_variance_plain()
  {
    double actualMean = 6.0;

    double runningSumSq = 0.0;
    for (unsigned int i = 0; i < 13; i++) {
      double diff = (*sequence)[i] - actualMean;
      runningSumSq += diff * diff;
    }

    double actualVar = runningSumSq / 12.0;  // Sample variance, not population

    double unifiedVariancePlain = sequence->unifiedSampleVariancePlain(false);

    CPPUNIT_ASSERT_EQUAL(actualVar, unifiedVariancePlain);
  }

  void test_set_gaussian()
  {
    unsigned int N = 1000000;
    QUESO::ScalarSequence<double> gaussian_sequence(*env, N, "");
    gaussian_sequence.setGaussian(5.0, 0.1);

    // Check the mean is close to 5.0
    double mean = gaussian_sequence.subMeanPlain();

    // Just use 0.001 delta to account for statistical error.  Not rigorous.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, mean, 0.001);

    // Check variance is close to 1.0
    double variance = gaussian_sequence.subSampleVariancePlain();

    // Just use 0.01 delta to account for statistical error.  Not rigorous.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1 * 0.1, variance, 0.001);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::ScalarSequence<double> >::Type sequence;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ScalarSequenceTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
