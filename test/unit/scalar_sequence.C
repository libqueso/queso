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
  CPPUNIT_TEST(test_set_uniform);
  CPPUNIT_TEST(test_sub_uniformly_sampled_cdf);
  CPPUNIT_TEST(test_unified_uniformly_sampled_cdf);
  CPPUNIT_TEST(test_compare_hists);
  CPPUNIT_TEST(test_sub_basic_cdf);
  CPPUNIT_TEST(test_population_variance);
  CPPUNIT_TEST(test_auto_covariance);
  CPPUNIT_TEST(test_sort);
  CPPUNIT_TEST(test_interquantile_range);
  CPPUNIT_TEST(test_sub_weight_cdf);
  CPPUNIT_TEST(test_scale_kde);
  CPPUNIT_TEST(test_gaussian_kde);
  CPPUNIT_TEST(test_positions_of_maximum);
  CPPUNIT_TEST(test_read);
  CPPUNIT_TEST(test_brooks_gelman);
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
    double subStd = sequence->subSampleStd(0, sequence->subSequenceSize(), actualMean);

    CPPUNIT_ASSERT_EQUAL(actualVar, subVariancePlain);
    CPPUNIT_ASSERT_EQUAL(std::sqrt(actualVar), subStd);
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
    double unifiedStd = sequence->unifiedSampleStd(false, 0, sequence->subSequenceSize(), actualMean);

    CPPUNIT_ASSERT_EQUAL(actualVar, unifiedVariancePlain);
    CPPUNIT_ASSERT_EQUAL(std::sqrt(actualVar), unifiedStd);
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

  void test_set_uniform()
  {
    unsigned int N = 1000000;
    QUESO::ScalarSequence<double> uniform_sequence(*env, N, "");

    double a = 3.0;
    double b = 7.0;
    uniform_sequence.setUniform(a, b);

    double mean = uniform_sequence.subMeanPlain();

    // Check mean is close to 5.0
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, mean, 0.01);
  }

  void test_sub_uniformly_sampled_cdf()
  {
    double min, max;
    std::vector<double> cdf;
    sequence->subUniformlySampledCdf(14, min, max, cdf);

    for (unsigned int i = 0; i < cdf.size(); i++) {
      unsigned int cdf_val = cdf[i] * 13;
      CPPUNIT_ASSERT_EQUAL(i, cdf_val);  // Expect linearly increasing ints
    }
  }

  void test_unified_uniformly_sampled_cdf()
  {
    double min, max;
    std::vector<double> cdf;
    sequence->unifiedUniformlySampledCdf(false, 14, min, max, cdf);

    for (unsigned int i = 0; i < cdf.size(); i++) {
      unsigned int cdf_val = cdf[i] * 13;
      CPPUNIT_ASSERT_EQUAL(i, cdf_val);  // Expect linearly increasing ints
    }
  }

  void test_compare_hists()
  {
    unsigned int num_points = 14;
    std::vector<unsigned int> bins1(num_points, 0);
    std::vector<unsigned int> bins2(num_points, 0);
    std::vector<double> centers(num_points, 0);

    QUESO::ScopedPtr<QUESO::UniformOneDGrid<double> >::Type grid;
    QUESO::UniformOneDGrid<double> * grid_rawptr = grid.get();

    // Compute histogram from subBasicHistogram
    double min, max;
    sequence->subMinMaxExtra(0, sequence->subSequenceSize(), min, max);
    sequence->subBasicHistogram(0, min, max, grid_rawptr, bins1);

    // Compute histogram from subHistogram
    sequence->subHistogram(0, min, max, centers, bins2);

    for (unsigned int i = 0; i < bins1.size(); i++) {
      CPPUNIT_ASSERT_EQUAL(bins1[i], bins2[i]);
    }

    // Test the unified version
    sequence->unifiedHistogram(false, 0, min, max, centers, bins2);
    for (unsigned int i = 0; i < bins1.size(); i++) {
      CPPUNIT_ASSERT_EQUAL(bins1[i], bins2[i]);
    }

    // Now test the weighted versions.  Let's overwrite bins2
    grid.reset(0);
    grid_rawptr = grid.get();
    sequence->subWeightHistogram(0, min, max, grid_rawptr, bins2);

    for (unsigned int i = 0; i < bins2.size(); i++) {
      CPPUNIT_ASSERT_EQUAL(bins1[i], bins2[i]);
    }

    sequence->subWeightHistogram(0, min, max, centers, bins2);
    for (unsigned int i = 0; i < bins2.size(); i++) {
      CPPUNIT_ASSERT_EQUAL(bins1[i], bins2[i]);
    }
  }

  void test_sub_basic_cdf()
  {
    std::vector<double> cdf;
    QUESO::UniformOneDGrid<double> * grid;

    sequence->subBasicCdf(14, grid, cdf);

    for (unsigned int i = 0; i < cdf.size(); i++) {
      unsigned int cdf_val = cdf[i] * 13;
      CPPUNIT_ASSERT_EQUAL(i, cdf_val);  // Expect linearly increasing ints
    }
  }

  void test_population_variance()
  {
    double actualMean = 6.0;
    double actualPopVar = 14.0;

    double subPopVariance = sequence->subPopulationVariance(0, sequence->subSequenceSize(), actualMean);
    double unifiedPopVariance = sequence->unifiedPopulationVariance(false, 0, sequence->subSequenceSize(), actualMean);

    CPPUNIT_ASSERT_EQUAL(actualPopVar, subPopVariance);
    CPPUNIT_ASSERT_EQUAL(actualPopVar, unifiedPopVariance);
  }

  void test_auto_covariance()
  {
    double actualMean = 6.0;
    double actualPopVar = 14.0;
    unsigned int numPos = sequence->subSequenceSize();

    // Compute auto covariance with zero lag.  Should be population var?
    double autoCov = sequence->autoCovariance(0, numPos, actualMean, 0);
    CPPUNIT_ASSERT_EQUAL(actualPopVar, autoCov);

    // This is normalised by autoCov, so result should be 1.0 with lag 0
    double autoCorrViaDef = sequence->autoCorrViaDef(0, numPos, 0);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrViaDef);

    std::vector<double> autoCorrsViaFft;
    sequence->autoCorrViaFft(0, numPos, 0, autoCorrsViaFft);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrsViaFft[0]);

    sequence->autoCorrViaFft(0, numPos, 1, autoCorrsViaFft[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrsViaFft[0]);
  }

  void test_sort()
  {
    QUESO::ScalarSequence<double> reversed_sequence(*env, 13, "");

    unsigned int size = sequence->subSequenceSize();
    for (unsigned int i = 0; i < size; i++) {
      reversed_sequence[i] = (*sequence)[size-1-i];
    }

    QUESO::ScalarSequence<double> sorted_sequence(*env, 13, "");
    reversed_sequence.subSort(0, sorted_sequence);

    for (unsigned int i = 0; i < size; i++) {
      CPPUNIT_ASSERT_EQUAL((*sequence)[i], sorted_sequence[i]);
    }

    reversed_sequence.unifiedSort(false, 0, sorted_sequence);
    for (unsigned int i = 0; i < size; i++) {
      CPPUNIT_ASSERT_EQUAL((*sequence)[i], sorted_sequence[i]);
    }
  }

  void test_interquantile_range()
  {
    double actual_iqr = 7.0;

    double iqr = sequence->subInterQuantileRange(0);
    CPPUNIT_ASSERT_EQUAL(actual_iqr, iqr);

    iqr = sequence->unifiedInterQuantileRange(false, 0);
    CPPUNIT_ASSERT_EQUAL(actual_iqr, iqr);
  }

  void test_sub_weight_cdf()
  {
    unsigned int num_points = 14;
    std::vector<double> cdf1;
    std::vector<double> cdf2;
    std::vector<double> centers(num_points, 0);
    QUESO::UniformOneDGrid<double> * grid;

    sequence->subWeightCdf(num_points, grid, cdf1);
    sequence->subWeightCdf(num_points, centers, cdf2);

    for (unsigned int i = 0; i < cdf1.size(); i++) {
      unsigned int cdf_val1 = cdf1[i] * 13;
      unsigned int cdf_val2 = cdf2[i] * 13;
      CPPUNIT_ASSERT_EQUAL(i, cdf_val1);  // Expect linearly increasing ints
      CPPUNIT_ASSERT_EQUAL(i, cdf_val2);  // Expect linearly increasing ints
    }
  }

  void test_scale_kde()
  {
    double TOL = 1e-12;
    double actualIrq = 7.0;

    // Regression solution
    double actualKDEScale = 2.471509395452762269940194528317;

    double KDEScale;

    KDEScale = sequence->subScaleForKde(0, actualIrq, 1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDEScale, KDEScale, TOL);

    KDEScale = sequence->unifiedScaleForKde(false, 0, actualIrq, 1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDEScale, KDEScale, TOL);
  }

  void test_gaussian_kde()
  {
    double TOL = 1e-12;

    QUESO::ScalarSequence<double> point(*env, 1, "");
    point[0] = -0.1;

    std::vector<double> positions(1, -0.1);
    std::vector<double> density(1, 0.0);
    point.subGaussian1dKde(0, 1.0, positions, density);

    double actualKDE = 1.0 / std::sqrt(2.0 * M_PI);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDE, density[0], TOL);

    point.unifiedGaussian1dKde(false, 0, 1.0, positions, density);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDE, density[0], TOL);
  }

  void test_positions_of_maximum()
  {
    QUESO::ScalarSequence<double> response(*env, 13, "");
    response[0] = 1.0;
    response[1] = 2.0;
    response[2] = 3.0;
    response[3] = 3.0;
    response[4] = 2.0;
    response[5] = 9.0;
    response[6] = 9.0;
    response[7] = 8.0;
    response[8] = 7.0;
    response[9] = 3.0;
    response[10] = 2.0;
    response[11] = -1.0;
    response[12] = -10.0;

    QUESO::ScalarSequence<double> positions(*env, 0, "");
    double max = sequence->subPositionsOfMaximum(response, positions);
    CPPUNIT_ASSERT_EQUAL(9.0, max);
    CPPUNIT_ASSERT_EQUAL(5.0, positions[0]);
    CPPUNIT_ASSERT_EQUAL(6.0, positions[1]);

    max = sequence->unifiedPositionsOfMaximum(response, positions);
    CPPUNIT_ASSERT_EQUAL(9.0, max);
    CPPUNIT_ASSERT_EQUAL(5.0, positions[0]);
    CPPUNIT_ASSERT_EQUAL(6.0, positions[1]);
  }

  void test_read()
  {
    QUESO::ScalarSequence<double> read_sequence(*env, 5, "");

    // Filename logic
    std::string sequenceFileName = "unit/read_sequence";
    const char * test_srcdir = std::getenv("srcdir");
    if (test_srcdir) {
      sequenceFileName = test_srcdir + ('/' + sequenceFileName);
    }

    read_sequence.unifiedReadContents(sequenceFileName, "m", 5);

    for (unsigned int i = 0; i < 5; i++) {
      double val = i + 1;
      CPPUNIT_ASSERT_EQUAL(val, read_sequence[i]);
    }
  }

  void test_brooks_gelman()
  {
    int returnVal;
    try {
      returnVal = 1;
      sequence->brooksGelmanConvMeasure(true, 0, 0);
    }
    catch (...) {
      returnVal = 0;
    }

    // We expect the exception to be caught because brooks gelman isn't
    // implemented
    CPPUNIT_ASSERT_EQUAL(0, returnVal);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::ScalarSequence<double> >::Type sequence;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ScalarSequenceTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
