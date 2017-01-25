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
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/SequenceOfVectors.h>
#include <queso/BoxSubset.h>

namespace QUESOTesting
{

class BaseVectorSequenceTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(BaseVectorSequenceTest);
  CPPUNIT_TEST(test_size);
  CPPUNIT_TEST(test_min_max);
  CPPUNIT_TEST(test_median);
  CPPUNIT_TEST(test_sample_var);
  CPPUNIT_TEST(test_box);
  CPPUNIT_TEST(test_gaussian);
  CPPUNIT_TEST(test_uniform);
  CPPUNIT_TEST(test_compute_filter_params);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));
    space.reset(new QUESO::VectorSpace<>(*env, "", 2, NULL));
    sequence.reset(new QUESO::SequenceOfVectors<>(*space, 13, ""));

    // Fill up the sequence with some not-so-random stuff
    QUESO::GslVector v(space->zeroVector());
    for (unsigned int i = 0; i < sequence->subSequenceSize(); i++) {
      v[0] = i;
      v[1] = i + 1;
      sequence->setPositionValues(i, v);
    }
  }

  void test_size()
  {
    unsigned int dim = 2;
    CPPUNIT_ASSERT_EQUAL(dim, sequence->vectorSizeGlobal());
  }

  void test_min_max()
  {
    QUESO::GslVector min(sequence->subMinPlain());
    CPPUNIT_ASSERT_EQUAL(0.0, min[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, min[1]);
    sequence->deleteStoredVectors();

    QUESO::GslVector max(sequence->subMaxPlain());
    CPPUNIT_ASSERT_EQUAL(12.0, max[0]);
    CPPUNIT_ASSERT_EQUAL(13.0, max[1]);
    sequence->deleteStoredVectors();

    QUESO::GslVector max2(sequence->unifiedMaxPlain());
    CPPUNIT_ASSERT_EQUAL(12.0, max2[0]);
    CPPUNIT_ASSERT_EQUAL(13.0, max2[1]);
    sequence->deleteStoredVectors();
  }

  void test_median()
  {
    QUESO::GslVector median(sequence->subMedianPlain());
    CPPUNIT_ASSERT_EQUAL(6.0, median[0]);
    CPPUNIT_ASSERT_EQUAL(7.0, median[1]);
    sequence->deleteStoredVectors();

    QUESO::GslVector median2(sequence->unifiedMedianPlain());
    CPPUNIT_ASSERT_EQUAL(6.0, median2[0]);
    CPPUNIT_ASSERT_EQUAL(7.0, median2[1]);
    sequence->deleteStoredVectors();
  }

  void test_sample_var()
  {
    double actualVar = 182.0 / 12.0;

    QUESO::GslVector var(sequence->subSampleVariancePlain());
    CPPUNIT_ASSERT_EQUAL(actualVar, var[0]);
    CPPUNIT_ASSERT_EQUAL(actualVar, var[1]);
    sequence->deleteStoredVectors();
  }

  void test_box()
  {
    CPPUNIT_ASSERT_EQUAL(0.0, sequence->subBoxPlain().minValues()[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, sequence->subBoxPlain().minValues()[1]);
    CPPUNIT_ASSERT_EQUAL(12.0, sequence->subBoxPlain().maxValues()[0]);
    CPPUNIT_ASSERT_EQUAL(13.0, sequence->subBoxPlain().maxValues()[1]);
  }

  void test_gaussian()
  {
    QUESO::GslVector mean(space->zeroVector());
    mean.cwSet(5.0);
    QUESO::GslVector std(space->zeroVector());
    std.cwSet(0.1);

    unsigned int N = 1000000;
    QUESO::SequenceOfVectors<> gaussian_seq(*space, N, "");
    for (unsigned int i = 0; i < N; i++) {
      gaussian_seq.setPositionValues(i, space->zeroVector());
    }

    gaussian_seq.setGaussian(mean, std);

    // Check the mean is close to (5, 5)
    QUESO::GslVector sampleMean(gaussian_seq.subMeanPlain());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mean[0], sampleMean[0], 0.001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mean[1], sampleMean[1], 0.001);

    // Check the var is close to 0.1
    QUESO::GslVector sampleVar(gaussian_seq.subSampleVariancePlain());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std[0]* std[0], sampleVar[0], 0.001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std[1]* std[1], sampleVar[1], 0.001);
  }

  void test_uniform()
  {
    QUESO::GslVector mins(space->zeroVector());
    mins.cwSet(3.0);
    QUESO::GslVector maxs(space->zeroVector());
    maxs.cwSet(7.0);

    unsigned int N = 1000000;
    QUESO::SequenceOfVectors<> uniform_seq(*space, N, "");
    for (unsigned int i = 0; i < N; i++) {
      uniform_seq.setPositionValues(i, space->zeroVector());
    }

    uniform_seq.setUniform(mins, maxs);

    // Check the mean is close to (5, 5)
    QUESO::GslVector sampleMean(uniform_seq.subMeanPlain());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, sampleMean[0], 0.01);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, sampleMean[1], 0.01);
  }

  void test_compute_filter_params()
  {
    unsigned int initialPos, spacing;

    // What is this even supposed to do?
    sequence->computeFilterParams(NULL, initialPos, spacing);

    unsigned int expected = 1;
    CPPUNIT_ASSERT_EQUAL(expected, spacing);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::VectorSpace<> >::Type space;
  typename QUESO::ScopedPtr<QUESO::SequenceOfVectors<> >::Type sequence;
};

CPPUNIT_TEST_SUITE_REGISTRATION(BaseVectorSequenceTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
