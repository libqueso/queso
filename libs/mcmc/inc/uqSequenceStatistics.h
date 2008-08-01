/* uq/libs/mcmc/inc/uqSequenceStatistics.h
 *
 * Copyright (C) 2008 The PECOS Team, http://pecos.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_SEQUENCE_STATISTICS_H__
#define __UQ_SEQUENCE_STATISTICS_H__

#include <uqMatrixOfStuff.h>
#include <iostream>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

template <class V>
void
uqSequenceMean(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  unsigned int                 numPos,
  V&                           mean)
{
  bool bRC = ((0                     <= initialPos         ) &&
              (0                     <  numPos             ) &&
              ((initialPos+numPos-1) <= (sequence.size()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceMean<V>()",
                      "invalid initial position or number of positions");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceMean<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  mean.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < mean.size(); ++j) {
      mean[j] += (*sequence[i])[j]/doubleLoopSize;
    }
  }

  return;
}

template <class V>
void
uqSequenceSampleVariance(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  unsigned int                 numPos,
  const V&                     mean,
  V&                           std)
{
  bool bRC = ((0                     <= initialPos         ) &&
              (0                     <  numPos             ) &&
              ((initialPos+numPos-1) <= (sequence.size()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceSampleVariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (sequence[0]->size() == std.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceSampleVariance<V>()",
                      "incompatible sizes between std vector and vectors in sequence");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceSampleVariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  std.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < std.size(); ++j) {
      double diff = (*sequence[i])[j] - mean[j];
      std[j] += diff*diff;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < std.size(); ++j) {
    std[j] /= (doubleLoopSize - 1.);
  }

  return;
}

template <class V>
void
uqSequencePopulationVariance(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  unsigned int                 numPos,
  const V&                     mean,
  V&                           std)
{
  bool bRC = ((0                     <= initialPos         ) &&
              (0                     <  numPos             ) &&
              ((initialPos+numPos-1) <= (sequence.size()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequencePopulationVariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (sequence[0]->size() == std.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequencePopulationVariance<V>()",
                      "incompatible sizes between std vector and vectors in sequence");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequencePopulationVariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  std.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < std.size(); ++j) {
      double diff = (*sequence[i])[j] - mean[j];
      std[j] += diff*diff;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < std.size(); ++j) {
    std[j] /= doubleLoopSize;
  }

  return;
}

template <class V>
void
uqSequenceAutoCovariance(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  unsigned int                 numPos,
  const V&                     mean,
  unsigned int                 lag,
  V&                           autoCov)
{
  bool bRC = ((0                     <= initialPos         ) &&
              (0                     <  numPos             ) &&
              ((initialPos+numPos-1) <= (sequence.size()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceAutoCovariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (numPos > lag);
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceAutoCovariance<V>()",
                      "lag is too large");

  bRC = (sequence[0]->size() == autoCov.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceAutoCovariance<V>()",
                      "incompatible sizes between autoCov vector and vectors in sequence");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqSequenceAutoCovariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  autoCov.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < autoCov.size(); ++j) {
      double diff1 = (*sequence[i    ])[j] - mean[j];
      double diff2 = (*sequence[i+lag])[j] - mean[j];
      autoCov[j] += diff1*diff2;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < autoCov.size(); ++j) {
    autoCov[j] /= doubleLoopSize;
  }

  return;
}

template <class V>
void
uqSequenceAutoCorrelations(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& lags,
  uqMatrixOfStuff<V>&              matrixOfAutoCorrs) // [numOfPos x numOfLags] matrix
{
  V subChainMean              (*(sequence[0]));
  V subChainAutoCovarianceLag0(*(sequence[0]));

  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    uqSequenceMean(sequence,
                   initialPositions[initialPosId],
                   sequence.size()-initialPositions[initialPosId],
                   subChainMean);
    uqSequenceAutoCovariance(sequence,
                             initialPositions[initialPosId],
                             sequence.size()-initialPositions[initialPosId],
                             subChainMean,
                             0, // lag
                             subChainAutoCovarianceLag0);
    for (unsigned int lagId = 0; lagId < lags.size(); lagId++) {
      uqSequenceAutoCovariance(sequence,
                               initialPositions[initialPosId],
                               sequence.size()-initialPositions[initialPosId],
                               subChainMean,
                               lags[lagId], // lag
                               matrixOfAutoCorrs(initialPosId,lagId));
      matrixOfAutoCorrs(initialPosId,lagId) /= subChainAutoCovarianceLag0; 
    }
  }

  return;
}

template <class V>
void
uqSequenceBMM(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& batchLengths,
  uqMatrixOfStuff<V>&              matrixOfBMM) // [numOfPos x numOfLengths] matrix
{
  V meanOfBatchMeans   (*(sequence[0]));
  V covLag0OfBatchMeans(*(sequence[0]));
  V covLag1OfBatchMeans(*(sequence[0]));

  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    for (unsigned int batchLengthId = 0; batchLengthId < batchLengths.size(); batchLengthId++) {
      unsigned int batchLength = batchLengths[batchLengthId];
      unsigned int numberOfBatches = (sequence.size() - initialPositions[initialPosId])/batchLength;

      V* tmpVector = new V(*(sequence[0])); // In order to contour the fact that 'batchMeans' is a vector of 'const V*', but needs to be set first
      std::vector<const V* > batchMeans(numberOfBatches,NULL);
      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        uqSequenceMean(sequence,
                       initialPositions[initialPosId] + batchId*batchLength,
                       batchLength,
                       *tmpVector);
        batchMeans[batchId] = new V(*tmpVector);
      }
      delete tmpVector;

      uqSequenceMean(batchMeans,
                     0,
                     batchMeans.size(),
                     meanOfBatchMeans);

      uqSequenceAutoCovariance(batchMeans,
                               0,
                               batchMeans.size(),
                               meanOfBatchMeans,
                               0, // lag
                               covLag0OfBatchMeans);

      uqSequenceAutoCovariance(batchMeans,
                               0,
                               batchMeans.size(),
                               meanOfBatchMeans,
                               1, // lag
                               covLag0OfBatchMeans);

      uqSequenceSampleVariance(batchMeans,
                               0,
                               batchMeans.size(),
                               meanOfBatchMeans,
                               matrixOfBMM(initialPosId,batchLengthId));
      matrixOfBMM(initialPosId,batchLengthId) /= (double) batchMeans.size();

      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        if (batchMeans[batchId] != NULL) delete batchMeans[batchId];
      }
    }
  }

  return;
}

template <class V>
void
uqSequencePSD(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& blockLengths,
  double                           hopSizeRatio,
  uqMatrixOfStuff<V>&              matrixOfPSDAtZero) // [numOfPos x numOfLengths] matrix
{
  unsigned int numParams = sequence[0]->size();

  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    unsigned int dataSize = sequence.size() - initialPositions[initialPosId];
    std::vector<double> data(dataSize,0.);
    std::vector<double> avgData(dataSize,0.);

    gsl_fft_real_workspace*        work;
    gsl_fft_real_wavetable*        wvtable;
    //gsl_fft_halfcomplex_wavetable* hc;

    work    = gsl_fft_real_workspace_alloc       (dataSize);
    wvtable = gsl_fft_real_wavetable_alloc       (dataSize);
    //hc      = gsl_fft_halfcomplex_wavetable_alloc(dataSize);

    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int blockLengthId = 0; blockLengthId < blockLengths.size(); blockLengthId++) {
        unsigned int blockLength = blockLengths[blockLengthId];
        unsigned int numberOfBlocks = (sequence.size() - initialPositions[initialPosId])/blockLength;
        unsigned int hopSize = (unsigned int) ( ((double) blockLength) * hopSizeRatio );
        hopSize = 999;

        for (unsigned int j = 0; j < dataSize; ++j) {
          avgData[j] = 0.;
        }
        for (unsigned int blockId = 0; blockId < numberOfBlocks; blockId++) {
          for (unsigned int j = 0; j < dataSize; ++j) {
            data[j] = (*(sequence[initialPositions[initialPosId]+j]))[i];
          }

          gsl_fft_real_transform(&data[0],1,dataSize,wvtable,work);

          for (unsigned int j = 0; j < dataSize; ++j) {
            avgData[j] += data[j]/(double) numberOfBlocks;
          }

          //double sumOfAllChainTerms = 0.;
          //for (unsigned int j = 0; j < dataSize; ++j) {
          //  sumOfAllChainTerms += (*(sequence[initialPositions[initialPosId]+j]))[i];
          //}
          //std::cout << "After FFT, for parameter i = "    << i
          //          << ", sumOfAllChainTerms = "          << sumOfAllChainTerms
          //          << ", sumOfAllChainTerms - dft[0] = " << sumOfAllChainTerms - data[0]
          //          << std::endl;
        }

        matrixOfPSDAtZero(initialPosId,blockLengthId)[i] = avgData[0];
      }
    }

    //gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_wavetable_free       (wvtable);
    gsl_fft_real_workspace_free       (work);
  }

  return;
}

template <class V>
void
uqSequenceGeweke(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  double                           ratioNa,
  double                           ratioNb,
  std::vector<double>&             vectorOfGeweke) // [numPos x 1] vector
{
  return;
}
#endif // __UQ_SEQUENCE_STATISTICS_H__

