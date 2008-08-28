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

#include <uq2dArrayOfStuff.h>
#include <iostream>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

template <class V>
void
uqVectorSequenceMean(
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
                      "uqVectorSequenceMean<V>()",
                      "invalid initial position or number of positions");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceMean<V>()",
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
uqVectorSequenceSampleVariance(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  unsigned int                 numPos,
  const V&                     mean,
  V&                           sampleVariance)
{
  bool bRC = ((0                     <= initialPos         ) &&
              (0                     <  numPos             ) &&
              ((initialPos+numPos-1) <= (sequence.size()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceSampleVariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (sequence[0]->size() == sampleVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceSampleVariance<V>()",
                      "incompatible sizes between sampleVariance vector and vectors in sequence");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceSampleVariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  sampleVariance.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < sampleVariance.size(); ++j) {
      double diff = (*sequence[i])[j] - mean[j];
      sampleVariance[j] += diff*diff;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < sampleVariance.size(); ++j) {
    sampleVariance[j] /= (doubleLoopSize - 1.);
  }

  return;
}

template <class V>
void
uqVectorSequencePopulationVariance(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  unsigned int                 numPos,
  const V&                     mean,
  V&                           populVariance)
{
  bool bRC = ((0                     <= initialPos         ) &&
              (0                     <  numPos             ) &&
              ((initialPos+numPos-1) <= (sequence.size()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequencePopulationVariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (sequence[0]->size() == populVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequencePopulationVariance<V>()",
                      "incompatible sizes between populVariance vector and vectors in sequence");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequencePopulationVariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  populVariance.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < populVariance.size(); ++j) {
      double diff = (*sequence[i])[j] - mean[j];
      populVariance[j] += diff*diff;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < populVariance.size(); ++j) {
    populVariance[j] /= doubleLoopSize;
  }

  return;
}

template <class V>
void
uqVectorSequenceAutoCovariance(
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
                      "uqVectorSequenceAutoCovariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (numPos > lag);
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "lag is too large");

  bRC = (sequence[0]->size() == autoCov.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "incompatible sizes between autoCov vector and vectors in sequence");

  bRC = (sequence[0]->size() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      sequence[0]->env().rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
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
uqVectorSequenceAutoCorrelations(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& lags,
  uq2dArrayOfStuff<V>&             _2dArrayOfAutoCorrs) // [numOfPos x numOfLags] matrix
{
  V subChainMean              (*(sequence[0]));
  V subChainAutoCovarianceLag0(*(sequence[0]));

  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    uqVectorSequenceMean(sequence,
                         initialPositions[initialPosId],
                         sequence.size()-initialPositions[initialPosId],
                         subChainMean);
    uqVectorSequenceAutoCovariance(sequence,
                                   initialPositions[initialPosId],
                                   sequence.size()-initialPositions[initialPosId],
                                   subChainMean,
                                   0, // lag
                                   subChainAutoCovarianceLag0);
    for (unsigned int lagId = 0; lagId < lags.size(); lagId++) {
      uqVectorSequenceAutoCovariance(sequence,
                                     initialPositions[initialPosId],
                                     sequence.size()-initialPositions[initialPosId],
                                     subChainMean,
                                     lags[lagId], // lag
                                     _2dArrayOfAutoCorrs(initialPosId,lagId));
      _2dArrayOfAutoCorrs(initialPosId,lagId) /= subChainAutoCovarianceLag0; 
    }
  }

  return;
}

template <class V>
void
uqVectorSequenceBMM(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& batchLengths,
  uq2dArrayOfStuff<V>&             _2dArrayOfBMM) // [numOfPos x numOfLengths] matrix
{
  V meanOfBatchMeans   (*(sequence[0]));
  V covLag0OfBatchMeans(*(sequence[0]));
  V covLag1OfBatchMeans(*(sequence[0]));

  V* tmpVector = new V(*(sequence[0])); // In order to contour the fact that 'batchMeans' is a vector of 'const V*', but needs to be set first
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    for (unsigned int batchLengthId = 0; batchLengthId < batchLengths.size(); batchLengthId++) {
      unsigned int batchLength = batchLengths[batchLengthId];
      unsigned int numberOfBatches = (sequence.size() - initialPositions[initialPosId])/batchLength;

      std::vector<const V* > batchMeans(numberOfBatches,NULL);
      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        uqVectorSequenceMean(sequence,
                             initialPositions[initialPosId] + batchId*batchLength,
                             batchLength,
                             *tmpVector);
        batchMeans[batchId] = new V(*tmpVector);
      }

      uqVectorSequenceMean(batchMeans,
                           0,
                           batchMeans.size(),
                           meanOfBatchMeans);

      uqVectorSequenceAutoCovariance(batchMeans,
                                     0,
                                     batchMeans.size(),
                                     meanOfBatchMeans,
                                     0, // lag
                                     covLag0OfBatchMeans);

      uqVectorSequenceAutoCovariance(batchMeans,
                                     0,
                                     batchMeans.size(),
                                     meanOfBatchMeans,
                                     1, // lag
                                     covLag0OfBatchMeans);

      uqVectorSequenceSampleVariance(batchMeans,
                                     0,
                                     batchMeans.size(),
                                     meanOfBatchMeans,
                                     _2dArrayOfBMM(initialPosId,batchLengthId));
      _2dArrayOfBMM(initialPosId,batchLengthId) /= (double) batchMeans.size(); // CHECK
      //_2dArrayOfBMM(initialPosId,batchLengthId) *= (double) (sequence.size() - initialPositions[initialPosId]); // CHECK

      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        if (batchMeans[batchId] != NULL) delete batchMeans[batchId];
      }
    }
  }
  delete tmpVector;

  return;
}

void
uqScalarSequencePSD(
  const std::vector<double>& data,
  unsigned int               numBlocks,
  double                     hopSizeRatio,
  std::vector<double>&       psdData) // [dataSize x 1] vector
{
  psdData.clear();
  unsigned int dataSize = data.size();

  double tmp = ((double) dataSize)/(( ((double) numBlocks) - 1. )*hopSizeRatio + 1.);
  unsigned int blockSize = (unsigned int) tmp;
  unsigned int hopSize   = (unsigned int) ( ((double) blockSize) * hopSizeRatio );
  tmp = ((double) dataSize) - ( ((double) numBlocks) - 1.) * ((double) hopSize) - ((double) blockSize);
#if 1
  unsigned int numberOfDiscardedDataElements = (unsigned int) tmp;
  std::cout << "N = "         << dataSize
            << ", #Blocks = " << numBlocks
            << ", R = "       << hopSize
            << ", B = "       << blockSize
            << ", overlap = " << blockSize - hopSize
            << ", [(#Blocks - 1) * R + B] = "       << (numBlocks-1)*hopSize + blockSize
            << ", numberOfDiscardedDataElements = " << numberOfDiscardedDataElements
            << ", tmp = "                           << tmp
            << std::endl;
#endif
  UQ_FATAL_TEST_MACRO(tmp < 0.,
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarSequencePSD()",
                      "eventual extra space for last block should not be negative");

  tmp = log((double) blockSize)/log(2.);
  double fractionalPart = tmp - ((double) ((unsigned int) tmp));
  if (fractionalPart > 0.) tmp += (1. - fractionalPart);
  unsigned int fftSize = (unsigned int) pow(2.,tmp);
  //std::cout << "fractionalPart = " << fractionalPart
  //          << ", B = "            << blockSize
  //          << ", fftSize = "      << fftSize
  //          << std::endl;

  gsl_fft_real_workspace*        wkSpace;
  gsl_fft_real_wavetable*        wvTable;
  //gsl_fft_halfcomplex_wavetable* hc;

  wkSpace = gsl_fft_real_workspace_alloc       (fftSize);
  wvTable = gsl_fft_real_wavetable_alloc       (fftSize);
  //hc      = gsl_fft_halfcomplex_wavetable_alloc(fftSize);

  unsigned int halfFFTSize = fftSize/2;
  psdData.resize(1+halfFFTSize,0.);
  for (unsigned int blockId = 0; blockId < numBlocks; blockId++) {
    // Padding
    std::vector<double> blockData(fftSize,0.);

    // Fill block using window on full data
    for (unsigned int j = 0; j < blockSize; ++j) {
      unsigned int dataPos = j+blockId*hopSize;
      UQ_FATAL_TEST_MACRO(dataPos >= dataSize,
                          UQ_UNAVAILABLE_RANK,
                          "uqScalarSequencePSD()",
                          "too large position to be accessed in data");
      blockData[j] = uqMiscHammingWindow(fftSize-1,j) * data[dataPos];
    }

    //double sumOfAllTerms = 0.;
    //for (unsigned int j = 0; j < fftSize; ++j) {
    //  sumOfAllTerms += blockData[j];
    //}
    gsl_fft_real_transform(&blockData[0],1,fftSize,wvTable,wkSpace);
    //std::cout << "After FFT"
    //          << ", sumOfAllTerms = "          << sumOfAllTerms
    //          << ", sumOfAllTerms - dft[0] = " << sumOfAllTerms - blockData[0]
    //          << std::endl;

    // Normalized spectral density: power per radians per sample
    double factor = 1./((double) numBlocks*blockSize); // /M_PI; // CHECK
    double realPartOfFFT = 0.;
    double imagPartOfFFT = 0.;
    for (unsigned int j = 0; j < psdData.size(); ++j) {
      if (j == 0) {
        realPartOfFFT = blockData[j];
        imagPartOfFFT = 0.;
      }
      else if (j < halfFFTSize) {
        realPartOfFFT = blockData[j];
        imagPartOfFFT = blockData[fftSize-j];
      }
      else if (j == halfFFTSize) {
        realPartOfFFT = blockData[j];
        imagPartOfFFT = 0.;
      }
      else {
        realPartOfFFT =  blockData[fftSize-j];
        imagPartOfFFT = -blockData[j];
      }
      psdData[j] += (realPartOfFFT*realPartOfFFT + imagPartOfFFT*imagPartOfFFT) * factor;
    }
  }

  //gsl_fft_halfcomplex_wavetable_free(hc);
  gsl_fft_real_wavetable_free       (wvTable);
  gsl_fft_real_workspace_free       (wkSpace);

  return;
}

template <class V>
void
uqVectorSequencePSDAtZero(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& numsOfBlocks,
  double                           hopSizeRatio,
  uq2dArrayOfStuff<V>&             _2dArrayOfPSDAtZero) // [numOfPos x numOfBlocks] matrix
{
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    unsigned int dataSize = sequence.size() - initialPositions[initialPosId];
    std::vector<double> data(dataSize,0.);

    unsigned int numParams = sequence[0]->size();
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSize; ++j) {
	data[j] = (*(sequence[initialPositions[initialPosId]+j]))[i];
      }
      for (unsigned int numsOfBlocksId = 0; numsOfBlocksId < numsOfBlocks.size(); numsOfBlocksId++) {
        unsigned int numBlocks = numsOfBlocks[numsOfBlocksId];
        std::vector<double> psdData(0,0.); // size will be determined by 'uqScalarSequencePSD()'
        uqScalarSequencePSD(data,
                            numBlocks,
                            hopSizeRatio,
                            psdData);
        _2dArrayOfPSDAtZero(initialPosId,numsOfBlocksId)[i] = psdData[0];
	//std::cout << "psdData[0] = " << psdData[0] << std::endl;
        std::cout << "psdData = zeros(" << psdData.size() << ",1);" << std::endl;
        for (unsigned j = 0; j < psdData.size(); ++j) {
    	  std::cout << "psdData(" << j+1 << ") = " << psdData[j] << ";" << std::endl;
        }
      } // for 'numsOfBlocksId'
    } // for 'i'
  }

  return;
}

template <class V>
void
uqVectorSequenceGeweke(
  const std::vector<const V*>&     sequence,
  const std::vector<unsigned int>& initialPositions,
  double                           ratioNa,
  double                           ratioNb,
  std::vector<V*>&                 vectorOfGeweke)
{
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    unsigned int fullDataSize = sequence.size() - initialPositions[initialPosId];
    unsigned int dataSizeA    = (unsigned int) (((double) fullDataSize) * ratioNa);
    unsigned int dataSizeB    = (unsigned int) (((double) fullDataSize) * ratioNb);
    unsigned int initialPosA  = initialPositions[initialPosId];
    unsigned int initialPosB  = sequence.size() - dataSizeB;

    V meanA(*(sequence[0]));
    uqVectorSequenceMean(sequence,
                         initialPosA,
                         dataSizeA,
                         meanA);

    V meanB(*(sequence[0]));
    uqVectorSequenceMean(sequence,
                         initialPosB,
                         dataSizeB,
                         meanB);

    unsigned int numParams = sequence[0]->size();

    V psdAtZeroA(*(sequence[0]));
    std::vector<double> dataA(dataSizeA,0.);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeA; ++j) {
	dataA[j] = (*(sequence[initialPosA+j]))[i];
      }
      std::vector<double> psdData(0,0.);
      uqScalarSequencePSD(dataA,
                          8,  // numBlocks
                          .5, // hopSizeRatio
                          psdData);
      psdAtZeroA[i] = psdData[0];
    } // for 'i'

    V psdAtZeroB(*(sequence[0]));
    std::vector<double> dataB(dataSizeB,0.);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeB; ++j) {
	dataB[j] = (*(sequence[initialPosB+j]))[i];
      }
      std::vector<double> psdData(0,0.);
      uqScalarSequencePSD(dataB,
                          8,  // numBlocks
                          .5, // hopSizeRatio
                          psdData);
      psdAtZeroB[i] = psdData[0];
    } // for 'i'

    vectorOfGeweke[initialPosId] = new V(*(sequence[0]));

    double doubleDataSizeA = (double) dataSizeA;
    double doubleDataSizeB = (double) dataSizeB;
    for (unsigned int i = 0; i < numParams; ++i) {
      (*(vectorOfGeweke[initialPosId]))[i] = (meanA[i] - meanB[i])/sqrt(psdAtZeroA[i]/doubleDataSizeA + psdAtZeroB[i]/doubleDataSizeB);
    }
  }

  return;
}

template <class V>
void
uqVectorSequenceMinMax(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPos,
  V&                           mins,
  V&                           maxs)
{
  unsigned int dataSize = sequence.size() - initialPos;
  unsigned int numParams = sequence[0]->size();
  std::vector<double> data(dataSize,0.);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPos+j]))[i];
    }
    std::vector<double>::iterator pos;
    pos = std::min_element(data.begin(), data.end());
    mins[i] = *pos;
    pos = std::max_element(data.begin(), data.end());
    maxs[i] = *pos;
  }

  return;
}

void
uqScalarSequenceHistogram(
  const std::vector<double>& data,
  double                     minHorizontalValue,
  double                     maxHorizontalValue,
  std::vector<double>&       centers,
  std::vector<double>&       bins)
{
  UQ_FATAL_TEST_MACRO(centers.size() != bins.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarSequenceHistogram()",
                      "vectors 'centers' and 'bins' have different sizes");

  UQ_FATAL_TEST_MACRO(bins.size() < 3,
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarSequenceHistogram()",
                      "number of 'bins' is too small: should be at least 3");

  for (unsigned int j = 0; j < bins.size(); ++j) {
    centers[j] = 0;
    bins[j] = 0;
  }

  double horizontalDelta = (maxHorizontalValue - minHorizontalValue)/(((double) bins.size()) - 2.);

  double minCenter = minHorizontalValue - horizontalDelta/2.;
  double maxCenter = maxHorizontalValue + horizontalDelta/2.;
  for (unsigned int j = 0; j < centers.size(); ++j) {
    double factor = ((double) j)/(((double) centers.size()) - 1.);
    centers[j] = (1. - factor) * minCenter + factor * maxCenter;
  }

  unsigned int dataSize = data.size();
  for (unsigned int j = 0; j < dataSize; ++j) {
    double value = data[j];
    if (value < minHorizontalValue) {
      bins[0]++;
    }
    else if (value >= maxHorizontalValue) {
      bins[bins.size()-1]++;
    }
    else {
      unsigned int index = 1 + (unsigned int) ((value - minHorizontalValue)/horizontalDelta);
      bins[index]++;
    }
  }

  return;
}

template <class V>
void
uqVectorSequenceHistogram(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPosition,
  unsigned int                 spacing,
  const V&                     minHorizontalValues,
  const V&                     maxHorizontalValues,
  std::vector<V*>&             centersForAllBins,
  std::vector<V*>&             binsForAllParams)
{
  UQ_FATAL_TEST_MACRO(centersForAllBins.size() != binsForAllParams.size(),
                      sequence[0]->env().rank(),
                      "uqVectorSequenceHistogram<V>()",
                      "vectors 'centers' and 'bins' have different sizes");

  for (unsigned int j = 0; j < binsForAllParams.size(); ++j) {
    centersForAllBins[j] = new V(*(sequence[0]));
    binsForAllParams [j] = new V(*(sequence[0]));
  }

  unsigned int dataSize = sequence.size() - initialPosition;
  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    std::vector<double> data(dataSize,0.);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPosition+j]))[i];
    }

    std::vector<double> centers(centersForAllBins.size(),0.);
    std::vector<double> bins   (binsForAllParams.size(), 0.);
    uqScalarSequenceHistogram(data,
                              minHorizontalValues[i],
                              maxHorizontalValues[i],
                              centers,
                              bins);

    for (unsigned int j = 0; j < bins.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(binsForAllParams [j]))[i] = bins[j];
    }
  }

  return;
}

template <class V>
void
uqVectorSequenceSort(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPosition,
  std::vector<V*>&             sortedSequence)
{
  UQ_FATAL_TEST_MACRO((sequence.size() - initialPosition) != sortedSequence.size(),
                      sequence[0]->env().rank(),
                      "uqVectorSequenceSort<V>()",
                      "incompatible sizes between vectors 'sequence' and 'sortedSequence'");

  for (unsigned int j = 0; j < sortedSequence.size(); ++j) {
    sortedSequence[j] = new V(*(sequence[0]));
  }

  unsigned int dataSize = sequence.size() - initialPosition;
  unsigned int numParams = sequence[0]->size();
  std::vector<double> data(dataSize,0.);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPosition+j]))[i];
    }

    std::sort(data.begin(), data.end());

    for (unsigned int j = 0; j < dataSize; ++j) {
      (*(sortedSequence[j]))[i] = data[j];
    }
  }

  return;
}

template <class V>
void
uqVectorSequenceInterQuantileRange(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPosition,
  unsigned int                 spacing,
  V&                           iqrs)
{
  unsigned int dataSize = sequence.size() - initialPosition;

  std::vector<V*> sortedSequence(dataSize,NULL);
  uqVectorSequenceSort(sequence,
                       initialPosition,
                       sortedSequence);

  unsigned int pos1 = (unsigned int) ( (((double) dataSize) + 1.)*1./4. - 1. );
  unsigned int pos3 = (unsigned int) ( (((double) dataSize) + 1.)*3./4. - 1. );

  double fraction1 = (((double) dataSize) + 1.)*1./4. - 1. - ((double) pos1);
  double fraction3 = (((double) dataSize) + 1.)*3./4. - 1. - ((double) pos3);

  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    double value1 = (1.-fraction1) * (*sortedSequence[pos1])[i] + fraction1 * (*sortedSequence[pos1+1])[i];
    double value3 = (1.-fraction3) * (*sortedSequence[pos3])[i] + fraction3 * (*sortedSequence[pos3+1])[i];
    iqrs[i] = value3 - value1;
  }

  return;
}

template <class V>
void
uqVectorSequenceScalesForKDE(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPosition,
  unsigned int                 spacing,
  const V&                     iqrs,
  V&                           scales)
{
  unsigned int dataSize = sequence.size() - initialPosition;

  V mean(*(sequence[0]));
  uqVectorSequenceMean(sequence,
                       initialPosition,
                       dataSize,
                       mean);

  V sampleVariance(*(sequence[0]));
  uqVectorSequenceSampleVariance(sequence,
                                 initialPosition,
                                 dataSize,
                                 mean,
                                 sampleVariance);

  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    if (iqrs[i] <= 0.) {
      scales[i] = 1.06*sqrt(sampleVariance[i])/pow(dataSize,1./5.);
    }
    else {
      scales[i] = 1.06*std::min(sqrt(sampleVariance[i]),iqrs[i]/1.34)/pow(dataSize,1./5.);
    }
  }

  return;
}

template <class V>
void
uqVectorSequenceGaussianKDE(
  const std::vector<const V*>& sequence,
  unsigned int                 initialPosition,
  unsigned int                 spacing,
  const std::vector<V*>&       evaluationPositions,
  const V&                     scales,
  std::vector<V*>&             densityValues)
{
  unsigned int dataSize = sequence.size() - initialPosition;
  unsigned int numEstimationsPerParam = evaluationPositions.size();

  for (unsigned int j = 0; j < numEstimationsPerParam; ++j) {
    densityValues[j] = new V(*(sequence[0]));
  }

  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    double scaleInv = 1./scales[i];
    for (unsigned int j = 0; j < numEstimationsPerParam; ++j) {
      double x = (*(evaluationPositions[j]))[i];
      double value = 0.;
      for (unsigned int k = 0; k < dataSize; ++k) {
        double xk = (*(sequence[initialPosition+k]))[i];
        value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
      }
      (*(densityValues[j]))[i] = scaleInv * (value/(double) numEstimationsPerParam);
    }
  }

  return;
}
#endif // __UQ_SEQUENCE_STATISTICS_H__

