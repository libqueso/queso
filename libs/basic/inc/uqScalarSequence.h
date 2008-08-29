/* uq/libs/basic/inc/uqScalarSequence.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
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

#ifndef __UQ_SCALAR_SEQUENCE_H__
#define __UQ_SCALAR_SEQUENCE_H__

#include <uqFft.h>
#include <uqEnvironment.h>
#include <uqDefines.h>
#include <vector>
#include <complex>
#include <gsl/gsl_randist.h>

template <class T>
class uqScalarSequenceClass
{
public:
  typedef typename std::vector<T>::iterator       seqScalarPositionIteratorTypedef;
  typedef typename std::vector<T>::const_iterator seqScalarPositionConstIteratorTypedef;
  uqScalarSequenceClass(const uqEnvironmentClass& env, unsigned int sequenceSize);
 ~uqScalarSequenceClass();

  const unsigned int sequenceSize      () const;
        void         resizeSequence    (unsigned int newSequenceSize);
        void         resetValues       ();
        void         erasePositions    (unsigned int posBegin, unsigned int posEnd);
  const T&           operator[]        (unsigned int positionId) const;
        T&           operator[]        (unsigned int positionId);
        void         setGaussian       (gsl_rng* rng, const T& mean, const T& stdDev);
        T            mean              (unsigned int initialPos,
                                        unsigned int numPos) const;
        T            sampleVariance    (unsigned int initialPos,
                                        unsigned int numPos,
                                        const T&     mean) const;
        T            populationVariance(unsigned int initialPos,
                                        unsigned int numPos,
                                        const T&     mean) const;
        T            autoCovariance    (unsigned int initialPos,
                                        unsigned int numPos,
                                        const T&     mean,
                                        unsigned int lag) const;
        T            autoCorrelation   (unsigned int initialPosition,
                                        unsigned int lag) const;
        T            bmm               (unsigned int initialPosition,
                                        unsigned int batchLength) const;
        void         psd               (unsigned int               initialPosition,
                                        unsigned int               numBlocks,
                                        double                     hopSizeRatio,
                                        std::vector<double>&       psdData) const;
        T            geweke            (unsigned int               initialPosition,
                                        double                     ratioNa,
                                        double                     ratioNb) const;
        void         minMax            (unsigned int               initialPos,
                                        T&                         minValue,
                                        T&                         maxValue) const;
        void         histogram         (unsigned int               initialPosition,
                                        unsigned int               spacing,
                                        const T&                   minHorizontalValue,
                                        const T&                   maxHorizontalValue,
                                        std::vector<T>&            centers,
                                        std::vector<unsigned int>& bins) const;
#if 0
        void         sort              (unsigned int               initialPosition,
                                        uqScalarSequenceClass<T>&  sortedSequence) const;
        void         interQuantileRange(unsigned int               initialPosition,
                                        unsigned int               spacing,
                                        T&                         iqrValue) const;
        void         scalesForKDE      (unsigned int               initialPosition,
                                        unsigned int               spacing,
                                        const T&                   iqrsValue,
                                        double&                    scaleValue) const;
        void         gaussianKDE       (unsigned int               initialPosition,
                                        unsigned int               spacing,
                                        const std::vector<T>&      evaluationPositions,
                                        double                     scaleValue,
                                        std::vector<double>&       densityValues) const;
#endif
private:
  const uqEnvironmentClass& m_env;
  std::vector<T>            m_seq;
};

template <class T>
uqScalarSequenceClass<T>::uqScalarSequenceClass(
  const uqEnvironmentClass& env,
        unsigned int        sequenceSize)
  :
  m_env(env),
  m_seq(sequenceSize)
{
}

template <class T>
uqScalarSequenceClass<T>::~uqScalarSequenceClass()
{
}

template <class T>
const unsigned int
uqScalarSequenceClass<T>::sequenceSize() const
{
  return m_seq.size();
}

template <class T>
void
uqScalarSequenceClass<T>::resizeSequence(unsigned int newSequenceSize)
{
  m_seq.resize(newSequenceSize);
  std::vector<T>(m_seq).swap(m_seq);

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::resetValues()
{
  return;
}

template <class T>
void
uqScalarSequenceClass<T>::erasePositions(unsigned int posBegin, unsigned int posEnd)
{
  seqScalarPositionIteratorTypedef positionIteratorBegin = m_seq.begin();
  if (posBegin < m_seq.size()) std::advance(positionIteratorBegin,posBegin);
  else                         positionIteratorBegin = m_seq.end();

  seqScalarPositionIteratorTypedef positionIteratorEnd = m_seq.begin();
  if (posEnd < m_seq.size()) std::advance(positionIteratorEnd,posEnd);
  else                       positionIteratorEnd = m_seq.end();

  m_seq.erase(positionIteratorBegin,positionIteratorEnd);
  UQ_FATAL_TEST_MACRO((posBegin != m_seq.size()),
                      m_env.rank(),
                      "uqScalarSequences<V>::erasePositions()",
                      "posBegin != m_seq.size()");

  return;
}

template <class T>
const T&
uqScalarSequenceClass<T>::operator[](unsigned int positionId) const
{
  return m_seq[positionId];
}

template <class T>
T&
uqScalarSequenceClass<T>::operator[](unsigned int positionId)
{
  return m_seq[positionId];
}

template <class T>
void
uqScalarSequenceClass<T>::setGaussian(gsl_rng* rng, const T& mean, const T& stdDev)
{
  unsigned int maxJ = this->sequenceSize();
  if (mean == 0) {
    for (unsigned int j = 0; j < maxJ; ++j) {
      m_seq[j] = gsl_ran_gaussian(rng,stdDev);
    }
  }
  else {
    for (unsigned int j = 0; j < maxJ; ++j) {
      m_seq[j] = mean + gsl_ran_gaussian(rng,stdDev);
    }
  }

  return;
}

template <class T>
T
uqScalarSequenceClass<T>::mean(
  unsigned int initialPos,
  unsigned int numPos) const
{
  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;

  double result = 0.;
  for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
    result += m_seq[j];
  }

  return result/(double) loopSize;
}

template <class T>
T
uqScalarSequenceClass<T>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const T&     mean) const
{
  return;
}

template <class T>
T
uqScalarSequenceClass<T>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const T&     mean) const
{
  return;
}

template <class T>
T
uqScalarSequenceClass<T>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const T&     mean,
  unsigned int lag) const
{
  return;
}

template <class T>
T
uqScalarSequenceClass<T>::autoCorrelation(
  unsigned int initialPosition,
  unsigned int lag) const
{
  return;
}

template <class T>
T
uqScalarSequenceClass<T>::bmm(
  unsigned int initialPosition,
  unsigned int batchLength) const
{
  return;
}

template <class T>
void
uqScalarSequenceClass<T>::psd(
  unsigned int         initialPosition, // FIX IT
  unsigned int         numBlocks,
  double               hopSizeRatio,
  std::vector<double>& psdData) const
{
  psdData.clear();
  unsigned int dataSize = this->sequenceSize();

  double tmp = ((double) dataSize)/(( ((double) numBlocks) - 1. )*hopSizeRatio + 1.);
  unsigned int blockSize = (unsigned int) tmp;
  unsigned int hopSize   = (unsigned int) ( ((double) blockSize) * hopSizeRatio );
  tmp = ((double) dataSize) - ( ((double) numBlocks) - 1.) * ((double) hopSize) - ((double) blockSize);
#if 0
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
                      "uqScalarSequenceClass<T>::psd()",
                      "eventual extra space for last block should not be negative");

  tmp = log((double) blockSize)/log(2.);
  double fractionalPart = tmp - ((double) ((unsigned int) tmp));
  if (fractionalPart > 0.) tmp += (1. - fractionalPart);
  unsigned int fftSize = (unsigned int) pow(2.,tmp);
  //std::cout << "fractionalPart = " << fractionalPart
  //          << ", B = "            << blockSize
  //          << ", fftSize = "      << fftSize
  //          << std::endl;

  uqFftClass<T> fftObj(m_env,fftSize);
  std::vector<std::complex<double> > resultData(fftSize,std::complex<double>(0.,0.));

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
                          "uqScalarSequenceClass<T>::psd()",
                          "too large position to be accessed in data");
      blockData[j] = uqMiscHammingWindow(fftSize-1,j) * m_seq[dataPos];
    }

    fftObj.forward(blockData,resultData);

    // Normalized spectral density: power per radians per sample
    double factor = 1./((double) numBlocks*blockSize); // /M_PI; // CHECK
    for (unsigned int j = 0; j < psdData.size(); ++j) {
      psdData[j] += norm(resultData[j]) * factor;
    }
  }

  return;
}

template <class T>
T
uqScalarSequenceClass<T>::geweke(
  unsigned int initialPosition,
  double       ratioNa,
  double       ratioNb) const
{
  return;
}

template <class T>
void
uqScalarSequenceClass<T>::minMax(
  unsigned int initialPos,
  T&           minValue,
  T&           maxValue) const
{
  seqScalarPositionConstIteratorTypedef positionIterator = m_seq.begin();
  std::advance(positionIterator,initialPos);

  seqScalarPositionConstIteratorTypedef pos;
  pos = std::min_element(positionIterator, m_seq.end());
  minValue = *pos;
  pos = std::max_element(positionIterator, m_seq.end());
  maxValue = *pos;

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::histogram(
  unsigned int               initialPosition,
  unsigned int               spacing,
  const T&                   minHorizontalValue,
  const T&                   maxHorizontalValue,
  std::vector<T>&            centers,
  std::vector<unsigned int>& bins) const
{
  UQ_FATAL_TEST_MACRO(centers.size() != bins.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarSequenceClass<T>::histogram()",
                      "vectors 'centers' and 'bins' have different sizes");

  UQ_FATAL_TEST_MACRO(bins.size() < 3,
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarSequenceClass<T>::histogram()",
                      "number of 'bins' is too small: should be at least 3");

  for (unsigned int j = 0; j < bins.size(); ++j) {
    centers[j] = 0.;
    bins[j] = 0;
  }

  double horizontalDelta = (maxHorizontalValue - minHorizontalValue)/(((double) bins.size()) - 2.);

  double minCenter = minHorizontalValue - horizontalDelta/2.;
  double maxCenter = maxHorizontalValue + horizontalDelta/2.;
  for (unsigned int j = 0; j < centers.size(); ++j) {
    double factor = ((double) j)/(((double) centers.size()) - 1.);
    centers[j] = (1. - factor) * minCenter + factor * maxCenter;
  }

  unsigned int dataSize = this->sequenceSize();
  for (unsigned int j = 0; j < dataSize; ++j) {
    double value = m_seq[j];
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
#if 0
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
#if 0
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
#endif
#endif // __UQ_SCALAR_SEQUENCE_H__
