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
        void         resetValues       (unsigned int initialPos, unsigned int numPos);
        void         erasePositions    (unsigned int initialPos, unsigned int numPos);
  const T&           operator[]        (unsigned int posId) const;
        T&           operator[]        (unsigned int posId);
        void         setGaussian       (gsl_rng* rng, const T& mean, const T& stdDev);

        T            mean              (unsigned int               initialPos,
                                        unsigned int               numPos) const;
        T            sampleVariance    (unsigned int               initialPos,
                                        unsigned int               numPos,
                                        const T&                   meanValue) const;
        T            populationVariance(unsigned int               initialPos,
                                        unsigned int               numPos,
                                        const T&                   meanValue) const;
        T            autoCovariance    (unsigned int               initialPos,
                                        unsigned int               numPos,
                                        const T&                   meanValue,
                                        unsigned int               lag) const;
        T            autoCorrViaDef    (unsigned int               initialPos,
                                        unsigned int               numPos,
                                        unsigned int               lag) const;
        void         autoCorrViaFft    (unsigned int               initialPos,
                                        unsigned int               numPos,
                                        unsigned int               maxLag,
                                        std::vector<T>&            autoCorrs) const;
        T            bmm               (unsigned int               initialPos,
                                        unsigned int               batchLength) const;
        void         psd               (unsigned int               initialPos,
                                        unsigned int               numBlocks,
                                        double                     hopSizeRatio,
                                        std::vector<double>&       psdSequence) const;
        T            geweke            (unsigned int               initialPos,
                                        double                     ratioNa,
                                        double                     ratioNb) const;
        void         minMax            (unsigned int               initialPos,
                                        T&                         minValue,
                                        T&                         maxValue) const;
        void         histogram         (unsigned int               initialPos,
                                        const T&                   minHorizontalValue,
                                        const T&                   maxHorizontalValue,
                                        std::vector<T>&            centers,
                                        std::vector<unsigned int>& bins) const;
        void         sort              (unsigned int               initialPos,
                                        uqScalarSequenceClass<T>&  sortedSequence) const;
        void         sort              ();
        T            interQuantileRange(unsigned int               initialPos) const;
        T            scaleForKDE       (unsigned int               initialPos,
                                        const T&                   iqrsValue) const;
        void         gaussianKDE       (unsigned int               initialPos,
                                        const std::vector<T>&      evaluationPositions,
                                        double                     scaleValue,
                                        std::vector<double>&       densityValues) const;

private:
        void         extractScalarSeq  (unsigned int               initialPos,
                                        unsigned int               spacing,
                                        unsigned int               numPos,
                                        uqScalarSequenceClass<T>&  scalarSeq) const;
        void         extractRawData    (unsigned int               initialPos,
                                        unsigned int               spacing,
                                        unsigned int               numPos,
                                        std::vector<double>&       rawData) const;

  const uqEnvironmentClass& m_env;
  std::vector<T>            m_seq;
};

template <class T>
uqScalarSequenceClass<T>::uqScalarSequenceClass(
  const uqEnvironmentClass& env,
        unsigned int        sequenceSize)
  :
  m_env(env),
  m_seq(sequenceSize,0.)
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
  if (newSequenceSize != this->sequenceSize()) {
    m_seq.resize(newSequenceSize,0.);
    std::vector<T>(m_seq).swap(m_seq);
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::resetValues(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequences<T>::resetValues()",
                      "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    m_seq[initialPos+j] = 0.;
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::erasePositions(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequences<T>::erasePositions()",
                      "invalid input data");

  seqScalarPositionIteratorTypedef posIteratorBegin = m_seq.begin();
  if (initialPos < this->sequenceSize()) std::advance(posIteratorBegin,initialPos);
  else                                   posIteratorBegin = m_seq.end();

  unsigned int posEnd = initialPos + numPos - 1;
  seqScalarPositionIteratorTypedef posIteratorEnd = m_seq.begin();
  if (posEnd < this->sequenceSize()) std::advance(posIteratorEnd,posEnd);
  else                               posIteratorEnd = m_seq.end();

  unsigned int oldSequenceSize = this->sequenceSize();
  m_seq.erase(posIteratorBegin,posIteratorEnd);
  UQ_FATAL_TEST_MACRO((oldSequenceSize - numPos) != this->sequenceSize(),
                      m_env.rank(),
                      "uqScalarSequences<T>::erasePositions()",
                      "(oldSequenceSize - numPos) != this->sequenceSize()");

  return;
}

template <class T>
const T&
uqScalarSequenceClass<T>::operator[](unsigned int posId) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<T>::operator[] const",
                      "posId > sequenceSize()");

  return m_seq[posId];
}

template <class T>
T&
uqScalarSequenceClass<T>::operator[](unsigned int posId)
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<T>::operator[]",
                      "posId > sequenceSize()");

  return m_seq[posId];
}

template <class T>
void
uqScalarSequenceClass<T>::setGaussian(gsl_rng* rng, const T& meanValue, const T& stdDev)
{
  unsigned int maxJ = this->sequenceSize();
  if (meanValue == 0) {
    for (unsigned int j = 0; j < maxJ; ++j) {
      m_seq[j] = gsl_ran_gaussian(rng,stdDev);
    }
  }
  else {
    for (unsigned int j = 0; j < maxJ; ++j) {
      m_seq[j] = meanValue + gsl_ran_gaussian(rng,stdDev);
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
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::mean()",
                      "invalid input data");

  unsigned int finalPosPlus1 = initialPos + numPos;
  T tmpSum = 0.;
  for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
    tmpSum += m_seq[j];
  }

  return tmpSum/(T) numPos;
}

template <class T>
T
uqScalarSequenceClass<T>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const T&     meanValue) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::sampleVariance()",
                      "invalid input data");

  unsigned int finalPosPlus1 = initialPos + numPos;
  T diff;
  T samValue = 0.;
  for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
    diff = m_seq[j] - meanValue;
    samValue += diff*diff;
  }

  T doubleNumPos = (T) numPos;
  samValue /= (doubleNumPos - 1.);

  return samValue;
}

template <class T>
T
uqScalarSequenceClass<T>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const T&     meanValue) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::populationVariance()",
                      "invalid input data");

  unsigned int finalPosPlus1 = initialPos + numPos;
  T diff;
  T popValue = 0.;
  for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
    diff = m_seq[j] - meanValue;
    popValue += diff*diff;
  }

  popValue /= (T) numPos;

  return popValue;
}

template <class T>
T
uqScalarSequenceClass<T>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const T&     meanValue,
  unsigned int lag) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (lag                 <  numPos              )); // lag should not be too large
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::autoCovariance()",
                      "invalid input data");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  T diff1;
  T diff2;
  T covValue = 0.;
  for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
    diff1 = m_seq[j    ] - meanValue;
    diff2 = m_seq[j+lag] - meanValue;
    covValue += diff1*diff2;
  }

  covValue /= (T) loopSize;

  return covValue;
}

template <class T>
T
uqScalarSequenceClass<T>::autoCorrViaDef(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int lag) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (lag                 <  numPos              )); // lag should not be too large
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::autoCorrViaDef()",
                      "invalid input data");

  T meanValue = this->mean(initialPos,
                           numPos);

  T covValueZero = this->autoCovariance(initialPos,
                                        numPos,
                                        meanValue,
                                        0); // lag

  T corrValue = this->autoCovariance(initialPos,
                                     numPos,
                                     meanValue,
                                     lag);

  return corrValue/covValueZero;
}

template <class T>
void
uqScalarSequenceClass<T>::autoCorrViaFft(
  unsigned int    initialPos,
  unsigned int    numPos,
  unsigned int    maxLag,
  std::vector<T>& autoCorrs) const
{
  unsigned int fftSize = 0;
  {
    double tmp = log((double) maxLag)/log(2.);
    double fractionalPart = tmp - ((double) ((unsigned int) tmp));
    if (fractionalPart > 0.) tmp += (1. - fractionalPart);
    unsigned int fftSize1 = (unsigned int) pow(2.,tmp+1.); // Yes, tmp+1
    fftSize1 = fftSize1; // To remove warning

    tmp = log((double) numPos)/log(2.);
    fractionalPart = tmp - ((double) ((unsigned int) tmp));
    if (fractionalPart > 0.) tmp += (1. - fractionalPart);
    unsigned int fftSize2 = (unsigned int) pow(2.,tmp+1);

    fftSize = fftSize2;
  }

  std::vector<double> rawData(numPos,0.);
  std::vector<std::complex<double> > resultData(0,std::complex<double>(0.,0.));
  uqFftClass<T> fftObj(m_env);

  // Forward FFT
  this->extractRawData(initialPos,
                       1, // spacing
                       numPos,
                       rawData);
  rawData.resize(fftSize,0.);

  //if (m_env.rank() == 0) {
  //  std::cout << "In uqScalarSequenceClass<T>::autoCorrViaFft()"
  //            << ": about to call fftObj.forward()"
  //            << " with rawData.size() = " << rawData.size()
  //            << ", fftSize = "            << fftSize
  //            << ", resultData.size() = "  << resultData.size()
  //            << std::endl;
  //}
  fftObj.forward(rawData,fftSize,resultData);

  // Inverse FFT
  for (unsigned int j = 0; j < fftSize; ++j) {
    rawData[j] = std::norm(resultData[j]);
  }
  //if (m_env.rank() == 0) {
  //  std::cout << "In uqScalarSequenceClass<T>::autoCorrViaFft()"
  //            << ": about to call fftObj.inverse()"
  //            << " with rawData.size() = " << rawData.size()
  //            << ", fftSize = "            << fftSize
  //            << ", resultData.size() = "  << resultData.size()
  //            << std::endl;
  //}
  fftObj.inverse(rawData,fftSize,resultData);

  if (m_env.rank() == 0) {
    std::cout << "In uqScalarSequenceClass<T>::autoCorrViaFft()"
              << ": returned succesfully from fftObj.inverse()"
              << std::endl;
  }

  // Prepare return data
  autoCorrs.resize(maxLag+1,0.); // Yes, +1
  for (unsigned int j = 0; j < autoCorrs.size(); ++j) {
    autoCorrs[j] = resultData[j].real()/((double) (numPos-j));
  }

  return;
}

template <class T>
T
uqScalarSequenceClass<T>::bmm(
  unsigned int initialPos,
  unsigned int batchLength) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()            ) &&
              (batchLength         < (this->sequenceSize()-initialPos)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequences<T>::bmm()",
                      "invalid input data");

  unsigned int numberOfBatches = (this->sequenceSize() - initialPos)/batchLength;
  uqScalarSequenceClass<T> batchMeans(m_env,numberOfBatches);

  for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
    batchMeans[batchId] = this->mean(initialPos + batchId*batchLength,
                                     batchLength);
  }

  T meanOfBatchMeans = batchMeans.mean(0,
                                       batchMeans.sequenceSize());

  //T covLag0OfBatchMeans = batchMeans.autoCovariance(0,
  //                                                  batchMeans.sequenceSize(),
  //                                                  meanOfBatchMeans,
  //                                                  0); // lag

  T bmmValue = batchMeans.sampleVariance(0,
                                         batchMeans.sequenceSize(),
                                         meanOfBatchMeans);

  bmmValue /= (T) batchMeans.sequenceSize();           // CHECK
//bmmValue *= (T) (this->sequenceSize() - initialPos); // CHECK

  return bmmValue;
}

template <class T>
void
uqScalarSequenceClass<T>::psd(
  unsigned int         initialPos,
  unsigned int         numBlocks,
  double               hopSizeRatio,
  std::vector<double>& psdResult) const
{
  bool bRC = (initialPos < this->sequenceSize());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::psd()",
                      "invalid input data");

  unsigned int dataSize = this->sequenceSize() - initialPos;

  double tmp = ((double) dataSize)/(( ((double) numBlocks) - 1. )*hopSizeRatio + 1.);
  unsigned int blockSize = (unsigned int) tmp;
  unsigned int hopSize   = (unsigned int) ( ((double) blockSize) * hopSizeRatio );
  tmp = ((double) dataSize) - ( ((double) numBlocks) - 1.) * ((double) hopSize) - ((double) blockSize);
#if 0
  unsigned int numberOfDiscardedDataElements = (unsigned int) tmp;
  std::cout << "initialPos = " << initialPos
            << ", N = "        << dataSize
            << ", #Blocks = "  << numBlocks
            << ", R = "        << hopSize
            << ", B = "        << blockSize
            << ", overlap = "  << blockSize - hopSize
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

  uqFftClass<T> fftObj(m_env);
  std::vector<std::complex<double> > resultData(0,std::complex<double>(0.,0.));

  unsigned int halfFFTSize = fftSize/2;
  psdResult.clear();
  psdResult.resize(1+halfFFTSize,0.);
  for (unsigned int blockId = 0; blockId < numBlocks; blockId++) {
    // Padding
    std::vector<double> blockData(fftSize,0.);

    // Fill block using window on full data
    unsigned int initialDataPos = initialPos + blockId*hopSize;
    for (unsigned int j = 0; j < blockSize; ++j) {
      unsigned int dataPos = initialDataPos + j;
      UQ_FATAL_TEST_MACRO(dataPos >= dataSize,
                          UQ_UNAVAILABLE_RANK,
                          "uqScalarSequenceClass<T>::psd()",
                          "too large position to be accessed in data");
      blockData[j] = uqMiscHammingWindow(fftSize-1,j) * m_seq[dataPos];
    }

    fftObj.forward(blockData,fftSize,resultData);

    // Normalized spectral density: power per radians per sample
    double factor = 1./((double) numBlocks*blockSize); // /M_PI; // CHECK
    for (unsigned int j = 0; j < psdResult.size(); ++j) {
      psdResult[j] += norm(resultData[j]) * factor;
    }
  }

  return;
}

template <class T>
T
uqScalarSequenceClass<T>::geweke(
  unsigned int initialPos,
  double       ratioNa,
  double       ratioNb) const
{
  double doubleFullDataSize = (double) (this->sequenceSize() - initialPos);
  uqScalarSequenceClass<T> tmpSeq(m_env,0);
  std::vector<double> psdResult(0,0.);

  unsigned int dataSizeA       = (unsigned int) (doubleFullDataSize * ratioNa);
  double       doubleDataSizeA = (double) dataSizeA;
  unsigned int initialPosA     = initialPos;
  this->extractScalarSeq(initialPosA,
                         1,
                         dataSizeA,
                         tmpSeq);
  double meanA = tmpSeq.mean(0,
                             dataSizeA);
  tmpSeq.psd(0,
             8,  // numBlocks
             .5, // hopSizeRatio
             psdResult);
  double psdA = psdResult[0];

  unsigned int dataSizeB       = (unsigned int) (doubleFullDataSize * ratioNb);
  double       doubleDataSizeB = (double) dataSizeB;
  unsigned int initialPosB     = this->sequenceSize() - dataSizeB;
  this->extractScalarSeq(initialPosB,
                         1,
                         dataSizeB,
                         tmpSeq);
  double meanB = tmpSeq.mean(0,
                             dataSizeB);
  tmpSeq.psd(0,
             8,  // numBlocks
             .5, // hopSizeRatio
             psdResult);
  double psdB = psdResult[0];

#if 0
  if (m_env.rank() == 0) {
    std::cout << "In uqScalarSequenceClass<T>::geweke()"
              << ", before computation of gewCoef"
              << ": meanA = "           << meanA
              << ", psdA = "            << psdA
              << ", doubleDataSizeA = " << doubleDataSizeA
              << ", meanB = "           << meanB
              << ", psdB = "            << psdB
              << ", doubleDataSizeB = " << doubleDataSizeB
              << std::endl;
  }
#endif
  double gewCoef = (meanA - meanB)/sqrt(psdA/doubleDataSizeA + psdB/doubleDataSizeB);

  return gewCoef;
}

template <class T>
void
uqScalarSequenceClass<T>::minMax(
  unsigned int initialPos,
  T&           minValue,
  T&           maxValue) const
{
  seqScalarPositionConstIteratorTypedef posIterator = m_seq.begin();
  std::advance(posIterator,initialPos);

  seqScalarPositionConstIteratorTypedef pos;
  pos = std::min_element(posIterator, m_seq.end());
  minValue = *pos;
  pos = std::max_element(posIterator, m_seq.end());
  maxValue = *pos;

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::histogram(
  unsigned int               initialPos,
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

template <class T>
void
uqScalarSequenceClass<T>::sort(
  unsigned int              initialPos,
  uqScalarSequenceClass<T>& sortedSequence) const
{
  unsigned int numPos = this->sequenceSize() - initialPos;
  sortedSequence.resize(numPos);
  this->extractScalarSeq(initialPos,
                         1,
                         numPos,
                         sortedSequence);
  sortedSequence.sort();

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::sort()
{
  std::sort(m_seq.begin(), m_seq.end());
  return;
}

template <class T>
T
uqScalarSequenceClass<T>::interQuantileRange(unsigned int initialPos) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;

  uqScalarSequenceClass sortedSequence(m_env,0);
  this->sort(initialPos,
             sortedSequence);

  unsigned int pos1 = (unsigned int) ( (((double) dataSize) + 1.)*1./4. - 1. );
  unsigned int pos3 = (unsigned int) ( (((double) dataSize) + 1.)*3./4. - 1. );

  double fraction1 = (((double) dataSize) + 1.)*1./4. - 1. - ((double) pos1);
  double fraction3 = (((double) dataSize) + 1.)*3./4. - 1. - ((double) pos3);

  //std::cout << "In uqScalarSequenceClass::interQuantileRange()"
  //          << ", initialPos = "            << initialPos
  //          << ", this->sequenceSize() = "  << this->sequenceSize()
  //          << ", dataSize = "              << dataSize
  //          << ", sortedSequence.size() = " << sortedSequence.size()
  //          << ", pos1 = "                  << pos1
  //          << ", pos3 = "                  << pos3
  //          << std::endl;
  T value1 = (1.-fraction1) * sortedSequence[pos1] + fraction1 * sortedSequence[pos1+1];
  T value3 = (1.-fraction3) * sortedSequence[pos3] + fraction3 * sortedSequence[pos3+1];
  T iqrsValue = value3 - value1;

  return iqrsValue;
}

template <class T>
T
uqScalarSequenceClass<T>::scaleForKDE(
  unsigned int initialPos,
  const T&     iqrsValue) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;

  T meanValue = this->mean(initialPos,
                           dataSize);

  T samValue = this->sampleVariance(initialPos,
                                    dataSize,
                                    meanValue);

  T scaleValue;
  if (iqrsValue <= 0.) {
    scaleValue = 1.06*sqrt(samValue)/pow(dataSize,1./5.);
   }
  else {
    scaleValue = 1.06*std::min(sqrt(samValue),iqrsValue/1.34)/pow(dataSize,1./5.);
  }

  return scaleValue;
}

template <class T>
void
uqScalarSequenceClass<T>::gaussianKDE(
  unsigned int          initialPos,
  const std::vector<T>& evaluationPositions,
  double                scaleValue,
  std::vector<double>&  densityValues) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numEstimations = evaluationPositions.size();

  double scaleInv = 1./scaleValue;
  for (unsigned int j = 0; j < numEstimations; ++j) {
    double x = evaluationPositions[j];
    double value = 0.;
    for (unsigned int k = 0; k < dataSize; ++k) {
      double xk = m_seq[initialPos+k];
      value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
    }
    densityValues[j] = scaleInv * (value/(double) numEstimations);
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::extractScalarSeq(
  unsigned int              initialPos,
  unsigned int              spacing,
  unsigned int              numPos,
  uqScalarSequenceClass<T>& scalarSeq) const
{
  scalarSeq.resizeSequence(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = m_seq[initialPos+j        ];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = m_seq[initialPos+j*spacing];
    }
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::extractRawData(
  unsigned int         initialPos,
  unsigned int         spacing,
  unsigned int         numPos,
  std::vector<double>& rawData) const
{
  rawData.resize(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = m_seq[initialPos+j        ];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = m_seq[initialPos+j*spacing];
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
  std::vector<double>&       psdResult) // [dataSize x 1] vector
{
  psdResult.clear();
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
  psdResult.resize(1+halfFFTSize,0.);
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
    for (unsigned int j = 0; j < psdResult.size(); ++j) {
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
      psdResult[j] += (realPartOfFFT*realPartOfFFT + imagPartOfFFT*imagPartOfFFT) * factor;
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
