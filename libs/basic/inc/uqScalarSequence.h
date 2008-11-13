/* uq/libs/basic/inc/uqScalarSequence.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
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
#include <uqMiscellaneous.h>
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
  uqScalarSequenceClass(const uqBaseEnvironmentClass& env, unsigned int sequenceSize);
 ~uqScalarSequenceClass();

  const unsigned int sequenceSize       () const;
        void         resizeSequence     (unsigned int newSequenceSize);
        void         resetValues        (unsigned int initialPos, unsigned int numPos);
        void         erasePositions     (unsigned int initialPos, unsigned int numPos);
  const T&           operator[]         (unsigned int posId) const;
        T&           operator[]         (unsigned int posId);
        void         setGaussian        (const gsl_rng* rng, const T& mean, const T& stdDev);
        void         setUniform         (const gsl_rng* rng, const T& a,    const T& b     );
        void         uniformlySampledMdf(unsigned int               numIntervals,
                                         T&                         minDomainValue,
                                         T&                         maxDomainValue,
                                         std::vector<T>&            mdfValues) const;
        void         uniformlySampledCdf(unsigned int               numIntervals,
                                         T&                         minDomainValue,
                                         T&                         maxDomainValue,
                                         std::vector<T>&            cdfValues) const;

        T            mean               (unsigned int               initialPos,
                                         unsigned int               numPos) const;
        T            sampleVariance     (unsigned int               initialPos,
                                         unsigned int               numPos,
                                         const T&                   meanValue) const;
        T            populationVariance (unsigned int               initialPos,
                                         unsigned int               numPos,
                                         const T&                   meanValue) const;
        T            autoCovariance     (unsigned int               initialPos,
                                         unsigned int               numPos,
                                         const T&                   meanValue,
                                         unsigned int               lag) const;
        T            autoCorrViaDef     (unsigned int               initialPos,
                                         unsigned int               numPos,
                                         unsigned int               lag) const;
        void         autoCorrViaFft     (unsigned int               initialPos,
                                         unsigned int               numPos,
                                         unsigned int               maxLag,
                                         std::vector<T>&            autoCorrs) const;
        void         autoCorrViaFft     (unsigned int               initialPos,
                                         unsigned int               numPos,
                                         unsigned int               numSum,
                                         T&                         autoCorrsSum) const;
        T            bmm                (unsigned int               initialPos,
                                         unsigned int               batchLength) const;
        void         psd                (unsigned int               initialPos,
                                         unsigned int               numBlocks,
                                         double                     hopSizeRatio,
                                         std::vector<double>&       psdSequence) const;
        T            geweke             (unsigned int               initialPos,
                                         double                     ratioNa,
                                         double                     ratioNb) const;
        void         minMax             (unsigned int               initialPos,
                                         T&                         minValue,
                                         T&                         maxValue) const;
        void         histogram          (unsigned int               initialPos,
                                         const T&                   minHorizontalValue,
                                         const T&                   maxHorizontalValue,
                                         std::vector<T>&            centers,
                                         std::vector<unsigned int>& bins) const;
        void         sort               (unsigned int               initialPos,
                                         uqScalarSequenceClass<T>&  sortedSequence) const;
        void         sort               ();
        T            interQuantileRange (unsigned int               initialPos) const;
        T            scaleForKDE        (unsigned int               initialPos,
                                         const T&                   iqrValue) const;
        double       gaussianKDE        (T                          evaluationParam) const;
        void         gaussianKDE        (unsigned int               initialPos,
                                         double                     scaleValue,
                                         const std::vector<T>&      evaluationParams,
                                         std::vector<double>&       densityValues) const;

private:
        void         extractScalarSeq   (unsigned int               initialPos,
                                         unsigned int               spacing,
                                         unsigned int               numPos,
                                         uqScalarSequenceClass<T>&  scalarSeq) const;
        void         extractRawData     (unsigned int               initialPos,
                                         unsigned int               spacing,
                                         unsigned int               numPos,
                                         std::vector<double>&       rawData) const;

  const uqBaseEnvironmentClass& m_env;
  std::vector<T>            m_seq;
};

template <class T>
uqScalarSequenceClass<T>::uqScalarSequenceClass(
  const uqBaseEnvironmentClass& env,
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
uqScalarSequenceClass<T>::setGaussian(const gsl_rng* rng, const T& meanValue, const T& stdDev)
{
  unsigned int maxJ = this->sequenceSize();
  if (meanValue == 0.) {
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
void
uqScalarSequenceClass<T>::setUniform(const gsl_rng* rng, const T& a, const T& b)
{
  unsigned int maxJ = this->sequenceSize();
  if (a == 0.) {
    if (b == 1.) {
      for (unsigned int j = 0; j < maxJ; ++j) {
        m_seq[j] = gsl_rng_uniform(rng);
      }
    }
    else {
      for (unsigned int j = 0; j < maxJ; ++j) {
        m_seq[j] = b*gsl_rng_uniform(rng);
      }
    }
  }
  else {
    if ((b-a) == 1.) {
      for (unsigned int j = 0; j < maxJ; ++j) {
        m_seq[j] = a + gsl_rng_uniform(rng);
      }
    }
    else {
      for (unsigned int j = 0; j < maxJ; ++j) {
        m_seq[j] = a + (b-a)*gsl_rng_uniform(rng);
      }
    }
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::uniformlySampledMdf(
  unsigned int    numEvaluationPoints,
  T&              minDomainValue,
  T&              maxDomainValue,
  std::vector<T>& mdfValues) const
{
  T                         tmpMinValue;
  T                         tmpMaxValue;
  std::vector<T>            centers(numEvaluationPoints,0.);
  std::vector<unsigned int> bins   (numEvaluationPoints,0.);

  minMax(0, // initialPos
         tmpMinValue,
         tmpMaxValue);
  histogram(0, // initialPos,
            tmpMinValue,
            tmpMaxValue,
            centers,
            bins);

  minDomainValue = centers[0];
  maxDomainValue = centers[centers.size()-1];
  T binSize = (maxDomainValue - minDomainValue)/((double)(centers.size() - 1));

  unsigned int totalCount = 0;
  for (unsigned int i = 0; i < numEvaluationPoints; ++i) {
    totalCount += bins[i];
  }

  mdfValues.clear();
  mdfValues.resize(numEvaluationPoints);
  for (unsigned int i = 0; i < numEvaluationPoints; ++i) {
    mdfValues[i] = ((T) bins[i])/((T) totalCount)/binSize;
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::uniformlySampledCdf(
  unsigned int    numEvaluationPoints,
  T&              minDomainValue,
  T&              maxDomainValue,
  std::vector<T>& cdfValues) const
{
  T                         tmpMinValue;
  T                         tmpMaxValue;
  std::vector<T>            centers(numEvaluationPoints,0.);
  std::vector<unsigned int> bins   (numEvaluationPoints,0.);

  minMax(0, // initialPos
         tmpMinValue,
         tmpMaxValue);
  histogram(0, // initialPos,
            tmpMinValue,
            tmpMaxValue,
            centers,
            bins);

  minDomainValue = centers[0];
  maxDomainValue = centers[centers.size()-1];

  unsigned int totalCount = 0;
  for (unsigned int i = 0; i < numEvaluationPoints; ++i) {
    totalCount += bins[i];
  }

  cdfValues.clear();
  cdfValues.resize(numEvaluationPoints);
  unsigned int partialCount = 0;
  for (unsigned int i = 0; i < numEvaluationPoints; ++i) {
    partialCount += bins[i];
    cdfValues[i] = ((T) partialCount)/((T) totalCount);
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

  samValue /= (((T) numPos) - 1.);

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
  T meanValue = this->mean(initialPos,
                           numPos);
  for (unsigned int j = 0; j < numPos; ++j) {
    rawData[j] -= meanValue; // IMPORTANT
  }

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
  //if (m_env.rank() == 0) {
  //  std::cout << "In uqScalarSequenceClass<T>::autoCorrViaFft()"
  //            << ": returned succesfully from fftObj.inverse()"
  //            << std::endl;
  //}

  // Prepare return data
  autoCorrs.resize(maxLag+1,0.); // Yes, +1
  for (unsigned int j = 0; j < autoCorrs.size(); ++j) {
    double ratio = ((double) j)/((double) (numPos-1));
    autoCorrs[j] = ( resultData[j].real()/resultData[0].real() )*(1.-ratio);
  }

  return;
}

template <class T>
void
uqScalarSequenceClass<T>::autoCorrViaFft(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int numSum,
  T&           autoCorrsSum) const
{
  //if (m_env.rank() == 0) {
  //  std::cout << "Entering uqScalarSequenceClass<T>::autoCorrViaFft(), for sum"
  //            << ": initialPos = " << initialPos
  //            << ", numPos = "     << numPos
  //            << std::endl;
  //}

  double tmp = log((double) numPos)/log(2.);
  double fractionalPart = tmp - ((double) ((unsigned int) tmp));
  if (fractionalPart > 0.) tmp += (1. - fractionalPart);
  unsigned int fftSize = (unsigned int) pow(2.,tmp+1);

  std::vector<double> rawData(numPos,0.);
  std::vector<std::complex<double> > resultData(0,std::complex<double>(0.,0.));
  uqFftClass<T> fftObj(m_env);

  // Forward FFT
  this->extractRawData(initialPos,
                       1, // spacing
                       numPos,
                       rawData);
  T meanValue = this->mean(initialPos,
                           numPos);
  for (unsigned int j = 0; j < numPos; ++j) {
    rawData[j] -= meanValue; // IMPORTANT
  }
  rawData.resize(fftSize,0.);

  //if (m_env.rank() == 0) {
  //  std::cout << "In uqScalarSequenceClass<T>::autoCorrViaFft(), for sum"
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
  fftObj.inverse(rawData,fftSize,resultData);

  //if (m_env.rank() == 0) {
  //  std::cout << "In uqScalarSequenceClass<T>::autoCorrViaFft(), for sum"
  //            << ": computed auto covariance for lag 0 = " << resultData[0].real()/((double) (numPos))
  //            << ", computed resultData[0].imag() = "      << resultData[0].imag()
  //            << std::endl;
  //}

  // Prepare return data
  autoCorrsSum = 0.;
  for (unsigned int j = 0; j < numSum; ++j) { // Yes, begin at lag '0'
    double ratio = ((double) j)/((double) (numPos-1));
    autoCorrsSum += ( resultData[j].real()/resultData[0].real() )*(1.-ratio);
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
  bool bRC = ((initialPos         < this->sequenceSize()                        ) &&
              (hopSizeRatio       != 0.                                         ) &&
              (numBlocks          <          (this->sequenceSize() - initialPos)) &&
              (fabs(hopSizeRatio) < (double) (this->sequenceSize() - initialPos)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<T>::psd()",
                      "invalid input data");

  unsigned int dataSize = this->sequenceSize() - initialPos;

  T meanValue = this->mean(initialPos,
                           dataSize);

  // Determine hopSize and blockSize
  unsigned int hopSize = 0;
  unsigned int blockSize = 0;
  if (hopSizeRatio <= -1.) {
    double overlapSize = -hopSizeRatio;
    double tmp = ((double) (dataSize + (numBlocks - 1)*overlapSize))/((double) numBlocks);
    blockSize = (unsigned int) tmp;
  }
  else if (hopSizeRatio < 0.) {
    double tmp = ((double) dataSize)/(( ((double) numBlocks) - 1. )*(-hopSizeRatio) + 1.);
    blockSize = (unsigned int) tmp;
    hopSize = (unsigned int) ( ((double) blockSize) * (-hopSizeRatio) );
  }
  else if (hopSizeRatio == 0.) {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqScalarSequenceClass<T>::psd()",
                        "hopSizeRatio == 0");
  }
  else if (hopSizeRatio < 1.) {
    double tmp = ((double) dataSize)/(( ((double) numBlocks) - 1. )*hopSizeRatio + 1.);
    blockSize = (unsigned int) tmp;
    hopSize = (unsigned int) ( ((double) blockSize) * hopSizeRatio );
  }
  else {
    hopSize = (unsigned int) hopSizeRatio;
    blockSize = dataSize - (numBlocks - 1)*hopSize;
  }

  int numberOfDiscardedDataElements = dataSize - (numBlocks-1)*hopSize - blockSize;
  //if (m_env.rank() == 0) {
  //  std::cout << "initialPos = "       << initialPos
  //            << ", N = "              << dataSize
  //            << ", #Blocks = "        << numBlocks
  //            << ", R (hop size) = "   << hopSize
  //            << ", B (block size) = " << blockSize
  //            << ", overlap = "        << blockSize - hopSize
  //            << ", [(#Blocks - 1) * R + B] = "       << (numBlocks-1)*hopSize + blockSize
  //            << ", numberOfDiscardedDataElements = " << numberOfDiscardedDataElements
  //            << std::endl;
  //}
  UQ_FATAL_TEST_MACRO(numberOfDiscardedDataElements < 0.,
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarSequenceClass<T>::psd()",
                      "eventual extra space for last block should not be negative");

  double tmp = log((double) blockSize)/log(2.);
  double fractionalPart = tmp - ((double) ((unsigned int) tmp));
  if (fractionalPart > 0.) tmp += (1. - fractionalPart);
  unsigned int fftSize = (unsigned int) pow(2.,tmp);
  //std::cout << "fractionalPart = " << fractionalPart
  //          << ", B = "            << blockSize
  //          << ", fftSize = "      << fftSize
  //          << std::endl;

  double modificationScale = 0.;
  for (unsigned int j = 0; j < blockSize; ++j) {
    double tmpValue = uqMiscHammingWindow(blockSize-1,j);
    modificationScale += tmpValue*tmpValue;
  }
  modificationScale = 1./modificationScale;

  std::vector<double> blockData(blockSize,0.);
  uqFftClass<T> fftObj(m_env);
  std::vector<std::complex<double> > fftResult(0,std::complex<double>(0.,0.));

#if 0
  unsigned int halfFFTSize = fftSize/2;
  psdResult.clear();
  psdResult.resize(1+halfFFTSize,0.);
  double factor = 1./M_PI/((double) numBlocks); // /((double) blockSize);
#else
  psdResult.clear();
  psdResult.resize(fftSize,0.);
  double factor = 1./2./M_PI/((double) numBlocks); // /((double) blockSize);
#endif

  for (unsigned int blockId = 0; blockId < numBlocks; blockId++) {
    // Fill block using Hamming window
    unsigned int initialDataPos = initialPos + blockId*hopSize;
    for (unsigned int j = 0; j < blockSize; ++j) {
      unsigned int dataPos = initialDataPos + j;
      UQ_FATAL_TEST_MACRO(dataPos >= dataSize,
                          UQ_UNAVAILABLE_RANK,
                          "uqScalarSequenceClass<T>::psd()",
                          "too large position to be accessed in data");
      blockData[j] = uqMiscHammingWindow(blockSize-1,j) * ( m_seq[dataPos] - meanValue ); // IMPORTANT
    }

    fftObj.forward(blockData,fftSize,fftResult);

    //if (m_env.rank() == 0) {
    //  std::cout << "blockData.size() = "   << blockData.size()
    //            << ", fftSize = "          << fftSize
    //            << ", fftResult.size() = " << fftResult.size()
    //            << ", psdResult.size() = " << psdResult.size()
    //            << std::endl;
    //}
    // Normalized spectral density: power per radians per sample
    for (unsigned int j = 0; j < psdResult.size(); ++j) {
      psdResult[j] += norm(fftResult[j])*modificationScale;
    }
  }

  for (unsigned int j = 0; j < psdResult.size(); ++j) {
    psdResult[j] *= factor;
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
             (unsigned int) sqrt((double) dataSizeA),  // numBlocks
             .5, // hopSizeRatio
             psdResult);
  double psdA = psdResult[0];
  double varOfMeanA = 2.*M_PI*psdA/doubleDataSizeA;

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
             (unsigned int) sqrt((double) dataSizeB),  // numBlocks
             .5, // hopSizeRatio
             psdResult);
  double psdB = psdResult[0];
  double varOfMeanB = 2.*M_PI*psdB/doubleDataSizeB;

  if (m_env.rank() == 0) {
    std::cout << "In uqScalarSequenceClass<T>::geweke()"
              << ", before computation of gewCoef"
              << ":\n"
              << ", dataSizeA = "       << dataSizeA
              << ", numBlocksA = "      << (unsigned int) sqrt((double) dataSizeA)
              << ", meanA = "           << meanA
              << ", psdA = "            << psdA
              << ", varOfMeanA = "      << varOfMeanA
              << "\n"
              << ", dataSizeB = "       << dataSizeB
              << ", numBlocksB = "      << (unsigned int) sqrt((double) dataSizeB)
              << ", meanB = "           << meanB
              << ", psdB = "            << psdB
              << ", varOfMeanB = "      << varOfMeanB
              << std::endl;
  }
  double gewCoef = (meanA - meanB)/sqrt(varOfMeanA + varOfMeanB);

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
  sortedSequence.resizeSequence(numPos);
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
  T iqrValue = value3 - value1;

  return iqrValue;
}

template <class T>
T
uqScalarSequenceClass<T>::scaleForKDE(
  unsigned int initialPos,
  const T&     iqrValue) const
{
  bool bRC = (initialPos <  this->sequenceSize());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<V>::scalesForKDE()",
                      "invalid input data");

  unsigned int dataSize = this->sequenceSize() - initialPos;

  T meanValue = this->mean(initialPos,
                           dataSize);

  T samValue = this->sampleVariance(initialPos,
                                    dataSize,
                                    meanValue);

  T scaleValue;
  if (iqrValue <= 0.) {
    scaleValue = 1.06*sqrt(samValue)/pow(dataSize,1./5.);
   }
  else {
    scaleValue = 1.06*std::min(sqrt(samValue),iqrValue/1.34)/pow(dataSize,1./5.);
  }

  return scaleValue;
}

template <class T>
double
uqScalarSequenceClass<T>::gaussianKDE(T evaluationParam) const
{
  T iqrValue = this->interQuantileRange(0); // Use the whole chain

  T scaleValue = this->scaleForKDE(0, // Use the whole chain
                                   iqrValue);

  unsigned int dataSize = this->sequenceSize();

  double scaleInv = 1./scaleValue;
  double x = evaluationParam;
  double value = 0.;
  for (unsigned int k = 0; k < dataSize; ++k) {
    double xk = m_seq[k];
    value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
  }

  return scaleInv * (value/(double) dataSize);
}

template <class T>
void
uqScalarSequenceClass<T>::gaussianKDE(
  unsigned int          initialPos,
  double                scaleValue,
  const std::vector<T>& evaluationParams,
  std::vector<double>&  densityValues) const
{
  bool bRC = ((initialPos                 <  this->sequenceSize()   ) &&
              (0                          <  evaluationParams.size()) &&
              (evaluationParams.size() == densityValues.size()      ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqScalarSequenceClass<V>::gaussianKDE()",
                      "invalid input data");

  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numEvals = evaluationParams.size();

  double scaleInv = 1./scaleValue;
  for (unsigned int j = 0; j < numEvals; ++j) {
    double x = evaluationParams[j];
    double value = 0.;
    for (unsigned int k = 0; k < dataSize; ++k) {
      double xk = m_seq[initialPos+k];
      value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
    }
    densityValues[j] = scaleInv * (value/(double) dataSize);
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

#endif // __UQ_SCALAR_SEQUENCE_H__
