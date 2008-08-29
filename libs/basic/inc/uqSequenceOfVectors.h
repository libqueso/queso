/* uq/libs/basic/inc/uqSequenceOfVectors.h
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

#ifndef __UQ_SEQUENCE_OF_VECTORS_H__
#define __UQ_SEQUENCE_OF_VECTORS_H__

#include <EpetraExt_DistArray.h>
#include <uq2dArrayOfStuff.h>
#include <uqScalarSequence.h>

template <class V>
class uqSequenceOfVectorsClass
{
public:
  typedef typename std::vector<V*>::iterator seqVectorPositionIteratorTypedef;
  uqSequenceOfVectorsClass(unsigned int sequenceSize, const V& vectorExample);
 ~uqSequenceOfVectorsClass();

  const unsigned int sequenceSize      () const;
  const unsigned int vectorSize        () const;
        void         resizeSequence    (unsigned int newSequenceSize);
        void         resetValues       ();
        void         erasePositions    (unsigned int posBegin, unsigned int posEnd);
  const V*           operator[]        (unsigned int positionId) const;
  const V*&          operator[]        (unsigned int positionId);
        void         setGaussian       (gsl_rng* rng, const V& meanVec, const V& stdDevVec);
        void         mean              (unsigned int initialPos,
                                        unsigned int numPos,
                                        V&           mean) const;
        void         sampleVariance    (unsigned int initialPos,
                                        unsigned int numPos,
                                        const V&     mean,
                                        V&           sampleVariance) const;
        void         populationVariance(unsigned int initialPos,
                                        unsigned int numPos,
                                        const V&     mean,
                                        V&           populVariance) const;

        void         autoCovariance    (unsigned int initialPos,
                                        unsigned int numPos,
                                        const V&     mean,
                                        unsigned int lag,
                                        V&           autoCov) const;
        void         autoCorrelations  (const std::vector<unsigned int>& initialPositions,
                                        const std::vector<unsigned int>& lags,
                                        uq2dArrayOfStuff<V>&             _2dArrayOfAutoCorrs) const; // [numOfPos x numOfLags] matrix
        void         bmm               (const std::vector<unsigned int>& initialPositions,
                                        const std::vector<unsigned int>& batchLengths,
                                        uq2dArrayOfStuff<V>&             _2dArrayOfBMM) const; // [numOfPos x numOfLengths] matrix
        //void         fftAlloc          ();
        //void         fftForward        (unsigned int fftSize);
        //void         fftInverse        (unsigned int fftSize);
        //void         fftFree           ();
        void         psdAtZero         (const std::vector<unsigned int>& initialPositions,
                                        const std::vector<unsigned int>& numsOfBlocks,
                                        double                           hopSizeRatio,
                                        uq2dArrayOfStuff<V>&             _2dArrayOfPSDAtZero) const; // [numOfPos x numOfBlocks] matrix
        void         geweke            (const std::vector<unsigned int>& initialPositions,
                                        double                           ratioNa,
                                        double                           ratioNb,
                                        std::vector<V*>&                 vectorOfGeweke) const;
        void         minMax            (unsigned int                     initialPos,
                                        V&                               mins,
                                        V&                               maxs) const;
        void         histogram         (unsigned int                     initialPosition,
                                        unsigned int                     spacing,
                                        const V&                         minHorizontalValues,
                                        const V&                         maxHorizontalValues,
                                        std::vector<V*>&                 centersForAllBins,
                                        std::vector<V*>&                 binsForAllParams) const;
        void         sort              (unsigned int                     initialPosition,
                                        std::vector<V*>&                 sortedSequence) const;
        void         interQuantileRange(unsigned int                     initialPosition,
                                        unsigned int                     spacing,
                                        V&                               iqrs) const;
        void         scalesForKDE      (unsigned int                     initialPosition,
                                        unsigned int                     spacing,
                                        const V&                         iqrs,
                                        V&                               scales) const;
        void         gaussianKDE       (unsigned int                     initialPosition,
                                        unsigned int                     spacing,
                                        const std::vector<V*>&           evaluationPositions,
                                        const V&                         scales,
                                        std::vector<V*>&                 densityValues) const;
private:
  const uqEnvironmentClass&  m_env;
  V                          m_vectorExample;
  std::vector<const V*>      m_seq;
};

template <class V>
uqSequenceOfVectorsClass<V>::uqSequenceOfVectorsClass(
  unsigned int sequenceSize,
  const V&     vectorExample)
  :
  m_env          (vectorExample.env()),
  m_vectorExample(vectorExample),
  m_seq          (sequenceSize,NULL)
{

  //if (m_env.rank() == 0) std::cout << "Entering uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << std::endl;

  //if (m_env.rank() == 0) std::cout << "Leaving uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << std::endl;
}

template <class V>
uqSequenceOfVectorsClass<V>::~uqSequenceOfVectorsClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_seq.size(); ++i) {
    if (m_seq[i]) delete m_seq[i];
  }
}

template <class V>
const unsigned int
uqSequenceOfVectorsClass<V>::sequenceSize() const
{
  return m_seq.size();
}

template <class V>
const unsigned int
uqSequenceOfVectorsClass<V>::vectorSize() const
{
  return m_vectorExample.size();
}

template <class V>
void
uqSequenceOfVectorsClass<V>::resizeSequence(unsigned int newSequenceSize)
{
  m_seq.resize(newSequenceSize,NULL);
  std::vector<const V*>(m_seq).swap(m_seq);

 return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::resetValues()
{
  exit(1);
  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::erasePositions(unsigned int posBegin, unsigned int posEnd)
{
  if (posBegin < this->sequenceSize()) {
    exit(1);
#if 0
    for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
      seqVectorPositionIteratorTypedef positionIteratorBegin = m_seq.begin();
      if (posBegin < sequenceSize()) std::advance(positionIteratorBegin,posBegin);
      else                           positionIteratorBegin = m_seq.end();

      seqVectorPositionIteratorTypedef positionIteratorEnd = m_seq.begin();
      if (posEnd < sequenceSize()) std::advance(positionIteratorEnd,posEnd);
      else                         positionIteratorEnd = m_seq.end();

      m_seq.erase(positionIteratorBegin,positionIteratorEnd);
      UQ_FATAL_TEST_MACRO((posBegin != sequenceSize()),
                          m_env.rank(),
                          "uqSequenceOfVectors<V>::erasePositions()",
                          "posBegin != sequenceSize()");
    }
#endif
  }

  return;
}

template <class V>
const V*
uqSequenceOfVectorsClass<V>::operator[](unsigned int positionId) const
{
  return (const V*) (m_seq[positionId]);
}

template <class V>
const V*&
uqSequenceOfVectorsClass<V>::operator[](unsigned int positionId)
{
  return m_seq[positionId];
}

template <class V>
void
uqSequenceOfVectorsClass<V>::setGaussian(gsl_rng* rng, const V& meanVec, const V& stdDevVec)
{
  exit(1);
  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           mean) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::mean()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::mean()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  mean.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < mean.size(); ++j) {
      mean[j] += (*m_seq[i])[j]/doubleLoopSize;
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     mean,
  V&           sampleVariance) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == sampleVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "incompatible sizes between sampleVariance vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  sampleVariance.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < sampleVariance.size(); ++j) {
      double diff = (*m_seq[i])[j] - mean[j];
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
uqSequenceOfVectorsClass<V>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     mean,
  V&           populVariance) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == populVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "incompatible sizes between populVariance vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  populVariance.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < populVariance.size(); ++j) {
      double diff = (*m_seq[i])[j] - mean[j];
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
uqSequenceOfVectorsClass<V>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     mean,
  unsigned int lag,
  V&           autoCov) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectors<V>::autoCovariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (numPos > lag);
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectors<V>::autoCovariance<V>()",
                      "lag is too large");

  bRC = (this->vectorSize() == autoCov.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectors<V>::autoCovariance<V>()",
                      "incompatible sizes between autoCov vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectors<V>::autoCovariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  autoCov.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < autoCov.size(); ++j) {
      double diff1 = (*m_seq[i    ])[j] - mean[j];
      double diff2 = (*m_seq[i+lag])[j] - mean[j];
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
uqSequenceOfVectorsClass<V>::autoCorrelations(
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& lags,
  uq2dArrayOfStuff<V>&             _2dArrayOfAutoCorrs) const // [numOfPos x numOfLags] matrix
{
  V subChainMean              (m_vectorExample);
  V subChainAutoCovarianceLag0(m_vectorExample);

  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    this->mean(initialPositions[initialPosId],
               this->sequenceSize()-initialPositions[initialPosId],
               subChainMean);
    this->autoCovariance(initialPositions[initialPosId],
                         this->sequenceSize()-initialPositions[initialPosId],
                         subChainMean,
                         0, // lag
                         subChainAutoCovarianceLag0);
    for (unsigned int lagId = 0; lagId < lags.size(); lagId++) {
      this->autoCovariance(initialPositions[initialPosId],
                           this->sequenceSize()-initialPositions[initialPosId],
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
uqSequenceOfVectorsClass<V>::bmm(
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& batchLengths,
  uq2dArrayOfStuff<V>&             _2dArrayOfBMM) const // [numOfPos x numOfLengths] matrix
{
  V meanOfBatchMeans   (m_vectorExample);
  V covLag0OfBatchMeans(m_vectorExample);
  V covLag1OfBatchMeans(m_vectorExample);

  V tmpVector(m_vectorExample); // In order to contour the fact that 'batchMeans' is a vector of 'const V*', but needs to be set first
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    for (unsigned int batchLengthId = 0; batchLengthId < batchLengths.size(); batchLengthId++) {
      unsigned int batchLength = batchLengths[batchLengthId];
      unsigned int numberOfBatches = (this->sequenceSize() - initialPositions[initialPosId])/batchLength;

      uqSequenceOfVectorsClass batchMeans(numberOfBatches,m_vectorExample);
      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        this->mean(initialPositions[initialPosId] + batchId*batchLength,
                   batchLength,
                   tmpVector);
        batchMeans[batchId] = new V(tmpVector);
      }

      batchMeans.mean(0,
                      batchMeans.sequenceSize(),
                      meanOfBatchMeans);

      batchMeans.autoCovariance(0,
                                batchMeans.sequenceSize(),
                                meanOfBatchMeans,
                                0, // lag
                                covLag0OfBatchMeans);

      batchMeans.autoCovariance(0,
                                batchMeans.sequenceSize(),
                                meanOfBatchMeans,
                                1, // lag
                                covLag0OfBatchMeans);

      batchMeans.sampleVariance(0,
                                batchMeans.sequenceSize(),
                                meanOfBatchMeans,
                                _2dArrayOfBMM(initialPosId,batchLengthId));

      _2dArrayOfBMM(initialPosId,batchLengthId) /= (double) batchMeans.sequenceSize(); // CHECK
      //_2dArrayOfBMM(initialPosId,batchLengthId) *= (double) (this->sequenceSize() - initialPositions[initialPosId]); // CHECK
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::psdAtZero(
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& numsOfBlocks,
  double                           hopSizeRatio,
  uq2dArrayOfStuff<V>&             _2dArrayOfPSDAtZero) const // [numOfPos x numOfBlocks] matrix
{
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    unsigned int dataSize = this->sequenceSize() - initialPositions[initialPosId];
    uqScalarSequenceClass<double> data(m_env,dataSize);

    unsigned int numParams = vectorSize();
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSize; ++j) {
	data[j] = (*(m_seq[initialPositions[initialPosId]+j]))[i];
      }
      for (unsigned int numsOfBlocksId = 0; numsOfBlocksId < numsOfBlocks.size(); numsOfBlocksId++) {
        unsigned int numBlocks = numsOfBlocks[numsOfBlocksId];
        std::vector<double> psdData(0,0.); // size will be determined by 'uqScalarSequencePSD()'
        data.psd(0,
                 numBlocks,
                 hopSizeRatio,
                 psdData);
        _2dArrayOfPSDAtZero(initialPosId,numsOfBlocksId)[i] = psdData[0];

	std::cout << "psdData[0] = " << psdData[0] << std::endl;

        //std::cout << "psdData = zeros(" << psdData.size() << ",1);" << std::endl;
        //for (unsigned j = 0; j < psdData.size(); ++j) {
    	//  std::cout << "psdData(" << j+1 << ") = " << psdData[j] << ";" << std::endl;
        //}
      } // for 'numsOfBlocksId'
    } // for 'i'
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::geweke(
  const std::vector<unsigned int>& initialPositions,
  double                           ratioNa,
  double                           ratioNb,
  std::vector<V*>&                 vectorOfGeweke) const
{
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    unsigned int fullDataSize = this->sequenceSize() - initialPositions[initialPosId];
    unsigned int dataSizeA    = (unsigned int) (((double) fullDataSize) * ratioNa);
    unsigned int dataSizeB    = (unsigned int) (((double) fullDataSize) * ratioNb);
    unsigned int initialPosA  = initialPositions[initialPosId];
    unsigned int initialPosB  = this->sequenceSize() - dataSizeB;

    V meanA(m_vectorExample);
    this->mean(initialPosA,
               dataSizeA,
               meanA);

    V meanB(m_vectorExample);
    this->mean(initialPosB,
               dataSizeB,
               meanB);

    unsigned int numParams = vectorSize();

    V psdAtZeroA(m_vectorExample);
    uqScalarSequenceClass<double> dataA(m_env,dataSizeA);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeA; ++j) {
	dataA[j] = (*(m_seq[initialPosA+j]))[i];
      }
      std::vector<double> psdData(0,0.);
      dataA.psd(0,
                8,  // numBlocks
                .5, // hopSizeRatio
                psdData);
      psdAtZeroA[i] = psdData[0];
    } // for 'i'

    V psdAtZeroB(m_vectorExample);
    uqScalarSequenceClass<double> dataB(m_env,dataSizeB);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeB; ++j) {
	dataB[j] = (*(m_seq[initialPosB+j]))[i];
      }
      std::vector<double> psdData(0,0.);
      dataB.psd(0,
                8,  // numBlocks
                .5, // hopSizeRatio
                psdData);
      psdAtZeroB[i] = psdData[0];
    } // for 'i'

    vectorOfGeweke[initialPosId] = new V(m_vectorExample);

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
uqSequenceOfVectorsClass<V>::minMax(
  unsigned int initialPos,
  V&           mins,
  V&           maxs) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numParams = vectorSize();
  std::vector<double> data(dataSize,0.);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }
    std::vector<double>::iterator pos;
    pos = std::min_element(data.begin(), data.end());
    mins[i] = *pos;
    pos = std::max_element(data.begin(), data.end());
    maxs[i] = *pos;
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::histogram(
  unsigned int     initialPosition,
  unsigned int     spacing,
  const V&         minHorizontalValues,
  const V&         maxHorizontalValues,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& binsForAllParams) const
{
  UQ_FATAL_TEST_MACRO(centersForAllBins.size() != binsForAllParams.size(),
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::histogram()",
                      "vectors 'centers' and 'bins' have different sizes");

  for (unsigned int j = 0; j < binsForAllParams.size(); ++j) {
    centersForAllBins[j] = new V(m_vectorExample);
    binsForAllParams [j] = new V(m_vectorExample);
  }

  unsigned int dataSize = this->sequenceSize() - initialPosition;
  unsigned int numParams = vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPosition+j]))[i];
    }

    std::vector<double      > centers(centersForAllBins.size(),0.);
    std::vector<unsigned int> bins   (binsForAllParams.size(), 0 );
    data.histogram(0,
                   1,
                   minHorizontalValues[i],
                   maxHorizontalValues[i],
                   centers,
                   bins);

    for (unsigned int j = 0; j < bins.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(binsForAllParams [j]))[i] = (double) bins[j];
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::sort(
  unsigned int     initialPosition,
  std::vector<V*>& sortedSequence) const
{
  UQ_FATAL_TEST_MACRO((this->sequenceSize() - initialPosition) != sortedSequence.size(),
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sort()",
                      "incompatible sizes between vectors 'sequence' and 'sortedSequence'");

  for (unsigned int j = 0; j < sortedSequence.size(); ++j) {
    sortedSequence[j] = new V(m_vectorExample);
  }

  unsigned int dataSize = this->sequenceSize() - initialPosition;
  unsigned int numParams = vectorSize();
  std::vector<double> data(dataSize,0.);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPosition+j]))[i];
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
uqSequenceOfVectorsClass<V>::interQuantileRange(
  unsigned int initialPosition,
  unsigned int spacing,
  V&           iqrs) const
{
  unsigned int dataSize = this->sequenceSize() - initialPosition;

  std::vector<V*> sortedSequence(dataSize,NULL);
  this->sort(initialPosition,
             sortedSequence);

  unsigned int pos1 = (unsigned int) ( (((double) dataSize) + 1.)*1./4. - 1. );
  unsigned int pos3 = (unsigned int) ( (((double) dataSize) + 1.)*3./4. - 1. );

  double fraction1 = (((double) dataSize) + 1.)*1./4. - 1. - ((double) pos1);
  double fraction3 = (((double) dataSize) + 1.)*3./4. - 1. - ((double) pos3);

  unsigned int numParams = vectorSize();
  //std::cout << "In uqSeqOfVecs::iqr()"
  //          << ", initialPosition = " << initialPosition
  //          << ", spacing = " << spacing
  //          << ", this->sequenceSize() = " << this->sequenceSize()
  //          << ", dataSize = " << dataSize
  //          << ", sortedSequence.size() = " << sortedSequence.size()
  //          << ", pos1 = " << pos1
  //          << ", pos3 = " << pos3
  //          << ", numParams = " << numParams
  //          << ", iqrs.size() = " << iqrs.size()
  //          << std::endl;
  for (unsigned int i = 0; i < numParams; ++i) {
    double value1 = (1.-fraction1) * (*sortedSequence[pos1])[i] + fraction1 * (*sortedSequence[pos1+1])[i];
    double value3 = (1.-fraction3) * (*sortedSequence[pos3])[i] + fraction3 * (*sortedSequence[pos3+1])[i];
    iqrs[i] = value3 - value1;
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::scalesForKDE(
  unsigned int initialPosition,
  unsigned int spacing,
  const V&     iqrs,
  V&           scales) const
{
  unsigned int dataSize = this->sequenceSize() - initialPosition;

  V mean(m_vectorExample);
  this->mean(initialPosition,
             dataSize,
             mean);

  V sampleVariance(m_vectorExample);
  this->sampleVariance(initialPosition,
                       dataSize,
                       mean,
                       sampleVariance);

  unsigned int numParams = vectorSize();
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
uqSequenceOfVectorsClass<V>::gaussianKDE(
  unsigned int           initialPosition,
  unsigned int           spacing,
  const std::vector<V*>& evaluationPositions,
  const V&               scales,
  std::vector<V*>&       densityValues) const
{
  unsigned int dataSize = this->sequenceSize() - initialPosition;
  unsigned int numEstimationsPerParam = evaluationPositions.size();

  for (unsigned int j = 0; j < numEstimationsPerParam; ++j) {
    densityValues[j] = new V(m_vectorExample);
  }

  unsigned int numParams = vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    double scaleInv = 1./scales[i];
    for (unsigned int j = 0; j < numEstimationsPerParam; ++j) {
      double x = (*(evaluationPositions[j]))[i];
      double value = 0.;
      for (unsigned int k = 0; k < dataSize; ++k) {
        double xk = (*(m_seq[initialPosition+k]))[i];
        value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
      }
      (*(densityValues[j]))[i] = scaleInv * (value/(double) numEstimationsPerParam);
    }
  }

  return;
}
#endif // __UQ_SEQUENCE_OF_VECTORS_H__

