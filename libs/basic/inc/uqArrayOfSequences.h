/* uq/libs/basic/inc/uqArrayOfSequences.h
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

#ifndef __UQ_ARRAY_OF_SEQUENCES_H__
#define __UQ_ARRAY_OF_SEQUENCES_H__

#include <EpetraExt_DistArray.h>
#include <uq2dArrayOfStuff.h>
#include <uqScalarSequence.h>

template <class V>
class uqArrayOfSequencesClass
{
public:
  uqArrayOfSequencesClass(unsigned int sequenceSize, const V& vectorExample);
 ~uqArrayOfSequencesClass();

  const unsigned int sequenceSize      () const;
  const unsigned int vectorSize        () const;
        void         resizeSequence    (unsigned int newSequenceSize);
        void         resetValues       ();
        void         erasePositions    (unsigned int posBegin, unsigned int posEnd);
        void         getPositionValues (unsigned int positionId,       V& vector) const;
        void         setPositionValues (unsigned int positionId, const V& vector);
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
                                        uqArrayOfSequencesClass&         sortedSequence) const;
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
  const uqEnvironmentClass&                            m_env;
  V                                                    m_vectorExample;
  EpetraExt::DistArray<uqScalarSequenceClass<double>*> m_scalarSequences;
};

template <class V>
uqArrayOfSequencesClass<V>::uqArrayOfSequencesClass(
  unsigned int sequenceSize,
  const V&     vectorExample)
  :
  m_env            (vectorExample.env()),
  m_vectorExample  (vectorExample),
  m_scalarSequences(vectorExample.map(),1)
{

  //if (m_env.rank() == 0) std::cout << "Entering uqArrayOfSequencesClass<V>::constructor()"
  //                                 << std::endl;

  //if (m_env.rank() == 0) std::cout << "In uqArrayOfSequencesClass<V>::constructor()"
  //                                 << "\n sequenceSize = "                 << sequenceSize
  //                                 << "\n m_scalarSequences.MyLength() = " << m_scalarSequences.MyLength()
  //                                 << std::endl;

  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0) = new uqScalarSequenceClass<double>(m_env,sequenceSize);
  }

  //if (m_env.rank() == 0) std::cout << "Leaving uqArrayOfSequencesClass<V>::constructor()"
  //                                 << std::endl;
}

template <class V>
uqArrayOfSequencesClass<V>::~uqArrayOfSequencesClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    if (m_scalarSequences(i,0)) delete m_scalarSequences(i,0);
  }
}

template <class V>
const unsigned int
uqArrayOfSequencesClass<V>::sequenceSize() const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);

  return tmp->m_scalarSequences(0,0)->sequenceSize();
}

template <class V>
const unsigned int
uqArrayOfSequencesClass<V>::vectorSize() const
{
  return m_vectorExample.size();
}

template <class V>
void
uqArrayOfSequencesClass<V>::resizeSequence(unsigned int newSequenceSize)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0)->resizeSequence(newSequenceSize);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::resetValues()
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0)->resetValues();
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::erasePositions(unsigned int posBegin, unsigned int posEnd)
{
  if (posBegin < this->sequenceSize()) {
    for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
      uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
      seq.erasePositions(posBegin,posEnd);
    }
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::getPositionValues(unsigned int positionId, V& vector) const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    vector[i] = (*(tmp->m_scalarSequences(i,0)))[positionId];
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::setPositionValues(unsigned int positionId, const V& vector)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq[positionId] = vector[i];
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::setGaussian(gsl_rng* rng, const V& meanVec, const V& stdDevVec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq.setGaussian(rng,meanVec[i],stdDevVec[i]);
  }
  return;
}


template <class V>
void
uqArrayOfSequencesClass<V>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           mean) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::mean()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::mean()",
                      "incompatible sizes between mean vector and vectors in sequence");

  mean.cwSet(0.);
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < mean.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    mean[i] = seq.mean(initialPos, numPos);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::sampleVariance(
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
                      "uqArrayOfSequencesClass<V>::sampleVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == sampleVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::sampleVariance()",
                      "incompatible sizes between sampleVariance vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::sampleVariance()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  sampleVariance.cwSet(0.);

  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < sampleVariance.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    double                               tmpMean = mean[i];
    double                               result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff = seq[j] - tmpMean;
      result += diff*diff;
    }
    sampleVariance[i] = result/(doubleLoopSize - 1.);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::populationVariance(
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
                      "uqArrayOfSequencesClass<V>::populationVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == populVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::populationVariance()",
                      "incompatible sizes between populVariance vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::populationVariance()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  populVariance.cwSet(0.);

  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < populVariance.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    double                               tmpMean = mean[i];
    double                               result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff = seq[j] - tmpMean;
      result += diff*diff;
    }
    populVariance[i] = result/doubleLoopSize;
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::autoCovariance(
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
                      "uqVectorSequenceAutoCovariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (numPos > lag);
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "lag is too large");

  bRC = (this->vectorSize() == autoCov.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "incompatible sizes between autoCov vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  autoCov.cwSet(0.);

  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < autoCov.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    double meanValue = mean[i];
    double result = 0;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff1 = seq[j]     - meanValue;
      double diff2 = seq[j+lag] - meanValue;
      result += diff1*diff2;
    }
    autoCov[i] = result/doubleLoopSize;
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::autoCorrelations(
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
uqArrayOfSequencesClass<V>::bmm(
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& batchLengths,
  uq2dArrayOfStuff<V>&             _2dArrayOfBMM) const // [numOfPos x numOfLengths] matrix
{
#if 0
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
      //_2dArrayOfBMM(initialPosId,batchLengthId) /= (double) batchMeans.size(); // CHECK
      _2dArrayOfBMM(initialPosId,batchLengthId) *= (double) (sequence.size() - initialPositions[initialPosId]); // CHECK

      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        if (batchMeans[batchId] != NULL) delete batchMeans[batchId];
      }
    }
  }
  delete tmpVector;

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::psdAtZero(
  const std::vector<unsigned int>& initialPositions,
  const std::vector<unsigned int>& numsOfBlocks,
  double                           hopSizeRatio,
  uq2dArrayOfStuff<V>&             _2dArrayOfPSDAtZero) const // [numOfPos x numOfBlocks] matrix
{
#if 0
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    unsigned int dataSize = sequence.size() - initialPositions[initialPosId];
    uqScalarSequenceClass<double> data(dataSize,0.);

    unsigned int numParams = sequence[0]->size();
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSize; ++j) {
	data[j] = (*(sequence[initialPositions[initialPosId]+j]))[i];
      }
      for (unsigned int numsOfBlocksId = 0; numsOfBlocksId < numsOfBlocks.size(); numsOfBlocksId++) {
        unsigned int numBlocks = numsOfBlocks[numsOfBlocksId];
        std::vector<double> psdData(0,0.); // size will be determined by 'uqScalarSequencePSD()'
        data.psd(numBlocks,
                 hopSizeRatio,
                 psdData);
        _2dArrayOfPSDAtZero(initialPosId,numsOfBlocksId)[i] = psdData[0];
	std::cout << "psdData[0] = " << psdData[0] << std::endl;
      } // for 'numsOfBlocksId'
    } // for 'i'
  }
#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::geweke(
  const std::vector<unsigned int>& initialPositions,
  double                           ratioNa,
  double                           ratioNb,
  std::vector<V*>&                 vectorOfGeweke) const
{
#if 0
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
    uqScalarSequenceClass<double> dataA(dataSizeA,0.);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeA; ++j) {
	dataA[j] = (*(sequence[initialPosA+j]))[i];
      }
      std::vector<double> psdData(0,0.);
      dataA.psd(8,  // numBlocks
                .5, // hopSizeRatio
                psdData);
      psdAtZeroA[i] = psdData[0];
    } // for 'i'

    V psdAtZeroB(*(sequence[0]));
    uqScalarSequenceClass<double> dataB(dataSizeB,0.);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeB; ++j) {
	dataB[j] = (*(sequence[initialPosB+j]))[i];
      }
      std::vector<double> psdData(0,0.);
      dataB.psd(8,  // numBlocks
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

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::minMax(
  unsigned int initialPos,
  V&           mins,
  V&           maxs) const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);

  unsigned int numParams = vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    seq.minMax(initialPos,mins[i],maxs[i]);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::histogram(
  unsigned int     initialPosition,
  unsigned int     spacing,
  const V&         minHorizontalValues,
  const V&         maxHorizontalValues,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& binsForAllParams) const
{
#if 0
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
    uqScalarSequenceClass<double> data(dataSize,0.);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPosition+j]))[i];
    }

    std::vector<double> centers(centersForAllBins.size(),0.);
    std::vector<double> bins   (binsForAllParams.size(), 0.);
    data.histogram(minHorizontalValues[i],
                   maxHorizontalValues[i],
                   centers,
                   bins);

    for (unsigned int j = 0; j < bins.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(binsForAllParams [j]))[i] = bins[j];
    }
  }

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::sort(
  unsigned int             initialPosition,
  uqArrayOfSequencesClass& sortedSequence) const
{
#if 0
  UQ_FATAL_TEST_MACRO((sequence.size() - initialPosition) != sortedSequence.size(),
                      sequence[0]->env().rank(),
                      "uqVectorSequenceSort<V>()",
                      "incompatible sizes between vectors 'sequence' and 'sortedSequence'");

  for (unsigned int j = 0; j < sortedSequence.size(); ++j) {
    sortedSequence[j] = new V(*(sequence[0]));
  }

  unsigned int dataSize = sequence.size() - initialPosition;
  unsigned int numParams = sequence[0]->size();
  uqScalarSequenceClass<double> data(dataSize,0.);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPosition+j]))[i];
    }

    std::sort(data.begin(), data.end());

    for (unsigned int j = 0; j < dataSize; ++j) {
      (*(sortedSequence[j]))[i] = data[j];
    }
  }

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::interQuantileRange(
  unsigned int initialPosition,
  unsigned int spacing,
  V&           iqrs) const
{
#if 0
  unsigned int dataSize = sequence.size() - initialPosition;

  uqArrayOfSequencesClass sortedSequence(dataSize,m_vectorExample);
  this->sort(initialPosition,
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

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::scalesForKDE(
  unsigned int initialPosition,
  unsigned int spacing,
  const V&     iqrs,
  V&           scales) const
{
#if 0
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

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::gaussianKDE(
  unsigned int           initialPosition,
  unsigned int           spacing,
  const std::vector<V*>& evaluationPositions,
  const V&               scales,
  std::vector<V*>&       densityValues) const
{
#if 0
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

#endif
  return;
}
#endif // __UQ_ARRAY_OF_SEQUENCES_H__

