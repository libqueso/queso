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

#include <uqChain.h>
#include <EpetraExt_DistArray.h>

template <class V>
class uqArrayOfSequencesClass : public uqChainBaseClass<V>
{
public:


  uqArrayOfSequencesClass(unsigned int sequenceSize, const V& vectorExample);
 ~uqArrayOfSequencesClass();

  const unsigned int sequenceSize      () const;
        void         resizeSequence    (unsigned int newSequenceSize);
        void         resetValues       (unsigned int initialPos, unsigned int numPos);
        void         erasePositions    (unsigned int initialPos, unsigned int numPos);
        void         getPositionValues (unsigned int posId,       V& vec) const;
        void         setPositionValues (unsigned int posId, const V& vec);
        void         setGaussian       (const gsl_rng* rng, const V& meanVec, const V& stdDevVec);
        void         setUniform        (const gsl_rng* rng, const V& aVec,    const V& bVec     );

        void         mean              (unsigned int             initialPos,
                                        unsigned int             numPos,
                                        V&                       meanVec) const;
        void         sampleVariance    (unsigned int             initialPos,
                                        unsigned int             numPos,
                                        const V&                 meanVec,
                                        V&                       samVec) const;
        void         populationVariance(unsigned int             initialPos,
                                        unsigned int             numPos,
                                        const V&                 meanVec,
                                        V&                       popVec) const;
        void         autoCovariance    (unsigned int             initialPos,
                                        unsigned int             numPos,
                                        const V&                 meanVec,
                                        unsigned int             lag,
                                        V&                       covVec) const;

        void         autoCorrViaDef    (unsigned int             initialPos,
                                        unsigned int             numPos,
                                        unsigned int             lag,
                                        V&                       corrVec) const;
        void         autoCorrViaFft    (unsigned int                     initialPos,
                                        unsigned int                     numPos,
                                        const std::vector<unsigned int>& lags,
                                        std::vector<V*>&                 corrVecs) const;
        void         autoCorrViaFft    (unsigned int             initialPos,
                                        unsigned int             numPos,
                                        unsigned int             numSum,
                                        V&                       autoCorrsSumVec) const;
        void         bmm               (unsigned int             initialPos,
                                        unsigned int             batchLength,
                                        V&                       bmmVec) const;
        void         fftForward        (unsigned int                        initialPos,
                                        unsigned int                        fftSize,
                                        unsigned int                        paramId,
                                        std::vector<std::complex<double> >& resultData) const;
        //void         fftInverse        (unsigned int fftSize);
        void         psd               (unsigned int             initialPos,
                                        unsigned int             numBlocks,
                                        double                   hopSizeRatio,
                                        unsigned int             paramId,
                                        std::vector<double>&     psdResult) const;
        void         psdAtZero         (unsigned int             initialPos,
                                        unsigned int             numBlocks,
                                        double                   hopSizeRatio,
                                        V&                       psdVec) const;
        void         geweke            (unsigned int             initialPos,
                                        double                   ratioNa,
                                        double                   ratioNb,
                                        V&                       gewVec) const;
        void         minMax            (unsigned int             initialPos,
                                        V&                       minVec,
                                        V&                       maxVec) const;
        void         histogram         (unsigned int             initialPos,
                                        const V&                 minVec,
                                        const V&                 maxVec,
                                        std::vector<V*>&         centersForAllBins,
                                        std::vector<V*>&         binsForAllParams) const;
        void         interQuantileRange(unsigned int             initialPos,
                                        V&                       iqrs) const;
        void         scalesForKDE      (unsigned int             initialPos,
                                        const V&                 iqrs,
                                        V&                       scales) const;
        void         gaussianKDE       (unsigned int             initialPos,
                                        const V&                 scales,
                                        const std::vector<V*>&   evaluationPositions,
                                        std::vector<V*>&         densityValues) const;
        void         write             (const std::string&       name,
                                        std::ofstream&           ofs) const;
        void         filter            (const std::vector<unsigned int>& idsOfUniquePositions);
        void         filter            (unsigned int             initialPos,
                                        unsigned int             spacing);

private:
        void         extractScalarSeq  (unsigned int                   initialPos,
                                        unsigned int                   spacing,
                                        unsigned int                   numPos,
                                        unsigned int                   paramId,
                                        uqScalarSequenceClass<double>& scalarSeq) const;
        void         extractRawData    (unsigned int                   initialPos,
                                        unsigned int                   spacing,
                                        unsigned int                   numPos,
                                        unsigned int                   paramId,
                                        std::vector<double>&           rawData) const;

  EpetraExt::DistArray<uqScalarSequenceClass<double>*> m_scalarSequences;

  using uqChainBaseClass<V>::m_env;
  using uqChainBaseClass<V>::m_vectorExample;
  using uqChainBaseClass<V>::m_fftObj;
};

template <class V>
uqArrayOfSequencesClass<V>::uqArrayOfSequencesClass(
  unsigned int sequenceSize,
  const V&     vectorExample)
  :
  uqChainBaseClass<V>(sequenceSize,vectorExample),
  m_scalarSequences  (vectorExample.map(),1)
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
void
uqArrayOfSequencesClass<V>::resizeSequence(unsigned int newSequenceSize)
{
  if (newSequenceSize != this->sequenceSize()) {
    for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
      m_scalarSequences(i,0)->resizeSequence(newSequenceSize);
    }
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::resetValues(
  unsigned int initialPos,
  unsigned int numPos)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0)->resetValues(initialPos,numPos);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::erasePositions(
  unsigned int initialPos,
  unsigned int numPos)
{
  if (initialPos < this->sequenceSize()) {
    for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
      uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
      seq.erasePositions(initialPos,numPos);
    }
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::getPositionValues(unsigned int posId, V& vec) const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    vec[i] = (*(tmp->m_scalarSequences(i,0)))[posId];
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::setPositionValues(unsigned int posId, const V& vec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq[posId] = vec[i];
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::setGaussian(const gsl_rng* rng, const V& meanVec, const V& stdDevVec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq.setGaussian(rng,meanVec[i],stdDevVec[i]);
  }
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::setUniform(const gsl_rng* rng, const V& aVec, const V& bVec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq.setUniform(rng,aVec[i],bVec[i]);
  }
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           meanVec) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::mean()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::mean()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  meanVec.cwSet(0.);
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < meanVec.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    meanVec[i] = seq.mean(initialPos, numPos);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::sampleVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == samVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::sampleVariance()",
                      "incompatible sizes between samVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::sampleVariance()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  samVec.cwSet(0.);

  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < samVec.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    double                               tmpMean = meanVec[i];
    double                               result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff = seq[j] - tmpMean;
      result += diff*diff;
    }
    samVec[i] = result/(doubleLoopSize - 1.);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::populationVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == popVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::populationVariance()",
                      "incompatible sizes between popVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::populationVariance()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  popVec.cwSet(0.);

  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < popVec.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    double                               tmpMean = meanVec[i];
    double                               result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff = seq[j] - tmpMean;
      result += diff*diff;
    }
    popVec[i] = result/doubleLoopSize;
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  unsigned int lag,
  V&           covVec) const
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

  bRC = (this->vectorSize() == covVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "incompatible sizes between covVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceAutoCovariance<V>()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  covVec.cwSet(0.);

  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  for (unsigned int i = 0; i < covVec.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    double meanValue = meanVec[i];
    double result = 0;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff1 = seq[j]     - meanValue;
      double diff2 = seq[j+lag] - meanValue;
      result += diff1*diff2;
    }
    covVec[i] = result/doubleLoopSize;
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::autoCorrViaDef(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int lag,
  V&           corrVec) const
{
  V subChainMean              (m_vectorExample);
  V subChainAutoCovarianceLag0(m_vectorExample);

  this->mean(initialPos,
             numPos,
             subChainMean);
  this->autoCovariance(initialPos,
                       numPos,
                       subChainMean,
                       0, // lag
                       subChainAutoCovarianceLag0);

  this->autoCovariance(initialPos,
                       numPos,
                       subChainMean,
                       lag,
                       corrVec);
  corrVec /= subChainAutoCovarianceLag0; 

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::autoCorrViaFft(
  unsigned int                     initialPos,
  unsigned int                     numPos,
  const std::vector<unsigned int>& lags,
  std::vector<V*>&                 corrVecs) const
{
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::bmm(
  unsigned int initialPos,
  unsigned int batchLength,
  V&           bmmVec) const
{
#if 0
  V meanOfBatchMeans   (*(sequence[0]));
  V covLag0OfBatchMeans(*(sequence[0]));
  V covLag1OfBatchMeans(*(sequence[0]));

  V tmpVector(m_vectorExample); // In order to contour the fact that 'batchMeans' is a vector of 'const V*', but needs to be set first
  for (unsigned int initialPosId = 0; initialPosId < initialPositions.size(); initialPosId++) {
    for (unsigned int batchLengthId = 0; batchLengthId < batchLengths.size(); batchLengthId++) {
      unsigned int batchLength = batchLengths[batchLengthId];
      unsigned int numberOfBatches = (sequence.size() - initialPositions[initialPosId])/batchLength;

      std::vector<const V* > batchMeans(numberOfBatches,NULL);
      for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
        uqVectorSequenceMean(sequence,
                             initialPositions[initialPosId] + batchId*batchLength,
                             batchLength,
                             tmpVector);
        batchMeans[batchId] = new V(tmpVector);
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
#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::autoCorrViaFft(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int numSum,
  V&           autoCorrsSumVec) const
{
#if 0
  bool bRC = ((initialPos             <  this->sequenceSize()) &&
              (0                      <  numPos              ) &&
              ((initialPos+numPos)    <= this->sequenceSize()) &&
              (autoCorrsSumVec.size() == this->vectorSize()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqArrayOfSequencesClass<V>::autoCorrViaFft(), for sum",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    data.autoCorrViaFft(0,
                        numPos,
                        autoCorrsSumVec[i]);
  }
#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::fftForward(
  unsigned int                        initialPos,
  unsigned int                        fftSize,
  unsigned int                        paramId,
  std::vector<std::complex<double> >& resultData) const
{
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::psd(
  unsigned int         initialPos,
  unsigned int         numBlocks,
  double               hopSizeRatio,
  unsigned int         paramId,
  std::vector<double>& psdResult) const
{
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::psdAtZero(
  unsigned int initialPos,
  unsigned int numBlocks,
  double       hopSizeRatio,
  V&           psdVec) const
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
        std::vector<double> psdSequence(0,0.); // size will be determined by 'uqScalarSequencePSD()'
        data.psd(numBlocks,
                 hopSizeRatio,
                 psdSequence);
        _2dArrayOfPSDAtZero(initialPosId,numsOfBlocksId)[i] = psdSequence[0];
	std::cout << "psdSequence[0] = " << psdSequence[0] << std::endl;
      } // for 'numsOfBlocksId'
    } // for 'i'
  }
#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::geweke(
  unsigned int initialPos,
  double       ratioNa,
  double       ratioNb,
  V&           gewVec) const
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

    V psdVecA(*(sequence[0]));
    uqScalarSequenceClass<double> dataA(dataSizeA,0.);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeA; ++j) {
	dataA[j] = (*(sequence[initialPosA+j]))[i];
      }
      std::vector<double> psdSequence(0,0.);
      dataA.psd(8,  // numBlocks
                .5, // hopSizeRatio
                psdSequence);
      psdVecA[i] = psdSequence[0];
    } // for 'i'

    V psdVecB(*(sequence[0]));
    uqScalarSequenceClass<double> dataB(dataSizeB,0.);
    for (unsigned int i = 0; i < numParams; ++i) {
      for (unsigned int j = 0; j < dataSizeB; ++j) {
	dataB[j] = (*(sequence[initialPosB+j]))[i];
      }
      std::vector<double> psdSequence(0,0.);
      dataB.psd(8,  // numBlocks
                .5, // hopSizeRatio
                psdSequence);
      psdVecB[i] = psdSequence[0];
    } // for 'i'

    vectorOfGeweke[initialPosId] = new V(*(sequence[0]));

    double doubleDataSizeA = (double) dataSizeA;
    double doubleDataSizeB = (double) dataSizeB;
    for (unsigned int i = 0; i < numParams; ++i) {
      (*(vectorOfGeweke[initialPosId]))[i] = (meanA[i] - meanB[i])/sqrt(psdVecA[i]/doubleDataSizeA + psdVecB[i]/doubleDataSizeB);
    }
  }

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::minMax(
  unsigned int initialPos,
  V&           minVec,
  V&           maxVec) const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    seq.minMax(initialPos,minVec[i],maxVec[i]);
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::histogram(
  unsigned int     initialPos,
  const V&         minVec,
  const V&         maxVec,
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

  unsigned int dataSize = sequence.size() - initialPos;
  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(dataSize,0.);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPos+j]))[i];
    }

    std::vector<double> centers(centersForAllBins.size(),0.);
    std::vector<double> bins   (binsForAllParams.size(), 0.);
    data.histogram(minVec[i],
                   maxVec[i],
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
uqArrayOfSequencesClass<V>::interQuantileRange(
  unsigned int initialPos,
  V&           iqrs) const
{
#if 0
  unsigned int dataSize = sequence.size() - initialPos;

  uqArrayOfSequencesClass sortedSequence(dataSize,m_vectorExample);
  this->sort(initialPos,
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
  unsigned int initialPos,
  const V&     iqrs,
  V&           scales) const
{
#if 0
  unsigned int dataSize = sequence.size() - initialPos;

  V mean(*(sequence[0]));
  uqVectorSequenceMean(sequence,
                       initialPos,
                       dataSize,
                       mean);

  V samVec(*(sequence[0]));
  uqVectorSequenceSampleVariance(sequence,
                                 initialPos,
                                 dataSize,
                                 mean,
                                 samVec);

  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    if (iqrs[i] <= 0.) {
      scales[i] = 1.06*sqrt(samVec[i])/pow(dataSize,1./5.);
    }
    else {
      scales[i] = 1.06*std::min(sqrt(samVec[i]),iqrs[i]/1.34)/pow(dataSize,1./5.);
    }
  }

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::gaussianKDE(
  unsigned int           initialPos,
  const V&               scales,
  const std::vector<V*>& evaluationPositions,
  std::vector<V*>&       densityValues) const
{
#if 0
  unsigned int dataSize = sequence.size() - initialPos;
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
        double xk = (*(sequence[initialPos+k]))[i];
        value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
      }
      (*(densityValues[j]))[i] = scaleInv * (value/(double) numEstimationsPerParam);
    }
  }

#endif
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::filter(const std::vector<unsigned int>& idsOfUniquePositions)
{
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::filter(
  unsigned int initialPos,
  unsigned int spacing)
{
  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::write(
  const std::string& chainName,
  std::ofstream&     ofs) const
{
  // Write chain
  ofs << chainName << " = zeros(" << this->sequenceSize()
      << ","                      << this->vectorSize()
      << ");"
      << std::endl;
  ofs << chainName << " = [";

  V tmpVec(m_vectorExample);
  unsigned int chainSize = this->sequenceSize();
  for (unsigned int j = 0; j < chainSize; ++j) {
    this->getPositionValues(j,tmpVec);
    ofs << tmpVec
        << std::endl;
  }
  ofs << "];\n";

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::extractScalarSeq(
  unsigned int                   initialPos,
  unsigned int                   spacing,
  unsigned int                   numPos,
  unsigned int                   paramId,
  uqScalarSequenceClass<double>& scalarSeq) const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(paramId,0));

  scalarSeq.resizeSequence(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = seq[paramId];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = seq[paramId];
    }
  }

  return;
}

template <class V>
void
uqArrayOfSequencesClass<V>::extractRawData(
  unsigned int         initialPos,
  unsigned int         spacing,
  unsigned int         numPos,
  unsigned int         paramId,
  std::vector<double>& rawData) const
{
  uqArrayOfSequencesClass<V>* tmp = const_cast<uqArrayOfSequencesClass<V>*>(this);
  uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(paramId,0));

  rawData.resize(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = seq[paramId];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = seq[paramId];
    }
  }

  return;
}
#endif // __UQ_ARRAY_OF_SEQUENCES_H__

