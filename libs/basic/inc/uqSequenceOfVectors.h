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

        void         mean              (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        V&                        meanVec) const;
        void         sampleVariance    (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        const V&                  meanVec,
                                        V&                        samVec) const;
        void         populationVariance(unsigned int              initialPos,
                                        unsigned int              numPos,
                                        const V&                  meanVec,
                                        V&                        popVec) const;

        void         autoCovariance    (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        const V&                  meanVec,
                                        unsigned int              lag,
                                        V&                        covVec) const;
        void         autoCorrelations  (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        unsigned int              lag,
                                        V&                        corrVec) const;
        void         bmm               (unsigned int              initialPos,
                                        unsigned int              batchLength,
                                        V&                        bmmVec) const;
        //void         fftAlloc          ();
        //void         fftForward        (unsigned int fftSize);
        //void         fftInverse        (unsigned int fftSize);
        //void         fftFree           ();
        void         psdAtZero         (unsigned int              initialPos,
                                        unsigned int              numBlocks,
                                        double                    hopSizeRatio,
                                        V&                        psdVec) const;
        void         geweke            (unsigned int              initialPos,
                                        double                    ratioNa,
                                        double                    ratioNb,
                                        V&                        gewVec) const;
        void         minMax            (unsigned int              initialPos,
                                        V&                        minVec,
                                        V&                        maxVec) const;
        void         histogram         (unsigned int              initialPos,
                                        unsigned int              spacing,
                                        const V&                  minVec,
                                        const V&                  maxVec,
                                        std::vector<V*>&          centersForAllBins,
                                        std::vector<V*>&          binsForAllParams) const;
        void         sort              (unsigned int              initialPos,
                                        std::vector<V*>&          sortedSequence) const; // Instead of uqSequenceOfVectorsClass&
        void         interQuantileRange(unsigned int              initialPos,
                                        unsigned int              spacing,
                                        V&                        iqrs) const;
        void         scalesForKDE      (unsigned int              initialPos,
                                        unsigned int              spacing,
                                        const V&                  iqrs,
                                        V&                        scales) const;
        void         gaussianKDE       (unsigned int              initialPos,
                                        unsigned int              spacing,
                                        const std::vector<V*>&    evaluationPositions,
                                        const V&                  scales,
                                        std::vector<V*>&          densityValues) const;
private:
        void         extractScalarSeq  (unsigned int                   initialPos,
                                        unsigned int                   spacing,
                                        unsigned int                   numPos,
                                        unsigned int                   paramId,
                                        uqScalarSequenceClass<double>& scalarSeq) const;
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
  V&           meanVec) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::mean()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::mean()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  meanVec.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < meanVec.size(); ++j) {
      meanVec[j] += (*m_seq[i])[j]/doubleLoopSize;
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == samVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "incompatible sizes between samVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  samVec.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < samVec.size(); ++j) {
      double diff = (*m_seq[i])[j] - meanVec[j];
      samVec[j] += diff*diff;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < samVec.size(); ++j) {
    samVec[j] /= (doubleLoopSize - 1.);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == popVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "incompatible sizes between popVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  popVec.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < popVec.size(); ++j) {
      double diff = (*m_seq[i])[j] - meanVec[j];
      popVec[j] += diff*diff;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < popVec.size(); ++j) {
    popVec[j] /= doubleLoopSize;
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  unsigned int lag,
  V&           covVec) const
{
  bool bRC = ((0                     <  numPos                 ) &&
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

  bRC = (this->vectorSize() == covVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectors<V>::autoCovariance<V>()",
                      "incompatible sizes between covVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectors<V>::autoCovariance<V>()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  covVec.cwSet(0.);
  for (unsigned int i = initialPos; i < finalPosPlus1; ++i) {
    for (unsigned int j = 0; j < covVec.size(); ++j) {
      double diff1 = (*m_seq[i    ])[j] - meanVec[j];
      double diff2 = (*m_seq[i+lag])[j] - meanVec[j];
      covVec[j] += diff1*diff2;
    }
  }

  double doubleLoopSize = (double) loopSize;
  for (unsigned int j = 0; j < covVec.size(); ++j) {
    covVec[j] /= doubleLoopSize;
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::autoCorrelations(
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
uqSequenceOfVectorsClass<V>::bmm(
  unsigned int initialPos,
  unsigned int batchLength,
  V&           bmmVec) const
{
  V meanOfBatchMeans   (m_vectorExample);
  V covLag0OfBatchMeans(m_vectorExample);
  V covLag1OfBatchMeans(m_vectorExample);

  unsigned int numberOfBatches = (this->sequenceSize() - initialPos)/batchLength;
  uqSequenceOfVectorsClass batchMeans(numberOfBatches,m_vectorExample);

  V tmpVector(m_vectorExample); // In order to contour the fact that 'batchMeans' is a vector of 'const V*', but needs to be set first
  for (unsigned int batchId = 0; batchId < numberOfBatches; batchId++) {
    this->mean(initialPos + batchId*batchLength,
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
                            bmmVec);

  bmmVec /= (double) batchMeans.sequenceSize(); // CHECK
//bmmVec *= (double) (this->sequenceSize() - initialPos); // CHECK

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::psdAtZero(
  unsigned int initialPos,
  unsigned int numBlocks,
  double       hopSizeRatio,
  V&           psdVec) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,dataSize);
  std::vector<double> psdSequence(0,0.); // size will be determined by 'uqScalarSequencePSD()'

  unsigned int numParams = vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }
    data.psd(0,
             numBlocks,
             hopSizeRatio,
             psdSequence);
    psdVec[i] = psdSequence[0];

    std::cout << "psdSequence[0] = " << psdSequence[0] << std::endl;

    //std::cout << "psdSequence = zeros(" << psdSequence.size() << ",1);" << std::endl;
    //for (unsigned j = 0; j < psdSequence.size(); ++j) {
    //  std::cout << "psdSequence(" << j+1 << ") = " << psdSequence[j] << ";" << std::endl;
    //}
  } // for 'i'

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::geweke(
  unsigned int initialPos,
  double       ratioNa,
  double       ratioNb,
  V&           gewVec) const
{
#if 1
  unsigned int numPos = this->sequenceSize() - initialPos;
  unsigned int numParams = vectorSize();
  uqScalarSequenceClass<double> data(m_env,numPos);

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    gewVec[i] = data.geweke(0,
                            ratioNa,
                            ratioNb);
  }
#else
  double       doubleFullDataSize = (double) (this->sequenceSize() - initialPos);
  unsigned int dataSizeA          = (unsigned int) (doubleFullDataSize * ratioNa);
  unsigned int dataSizeB          = (unsigned int) (doubleFullDataSize * ratioNb);
  unsigned int initialPosA        = initialPos;
  unsigned int initialPosB        = this->sequenceSize() - dataSizeB;

  V meanA(m_vectorExample);
  this->mean(initialPosA,
             dataSizeA,
             meanA);

  V meanB(m_vectorExample);
  this->mean(initialPosB,
             dataSizeB,
             meanB);

  unsigned int numParams = vectorSize();

  std::vector<double> psdSequence(0,0.);
  V psdVecA(m_vectorExample);
  uqScalarSequenceClass<double> dataA(m_env,dataSizeA);
  for (unsigned int i = 0; i < numParams; ++i) {
   for (unsigned int j = 0; j < dataSizeA; ++j) {
      dataA[j] = (*(m_seq[initialPosA+j]))[i];
    }
    dataA.psd(0,
              8,  // numBlocks
              .5, // hopSizeRatio
              psdSequence);
    psdVecA[i] = psdSequence[0];
  } // for 'i'

  V psdVecB(m_vectorExample);
  uqScalarSequenceClass<double> dataB(m_env,dataSizeB);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSizeB; ++j) {
      dataB[j] = (*(m_seq[initialPosB+j]))[i];
    }
    dataB.psd(0,
              8,  // numBlocks
              .5, // hopSizeRatio
              psdSequence);
    psdVecB[i] = psdSequence[0];
  } // for 'i'

  double doubleDataSizeA = (double) dataSizeA;
  double doubleDataSizeB = (double) dataSizeB;
#if 0
  if (m_env.rank() == 0) {
    std::cout << "In uqSequenceOfVectorsClass<T>::geweke()"
              << ", before computation of gewCoef"
              << ": meanA "             << meanA[0]
              << ", psdA = "            << psdVecA[0]
              << ", doubleDataSizeA = " << doubleDataSizeA
              << ", meanB = "           << meanB[0]
              << ", psdB = "            << psdVecB[0]
              << ", doubleDataSizeB = " << doubleDataSizeB
              << std::endl;
  }
#endif
  for (unsigned int i = 0; i < numParams; ++i) {
    gewVec[i] = (meanA[i] - meanB[i])/sqrt(psdVecA[i]/doubleDataSizeA + psdVecB[i]/doubleDataSizeB);
  }
#endif
  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::minMax(
  unsigned int initialPos,
  V&           minVec,
  V&           maxVec) const
{
  unsigned int numPos = this->sequenceSize() - initialPos;
  unsigned int numParams = vectorSize();
  uqScalarSequenceClass<double> data(m_env,numPos);

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.minMax(0,minVec[i],maxVec[i]);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::histogram(
  unsigned int     initialPos,
  unsigned int     spacing,
  const V&         minVec,
  const V&         maxVec,
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

  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numParams = vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > centers(centersForAllBins.size(),0.);
    std::vector<unsigned int> bins   (binsForAllParams.size(), 0 );
    data.histogram(0,
                   1,
                   minVec[i],
                   maxVec[i],
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
  unsigned int     initialPos,
  std::vector<V*>& sortedSequence) const
{
  UQ_FATAL_TEST_MACRO((this->sequenceSize() - initialPos) != sortedSequence.size(),
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sort()",
                      "incompatible sizes between vectors 'sequence' and 'sortedSequence'");

  for (unsigned int j = 0; j < sortedSequence.size(); ++j) {
    sortedSequence[j] = new V(m_vectorExample);
  }

  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numParams = vectorSize();
  std::vector<double> data(dataSize,0.);
  for (unsigned int i = 0; i < numParams; ++i) {
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
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
  unsigned int initialPos,
  unsigned int spacing,
  V&           iqrs) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;

  std::vector<V*> sortedSequence(dataSize,NULL);
  this->sort(initialPos,
             sortedSequence);

  unsigned int pos1 = (unsigned int) ( (((double) dataSize) + 1.)*1./4. - 1. );
  unsigned int pos3 = (unsigned int) ( (((double) dataSize) + 1.)*3./4. - 1. );

  double fraction1 = (((double) dataSize) + 1.)*1./4. - 1. - ((double) pos1);
  double fraction3 = (((double) dataSize) + 1.)*3./4. - 1. - ((double) pos3);

  unsigned int numParams = vectorSize();
  //std::cout << "In uqSeqOfVecs::iqr()"
  //          << ", initialPos = " << initialPos
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
  unsigned int initialPos,
  unsigned int spacing,
  const V&     iqrs,
  V&           scales) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;

  V meanVec(m_vectorExample);
  this->mean(initialPos,
             dataSize,
             meanVec);

  V samVec(m_vectorExample);
  this->sampleVariance(initialPos,
                       dataSize,
                       meanVec,
                       samVec);

  unsigned int numParams = vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    if (iqrs[i] <= 0.) {
      scales[i] = 1.06*sqrt(samVec[i])/pow(dataSize,1./5.);
    }
    else {
      scales[i] = 1.06*std::min(sqrt(samVec[i]),iqrs[i]/1.34)/pow(dataSize,1./5.);
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::gaussianKDE(
  unsigned int           initialPos,
  unsigned int           spacing,
  const std::vector<V*>& evaluationPositions,
  const V&               scales,
  std::vector<V*>&       densityValues) const
{
  unsigned int dataSize = this->sequenceSize() - initialPos;
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
        double xk = (*(m_seq[initialPos+k]))[i];
        value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
      }
      (*(densityValues[j]))[i] = scaleInv * (value/(double) numEstimationsPerParam);
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::extractScalarSeq(
  unsigned int                   initialPos,
  unsigned int                   spacing,
  unsigned int                   numPos,
  unsigned int                   paramId,
  uqScalarSequenceClass<double>& scalarSeq) const
{
  scalarSeq.resizeSequence(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = (*(m_seq[initialPos+j        ]))[paramId];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = (*(m_seq[initialPos+j*spacing]))[paramId];
    }
  }

  return;
}
#endif // __UQ_SEQUENCE_OF_VECTORS_H__

