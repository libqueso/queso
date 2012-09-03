//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_ARRAY_OF_SEQUENCES_H__
#define __UQ_ARRAY_OF_SEQUENCES_H__

#include <uqVectorSequence.h>

template <class V, class M>
class uqArrayOfSequencesClass : public uqBaseVectorSequenceClass<V,M>
{
public:


  uqArrayOfSequencesClass(const uqVectorSpaceClass<V,M>& vectorSpace,
                          unsigned int                   subSequenceSize,
                          const std::string&             name);
 ~uqArrayOfSequencesClass();

        unsigned int subSequenceSize      () const;
        void         resizeSequence       (unsigned int newSubSequenceSize);
        void         resetValues          (unsigned int initialPos, unsigned int numPos);
        void         erasePositions       (unsigned int initialPos, unsigned int numPos);
        void         getPositionValues    (unsigned int posId,       V& vec) const;
        void         setPositionValues    (unsigned int posId, const V& vec);
        void         setGaussian          (const gsl_rng* rng, const V& meanVec, const V& stdDevVec);
        void         setUniform           (const gsl_rng* rng, const V& aVec,    const V& bVec     );
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
        void         uniformlySampledMdf  (const V&                       numEvaluationPointsVec,
                                           uqArrayOfOneDGridsClass <V,M>& mdfGrids,
                                           uqArrayOfOneDTablesClass<V,M>& mdfValues) const;
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
        void         uniformlySampledCdf  (const V&                       numEvaluationPointsVec,
                                           uqArrayOfOneDGridsClass <V,M>& cdfGrids,
                                           uqArrayOfOneDTablesClass<V,M>& cdfValues) const;
#endif
        void         mean                 (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           V&                       meanVec) const;
        void         unifiedMean          (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           V&                       unifiedMeanVec) const;
        void         sampleVariance       (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           const V&                 meanVec,
                                           V&                       samVec) const;
        void         unifiedSampleVariance(unsigned int             initialPos,
                                           unsigned int             numPos,
                                           const V&                 meanVec,
                                           V&                       samVec) const;
        void         populationVariance   (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           const V&                 meanVec,
                                           V&                       popVec) const;
        void         autoCovariance       (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           const V&                 meanVec,
                                           unsigned int             lag,
                                           V&                       covVec) const;

        void         autoCorrViaDef       (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           unsigned int             lag,
                                           V&                       corrVec) const;
        void         autoCorrViaFft       (unsigned int                     initialPos,
                                           unsigned int                     numPos,
                                           const std::vector<unsigned int>& lags,
                                           std::vector<V*>&                 corrVecs) const;
        void         autoCorrViaFft       (unsigned int             initialPos,
                                           unsigned int             numPos,
                                           unsigned int             numSum,
                                           V&                       autoCorrsSumVec) const;
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
        void         bmm                  (unsigned int             initialPos,
                                           unsigned int             batchLength,
                                           V&                       bmmVec) const;
        void         fftForward           (unsigned int                        initialPos,
                                           unsigned int                        fftSize,
                                           unsigned int                        paramId,
                                           std::vector<std::complex<double> >& resultData) const;
      //void         fftInverse           (unsigned int fftSize);
        void         psd                  (unsigned int             initialPos,
                                           unsigned int             numBlocks,
                                           double                   hopSizeRatio,
                                           unsigned int             paramId,
                                           std::vector<double>&     psdResult) const;
        void         psdAtZero            (unsigned int             initialPos,
                                           unsigned int             numBlocks,
                                           double                   hopSizeRatio,
                                           V&                       psdVec) const;
        void         geweke               (unsigned int             initialPos,
                                           double                   ratioNa,
                                           double                   ratioNb,
                                           V&                       gewVec) const;
#endif
        void         minMax               (unsigned int             initialPos,
                                           V&                       minVec,
                                           V&                       maxVec) const;
        void         histogram            (unsigned int             initialPos,
                                           const V&                 minVec,
                                           const V&                 maxVec,
                                           std::vector<V*>&         centersForAllBins,
                                           std::vector<V*>&         quanttsForAllBins) const;
        void         interQuantileRange   (unsigned int             initialPos,
                                           V&                       iqrs) const;
        void         scalesForKDE         (unsigned int             initialPos,
                                           const V&                 iqrs,
                                           unsigned int             kdeDimension,
                                           V&                       scales) const;
        void         gaussianKDE          (const V&                 evaluationParamVec,
                                           V&                       densityVec) const;
        void         gaussianKDE          (unsigned int             initialPos,
                                           const V&                 scales,
                                           const std::vector<V*>&   evaluationParamVecs,
                                           std::vector<V*>&         densityVecs) const;
        void         writeContents        (std::ofstream&                   ofsvar) const;
        void         unifiedWriteContents (std::ofstream&                   ofsvar) const;
        void         unifiedWriteContents (const std::string&               fileName,
                                           const std::string&               fileType) const;
        void         unifiedReadContents  (const std::string&               fileName,
                                           const std::string&               fileType,
                                           const unsigned int               subSequenceSize);
        void         select               (const std::vector<unsigned int>& idsOfUniquePositions);
        void         filter               (unsigned int                     initialPos,
                                           unsigned int                     spacing);

        void         extractScalarSeq     (unsigned int                   initialPos,
                                           unsigned int                   spacing,
                                           unsigned int                   numPos,
                                           unsigned int                   paramId,
                                           uqScalarSequenceClass<double>& scalarSeq) const;
private:
        void         extractRawData       (unsigned int                   initialPos,
                                           unsigned int                   spacing,
                                           unsigned int                   numPos,
                                           unsigned int                   paramId,
                                           std::vector<double>&           rawData) const;

	uqDistArrayClass<uqScalarSequenceClass<double>*> m_scalarSequences;

  using uqBaseVectorSequenceClass<V,M>::m_env;
  using uqBaseVectorSequenceClass<V,M>::m_vectorSpace;
  using uqBaseVectorSequenceClass<V,M>::m_name;
  using uqBaseVectorSequenceClass<V,M>::m_fftObj;
};

template <class V, class M>
uqArrayOfSequencesClass<V,M>::uqArrayOfSequencesClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   subSequenceSize,
  const std::string&             name)
  :
  uqBaseVectorSequenceClass<V,M>(vectorSpace,subSequenceSize,name),
  m_scalarSequences             (m_vectorSpace.map(),1)
{

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering uqArrayOfSequencesClass<V,M>::constructor()"
  //                         << std::endl;
  //}

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "In uqArrayOfSequencesClass<V,M>::constructor()"
  //                         << "\n subSequenceSize = "              << subSequenceSize
  //                         << "\n m_scalarSequences.MyLength() = " << m_scalarSequences.MyLength()
  //                         << std::endl;
  //}

  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0) = new uqScalarSequenceClass<double>(m_env,subSequenceSize);
  }

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving uqArrayOfSequencesClass<V,M>::constructor()"
  //                         << std::endl;
  //}
}

template <class V, class M>
uqArrayOfSequencesClass<V,M>::~uqArrayOfSequencesClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    if (m_scalarSequences(i,0)) delete m_scalarSequences(i,0);
  }
}

template <class V, class M>
unsigned int
uqArrayOfSequencesClass<V,M>::subSequenceSize() const
{
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);

  return tmp->m_scalarSequences(0,0)->subSequenceSize();
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::resizeSequence(unsigned int newSubSequenceSize)
{
  if (newSubSequenceSize != this->subSequenceSize()) {
    for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
      m_scalarSequences(i,0)->resizeSequence(newSubSequenceSize);
    }
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::resetValues(
  unsigned int initialPos,
  unsigned int numPos)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0)->resetValues(initialPos,numPos);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::erasePositions(
  unsigned int initialPos,
  unsigned int numPos)
{
  if (initialPos < this->subSequenceSize()) {
    for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
      uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
      seq.erasePositions(initialPos,numPos);
    }
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::getPositionValues(unsigned int posId, V& vec) const
{
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    vec[i] = (*(tmp->m_scalarSequences(i,0)))[posId];
  }

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::setPositionValues(unsigned int posId, const V& vec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq[posId] = vec[i];
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::setGaussian(const gsl_rng* rng, const V& meanVec, const V& stdDevVec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq.setGaussian(rng,meanVec[i],stdDevVec[i]);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::setUniform(const gsl_rng* rng, const V& aVec, const V& bVec)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    uqScalarSequenceClass<double>& seq = *(m_scalarSequences(i,0));
    seq.setUniform(rng,aVec[i],bVec[i]);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::uniformlySampledCdf(
  const V&                       numEvaluationPointsVec,
  uqArrayOfOneDGridsClass <V,M>& cdfGrids,
  uqArrayOfOneDTablesClass<V,M>& cdfValues) const
{
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
  V minCdfValues(m_vectorSpace.zeroVector());
  V maxCdfValues(m_vectorSpace.zeroVector());
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));

    unsigned int numEvaluationPoints = (unsigned int) numEvaluationPointsVec[i];
    std::vector<double> aCdf(0);
    seq.subUniformlySampledCdf(numEvaluationPoints,
                            minCdfValues[i],
                            maxCdfValues[i],
                            aCdf);
    cdfValues.setOneDTable(i,aCdf);
  }
  cdfGrids.setUniformGrids(numEvaluationPointsVec,
                           minCdfValues,
                           maxCdfValues);

  return;
}
#endif
template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           meanVec) const
{
  bool bRC = ((0                     <= initialPos                 ) &&
              (0                     <  numPos                     ) &&
              ((initialPos+numPos-1) <= (this->subSequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::mean()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::mean()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  meanVec.cwSet(0.);
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
  for (unsigned int i = 0; i < meanVec.size(); ++i) {
    const uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    meanVec[i] = seq.subMeanExtra(initialPos, numPos);
  }

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((0                     <= initialPos                 ) &&
              (0                     <  numPos                     ) &&
              ((initialPos+numPos-1) <= (this->subSequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::sampleVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == samVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::sampleVariance()",
                      "incompatible sizes between samVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::sampleVariance()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  samVec.cwSet(0.);

  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((0                     <= initialPos                 ) &&
              (0                     <  numPos                     ) &&
              ((initialPos+numPos-1) <= (this->subSequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::populationVariance()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == popVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::populationVariance()",
                      "incompatible sizes between popVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::populationVariance()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  popVec.cwSet(0.);

  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  unsigned int lag,
  V&           covVec) const
{
  bool bRC = ((0                     <= initialPos                 ) &&
              (0                     <  numPos                     ) &&
              ((initialPos+numPos-1) <= (this->subSequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqVectorSequenceAutoCovariance<V,M>()",
                      "invalid initial position or number of positions");

  bRC = (numPos > lag);
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqVectorSequenceAutoCovariance<V,M>()",
                      "lag is too large");

  bRC = (this->vectorSize() == covVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqVectorSequenceAutoCovariance<V,M>()",
                      "incompatible sizes between covVec vector and vectors in sequence");

  bRC = (this->vectorSize() == meanVec.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqVectorSequenceAutoCovariance<V,M>()",
                      "incompatible sizes between meanVec vector and vectors in sequence");

  unsigned int loopSize      = numPos - lag;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  covVec.cwSet(0.);

  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::autoCorrViaDef(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int lag,
  V&           corrVec) const
{
  V subChainMean              (m_vectorSpace.zeroVector());
  V subChainAutoCovarianceLag0(m_vectorSpace.zeroVector());

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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::autoCorrViaFft(
  unsigned int                     initialPos,
  unsigned int                     numPos,
  const std::vector<unsigned int>& lags,
  std::vector<V*>&                 corrVecs) const
{
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::autoCorrViaFft(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int numSum,
  V&           autoCorrsSumVec) const
{
#if 0
  bool bRC = ((initialPos             <  this->subSequenceSize()) &&
              (0                      <  numPos                 ) &&
              ((initialPos+numPos)    <= this->subSequenceSize()) &&
              (autoCorrsSumVec.size() == this->vectorSize()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::autoCorrViaFft(), for sum",
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

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::bmm(
  unsigned int initialPos,
  unsigned int batchLength,
  V&           bmmVec) const
{
#if 0
  V meanOfBatchMeans   (*(sequence[0]));
  V covLag0OfBatchMeans(*(sequence[0]));
  V covLag1OfBatchMeans(*(sequence[0]));

  V tmpVector(m_vectorSpace.zeroVector()); // In order to contour the fact that 'batchMeans' is a vector of 'const V*', but needs to be set first
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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::fftForward(
  unsigned int                        initialPos,
  unsigned int                        fftSize,
  unsigned int                        paramId,
  std::vector<std::complex<double> >& resultData) const
{
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::psd(
  unsigned int         initialPos,
  unsigned int         numBlocks,
  double               hopSizeRatio,
  unsigned int         paramId,
  std::vector<double>& psdResult) const
{
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::psdAtZero(
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
        //if (m_env.subDisplayFile()) {
	//  *m_env.subDisplayFile() << "psdSequence[0] = " << psdSequence[0] << std::endl;
        //}
      } // for 'numsOfBlocksId'
    } // for 'i'
  }
#endif
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::geweke(
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
      (*(vectorOfGeweke[initialPosId]))[i] = (meanA[i] - meanB[i])/std::sqrt(psdVecA[i]/doubleDataSizeA + psdVecB[i]/doubleDataSizeB);
    }
  }

#endif
  return;
}
#endif // #ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::minMax(
  unsigned int initialPos,
  V&           minVec,
  V&           maxVec) const
{
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double>& seq = *(tmp->m_scalarSequences(i,0));
    seq.subMinMaxExtra(initialPos,minVec[i],maxVec[i]);
  }

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::histogram(
  unsigned int     initialPos,
  const V&         minVec,
  const V&         maxVec,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& quanttsForAllBins) const
{
#if 0
  UQ_FATAL_TEST_MACRO(centersForAllBins.size() != quanttsForAllBins.size(),
                      sequence[0]->env().worldRank(),
                      "uqVectorSequenceHistogram<V,M>()",
                      "vectors 'centers' and 'quantts' have different sizes");

  for (unsigned int j = 0; j < quanttsForAllBins.size(); ++j) {
    centersForAllBins[j] = new V(*(sequence[0]));
    quanttsForAllBins [j] = new V(*(sequence[0]));
  }

  unsigned int dataSize = sequence.size() - initialPos;
  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(dataSize,0.);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(sequence[initialPos+j]))[i];
    }

    std::vector<double> centers(centersForAllBins.size(),0.);
    std::vector<double> quantts(quanttsForAllBins.size(),0.);
    data.histogram(minVec[i],
                   maxVec[i],
                   centers,
                   quantts);

    for (unsigned int j = 0; j < quantts.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(quanttsForAllBins[j]))[i] = quantts[j];
    }
  }

#endif
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::interQuantileRange(
  unsigned int initialPos,
  V&           iqrs) const
{
#if 0
  unsigned int dataSize = sequence.size() - initialPos;

  uqArrayOfSequencesClass sortedSequence(dataSize,m_vectorSpace.zeroVector());
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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::scalesForKDE(
  unsigned int initialPos,
  const V&     iqrs,
  unsigned int kdeDimension,
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
      scales[i] = 1.06*std::sqrt(samVec[i])/std::pow(dataSize,1./5.);
    }
    else {
      scales[i] = 1.06*std::min(std::sqrt(samVec[i]),iqrs[i]/1.34)/std::pow(dataSize,1./5.);
    }
  }

#endif
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::gaussianKDE(
  const V& evaluationParamVec,
  V&       densityVec) const
{
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::gaussianKDE(
  unsigned int           initialPos,
  const V&               scales,
  const std::vector<V*>& evaluationParamVecs,
  std::vector<V*>&       densityVecs) const
{
#if 0
  unsigned int dataSize = sequence.size() - initialPos;
  unsigned int numEstimationsPerParam = evaluationParamVecs.size();

  for (unsigned int j = 0; j < numEstimationsPerParam; ++j) {
    densityVecs[j] = new V(*(sequence[0]));
  }

  unsigned int numParams = sequence[0]->size();
  for (unsigned int i = 0; i < numParams; ++i) {
    double scaleInv = 1./scales[i];
    for (unsigned int j = 0; j < numEstimationsPerParam; ++j) {
      double x = (*(evaluationParamVecs[j]))[i];
      double value = 0.;
      for (unsigned int k = 0; k < dataSize; ++k) {
        double xk = (*(sequence[initialPos+k]))[i];
        value += uqMiscGaussianDensity((x-xk)*scaleInv,0.,1.);
      }
      (*(densityVecs[j]))[i] = scaleInv * (value/(double) numEstimationsPerParam);
    }
  }

#endif
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::select(const std::vector<unsigned int>& idsOfUniquePositions)
{
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::filter(
  unsigned int initialPos,
  unsigned int spacing)
{
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::writeContents(std::ofstream& ofsvar) const
{
  // Write chain
  ofsvar << m_name << "_sub" << m_env.subIdString() << " = zeros(" << this->subSequenceSize()
         << ","                                                    << this->vectorSize()
         << ");"
         << std::endl;
  ofsvar << m_name << "_sub" << m_env.subIdString() << " = [";

  V tmpVec(m_vectorSpace.zeroVector());
  unsigned int chainSize = this->subSequenceSize();
  for (unsigned int j = 0; j < chainSize; ++j) {
    this->getPositionValues(j,tmpVec);
    ofsvar << tmpVec
           << std::endl;
  }
  ofsvar << "];\n";

  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::unifiedWriteContents(std::ofstream& ofsvar) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::unifiedWriteContents(1)",
                      "not implemented yet");
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::unifiedWriteContents(
  const std::string& fileName,
  const std::string& fileType) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::unifiedWriteContents(2)",
                      "not implemented yet");
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::unifiedReadContents(
  const std::string& fileName,
  const std::string& fileType,
  const unsigned int subSequenceSize)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqArrayOfSequencesClass<V,M>::unifiedReadContents()",
                      "not implemented yet");
  return;
}

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::extractScalarSeq(
  unsigned int                   initialPos,
  unsigned int                   spacing,
  unsigned int                   numPos,
  unsigned int                   paramId,
  uqScalarSequenceClass<double>& scalarSeq) const
{
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
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

template <class V, class M>
void
uqArrayOfSequencesClass<V,M>::extractRawData(
  unsigned int         initialPos,
  unsigned int         spacing,
  unsigned int         numPos,
  unsigned int         paramId,
  std::vector<double>& rawData) const
{
  uqArrayOfSequencesClass<V,M>* tmp = const_cast<uqArrayOfSequencesClass<V,M>*>(this);
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

