/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_SEQUENCE_OF_VECTORS_H__
#define __UQ_SEQUENCE_OF_VECTORS_H__

#include <uqVectorSequence.h>
#define UQ_SEQ_VEC_USES_SCALAR_SEQ_CODE
#undef UQ_SEQ_VEC_USES_OPERATOR

template <class V, class M>
class uqSequenceOfVectorsClass : public uqBaseVectorSequenceClass<V,M>
{
public:
  typedef typename std::vector<const V*>::const_iterator seqVectorPositionConstIteratorTypedef;
  typedef typename std::vector<const V*>::iterator       seqVectorPositionIteratorTypedef;
  uqSequenceOfVectorsClass(const uqVectorSpaceClass<V,M>& vectorSpace,
                           unsigned int                   sequenceSize,
                           const std::string&             name);
 ~uqSequenceOfVectorsClass();

        unsigned int sequenceSize              () const;
        void         resizeSequence            (unsigned int newSequenceSize);
        void         resetValues               (unsigned int initialPos, unsigned int numPos);
        void         erasePositions            (unsigned int initialPos, unsigned int numPos);
#ifdef UQ_SEQ_VEC_USES_OPERATOR
	const V*     operator[]                (unsigned int posId) const;
	const V*&    operator[]                (unsigned int posId);
#endif
        void         getPositionValues         (unsigned int posId,       V& vec) const;
        void         setPositionValues         (unsigned int posId, const V& vec);
        void         setGaussian               (const gsl_rng* rng, const V& meanVec, const V& stdDevVec);
        void         setUniform                (const gsl_rng* rng, const V& aVec,    const V& bVec     );
        void         uniformlySampledMdf       (const V&                            numEvaluationPointsVec,
                                                uqArrayOfOneDGridsClass <V,M>&      mdfGrids,
                                                uqArrayOfOneDTablesClass<V,M>&      mdfValues) const;
        void         uniformlySampledCdf       (const V&                            numEvaluationPointsVec,
                                                uqArrayOfOneDGridsClass <V,M>&      cdfGrids,
                                                uqArrayOfOneDTablesClass<V,M>&      cdfValues) const;
        void         unifiedUniformlySampledCdf(const V&                            numEvaluationPointsVec,
                                                uqArrayOfOneDGridsClass <V,M>&      unifiedCdfGrids,
                                                uqArrayOfOneDTablesClass<V,M>&      unifiedCdfValues) const;

        void         mean                      (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                V&                                  meanVec) const;
        void         sampleVariance            (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                const V&                            meanVec,
                                                V&                                  samVec) const;
        void         populationVariance        (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                const V&                            meanVec,
                                                V&                                  popVec) const;
        void         autoCovariance            (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                const V&                            meanVec,
                                                unsigned int                        lag,
                                                V&                                  covVec) const;

        void         autoCorrViaDef            (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                unsigned int                        lag,
                                                V&                                  corrVec) const;
        void         autoCorrViaFft            (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                const std::vector<unsigned int>&    lags,
                                                std::vector<V*>&                    corrVecs) const;
        void         autoCorrViaFft            (unsigned int                        initialPos,
                                                unsigned int                        numPos,
                                                unsigned int                        numSum,
                                                V&                                  autoCorrsSumVec) const;
        void         bmm                       (unsigned int                        initialPos,
                                                unsigned int                        batchLength,
                                                V&                                  bmmVec) const;
        void         fftForward                (unsigned int                        initialPos,
                                                unsigned int                        fftSize,
                                                unsigned int                        paramId,
                                                std::vector<std::complex<double> >& fftResult) const;
      //void         fftInverse                (unsigned int fftSize);
        void         psd                       (unsigned int                        initialPos,
                                                unsigned int                        numBlocks,
                                                double                              hopSizeRatio,
                                                unsigned int                        paramId,
                                                std::vector<double>&                psdResult) const;
        void         psdAtZero                 (unsigned int                        initialPos,
                                                unsigned int                        numBlocks,
                                                double                              hopSizeRatio,
                                                V&                                  psdVec) const;
        void         geweke                    (unsigned int                        initialPos,
                                                double                              ratioNa,
                                                double                              ratioNb,
                                                V&                                  gewVec) const;
        void         minMax                    (unsigned int                        initialPos,
                                                V&                                  minVec,
                                                V&                                  maxVec) const;
        void         unifiedMinMax             (unsigned int                        initialPos,
                                                V&                                  unifiedMinVec,
                                                V&                                  unifiedMaxVec) const;
        void         histogram                 (unsigned int                        initialPos,
                                                const V&                            minVec,
                                                const V&                            maxVec,
                                                std::vector<V*>&                    centersForAllBins,
                                                std::vector<V*>&                    binsForAllParams) const;
        void         unifiedHistogram          (unsigned int                        initialPos,
                                                const V&                            unifiedMinVec,
                                                const V&                            unifiedMaxVec,
                                                std::vector<V*>&                    unifiedCentersForAllBins,
                                                std::vector<V*>&                    unifiedBinsForAllParams) const;
        void         interQuantileRange        (unsigned int                        initialPos,
                                                V&                                  iqrVec) const;
        void         unifiedInterQuantileRange (unsigned int                        initialPos,
                                                V&                                  unifiedIqrVec) const;
        void         scalesForKDE              (unsigned int                        initialPos,
                                                const V&                            iqrVec,
                                                V&                                  scaleVec) const;
        void         unifiedScalesForKDE       (unsigned int                        initialPos,
                                                const V&                            unifiedIqrVec,
                                                V&                                  unifiedScaleVec) const;
        void         gaussianKDE               (const V&                            evalParamVec,
                                                V&                                  densityVec) const;
        void         gaussianKDE               (unsigned int                        initialPos,
                                                const V&                            scaleVec,
                                                const std::vector<V*>&              evalParamVecs,
                                                std::vector<V*>&                    densityVecs) const;
        void         unifiedGaussianKDE        (unsigned int                        initialPos,
                                                const V&                            unifiedScaleVec,
                                                const std::vector<V*>&              unifiedEvalParamVecs,
                                                std::vector<V*>&                    unifiedDensityVecs) const;
        void         printContents             (std::ofstream&                      ofsvar) const;
        void         select                    (const std::vector<unsigned int>&    idsOfUniquePositions);
        void         filter                    (unsigned int                        initialPos,
                                                unsigned int                        spacing);

private:
        void         extractScalarSeq          (unsigned int                        initialPos,
                                                unsigned int                        spacing,
                                                unsigned int                        numPos,
                                                unsigned int                        paramId,
                                                uqScalarSequenceClass<double>&      scalarSeq) const;
        void         extractRawData            (unsigned int                        initialPos,
                                                unsigned int                        spacing,
                                                unsigned int                        numPos,
                                                unsigned int                        paramId,
                                                std::vector<double>&                rawData) const;

  std::vector<const V*> m_seq;

  using uqBaseVectorSequenceClass<V,M>::m_env;
  using uqBaseVectorSequenceClass<V,M>::m_vectorSpace;
  using uqBaseVectorSequenceClass<V,M>::m_name;
  using uqBaseVectorSequenceClass<V,M>::m_fftObj;
};

template <class V, class M>
uqSequenceOfVectorsClass<V,M>::uqSequenceOfVectorsClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   sequenceSize,
  const std::string&             name)
  :
  uqBaseVectorSequenceClass<V,M>(vectorSpace,sequenceSize,name),
  m_seq                         (sequenceSize,NULL)
{

  //if (m_env.subScreenFile()) {
  //  *m_env.subScreenFile() << "Entering uqSequenceOfVectorsClass<V,M>::constructor()"
  //                         << std::endl;
  //}

  //if (m_env.subScreenFile()) {
  //  *m_env.subScreenFile() << "Leaving uqSequenceOfVectorsClass<V,M>::constructor()"
  //                         << std::endl;
  //}
}

template <class V, class M>
uqSequenceOfVectorsClass<V,M>::~uqSequenceOfVectorsClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_seq.size(); ++i) {
    if (m_seq[i]) delete m_seq[i];
  }
}

template <class V, class M>
unsigned int
uqSequenceOfVectorsClass<V,M>::sequenceSize() const
{
  return m_seq.size();
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::resizeSequence(unsigned int newSequenceSize)
{
  if (newSequenceSize != this->sequenceSize()) {
    m_seq.resize(newSequenceSize,NULL);
    std::vector<const V*>(m_seq).swap(m_seq);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

 return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::resetValues(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  if ((bRC == false) && (m_env.subScreenFile())) {
    *m_env.subScreenFile() << "In uqSequenceOfVectorsClass<V,M>::resetValues()"
                           << ", initialPos = "           << initialPos
                           << ", this->sequenceSize() = " << this->sequenceSize()
                           << ", numPos = "               << numPos
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::resetValues()",
                      "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) {
      delete m_seq[initialPos+j];
      m_seq[initialPos+j] = NULL;
    }
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::erasePositions(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::erasePositions()",
                      "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) delete m_seq[initialPos+j];
  }

  seqVectorPositionIteratorTypedef posIteratorBegin = m_seq.begin();
  if (initialPos < this->sequenceSize()) std::advance(posIteratorBegin,initialPos);
  else                                   posIteratorBegin = m_seq.end();

  unsigned int posEnd = initialPos + numPos - 1;
  seqVectorPositionIteratorTypedef posIteratorEnd = m_seq.begin();
  if (posEnd < this->sequenceSize()) std::advance(posIteratorEnd,posEnd);
  else                               posIteratorEnd = m_seq.end();

  unsigned int oldSequenceSize = this->sequenceSize();
  m_seq.erase(posIteratorBegin,posIteratorEnd);
  UQ_FATAL_TEST_MACRO((oldSequenceSize - numPos) != this->sequenceSize(),
                      m_env.rank(),
                      "uqSequenceOfVectors::erasePositions()",
                      "(oldSequenceSize - numPos) != this->sequenceSize()");

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

#ifdef UQ_SEQ_VEC_USES_OPERATOR
template <class V, class M>
const V*
uqSequenceOfVectorsClass<V,M>::operator[](unsigned int posId) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V,M>::operator[] const",
                      "posId > sequenceSize()");

  return (const V*) (m_seq[posId]);
}

template <class V, class M>
const V*&
uqSequenceOfVectorsClass<V,M>::operator[](unsigned int posId)
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V,M>::operator[] const",
                      "posId > sequenceSize()");

  return m_seq[posId];
}
#endif

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::getPositionValues(unsigned int posId, V& vec) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V,M>::getPositionValues()",
                      "posId > sequenceSize()");

  vec = *(const_cast<V*>(m_seq[posId]));

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::setPositionValues(unsigned int posId, const V& vec)
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V,M>::setPositionValues()",
                      "posId > sequenceSize()");

  if (m_seq[posId] != NULL) delete m_seq[posId];
  m_seq[posId] = new V(vec);

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::setGaussian(const gsl_rng* rng, const V& meanVec, const V& stdDevVec)
{
  V gaussianVector(m_vectorSpace.zeroVector());
  for (unsigned int j = 0; j < this->sequenceSize(); ++j) {
    gaussianVector.cwSetGaussian(m_env.rng(),meanVec,stdDevVec);
    if (m_seq[j] != NULL) delete m_seq[j];
    m_seq[j] = new V(gaussianVector);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}


template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::setUniform(const gsl_rng* rng, const V& aVec, const V& bVec)
{
  V uniformVector(m_vectorSpace.zeroVector());
  for (unsigned int j = 0; j < this->sequenceSize(); ++j) {
    uniformVector.cwSetUniform(m_env.rng(),aVec,bVec);
    if (m_seq[j] != NULL) delete m_seq[j];
    m_seq[j] = new V(uniformVector);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::uniformlySampledMdf(
  const V&                       numEvaluationPointsVec,
  uqArrayOfOneDGridsClass <V,M>& mdfGrids,
  uqArrayOfOneDTablesClass<V,M>& mdfValues) const
{
  V minDomainValues(m_vectorSpace.zeroVector());
  V maxDomainValues(m_vectorSpace.zeroVector());

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,              // initialPos
                           1,              // spacing
                           sequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aMdf(0);
    data.uniformlySampledMdf((unsigned int) numEvaluationPointsVec[i],
                             minDomainValues[i],
                             maxDomainValues[i],
                             aMdf);
    mdfValues.setOneDTable(i,aMdf);
  }

  mdfGrids.setUniformGrids(numEvaluationPointsVec,
                           minDomainValues,
                           maxDomainValues);

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::uniformlySampledCdf(
  const V&                       numEvaluationPointsVec,
  uqArrayOfOneDGridsClass <V,M>& cdfGrids,
  uqArrayOfOneDTablesClass<V,M>& cdfValues) const
{
  V minDomainValues(m_vectorSpace.zeroVector());
  V maxDomainValues(m_vectorSpace.zeroVector());

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,              // initialPos
                           1,              // spacing
                           sequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aCdf(0);
    data.uniformlySampledCdf((unsigned int) numEvaluationPointsVec[i],
                             minDomainValues[i],
                             maxDomainValues[i],
                             aCdf);
    cdfValues.setOneDTable(i,aCdf);
  }

  cdfGrids.setUniformGrids(numEvaluationPointsVec,
                           minDomainValues,
                           maxDomainValues);

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedUniformlySampledCdf(
  const V&                       numEvaluationPointsVec,
  uqArrayOfOneDGridsClass <V,M>& unifiedCdfGrids,
  uqArrayOfOneDTablesClass<V,M>& unifiedCdfValues) const
{
  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 10)) {
    *m_env.subScreenFile() << "Entering uqSequenceOfVectorsClass<V,M>::unifiedUniformlySampledCdf()"
                           << std::endl;
  }

  V unifiedMinDomainValues(m_vectorSpace.zeroVector());
  V unifiedMaxDomainValues(m_vectorSpace.zeroVector());

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,              // initialPos
                           1,              // spacing
                           sequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aCdf(0);
    data.unifiedUniformlySampledCdf(m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1,
                                    (unsigned int) numEvaluationPointsVec[i],
                                    unifiedMinDomainValues[i],
                                    unifiedMaxDomainValues[i],
                                    aCdf);
    unifiedCdfValues.setOneDTable(i,aCdf);
  }

  unifiedCdfGrids.setUniformGrids(numEvaluationPointsVec,
                                  unifiedMinDomainValues,
                                  unifiedMaxDomainValues);

  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 10)) {
    *m_env.subScreenFile() << "Leaving uqSequenceOfVectorsClass<V,M>::unifiedUniformlySampledCdf()"
                           << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           meanVec) const
{
  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 5)) {
    *m_env.subScreenFile() << "Entering uqSequenceOfVectorsClass<V,M>::mean()"
                           << ": initialPos = "         << initialPos
                           << ", numPos = "             << numPos
                           << ", full sequence size = " << this->sequenceSize()
                           << std::endl;
  }

  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ));
  if ((bRC == false) && (m_env.subScreenFile())) {
    *m_env.subScreenFile() << "In uqSequenceOfVectorsClass<V,M>::mean()"
                           << ", initialPos = "           << initialPos
                           << ", this->sequenceSize() = " << this->sequenceSize()
                           << ", numPos = "               << numPos
                           << ", this->vectorSize() = "   << this->vectorSize()
                           << ", meanVec.size() = "       << meanVec.size()
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::mean()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    meanVec[i] = data.mean(0,
                           numPos);
  }

  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 5)) {
    *m_env.subScreenFile() << "Leaving uqSequenceOfVectorsClass<V,M>::mean()"
                           << ": initialPos = "         << initialPos
                           << ", numPos = "             << numPos
                           << ", full sequence size = " << this->sequenceSize()
                           << ", meanVec = "            << meanVec
                           << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ) &&
              (this->vectorSize()  == samVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::sampleVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    samVec[i] = data.sampleVariance(0,
                                    numPos,
                                    meanVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ) &&
              (this->vectorSize()  == popVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::populationVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    popVec[i] = data.populationVariance(0,
                                        numPos,
                                        meanVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  unsigned int lag,
  V&           covVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ) &&
              (lag                 <  numPos              ) && // lag should not be too large
              (this->vectorSize()  == covVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCovariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    covVec[i] = data.autoCovariance(0,
                                    numPos,
                                    meanVec[i],
                                    lag);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::autoCorrViaDef(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int lag,
  V&           corrVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (lag                 <  numPos              ) && // lag should not be too large
              (this->vectorSize()  == corrVec.size()      ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCorrViaDef()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    corrVec[i] = data.autoCorrViaDef(0,
                                     numPos,
                                     lag);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::autoCorrViaFft(
  unsigned int                     initialPos,
  unsigned int                     numPos,
  const std::vector<unsigned int>& lags,
  std::vector<V*>&                 corrVecs) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (0                   < lags.size()          ) &&
              (lags[lags.size()-1] <  numPos              )); // lag should not be too large
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCorrViaFft()",
                      "invalid input data");

  for (unsigned int j = lags.size(); j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }
  corrVecs.resize(lags.size(),NULL);
  for (unsigned int j = 0;           j < corrVecs.size(); ++j) {
    if (corrVecs[j] == NULL) corrVecs[j] = new V(m_vectorSpace.zeroVector());
  }

  uqScalarSequenceClass<double> data(m_env,0);
  unsigned int maxLag = lags[lags.size()-1];
  std::vector<double> autoCorrs(maxLag+1,0.); // Yes, +1

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    //if (m_env.subScreenFile()) {
    //  *m_env.subScreenFile() << "In uqSequenceOfVectorsClass<V,M>::autoCorrViaFft()"
    //                         << ": about to call data.autoCorrViaFft() for paramId = " << i
    //                         << ", with numPos = "      << numPos
    //                         << ", maxLag = "           << maxLag
    //                         << ", autoCorrs.size() = " << autoCorrs.size()
    //                         << std::endl;
    //}
    data.autoCorrViaFft(0,
                        numPos,
                        maxLag,
                        autoCorrs);

    for (unsigned int j = 0; j < lags.size(); ++j) {
      (*(corrVecs[j]))[i] = autoCorrs[lags[j]];
    }
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::autoCorrViaFft(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int numSum,
  V&           autoCorrsSumVec) const
{
  bool bRC = ((initialPos             <  this->sequenceSize()) &&
              (0                      <  numPos              ) &&
              ((initialPos+numPos)    <= this->sequenceSize()) &&
              (0                      <  numSum              ) &&
              (numSum                 <= numPos              ) &&
              (autoCorrsSumVec.size() == this->vectorSize()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCorrViaFft(), for sum",
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
                        numSum,
                        autoCorrsSumVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::bmm(
  unsigned int initialPos,
  unsigned int batchLength,
  V&           bmmVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()            ) &&
              (batchLength         < (this->sequenceSize()-initialPos)) &&
              (this->vectorSize()  == bmmVec.size()                   ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::bmm()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->sequenceSize()-initialPos,
                           i,
                           data);
    bmmVec[i] = data.bmm(0,
                         batchLength);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::fftForward(
  unsigned int                        initialPos,
  unsigned int                        fftSize,
  unsigned int                        paramId,
  std::vector<std::complex<double> >& fftResult) const
{
  bool bRC = ((initialPos           <  this->sequenceSize()) &&
              (paramId              <  this->vectorSize()  ) &&
              (0                    <  fftSize             ) &&
              ((initialPos+fftSize) <= this->sequenceSize()) &&
              (fftSize              <  this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::fftForward()",
                      "invalid input data");

  std::vector<double> rawData(fftSize,0.);
  this->extractRawData(initialPos,
                       1, // spacing
                       fftSize,
                       paramId,
                       rawData);

  m_fftObj->forward(rawData,fftSize,fftResult);

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::psd(
  unsigned int         initialPos,
  unsigned int         numBlocks,
  double               hopSizeRatio,
  unsigned int         paramId,
  std::vector<double>& psdResult) const
{
  bool bRC = ((initialPos < this->sequenceSize()) &&
              (paramId    < this->vectorSize()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::psd()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  this->extractScalarSeq(initialPos,
                         1, // spacing
                         this->sequenceSize()-initialPos,
                         paramId,
                         data);
  data.psd(0,
           numBlocks,
           hopSizeRatio,
           psdResult);

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::psdAtZero(
  unsigned int initialPos,
  unsigned int numBlocks,
  double       hopSizeRatio,
  V&           psdVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == psdVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::psdAtZero()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);
  std::vector<double> psdResult(0,0.); // size will be determined by 'data.psd()'

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->sequenceSize()-initialPos,
                           i,
                           data);
    data.psd(0,
             numBlocks,
             hopSizeRatio,
             psdResult);
    psdVec[i] = psdResult[0];
    //*m_env.subScreenFile() << "psdResult[0] = " << psdResult[0] << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::geweke(
  unsigned int initialPos,
  double       ratioNa,
  double       ratioNb,
  V&           gewVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == gewVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::geweke()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
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

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::minMax(
  unsigned int initialPos,
  V&           minVec,
  V&           maxVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == minVec.size()       ) &&
              (this->vectorSize() == maxVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::minMax()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  unsigned int numParams = this->vectorSize();
  uqScalarSequenceClass<double> data(m_env,0);

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

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedMinMax(
  unsigned int initialPos,
  V&           unifiedMinVec,
  V&           unifiedMaxVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == unifiedMinVec.size()) &&
              (this->vectorSize() == unifiedMaxVec.size()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedMinMax()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  unsigned int numParams = this->vectorSize();
  uqScalarSequenceClass<double> data(m_env,0);

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.unifiedMinMax(m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1,
                       0,
                       unifiedMinVec[i],
                       unifiedMaxVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::histogram(
  unsigned int     initialPos,
  const V&         minVec,
  const V&         maxVec,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& binsForAllParams) const
{
  bool bRC = ((initialPos               <  this->sequenceSize()    ) &&
              (this->vectorSize()       == minVec.size()           ) &&
              (this->vectorSize()       == maxVec.size()           ) &&
              (0                        <  centersForAllBins.size()) &&
              (centersForAllBins.size() == binsForAllParams.size() ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::histogram()",
                      "invalid input data");

  for (unsigned int j = 0; j < binsForAllParams.size(); ++j) {
    centersForAllBins[j] = new V(m_vectorSpace.zeroVector());
    binsForAllParams [j] = new V(m_vectorSpace.zeroVector());
  }

  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > centers(centersForAllBins.size(),0.);
    std::vector<unsigned int> bins   (binsForAllParams.size(), 0 );
    data.histogram(0,
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

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedHistogram(
  unsigned int     initialPos,
  const V&         unifiedMinVec,
  const V&         unifiedMaxVec,
  std::vector<V*>& unifiedCentersForAllBins,
  std::vector<V*>& unifiedBinsForAllParams) const
{
  bool bRC = ((initialPos                      <  this->sequenceSize()           ) &&
              (this->vectorSize()              == unifiedMinVec.size()           ) &&
              (this->vectorSize()              == unifiedMaxVec.size()           ) &&
              (0                               <  unifiedCentersForAllBins.size()) &&
              (unifiedCentersForAllBins.size() == unifiedBinsForAllParams.size() ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedHistogram()",
                      "invalid input data");

  for (unsigned int j = 0; j < unifiedBinsForAllParams.size(); ++j) {
    unifiedCentersForAllBins[j] = new V(m_vectorSpace.zeroVector());
    unifiedBinsForAllParams [j] = new V(m_vectorSpace.zeroVector());
  }

  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > unifiedCenters(unifiedCentersForAllBins.size(),0.);
    std::vector<unsigned int> unifiedBins   (unifiedBinsForAllParams.size(), 0 );
    data.unifiedHistogram(m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1,
                          0,
                          unifiedMinVec[i],
                          unifiedMaxVec[i],
                          unifiedCenters,
                          unifiedBins);

    for (unsigned int j = 0; j < unifiedBins.size(); ++j) {
      (*(unifiedCentersForAllBins[j]))[i] = unifiedCenters[j];
      (*(unifiedBinsForAllParams [j]))[i] = (double) unifiedBins[j];
    }
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::interQuantileRange(
  unsigned int initialPos,
  V&           iqrVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == iqrVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::interQuantileRange()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    iqrVec[i] = data.interQuantileRange(0);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedInterQuantileRange(
  unsigned int initialPos,
  V&           unifiedIqrVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == unifiedIqrVec.size()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedInterQuantileRange()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedIqrVec[i] = data.unifiedInterQuantileRange(m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1,
                                                      0);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::scalesForKDE(
  unsigned int initialPos,
  const V&     iqrVec,
  V&           scaleVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == iqrVec.size()       ) &&
              (this->vectorSize() == scaleVec.size()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::scalesForKDE()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    scaleVec[i] = data.scaleForKDE(0,
                                   iqrVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedScalesForKDE(
  unsigned int initialPos,
  const V&     unifiedIqrVec,
  V&           unifiedScaleVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()  ) &&
              (this->vectorSize() == unifiedIqrVec.size()  ) &&
              (this->vectorSize() == unifiedScaleVec.size()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedScalesForKDE()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedScaleVec[i] = data.unifiedScaleForKDE(m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1,
                                                 0,
                                                 unifiedIqrVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::gaussianKDE(
  const V& evalParamVec,
        V& densityVec) const
{
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0, // Use the whole chain
                           1, // spacing
                           this->sequenceSize(),
                           i,
                           data);

    densityVec[i] = data.gaussianKDE(evalParamVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::gaussianKDE(
  unsigned int           initialPos,
  const V&               scaleVec,
  const std::vector<V*>& evalParamVecs,
  std::vector<V*>&       densityVecs) const
{
  bool bRC = ((initialPos           <  this->sequenceSize()) &&
              (this->vectorSize()   == scaleVec.size()     ) &&
              (0                    <  evalParamVecs.size()) &&
              (evalParamVecs.size() == densityVecs.size()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::gaussianKDE()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numEvals = evalParamVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    densityVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> evalParams(numEvals,0.);
  std::vector<double> densities  (numEvals,0.);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      evalParams[j] = (*evalParamVecs[j])[i];
    }

    data.gaussianKDE(0,
                     scaleVec[i],
                     evalParams,
                     densities);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*densityVecs[j])[i] = densities[j];
    }
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedGaussianKDE(
  unsigned int           initialPos,
  const V&               unifiedScaleVec,
  const std::vector<V*>& unifiedEvalParamVecs,
  std::vector<V*>&       unifiedDensityVecs) const
{
  bool bRC = ((initialPos                  <  this->sequenceSize()       ) &&
              (this->vectorSize()          == unifiedScaleVec.size()     ) &&
              (0                           <  unifiedEvalParamVecs.size()) &&
              (unifiedEvalParamVecs.size() == unifiedDensityVecs.size()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V,M>::gaussianKDE()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numEvals = unifiedEvalParamVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    unifiedDensityVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> unifiedEvalParams(numEvals,0.);
  std::vector<double> unifiedDensities (numEvals,0.);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      unifiedEvalParams[j] = (*unifiedEvalParamVecs[j])[i];
    }

    data.unifiedGaussianKDE(m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1,
                            0,
                            unifiedScaleVec[i],
                            unifiedEvalParams,
                            unifiedDensities);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*unifiedDensityVecs[j])[i] = unifiedDensities[j];
    }
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::select(const std::vector<unsigned int>& idsOfUniquePositions)
{
  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::filter(
  unsigned int initialPos,
  unsigned int spacing)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Entering uqSequenceOfVectorsClass<V,M>::filter()"
                           << ": initialPos = "   << initialPos
                           << ", spacing = "      << spacing
                           << ", sequenceSize = " << this->sequenceSize()
                           << std::endl;
  }

  unsigned int i = 0;
  unsigned int j = initialPos;
  unsigned int originalSequenceSize = this->sequenceSize();
  while (j < originalSequenceSize) {
    if (i != j) {
      //*m_env.subScreenFile() << i << "--" << j << " ";
      delete m_seq[i];
      m_seq[i] = new V(*(m_seq[j]));
    }
    i++;
    j += spacing;
  }

  this->resetValues(i,originalSequenceSize-i);
  this->resizeSequence(i);

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Leaving uqSequenceOfVectorsClass<V,M>::filter()"
                           << ": initialPos = "   << initialPos
                           << ", spacing = "      << spacing
                           << ", sequenceSize = " << this->sequenceSize()
                           << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::printContents(std::ofstream& ofsvar) const
{
  ofsvar << m_name << "_subenv" << m_env.subIdString() << " = zeros(" << this->sequenceSize()
         << ","                                                       << this->vectorSize()
         << ");"
         << std::endl;
  ofsvar << m_name << "_subenv" << m_env.subIdString() << " = [";
  unsigned int chainSize = this->sequenceSize();
  for (unsigned int j = 0; j < chainSize; ++j) {
    bool savedVectorPrintState = m_seq[j]->getPrintHorizontally();
    m_seq[j]->setPrintHorizontally(true);
    ofsvar << *(m_seq[j])
           << std::endl;
    m_seq[j]->setPrintHorizontally(savedVectorPrintState);
  }
  ofsvar << "];\n";

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::extractScalarSeq(
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

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::extractRawData(
  unsigned int         initialPos,
  unsigned int         spacing,
  unsigned int         numPos,
  unsigned int         paramId,
  std::vector<double>& rawData) const
{
  rawData.resize(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = (*(m_seq[initialPos+j        ]))[paramId];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = (*(m_seq[initialPos+j*spacing]))[paramId];
    }
  }

  return;
}
#endif // __UQ_SEQUENCE_OF_VECTORS_H__

