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
                           unsigned int                   subSequenceSize,
                           const std::string&             name);
 ~uqSequenceOfVectorsClass();

  uqSequenceOfVectorsClass<V,M>& operator= (const uqSequenceOfVectorsClass<V,M>& rhs);

        unsigned int subSequenceSize            () const;
        void         resizeSequence             (unsigned int newSubSequenceSize);
        void         resetValues                (unsigned int initialPos, unsigned int numPos);
        void         erasePositions             (unsigned int initialPos, unsigned int numPos);
#ifdef UQ_SEQ_VEC_USES_OPERATOR
	const V*     operator[]                 (unsigned int posId) const;
	const V*&    operator[]                 (unsigned int posId);
#endif
        void         getPositionValues          (unsigned int posId,       V& vec) const;
        void         setPositionValues          (unsigned int posId, const V& vec);
        void         setGaussian                (const gsl_rng* rng, const V& meanVec, const V& stdDevVec);
        void         setUniform                 (const gsl_rng* rng, const V& aVec,    const V& bVec     );
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
        void         subUniformlySampledMdf     (const V&                             numEvaluationPointsVec,
                                                 uqArrayOfOneDGridsClass <V,M>&       mdfGrids,
                                                 uqArrayOfOneDTablesClass<V,M>&       mdfValues) const;
#endif
        void         subUniformlySampledCdf     (const V&                             numEvaluationPointsVec,
                                                 uqArrayOfOneDGridsClass <V,M>&       cdfGrids,
                                                 uqArrayOfOneDTablesClass<V,M>&       cdfValues) const;
        void         unifiedUniformlySampledCdf (const V&                             numEvaluationPointsVec,
                                                 uqArrayOfOneDGridsClass <V,M>&       unifiedCdfGrids,
                                                 uqArrayOfOneDTablesClass<V,M>&       unifiedCdfValues) const;

        void         subMean                    (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 V&                                   meanVec) const;
        void         unifiedMean                (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 V&                                   unifiedMeanVec) const;
        void         subSampleVariance          (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 const V&                             meanVec,
                                                 V&                                   samVec) const;
        void         unifiedSampleVariance      (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 const V&                             unifiedMeanVec,
                                                 V&                                   unifiedSamVec) const;
        void         subPopulationVariance      (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 const V&                             meanVec,
                                                 V&                                   popVec) const;
        void         autoCovariance             (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 const V&                             meanVec,
                                                 unsigned int                         lag,
                                                 V&                                   covVec) const;
 
        void         autoCorrViaDef             (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 unsigned int                         lag,
                                                 V&                                   corrVec) const;
        void         autoCorrViaFft             (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 const std::vector<unsigned int>&     lags,
                                                 std::vector<V*>&                     corrVecs) const;
        void         autoCorrViaFft             (unsigned int                         initialPos,
                                                 unsigned int                         numPos,
                                                 unsigned int                         numSum,
                                                 V&                                   autoCorrsSumVec) const;
        void         bmm                        (unsigned int                         initialPos,
                                                 unsigned int                         batchLength,
                                                 V&                                   bmmVec) const;
        void         fftForward                 (unsigned int                         initialPos,
                                                 unsigned int                         fftSize,
                                                 unsigned int                         paramId,
                                                 std::vector<std::complex<double> >&  fftResult) const;
      //void         fftInverse                 (unsigned int fftSize);
        void         psd                        (unsigned int                         initialPos,
                                                 unsigned int                         numBlocks,
                                                 double                               hopSizeRatio,
                                                 unsigned int                         paramId,
                                                 std::vector<double>&                 psdResult) const;
        void         psdAtZero                  (unsigned int                         initialPos,
                                                 unsigned int                         numBlocks,
                                                 double                               hopSizeRatio,
                                                 V&                                   psdVec) const;
        void         geweke                     (unsigned int                         initialPos,
                                                 double                               ratioNa,
                                                 double                               ratioNb,
                                                 V&                                   gewVec) const;
        void         meanStacc                  (unsigned int                         initialPos,
                                                 V&                                   meanStaccVec) const;
        void         subMinMax                  (unsigned int                         initialPos,
                                                 V&                                   minVec,
                                                 V&                                   maxVec) const;
        void         unifiedMinMax              (unsigned int                         initialPos,
                                                 V&                                   unifiedMinVec,
                                                 V&                                   unifiedMaxVec) const;
        void         subHistogram               (unsigned int                         initialPos,
                                                 const V&                             minVec,
                                                 const V&                             maxVec,
                                                 std::vector<V*>&                     centersForAllBins,
                                                 std::vector<V*>&                     quanttsForAllBins) const;
        void         unifiedHistogram           (unsigned int                         initialPos,
                                                 const V&                             unifiedMinVec,
                                                 const V&                             unifiedMaxVec,
                                                 std::vector<V*>&                     unifiedCentersForAllBins,
                                                 std::vector<V*>&                     unifiedQuanttsForAllBins) const;
        void         subCdfStacc                (unsigned int                         initialPos,
                                                 const std::vector<V*>&               evalPositionsVecs,
                                                 std::vector<V*>&                     cdfStaccVecs) const;
        void         subInterQuantileRange      (unsigned int                         initialPos,
                                                 V&                                   iqrVec) const;
        void         unifiedInterQuantileRange  (unsigned int                         initialPos,
                                                 V&                                   unifiedIqrVec) const;
        void         subScalesForKDE            (unsigned int                         initialPos,
                                                 const V&                             iqrVec,
                                                 V&                                   scaleVec) const;
        void         unifiedScalesForKDE        (unsigned int                         initialPos,
                                                 const V&                             unifiedIqrVec,
                                                 V&                                   unifiedScaleVec) const;
      //void         sabGaussianKDE             (const V&                             evalParamVec,
      //                                         V&                                   densityVec) const;
        void         subGaussianKDE             (unsigned int                         initialPos,
                                                 const V&                             scaleVec,
                                                 const std::vector<V*>&               evalParamVecs,
                                                 std::vector<V*>&                     densityVecs) const;
        void         unifiedGaussianKDE         (unsigned int                         initialPos,
                                                 const V&                             unifiedScaleVec,
                                                 const std::vector<V*>&               unifiedEvalParamVecs,
                                                 std::vector<V*>&                     unifiedDensityVecs) const;
      //void         subWriteContents           (std::ofstream&                       ofsvar) const;
      //void         unifiedWriteContents       (std::ofstream&                       ofsvar) const;
        void         subWriteContents           (const std::string&                   fileName,
                                                 const std::set<unsigned int>&        allowedSubEnvIds) const;
        void         unifiedWriteContents       (const std::string&                   fileName) const;
        void         unifiedReadContents        (const std::string&                   fileName,
                                                 const unsigned int                   subSequenceSize);
        void         select                     (const std::vector<unsigned int>&     idsOfUniquePositions);
        void         filter                     (unsigned int                         initialPos,
                                                 unsigned int                         spacing);

        double       estimateConvBrooksGelman   (unsigned int                         initialPos,
                                                 unsigned int                         numPos) const;

        void         append                     (const uqSequenceOfVectorsClass<V,M>& src,
                                                 unsigned int                         initialPos,
                                                 unsigned int                         numPos);
private:
        void         copy                       (const uqSequenceOfVectorsClass<V,M>& src);
        void         extractScalarSeq           (unsigned int                         initialPos,
                                                 unsigned int                         spacing,
                                                 unsigned int                         numPos,
                                                 unsigned int                         paramId,
                                                 uqScalarSequenceClass<double>&       scalarSeq) const;
        void         extractRawData             (unsigned int                         initialPos,
                                                 unsigned int                         spacing,
                                                 unsigned int                         numPos,
                                                 unsigned int                         paramId,
                                                 std::vector<double>&                 rawData) const;

  std::vector<const V*> m_seq;

  using uqBaseVectorSequenceClass<V,M>::m_env;
  using uqBaseVectorSequenceClass<V,M>::m_vectorSpace;
  using uqBaseVectorSequenceClass<V,M>::m_name;
  using uqBaseVectorSequenceClass<V,M>::m_fftObj;
};

template <class V, class M>
uqSequenceOfVectorsClass<V,M>::uqSequenceOfVectorsClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   subSequenceSize,
  const std::string&             name)
  :
  uqBaseVectorSequenceClass<V,M>(vectorSpace,subSequenceSize,name),
  m_seq                         (subSequenceSize,NULL)
{

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::constructor()"
  //                         << std::endl;
  //}

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::constructor()"
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
uqSequenceOfVectorsClass<V,M>&
uqSequenceOfVectorsClass<V,M>::operator= (const uqSequenceOfVectorsClass<V,M>& rhs)
{
  this->copy(rhs);
  return *this;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::copy(const uqSequenceOfVectorsClass<V,M>& src)
{
  uqBaseVectorSequenceClass<V,M>::copy(src);
  for (unsigned int i = 0; i < (unsigned int) m_seq.size(); ++i) {
    if (m_seq[i]) delete m_seq[i];
  }
  m_seq.resize(src.subSequenceSize(),NULL);
  for (unsigned int i = 0; i < m_seq.size(); ++i) {
    m_seq[i] = new V(*(src.m_seq[i]));
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::append(
  const uqSequenceOfVectorsClass<V,M>& src,
  unsigned int                         initialPos,
  unsigned int                         numPos)
{
  UQ_FATAL_TEST_MACRO((src.subSequenceSize() < (initialPos+1)),
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<T>::append()",
                      "initialPos is too big");

  UQ_FATAL_TEST_MACRO((src.subSequenceSize() < (initialPos+numPos)),
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<T>::append()",
                      "numPos is too big");

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();
  unsigned int currentSize = this->subSequenceSize();
  m_seq.resize(currentSize+numPos,NULL);
  for (unsigned int i = 0; i < numPos; ++i) {
    m_seq[currentSize+i] = new V(*(src.m_seq[initialPos+i]));
  }

  return;
}

template <class V, class M>
unsigned int
uqSequenceOfVectorsClass<V,M>::subSequenceSize() const
{
  return m_seq.size();
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::resizeSequence(unsigned int newSubSequenceSize)
{
  if (newSubSequenceSize != this->subSequenceSize()) {
    if (newSubSequenceSize < this->subSequenceSize()) {
      this->resetValues(newSubSequenceSize,this->subSequenceSize()-newSubSequenceSize);
    }
    m_seq.resize(newSubSequenceSize,NULL);
    std::vector<const V*>(m_seq).swap(m_seq);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

 return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::resetValues(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqSequenceOfVectorsClass<V,M>::resetValues()"
                           << ", initialPos = "              << initialPos
                           << ", this->subSequenceSize() = " << this->subSequenceSize()
                           << ", numPos = "                  << numPos
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
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
  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::erasePositions()",
                      "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) delete m_seq[initialPos+j];
  }

  seqVectorPositionIteratorTypedef posIteratorBegin = m_seq.begin();
  if (initialPos < this->subSequenceSize()) std::advance(posIteratorBegin,initialPos);
  else                                      posIteratorBegin = m_seq.end();

  unsigned int posEnd = initialPos + numPos - 1;
  seqVectorPositionIteratorTypedef posIteratorEnd = m_seq.begin();
  if (posEnd < this->subSequenceSize()) std::advance(posIteratorEnd,posEnd);
  else                                  posIteratorEnd = m_seq.end();

  unsigned int oldSubSequenceSize = this->subSequenceSize();
  m_seq.erase(posIteratorBegin,posIteratorEnd);
  UQ_FATAL_TEST_MACRO((oldSubSequenceSize - numPos) != this->subSequenceSize(),
                      m_env.fullRank(),
                      "uqSequenceOfVectors::erasePositions()",
                      "(oldSubSequenceSize - numPos) != this->subSequenceSize()");

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

#ifdef UQ_SEQ_VEC_USES_OPERATOR
template <class V, class M>
const V*
uqSequenceOfVectorsClass<V,M>::operator[](unsigned int posId) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->subSequenceSize()),
                      m_env.fullRank(),
                      "uqSequenceOfVectorss<V,M>::operator[] const",
                      "posId > subSequenceSize()");

  return (const V*) (m_seq[posId]);
}

template <class V, class M>
const V*&
uqSequenceOfVectorsClass<V,M>::operator[](unsigned int posId)
{
  UQ_FATAL_TEST_MACRO((posId >= this->subSequenceSize()),
                      m_env.fullRank(),
                      "uqSequenceOfVectorss<V,M>::operator[] const",
                      "posId > subSequenceSize()");

  return m_seq[posId];
}
#endif

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::getPositionValues(unsigned int posId, V& vec) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->subSequenceSize()),
                      m_env.fullRank(),
                      "uqSequenceOfVectorss<V,M>::getPositionValues()",
                      "posId > subSequenceSize()");

  vec = *(const_cast<V*>(m_seq[posId]));

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::setPositionValues(unsigned int posId, const V& vec)
{
  UQ_FATAL_TEST_MACRO((posId >= this->subSequenceSize()),
                      m_env.fullRank(),
                      "uqSequenceOfVectorss<V,M>::setPositionValues()",
                      "posId > subSequenceSize()");

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
  for (unsigned int j = 0; j < this->subSequenceSize(); ++j) {
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
  for (unsigned int j = 0; j < this->subSequenceSize(); ++j) {
    uniformVector.cwSetUniform(m_env.rng(),aVec,bVec);
    if (m_seq[j] != NULL) delete m_seq[j];
    m_seq[j] = new V(uniformVector);
  }

  uqBaseVectorSequenceClass<V,M>::deleteStoredVectors();

  return;
}

#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subUniformlySampledMdf(
  const V&                       numEvaluationPointsVec,
  uqArrayOfOneDGridsClass <V,M>& mdfGrids,
  uqArrayOfOneDTablesClass<V,M>& mdfValues) const
{
  V minDomainValues(m_vectorSpace.zeroVector());
  V maxDomainValues(m_vectorSpace.zeroVector());

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,                 // initialPos
                           1,                 // spacing
                           subSequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aMdf(0);
    data.subUniformlySampledMdf((unsigned int) numEvaluationPointsVec[i],
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
#endif
template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subUniformlySampledCdf(
  const V&                       numEvaluationPointsVec,
  uqArrayOfOneDGridsClass <V,M>& cdfGrids,
  uqArrayOfOneDTablesClass<V,M>& cdfValues) const
{
  V minDomainValues(m_vectorSpace.zeroVector());
  V maxDomainValues(m_vectorSpace.zeroVector());

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,                 // initialPos
                           1,                 // spacing
                           subSequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aCdf(0);
    data.subUniformlySampledCdf((unsigned int) numEvaluationPointsVec[i],
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
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::unifiedUniformlySampledCdf()"
                           << std::endl;
  }

  V unifiedMinDomainValues(m_vectorSpace.zeroVector());
  V unifiedMaxDomainValues(m_vectorSpace.zeroVector());

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,                 // initialPos
                           1,                 // spacing
                           subSequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aCdf(0);
    data.unifiedUniformlySampledCdf(m_vectorSpace.numOfProcsForStorage() == 1,
                                    (unsigned int) numEvaluationPointsVec[i],
                                    unifiedMinDomainValues[i],
                                    unifiedMaxDomainValues[i],
                                    aCdf);
    unifiedCdfValues.setOneDTable(i,aCdf);
  }

  unifiedCdfGrids.setUniformGrids(numEvaluationPointsVec,
                                  unifiedMinDomainValues,
                                  unifiedMaxDomainValues);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::unifiedUniformlySampledCdf()"
                           << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subMean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           meanVec) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::subMean()"
                           << ": initialPos = "         << initialPos
                           << ", numPos = "             << numPos
                           << ", full sequence size = " << this->subSequenceSize()
                           << std::endl;
  }

  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()) &&
              (this->vectorSizeLocal()  == meanVec.sizeLocal()         ));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqSequenceOfVectorsClass<V,M>::subMean()"
                           << ", initialPos = "              << initialPos
                           << ", this->subSequenceSize() = " << this->subSequenceSize()
                           << ", numPos = "                  << numPos
                           << ", this->vectorSizeLocal() = " << this->vectorSizeLocal()
                           << ", meanVec.sizeLocal() = "     << meanVec.sizeLocal()
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subMean()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    meanVec[i] = data.subMean(0,
                              numPos);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::subMean()"
                           << ": initialPos = "         << initialPos
                           << ", numPos = "             << numPos
                           << ", full sequence size = " << this->subSequenceSize()
                           << ", meanVec = "            << meanVec
                           << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedMean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           unifiedMeanVec) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::unifiedMean()"
                           << ": initialPos = "         << initialPos
                           << ", numPos = "             << numPos
                           << ", full sequence size = " << this->subSequenceSize()
                           << std::endl;
  }

  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()) &&
              (this->vectorSizeLocal()  == unifiedMeanVec.sizeLocal()  ));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqSequenceOfVectorsClass<V,M>::unifiedMean()"
                           << ", initialPos = "                 << initialPos
                           << ", this->subSequenceSize() = "    << this->subSequenceSize()
                           << ", numPos = "                     << numPos
                           << ", this->vectorSizeLocal() = "    << this->vectorSizeLocal()
                           << ", unifiedMeanVec.sizeLocal() = " << unifiedMeanVec.sizeLocal()
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedMean()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedMeanVec[i] = data.unifiedMean(m_vectorSpace.numOfProcsForStorage() == 1,
                                         0,
                                         numPos);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::unifiedMean()"
                           << ": initialPos = "         << initialPos
                           << ", numPos = "             << numPos
                           << ", full sequence size = " << this->subSequenceSize()
                           << ", unifiedMeanVec = "     << unifiedMeanVec
                           << std::endl;
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subSampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (this->vectorSizeLocal() == samVec.sizeLocal()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subSampleVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    samVec[i] = data.subSampleVariance(0,
                                       numPos,
                                       meanVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedSampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     unifiedMeanVec,
  V&           unifiedSamVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()   ) &&
              (0                       <  numPos                    ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()   ) &&
              (this->vectorSizeLocal() == unifiedMeanVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedSamVec.sizeLocal() ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedSampleVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedSamVec[i] = data.unifiedSampleVariance(m_vectorSpace.numOfProcsForStorage() == 1,
                                                  0,
                                                  numPos,
                                                  unifiedMeanVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subPopulationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (this->vectorSizeLocal() == popVec.sizeLocal()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subPopulationVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    popVec[i] = data.subPopulationVariance(0,
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
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (lag                     <  numPos                 ) && // lag should not be too large
              (this->vectorSizeLocal() == covVec.sizeLocal()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCovariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
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
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (lag                     <  numPos                 ) && // lag should not be too large
              (this->vectorSizeLocal() == corrVec.sizeLocal()    ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCorrViaDef()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
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
  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()) &&
              (0                   <  lags.size()            ) &&
              (lags[lags.size()-1] <  numPos                 )); // lag should not be too large
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCorrViaFft()",
                      "invalid input data");

  for (unsigned int j = lags.size(); j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }
  corrVecs.resize(lags.size(),NULL);
  for (unsigned int j = 0;           j < corrVecs.size(); ++j) {
    if (corrVecs[j] == NULL) corrVecs[j] = new V(m_vectorSpace.zeroVector());
  }

  uqScalarSequenceClass<double> data(m_env,0,"");
  unsigned int maxLag = lags[lags.size()-1];
  std::vector<double> autoCorrs(maxLag+1,0.); // Yes, +1

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqSequenceOfVectorsClass<V,M>::autoCorrViaFft()"
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
  bool bRC = ((initialPos             <  this->subSequenceSize()) &&
              (0                      <  numPos                 ) &&
              ((initialPos+numPos)    <= this->subSequenceSize()) &&
              (0                      <  numSum                 ) &&
              (numSum                 <= numPos                 ) &&
              (autoCorrsSumVec.sizeLocal() == this->vectorSizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::autoCorrViaFft(), for sum",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
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
  bool bRC = ((initialPos              <  this->subSequenceSize()            ) &&
              (batchLength             < (this->subSequenceSize()-initialPos)) &&
              (this->vectorSizeLocal() == bmmVec.sizeLocal()                 ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::bmm()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->subSequenceSize()-initialPos,
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
  bool bRC = ((initialPos           <  this->subSequenceSize()) &&
              (paramId              <  this->vectorSizeLocal()) &&
              (0                    <  fftSize                ) &&
              ((initialPos+fftSize) <= this->subSequenceSize()) &&
              (fftSize              <  this->subSequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
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
  bool bRC = ((initialPos < this->subSequenceSize()) &&
              (paramId    < this->vectorSizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::psd()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");

  this->extractScalarSeq(initialPos,
                         1, // spacing
                         this->subSequenceSize()-initialPos,
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
  bool bRC = ((initialPos         <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == psdVec.sizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::psdAtZero()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0,"");
  std::vector<double> psdResult(0,0.); // size will be determined by 'data.psd()'

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->subSequenceSize()-initialPos,
                           i,
                           data);
    data.psd(0,
             numBlocks,
             hopSizeRatio,
             psdResult);
    psdVec[i] = psdResult[0];
    //*m_env.subDisplayFile() << "psdResult[0] = " << psdResult[0] << std::endl;
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
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == gewVec.sizeLocal()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::geweke()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
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
uqSequenceOfVectorsClass<V,M>::meanStacc(
  unsigned int initialPos,
  V&           meanStaccVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize() ) &&
              (this->vectorSizeLocal() == meanStaccVec.sizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::meanStacc()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    meanStaccVec[i] = data.meanStacc(0);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subMinMax(
  unsigned int initialPos,
  V&           minVec,
  V&           maxVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == minVec.sizeLocal()     ) &&
              (this->vectorSizeLocal() == maxVec.sizeLocal()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subMinMax()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  uqScalarSequenceClass<double> data(m_env,0,"");

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.subMinMax(0,minVec[i],maxVec[i]);
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
  bool bRC = ((initialPos              <  this->subSequenceSize()  ) &&
              (this->vectorSizeLocal() == unifiedMinVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedMaxVec.sizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedMinMax()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  uqScalarSequenceClass<double> data(m_env,0,"");

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.unifiedMinMax(m_vectorSpace.numOfProcsForStorage() == 1,
                       0,
                       unifiedMinVec[i],
                       unifiedMaxVec[i]);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subHistogram(
  unsigned int     initialPos,
  const V&         minVec,
  const V&         maxVec,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& quanttsForAllBins) const
{
  bool bRC = ((initialPos               <  this->subSequenceSize() ) &&
              (this->vectorSizeLocal()  == minVec.sizeLocal()      ) &&
              (this->vectorSizeLocal()  == maxVec.sizeLocal()      ) &&
              (0                        <  centersForAllBins.size()) &&
              (centersForAllBins.size() == quanttsForAllBins.size()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subHistogram()",
                      "invalid input data");

  for (unsigned int j = 0; j < quanttsForAllBins.size(); ++j) {
    centersForAllBins[j] = new V(m_vectorSpace.zeroVector());
    quanttsForAllBins [j] = new V(m_vectorSpace.zeroVector());
  }

  unsigned int dataSize = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize,"");
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > centers(centersForAllBins.size(),0.);
    std::vector<unsigned int> quantts(quanttsForAllBins.size(), 0 );
    data.subHistogram(0,
                      minVec[i],
                      maxVec[i],
                      centers,
                      quantts);

    for (unsigned int j = 0; j < quantts.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(quanttsForAllBins[j]))[i] = (double) quantts[j];
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
  std::vector<V*>& unifiedQuanttsForAllBins) const
{
  bool bRC = ((initialPos                      <  this->subSequenceSize()        ) &&
              (this->vectorSizeLocal()         == unifiedMinVec.sizeLocal()      ) &&
              (this->vectorSizeLocal()         == unifiedMaxVec.sizeLocal()      ) &&
              (0                               <  unifiedCentersForAllBins.size()) &&
              (unifiedCentersForAllBins.size() == unifiedQuanttsForAllBins.size()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedHistogram()",
                      "invalid input data");

  for (unsigned int j = 0; j < unifiedQuanttsForAllBins.size(); ++j) {
    unifiedCentersForAllBins[j] = new V(m_vectorSpace.zeroVector());
    unifiedQuanttsForAllBins [j] = new V(m_vectorSpace.zeroVector());
  }

  unsigned int dataSize = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize,"");
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > unifiedCenters(unifiedCentersForAllBins.size(),0.);
    std::vector<unsigned int> unifiedQuantts(unifiedQuanttsForAllBins.size(), 0 );
    data.unifiedHistogram(m_vectorSpace.numOfProcsForStorage() == 1,
                          0,
                          unifiedMinVec[i],
                          unifiedMaxVec[i],
                          unifiedCenters,
                          unifiedQuantts);

    for (unsigned int j = 0; j < unifiedQuantts.size(); ++j) {
      (*(unifiedCentersForAllBins[j]))[i] = unifiedCenters[j];
      (*(unifiedQuanttsForAllBins[j]))[i] = (double) unifiedQuantts[j];
    }
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subCdfStacc(
  unsigned int           initialPos,
  const std::vector<V*>& evalPositionsVecs,
  std::vector<V*>&       cdfStaccVecs) const
{
  bool bRC = ((initialPos               <  this->subSequenceSize() ) &&
              (0                        <  evalPositionsVecs.size()) &&
              (evalPositionsVecs.size() == cdfStaccVecs.size()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subCdfStacc()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numEvals = evalPositionsVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    cdfStaccVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> evalPositions(numEvals,0.);
  std::vector<double> cdfStaccs    (numEvals,0.);

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      evalPositions[j] = (*evalPositionsVecs[j])[i];
    }

    data.subCdfStacc(0,
                     evalPositions,
                     cdfStaccs);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*cdfStaccVecs[j])[i] = cdfStaccs[j];
    }
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subInterQuantileRange(
  unsigned int initialPos,
  V&           iqrVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == iqrVec.sizeLocal()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subInterQuantileRange()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    iqrVec[i] = data.subInterQuantileRange(0);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedInterQuantileRange(
  unsigned int initialPos,
  V&           unifiedIqrVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()  ) &&
              (this->vectorSizeLocal() == unifiedIqrVec.sizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedInterQuantileRange()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedIqrVec[i] = data.unifiedInterQuantileRange(m_vectorSpace.numOfProcsForStorage() == 1,
                                                      0);
  }

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subScalesForKDE(
  unsigned int initialPos,
  const V&     iqrVec,
  V&           scaleVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == iqrVec.sizeLocal()     ) &&
              (this->vectorSizeLocal() == scaleVec.sizeLocal()   ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subScalesForKDE()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    scaleVec[i] = data.subScaleForKDE(0,
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
  bool bRC = ((initialPos              <  this->subSequenceSize()    ) &&
              (this->vectorSizeLocal() == unifiedIqrVec.sizeLocal()  ) &&
              (this->vectorSizeLocal() == unifiedScaleVec.sizeLocal()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedScalesForKDE()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedScaleVec[i] = data.unifiedScaleForKDE(m_vectorSpace.numOfProcsForStorage() == 1,
                                                 0,
                                                 unifiedIqrVec[i]);
  }

  return;
}
#if 0
template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::sabGaussianKDE(
  const V& evalParamVec,
        V& densityVec) const
{
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0, // Use the whole chain
                           1, // spacing
                           this->subSequenceSize(),
                           i,
                           data);

    densityVec[i] = data.sabGaussianKDE(evalParamVec[i]);
  }

  return;
}
#endif
template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subGaussianKDE(
  unsigned int           initialPos,
  const V&               scaleVec,
  const std::vector<V*>& evalParamVecs,
  std::vector<V*>&       densityVecs) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == scaleVec.sizeLocal()   ) &&
              (0                       <  evalParamVecs.size()   ) &&
              (evalParamVecs.size()    == densityVecs.size()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subGaussianKDE()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numEvals = evalParamVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    densityVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> evalParams(numEvals,0.);
  std::vector<double> densities  (numEvals,0.);

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      evalParams[j] = (*evalParamVecs[j])[i];
    }

    data.subGaussianKDE(0,
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
  bool bRC = ((initialPos                  <  this->subSequenceSize()    ) &&
              (this->vectorSizeLocal()     == unifiedScaleVec.sizeLocal()) &&
              (0                           <  unifiedEvalParamVecs.size()) &&
              (unifiedEvalParamVecs.size() == unifiedDensityVecs.size()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedGaussianKDE()",
                      "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0,"");

  unsigned int numEvals = unifiedEvalParamVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    unifiedDensityVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> unifiedEvalParams(numEvals,0.);
  std::vector<double> unifiedDensities (numEvals,0.);

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      unifiedEvalParams[j] = (*unifiedEvalParamVecs[j])[i];
    }

    data.unifiedGaussianKDE(m_vectorSpace.numOfProcsForStorage() == 1,
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
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::filter()"
                           << ": initialPos = "      << initialPos
                           << ", spacing = "         << spacing
                           << ", subSequenceSize = " << this->subSequenceSize()
                           << std::endl;
  }

  unsigned int i = 0;
  unsigned int j = initialPos;
  unsigned int originalSubSequenceSize = this->subSequenceSize();
  while (j < originalSubSequenceSize) {
    if (i != j) {
      //*m_env.subDisplayFile() << i << "--" << j << " ";
      delete m_seq[i];
      m_seq[i] = new V(*(m_seq[j]));
    }
    i++;
    j += spacing;
  }

  this->resetValues(i,originalSubSequenceSize-i);
  this->resizeSequence(i);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::filter()"
                           << ": initialPos = "      << initialPos
                           << ", spacing = "         << spacing
                           << ", subSequenceSize = " << this->subSequenceSize()
                           << std::endl;
  }

  return;
}

template <class V, class M>
double
uqSequenceOfVectorsClass<V,M>::estimateConvBrooksGelman(
  unsigned int initialPos,
  unsigned int numPos) const
{
  // This method requires *at least* two sequences. Error if there is only one.
  UQ_FATAL_RC_MACRO( (m_env.numSubEnvironments() < 2),
		     UQ_UNAVAILABLE_RANK,
		     "uqSequenceOfVectorsClass<V,M>::estimateConvBrooksGelman()",
		     "At least two sequences required for Brooks-Gelman convergence test.");

  // TODO: Need special case for 1-dimensional parameter space.

  // Initialize with garbage to give the user a clue something is funky.
  double convMeasure = -1.0;

  // We only do the work on the subenvironment where the sequence data
  // is stashed.
  if( m_env.inter0Rank() >= 0 )
    {
      // Sanity Checks
      
      
      // REMEMBER: \psi is a *vector* of parameters    
      // Get quantities we will use several times
      V psi_j_dot = m_vectorSpace.zeroVector();
      V psi_dot_dot = m_vectorSpace.zeroVector();
      V work = m_vectorSpace.zeroVector();

      // m = number of chains > 1
      // n = number of steps for which we are computing the metric
      int m = m_env.numSubEnvironments();
      int n = numPos;
	
      this->subMean( initialPos, numPos , psi_j_dot );
      this->unifiedMean( initialPos, numPos, psi_dot_dot );

      /* W = \frac{1}{m*(n-1)}*\sum_{j=1}^m \sum{t=1}^n 
	 (\psi_{jt} - \overline{\psi_{j\cdot}})*(\psi_{jt} - \overline{\psi_{j\cdot}})^T
	 This corresponds to the "within-sequence" covariance matrix. */
      M* W_local = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );
      M* W = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );
      V  psi_j_t = m_vectorSpace.zeroVector();

      // Sum within the chain
      for( unsigned int t = initialPos; t < initialPos+numPos; ++t )
	{
	  psi_j_t = *(m_seq[t]);

	  work = psi_j_t - psi_j_dot;

	  (*W_local) += matrixProduct( work, work );
	}

      // Now do the sum over the chains
      // W will be available on all inter0 processors
      W_local->mpiSum( m_env.inter0Comm().Comm(), (*W) );

      (*W) = 1.0/(double(m)*(double(n)-1.0)) * (*W);

      // Need to delete pointers to temporary covariance matrices
      delete W_local;

      /* B/n = \frac{1}{m-1}\sum_{j=1}^m
	 (\overline{\psi_{j\cdot}} - \overline{\psi_{\cdot \cdot}})*
	 (\overline{\psi_{j\cdot}} - \overline{\psi_{\cdot \cdot}})^T
	 This corresponds to the "between-sequence" covariance matrix. */
      M* B_over_n_local = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );
      M* B_over_n = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );

      work = psi_j_dot - psi_dot_dot;
      (*B_over_n_local) = matrixProduct( work, work );

      B_over_n_local->mpiSum( m_env.inter0Comm().Comm(), (*B_over_n) );

      // Need to delete pointers to temporary covariance matrices
      delete B_over_n_local;

      (*B_over_n) = 1.0/(double(m)-1.0) * (*B_over_n); 

      /* R_p = (n-1)/n + (m+1)/m * \lambda
	 \lambda = largest eigenvalue of W^{-1}*B/n */
      M* A = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );

      W->invertMultiply( *B_over_n, *A );

      // Need to delete pointers to temporary covariance matrices
      delete W;
      delete B_over_n;

      double eigenValue;
      V eigenVector = m_vectorSpace.zeroVector(); 

      A->largestEigen( eigenValue, eigenVector );

      // Need to delete pointers to temporary covariance matrices
      delete A;

      // Now, finally compute the final convMeasure
      convMeasure = (double(n)-1.0)/double(n) + (double(m)+1.0)/double(m)*eigenValue;

    } // End of check on inter0Rank

  //TODO: Do we need a Barrier here?
  //TODO: Error checking on MPI.
  m_env.fullComm().Barrier();

  return convMeasure;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subWriteContents(
  const std::string&            fileName,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  bool okSituation = (m_env.subRank() >= 0);
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subWriteContents()",
                      "unexpected subRank");

  std::ofstream* ofsVar = NULL;
  m_env.openOutputFile(fileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                       allowedSubEnvIds,
                       false, // A 'true' causes problems when the user chooses (via options
                              // in the input file) to use just one file for all outputs.
                       ofsVar);

  if (ofsVar) {
    *ofsVar << m_name << "_sub" << m_env.subIdString() << " = zeros(" << this->subSequenceSize()
            << ","                                                    << this->vectorSizeLocal()
            << ");"
            << std::endl;
    *ofsVar << m_name << "_sub" << m_env.subIdString() << " = [";
    unsigned int chainSize = this->subSequenceSize();
    for (unsigned int j = 0; j < chainSize; ++j) {
      bool savedVectorPrintState = m_seq[j]->getPrintHorizontally();
      m_seq[j]->setPrintHorizontally(true);
      *ofsVar << *(m_seq[j])
              << std::endl;
      m_seq[j]->setPrintHorizontally(savedVectorPrintState);
    }
    *ofsVar << "];\n";
  }

  if (ofsVar) {
    ofsVar->close();
  }

  return;
}
#if 0
template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::subWriteContents(std::ofstream& ofsvar) const
{
  bool okSituation = (m_env.subRank() >= 0);
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::subWriteContents()",
                      "unexpected subRank");

  ofsvar << m_name << "_sub" << m_env.subIdString() << " = zeros(" << this->subSequenceSize()
         << ","                                                    << this->vectorSizeLocal()
         << ");"
         << std::endl;
  ofsvar << m_name << "_sub" << m_env.subIdString() << " = [";
  unsigned int chainSize = this->subSequenceSize();
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
uqSequenceOfVectorsClass<V,M>::unifiedWriteContents(std::ofstream& ofsvar) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::unifiedWriteContents(1)",
                      "not implemented yet");
  return;
}
#endif

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedWriteContents(const std::string& fileName) const
{
  m_env.fullComm().Barrier();
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::unifiedWriteContents()"
                            << ": fullRank "       << m_env.fullRank()
                            << ", subEnvironment " << m_env.subId()
                            << ", subRank "        << m_env.subRank()
                            << ", inter0Rank "     << m_env.inter0Rank()
                            << ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                            << ", fileName = "     << fileName
                            << std::endl;
  }

  if (m_env.inter0Rank() >= 0) {
    for (unsigned int r = 0; r < (unsigned int) m_env.inter0Comm().NumProc(); ++r) {
      if (m_env.inter0Rank() == (int) r) {
        // My turn
        std::ofstream* unifiedOfsVar = NULL;
        // bool writeOver = (r == 0);
        bool writeOver = false; // A 'true' causes problems when the user chooses (via options
                                // in the input file) to use just one file for all outputs.
        m_env.openUnifiedOutputFile(fileName,
                                    UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                                    writeOver,
                                    unifiedOfsVar);

        if (r == 0) {
          *unifiedOfsVar << m_name << "_unified" << " = zeros(" << this->subSequenceSize()*m_env.inter0Comm().NumProc()
                         << ","                                 << this->vectorSizeLocal()
                         << ");"
                         << std::endl;
          *unifiedOfsVar << m_name << "_unified" << " = [";
        }

        unsigned int chainSize = this->subSequenceSize();
        for (unsigned int j = 0; j < chainSize; ++j) {
          bool savedVectorPrintState = m_seq[j]->getPrintHorizontally();
          m_seq[j]->setPrintHorizontally(true);
          *unifiedOfsVar << *(m_seq[j])
                         << std::endl;
          m_seq[j]->setPrintHorizontally(savedVectorPrintState);
        }

        unifiedOfsVar->close();
      }
      m_env.inter0Comm().Barrier();
    }

    if (m_env.inter0Rank() == 0) {
      std::ofstream* unifiedOfsVar = NULL;
      m_env.openUnifiedOutputFile(fileName,
                                  UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                                  false, // Yes, 'writeOver = false' in order to close the array for matlab
                                  unifiedOfsVar);
      *unifiedOfsVar << "];\n";
      unifiedOfsVar->close();
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::unifiedWriteContents()"
                            << ", fileName = " << fileName
                            << std::endl;
  }
  m_env.fullComm().Barrier();

  return;
}

template <class V, class M>
void
uqSequenceOfVectorsClass<V,M>::unifiedReadContents(
  const std::string& fileName,
  const unsigned int subReadSize)
{
  double unifiedReadSize = subReadSize*m_env.inter0Comm().NumProc();

  m_env.fullComm().Barrier();
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqSequenceOfVectorsClass<V,M>::unifiedReadContents()"
                            << ": fullRank "                       << m_env.fullRank()
                            << ", subEnvironment "                 << m_env.subId()
                            << ", subRank "                        << m_env.subRank()
                            << ", inter0Rank "                     << m_env.inter0Rank()
                            << ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                            << ", fileName = "                     << fileName
                            << ", subReadSize = "              << subReadSize
                            << ", unifiedReadSize = "          << unifiedReadSize
                            << std::endl;
  }

  this->resizeSequence(subReadSize);

  if (m_env.inter0Rank() >= 0) {
    // In the logic below, the id of a line' begins with value 0 (zero)
    unsigned int idOfMyFirstLine = 1 + m_env.inter0Rank()*subReadSize;
    unsigned int idOfMyLastLine = (1 + m_env.inter0Rank())*subReadSize;
    unsigned int numParams = this->vectorSizeLocal();

    for (unsigned int r = 0; r < (unsigned int) m_env.inter0Comm().NumProc(); ++r) {
      if (m_env.inter0Rank() == (int) r) {
        // My turn
        std::ifstream* ifsvar = new std::ifstream((fileName+"."+UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT).c_str(), std::ofstream::in);
        UQ_FATAL_TEST_MACRO((ifsvar == NULL) || (ifsvar->is_open() == false),
                            m_env.fullRank(),
                            "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                            "file with fileName could not be found");

        if (m_env.inter0Rank() == 0) {
          // Read number of chain positions in the file by taking care of the first line,
          // which resembles something like 'variable_name = zeros(n_positions,m_params);'
	  std::string tmpString;

          // Read 'variable name' string
          *ifsvar >> tmpString;
	  //std::cout << "Just read '" << tmpString << "'" << std::endl;

          // Read '=' sign
          *ifsvar >> tmpString;
	  //std::cout << "Just read '" << tmpString << "'" << std::endl;
          UQ_FATAL_TEST_MACRO(tmpString != "=",
                              m_env.fullRank(),
                              "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                              "string should be the '=' sign");

          // Read     'zeros(n_positions,n_params)' string
          // Position  0123456
          *ifsvar >> tmpString;
	  //std::cout << "Just read '" << tmpString << "'" << std::endl;
          unsigned int posInTmpString = 6;

          // Isolate 'n_positions' in a string
          char nPositionsString[tmpString.size()-posInTmpString+1];
          unsigned int posInPositionsString = 0;
          do {
            UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                                m_env.fullRank(),
                                "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                                "symbol ',' not found in first line of file");
            nPositionsString[posInPositionsString++] = tmpString[posInTmpString++];
          } while (tmpString[posInTmpString] != ',');
          nPositionsString[posInPositionsString] = '\0';

          // Isolate 'n_params' in a string
          posInTmpString++; // Avoid reading ',' char
          char nParamsString[tmpString.size()-posInTmpString+1];
          unsigned int posInParamsString = 0;
          do {
            UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                                m_env.fullRank(),
                                "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                                "symbol ')' not found in first line of file");
            nParamsString[posInParamsString++] = tmpString[posInTmpString++];
          } while (tmpString[posInTmpString] != ')');
          nParamsString[posInParamsString] = '\0';

          // Convert 'n_positions' and 'n_params' strings to numbers
          unsigned int sizeOfChainInFile = (unsigned int) strtod(nPositionsString,NULL);
          unsigned int numParamsInFile   = (unsigned int) strtod(nParamsString,   NULL);
          if (m_env.subDisplayFile()) {
            *m_env.subDisplayFile() << "In uqSequenceOfVectorsClass<V,M>::unifiedReadContents()"
                                   << ": fullRank "            << m_env.fullRank()
                                   << ", sizeOfChainInFile = " << sizeOfChainInFile
                                   << ", numParamsInFile = "   << numParamsInFile
                                   << std::endl;
          }

          // Check if [size of chain in file] >= [requested unified sequence size]
          UQ_FATAL_TEST_MACRO(sizeOfChainInFile < unifiedReadSize,
                              m_env.fullRank(),
                              "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                              "size of chain in file is not big enough");

          // Check if [num params in file] == [num params in current chain]
          UQ_FATAL_TEST_MACRO(numParamsInFile != numParams,
                              m_env.fullRank(),
                              "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                              "number of parameters of chain in file is different than number of parameters in this chain object");
        }

        // Code common to any core in 'inter0Comm', including core of rank 0
        unsigned int maxCharsPerLine = 64*numParams; // Up to about 60 characters to represent each parameter value

        unsigned int lineId = 0;
        while (lineId < idOfMyFirstLine) {
          ifsvar->ignore(maxCharsPerLine,'\n');
          lineId++;
        };

        if (m_env.inter0Rank() == 0) {
          // Take care of initial part of the first data line,
          // which resembles something like 'variable_name = [value1 value2 ...'
	  std::string tmpString;

          // Read 'variable name' string
          *ifsvar >> tmpString;
	  //std::cout << "Core 0 just read '" << tmpString << "'" << std::endl;

          // Read '=' sign
          *ifsvar >> tmpString;
	  //std::cout << "Core 0 just read '" << tmpString << "'" << std::endl;
          UQ_FATAL_TEST_MACRO(tmpString != "=",
                              m_env.fullRank(),
                              "uqSequenceOfVectorsClass<V,M>::unifiedReadContents()",
                              "in core 0, string should be the '=' sign");

          // Take into account the ' [' portion
	  std::streampos tmpPos = ifsvar->tellg();
          ifsvar->seekg(tmpPos+(std::streampos)2);
        }

        V tmpVec(m_vectorSpace.zeroVector());
        while (lineId <= idOfMyLastLine) {
          for (unsigned int i = 0; i < numParams; ++i) {
            *ifsvar >> tmpVec[i];
          }
          this->setPositionValues(lineId - idOfMyFirstLine, tmpVec);
          lineId++;
        };

        ifsvar->close();
      }
      m_env.inter0Comm().Barrier();
    }
  }
  else {
    V tmpVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 1; i < subReadSize; ++i) {
      this->setPositionValues(i,tmpVec);
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqSequenceOfVectorsClass<V,M>::unifiedReadContents()"
                            << ", fileName = " << fileName
                            << std::endl;
  }
  m_env.fullComm().Barrier();

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

