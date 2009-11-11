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

#ifndef __UQ_VECTOR_SEQUENCE_H__
#define __UQ_VECTOR_SEQUENCE_H__

#include <uqVectorSubset.h>
#include <uqScalarSequence.h>
#include <uqSequenceStatisticalOptions.h>
#include <uqArrayOfOneDGrids.h>
#include <uqArrayOfOneDTables.h>
#include <uq2dArrayOfStuff.h>
#include <sys/time.h>
#include <fstream>

template <class V, class M>
class uqBaseVectorSequenceClass
{
public:
           uqBaseVectorSequenceClass(const uqVectorSpaceClass<V,M>& vectorSpace,
                                     unsigned int                   subSequenceSize,
                                     const std::string&             name);
  virtual ~uqBaseVectorSequenceClass();

  virtual  unsigned int             subSequenceSize            () const = 0;
           unsigned int             unifiedSequenceSize        () const;
           unsigned int             vectorSizeLocal            () const;
           unsigned int             vectorSizeGlobal           () const;
  const    uqVectorSpaceClass<V,M>& vectorSpace                () const;
  const    std::string&             name                       () const;
           void                     setName                    (const std::string& newName);
           void                     clear                      ();
  const    V&                       subMinValues               () const;
  const    V&                       unifiedMinValues           () const;
  const    V&                       subMaxValues               () const;
  const    V&                       unifiedMaxValues           () const;
  const    V&                       subMeanValues              () const;
  const    V&                       unifiedMeanValues          () const;
  const    V&                       subSampleVarianceValues    () const;
  const    V&                       unifiedSampleVarianceValues() const;
  const    uqBoxSubsetClass<V,M>&   subValuesBox               () const;
  const    uqBoxSubsetClass<V,M>&   unifiedValuesBox           () const;
           void                     deleteStoredVectors        ();

  virtual  void                     resizeSequence             (unsigned int newSubSequenceSize) = 0;
  virtual  void                     resetValues                (unsigned int initialPos, unsigned int numPos) = 0;
  virtual  void                     erasePositions             (unsigned int initialPos, unsigned int numPos) = 0;
  virtual  void                     getPositionValues          (unsigned int posId,       V& vec) const = 0;
  virtual  void                     setPositionValues          (unsigned int posId, const V& vec) = 0;
  virtual  void                     setGaussian                (const gsl_rng* rng, const V& meanVec, const V& stdDevVec) = 0;
  virtual  void                     setUniform                 (const gsl_rng* rng, const V& aVec,    const V& bVec     ) = 0;
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  virtual  void                     subUniformlySampledMdf     (const V&                       numEvaluationPointsVec,
                                                                uqArrayOfOneDGridsClass <V,M>& mdfGrids,
                                                                uqArrayOfOneDTablesClass<V,M>& mdfValues) const = 0;
#endif
  virtual  void                     subUniformlySampledCdf     (const V&                       numEvaluationPointsVec,
                                                                uqArrayOfOneDGridsClass <V,M>& cdfGrids,
                                                                uqArrayOfOneDTablesClass<V,M>& cdfValues) const = 0;
  virtual  void                     unifiedUniformlySampledCdf (const V&                       numEvaluationPointsVec,
                                                                uqArrayOfOneDGridsClass <V,M>& unifiedCdfGrids,
                                                                uqArrayOfOneDTablesClass<V,M>& unifieddfValues) const = 0;

  virtual  void                     subMean                    (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                V&                                       meanVec) const = 0;
  virtual  void                     unifiedMean                (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                V&                                       unifiedMeanVec) const = 0;
  virtual  void                     subSampleVariance          (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                const V&                                 meanVec,
                                                                V&                                       samVec) const = 0;
  virtual  void                     unifiedSampleVariance      (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                const V&                                 unifiedMeanVec,
                                                                V&                                       unifiedSamVec) const = 0;
  virtual  void                     subPopulationVariance      (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                const V&                                 meanVec,
                                                                V&                                       popVec) const = 0;
  virtual  void                     unifiedPopulationVariance  (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                const V&                                 unifiedMeanVec,
                                                                V&                                       unifiedPopVec) const = 0;

  virtual  void                     autoCovariance             (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                const V&                                 meanVec,
                                                                unsigned int                             lag,
                                                                V&                                       covVec) const = 0;
  virtual  void                     autoCorrViaDef             (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                unsigned int                             lag,
                                                                V&                                       corrVec) const = 0;
  virtual  void                     autoCorrViaFft             (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                const std::vector<unsigned int>&         lags,
                                                                std::vector<V*>&                         corrVecs) const = 0;
  virtual  void                     autoCorrViaFft             (unsigned int                             initialPos,
                                                                unsigned int                             numPos,
                                                                unsigned int                             numSum,
                                                                V&                                       autoCorrsSumVec) const = 0;
  virtual  void                     bmm                        (unsigned int                             initialPos,
                                                                unsigned int                             batchLength,
                                                                V&                                       bmmVec) const = 0;
  virtual  void                     fftForward                 (unsigned int                             initialPos,
                                                                unsigned int                             fftSize,
                                                                unsigned int                             paramId,
                                                                std::vector<std::complex<double> >&      fftResult) const = 0;
//virtual  void                     fftInverse                 (unsigned int                             fftSize) = 0;
  virtual  void                     psd                        (unsigned int                             initialPos,
                                                                unsigned int                             numBlocks,
                                                                double                                   hopSizeRatio,
                                                                unsigned int                             paramId,
                                                                std::vector<double>&                     psdResult) const = 0;
  virtual  void                     psdAtZero                  (unsigned int                             initialPos,
                                                                unsigned int                             numBlocks,
                                                                double                                   hopSizeRatio,
                                                                V&                                       psdVec) const = 0;
  virtual  void                     geweke                     (unsigned int                             initialPos,
                                                                double                                   ratioNa,
                                                                double                                   ratioNb,
                                                                V&                                       gewVec) const = 0;
  virtual  void                     meanStacc                  (unsigned int                             initialPos,
                                                                V&                                       meanStaccVec) const = 0;
  virtual  void                     subMinMax                  (unsigned int                             initialPos,
                                                                V&                                       minVec,
                                                                V&                                       maxVec) const = 0;
  virtual  void                     unifiedMinMax              (unsigned int                             initialPos,
                                                                V&                                       unifiedMinVec,
                                                                V&                                       unifiedMaxVec) const = 0;
  virtual  void                     subHistogram               (unsigned int                             initialPos,
                                                                const V&                                 minVec,
                                                                const V&                                 maxVec,
                                                                std::vector<V*>&                         centersForAllBins,
                                                                std::vector<V*>&                         quanttsForAllBins) const = 0;
  virtual  void                     unifiedHistogram           (unsigned int                             initialPos,
                                                                const V&                                 unifiedMinVec,
                                                                const V&                                 unifiedMaxVec,
                                                                std::vector<V*>&                         unifiedCentersForAllBins,
                                                                std::vector<V*>&                         unifiedQuanttsForAllBins) const = 0;
  virtual  void                     subCdfStacc                (unsigned int                             initialPos,
                                                                const std::vector<V*>&                   evalPositionsVecs,
                                                                std::vector<V*>&                         cdfStaccVecs) const = 0;
  virtual  void                     subInterQuantileRange      (unsigned int                             initialPos,
                                                                V&                                       iqrVec) const = 0;
  virtual  void                     unifiedInterQuantileRange  (unsigned int                             initialPos,
                                                                V&                                       unifiedIqrVec) const = 0;
  virtual  void                     subScalesForKDE            (unsigned int                             initialPos,
                                                                const V&                                 iqrVec,
                                                                unsigned int                             kdeDimension,
                                                                V&                                       scaleVec) const = 0;
  virtual  void                     unifiedScalesForKDE        (unsigned int                             initialPos,
                                                                const V&                                 unifiedIqrVec,
                                                                unsigned int                             kdeDimension,
                                                                V&                                       unifiedScaleVec) const = 0;
//virtual  void                     sabGaussianKDE             (const V&                                 evaluationParamVec,
//                                                              V&                                       densityVec) const = 0;
  virtual  void                     subGaussianKDE             (unsigned int                             initialPos,
                                                                const V&                                 scaleVec,
                                                                const std::vector<V*>&                   evaluationParamVecs,
                                                                std::vector<V*>&                         densityVecs) const = 0;
  virtual  void                     unifiedGaussianKDE         (unsigned int                             initialPos,
                                                                const V&                                 unifiedScaleVec,
                                                                const std::vector<V*>&                   unifiedEvaluationParamVecs,
                                                                std::vector<V*>&                         unifiedDensityVecs) const = 0;
//virtual  void                     subWriteContents           (std::ofstream&                           ofsvar) const = 0;
//virtual  void                     unifiedWriteContents       (std::ofstream&                           ofsvar) const = 0;
  virtual  void                     subWriteContents           (const std::string&                       fileName,
                                                                const std::set<unsigned int>&            allowedSubEnvIds) const = 0;
  virtual  void                     unifiedWriteContents       (const std::string&                       fileName) const = 0;
  virtual  void                     unifiedReadContents        (const std::string&                       fileName,
                                                                const unsigned int                       subSequenceSize) = 0;
  virtual  void                     select                     (const std::vector<unsigned int>&         idsOfUniquePositions) = 0;
  virtual  void                     filter                     (unsigned int                             initialPos,
                                                                unsigned int                             spacing) = 0;

           void                     computeStatistics          (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                std::ofstream*                           passedOfs);

           void                     computeFilterParams        (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                std::ofstream*                           passedOfs,
                                                                unsigned int   &                         initialPos,
                                                                unsigned int   &                         spacing);

  virtual  double                   estimateConvBrooksGelman   (unsigned int                             initialPos,
                                                                unsigned int                             numPos) const = 0;

protected:
           void                     copy                       (const uqBaseVectorSequenceClass<V,M>&    src);
           void                     computeMeanVars            (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                std::ofstream*                           passedOfs,
                                                                V*                                       subMeanPtr,
                                                                V*                                       subSampleVarPtr,
                                                                V*                                       subPopulVarPtr);
           void                     computeBMM                 (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                std::ofstream*                           passedOfs);
           void                     computeFFT                 (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                std::ofstream*                           passedOfs);
           void                     computePSD                 (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                std::ofstream*                           passedOfs);
           void                     computePSDAtZero           (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                std::ofstream*                           passedOfs);
           void                     computeGeweke              (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                std::ofstream*                           passedOfs);
           void                     computeMeanStacc           (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                std::ofstream*                           passedOfs);
           void                     computeAutoCorrViaDef      (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                const std::vector<unsigned int>&         lagsForCorrs,
                                                                std::ofstream*                           passedOfs);
           void                     computeAutoCorrViaFFT      (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                const std::vector<unsigned int>&         initialPosForStatistics,
                                                                const std::vector<unsigned int>&         lagsForCorrs,
                                                                std::ofstream*                           passedOfs);
           void                     computeHistCdfstaccKde     (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                std::ofstream*                           passedOfs);
           void                     computeCovCorrMatrices     (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                std::ofstream*                           passedOfs);

  virtual  void                     extractScalarSeq           (unsigned int                             initialPos,
                                                                unsigned int                             spacing,
                                                                unsigned int                             numPos,
                                                                unsigned int                             paramId,
                                                                uqScalarSequenceClass<double>&           scalarSeq) const = 0;
  virtual  void                     extractRawData             (unsigned int                             initialPos,
                                                                unsigned int                             spacing,
                                                                unsigned int                             numPos,
                                                                unsigned int                             paramId,
                                                                std::vector<double>&                     rawData) const = 0;

  const uqBaseEnvironmentClass&  m_env;
  const uqVectorSpaceClass<V,M>& m_vectorSpace;
  std::string                    m_name;

  mutable uqFftClass<double>*    m_fftObj;
  mutable V*                     m_subMinValues;
  mutable V*                     m_unifiedMinValues;
  mutable V*                     m_subMaxValues;
  mutable V*                     m_unifiedMaxValues;
  mutable V*                     m_subMeanValues;
  mutable V*                     m_unifiedMeanValues;
  mutable V*                     m_subSampleVarianceValues;
  mutable V*                     m_unifiedSampleVarianceValues;
  mutable uqBoxSubsetClass<V,M>* m_subValuesBox;
  mutable uqBoxSubsetClass<V,M>* m_unifiedValuesBox;
};

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::uqBaseVectorSequenceClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   subSequenceSize,
  const std::string&             name)
  :
  m_env                        (vectorSpace.env()),
  m_vectorSpace                (vectorSpace),
  m_name                       (name),
  m_fftObj                     (new uqFftClass<double>(m_env)),
  m_subMinValues               (NULL),
  m_unifiedMinValues           (NULL),
  m_subMaxValues               (NULL),
  m_unifiedMaxValues           (NULL),
  m_subMeanValues              (NULL),
  m_unifiedMeanValues          (NULL),
  m_subSampleVarianceValues    (NULL),
  m_unifiedSampleVarianceValues(NULL),
  m_subValuesBox               (NULL),
  m_unifiedValuesBox           (NULL)
{
}

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::~uqBaseVectorSequenceClass()
{
  //clear();
  deleteStoredVectors();
  if (m_fftObj != NULL) delete m_fftObj;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::copy(const uqBaseVectorSequenceClass<V,M>& src)
{
  // FIX ME: should check environments as well ???

  UQ_FATAL_TEST_MACRO(m_vectorSpace.dimLocal() != src.m_vectorSpace.dimLocal(),
                      m_env.fullRank(),
                      "uqSequenceOfVectorsClass<V,M>::copy()",
                      "incompatible vector space dimensions");

  m_name = src.m_name;
  deleteStoredVectors();

  return;
}

template <class V, class M>
unsigned int
uqBaseVectorSequenceClass<V,M>::unifiedSequenceSize() const
{
  unsigned int unifiedNumSamples = 0;

  bool useOnlyInter0Comm = (m_vectorSpace.numOfProcsForStorage() == 1);

  if (useOnlyInter0Comm) {
    if (m_env.inter0Rank() >= 0) {
      unsigned int subNumSamples = this->subSequenceSize();
      int mpiRC = MPI_Allreduce((void *) &subNumSamples, (void *) &unifiedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.fullRank(),
                          "uqBaseVectorSequenceClass<V,M>::unifiedSequenceSize()",
                          "failed MPI_Allreduce() for unifiedSequenceSize()");
    }
    else {
      // Node not in the 'inter0' communicator
      unifiedNumSamples = this->subSequenceSize();
    }
  }
  else {
    UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                        m_env.fullRank(),
                        "uqBaseVectorSequenceClass<V,M>::unifiedSequenceSize()",
                        "parallel vectors not supported yet");
  }

  return unifiedNumSamples;
}

template <class V, class M>
unsigned int
uqBaseVectorSequenceClass<V,M>::vectorSizeLocal() const
{
  return m_vectorSpace.dimLocal();
}

template <class V, class M>
unsigned int
uqBaseVectorSequenceClass<V,M>::vectorSizeGlobal() const
{
  return m_vectorSpace.dimGlobal();
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorSequenceClass<V,M>::vectorSpace() const
{
  return m_vectorSpace;
}

template <class V, class M>
const std::string&
uqBaseVectorSequenceClass<V,M>::name() const
{
  return m_name;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::setName(const std::string& newName)
{
  m_name = newName;
  return;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::clear()
{
  unsigned int numPos = this->subSequenceSize();
  if (numPos) {
    this->resetValues(0,numPos);
    this->resizeSequence(0);
  }

  return;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subMinValues() const
{
  if (m_subMinValues == NULL) {
    m_subMinValues = m_vectorSpace.newVector();
    if (m_subMaxValues == NULL) m_subMaxValues = m_vectorSpace.newVector();
    subMinMax(0,*m_subMinValues,*m_subMaxValues);
  }

  return *m_subMinValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMinValues() const
{
  if (m_unifiedMinValues == NULL) {
    m_unifiedMinValues = m_vectorSpace.newVector();
    if (m_unifiedMaxValues == NULL) m_unifiedMaxValues = m_vectorSpace.newVector();
    unifiedMinMax(0,*m_unifiedMinValues,*m_unifiedMaxValues);
  }

  return *m_unifiedMinValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subMaxValues() const
{
  if (m_subMaxValues == NULL) {
    if (m_subMinValues == NULL) m_subMinValues = m_vectorSpace.newVector();
    m_subMaxValues = m_vectorSpace.newVector();
    subMinMax(0,*m_subMinValues,*m_subMaxValues);
  }

  return *m_subMaxValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMaxValues() const
{
  if (m_unifiedMaxValues == NULL) {
    if (m_unifiedMinValues == NULL) m_unifiedMinValues = m_vectorSpace.newVector();
    m_unifiedMaxValues = m_vectorSpace.newVector();
    unifiedMinMax(0,*m_unifiedMinValues,*m_unifiedMaxValues);
  }

  return *m_unifiedMaxValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subMeanValues() const
{
  if (m_subMeanValues == NULL) {
    m_subMeanValues = m_vectorSpace.newVector();
    subMean(0,subSequenceSize(),*m_subMeanValues);
  }

  return *m_subMeanValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMeanValues() const
{
  if (m_unifiedMeanValues == NULL) {
    m_unifiedMeanValues = m_vectorSpace.newVector();
    unifiedMean(0,subSequenceSize(),*m_unifiedMeanValues);
  }

  return *m_unifiedMeanValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subSampleVarianceValues() const
{
  if (m_subSampleVarianceValues == NULL) {
    m_subSampleVarianceValues = m_vectorSpace.newVector();
    subSampleVariance(0,subSequenceSize(),subMeanValues(),*m_subSampleVarianceValues);
  }

  return *m_subSampleVarianceValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedSampleVarianceValues() const
{
  if (m_unifiedSampleVarianceValues == NULL) {
    m_unifiedSampleVarianceValues = m_vectorSpace.newVector();
    unifiedSampleVariance(0,subSequenceSize(),unifiedMeanValues(),*m_unifiedSampleVarianceValues);
  }

  return *m_unifiedSampleVarianceValues;
}

template <class V, class M>
const uqBoxSubsetClass<V,M>&
uqBaseVectorSequenceClass<V,M>::subValuesBox() const
{
  if (m_subValuesBox == NULL) {
    m_subValuesBox = new uqBoxSubsetClass<V,M>(m_name.c_str(),
                                               m_vectorSpace,
                                               this->subMinValues(),
                                               this->subMaxValues());
  }

  return *m_subValuesBox;
}

template <class V, class M>
const uqBoxSubsetClass<V,M>&
uqBaseVectorSequenceClass<V,M>::unifiedValuesBox() const
{
  if (m_unifiedValuesBox == NULL) {
    m_unifiedValuesBox = new uqBoxSubsetClass<V,M>(m_name.c_str(),
                                                   m_vectorSpace,
                                                   this->unifiedMinValues(),
                                                   this->unifiedMaxValues());
  }

  return *m_unifiedValuesBox;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::deleteStoredVectors()
{
  delete m_subMinValues;
  delete m_unifiedMinValues;
  delete m_subMaxValues;
  delete m_unifiedMaxValues;
  delete m_subMeanValues;
  delete m_unifiedMeanValues;
  delete m_subSampleVarianceValues;
  delete m_unifiedSampleVarianceValues;
  delete m_subValuesBox;
  delete m_unifiedValuesBox;

  m_subMinValues                = NULL;
  m_unifiedMinValues            = NULL;
  m_subMaxValues                = NULL;
  m_unifiedMaxValues            = NULL;
  m_subMeanValues               = NULL;
  m_unifiedMeanValues           = NULL;
  m_subSampleVarianceValues     = NULL;
  m_unifiedSampleVarianceValues = NULL;
  m_subValuesBox                = NULL;
  m_unifiedValuesBox            = NULL;

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeStatistics(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                           passedOfs)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing statistics for chain " << m_name << " ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  bool okSituation = ((passedOfs == NULL                            ) ||
                      ((passedOfs != NULL) && (m_env.subRank() >= 0)));
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.fullRank(),
                      "uqBaseVectorSequenceClass<V,M>::computeStatistics()",
                      "unexpected combination of file pointer and subRank");

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  // Set initial positions for the computation of chain statistics
  std::vector<unsigned int> initialPosForStatistics(statisticalOptions.initialDiscardedPortions().size(),0);
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    initialPosForStatistics[i] = (unsigned int) (statisticalOptions.initialDiscardedPortions()[i] * (double) (this->subSequenceSize()-1));
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeStatistics()"
                              << ": statisticalOptions.initialDiscardedPortions()[" << i << "] = " << statisticalOptions.initialDiscardedPortions()[i]
                              << ", initialPosForStatistics[" << i << "] = " << initialPosForStatistics[i]
                              << std::endl;
    }
  }
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeStatistics(): initial positions for statistics =";
    for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
      *m_env.subDisplayFile() << " " << initialPosForStatistics[i];
    }
    *m_env.subDisplayFile() << std::endl;
  }

  //****************************************************
  // Compute mean, sample std, population std
  //****************************************************
  this->computeMeanVars(statisticalOptions,
                        passedOfs,
                        NULL,
                        NULL,
                        NULL);

  //****************************************************
  // Compute variance of sample mean through the 'batch means method' (BMM)
  //****************************************************
  if ((statisticalOptions.bmmRun()               ) &&
      (initialPosForStatistics.size()         > 0) &&
      (statisticalOptions.bmmLengths().size() > 0)) { 
    this->computeBMM(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }

  //****************************************************
  // Compute FFT of chain, for one parameter only
  //****************************************************
  if ((statisticalOptions.fftCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    this->computeFFT(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain, for one parameter only
  //****************************************************
  if ((statisticalOptions.psdCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    this->computePSD(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain at zero frequency
  //****************************************************
  if ((statisticalOptions.psdAtZeroCompute()             ) &&
      (initialPosForStatistics.size()                 > 0) &&
      (statisticalOptions.psdAtZeroNumBlocks().size() > 0)) { 
    this->computePSDAtZero(statisticalOptions,
                           initialPosForStatistics,
                           passedOfs);
  }

  //****************************************************
  // Compute Geweke
  //****************************************************
  if ((statisticalOptions.gewekeCompute()) &&
      (initialPosForStatistics.size() > 0)) {
    this->computeGeweke(statisticalOptions,
                        initialPosForStatistics,
                        passedOfs);
  }

  //****************************************************
  // Compute mean statistical accuracy
  //****************************************************
  if ((statisticalOptions.meanStaccCompute()) &&
      (initialPosForStatistics.size() > 0   )) {
    this->computeMeanStacc(statisticalOptions,
                           initialPosForStatistics,
                           passedOfs);
  }

  // Set lags for the computation of chain autocorrelations
  std::vector<unsigned int> lagsForCorrs(statisticalOptions.autoCorrNumLags(),1);
  for (unsigned int i = 1; i < lagsForCorrs.size(); ++i) {
    lagsForCorrs[i] = statisticalOptions.autoCorrSecondLag() + (i-1)*statisticalOptions.autoCorrLagSpacing();
  }

  //****************************************************
  // Compute autocorrelation coefficients via definition
  //****************************************************
  if ((statisticalOptions.autoCorrComputeViaDef()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    this->computeAutoCorrViaDef(statisticalOptions,
                                initialPosForStatistics,
                                lagsForCorrs,
                                passedOfs);
  }

  //****************************************************
  // Compute autocorrelation coefficients via FFT
  //****************************************************
  if ((statisticalOptions.autoCorrComputeViaFft()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    this->computeAutoCorrViaFFT(statisticalOptions,
                                initialPosForStatistics,
                                lagsForCorrs,
                                passedOfs);
  }

  //****************************************************
  // Compute histogram and/or cdf stacc and/or Kde
  //****************************************************
  if ((statisticalOptions.histCompute    ()) ||
      (statisticalOptions.cdfStaccCompute()) ||
      (statisticalOptions.kdeCompute     ())) {
    this->computeHistCdfstaccKde(statisticalOptions,
                                 passedOfs);
  }

  //****************************************************
  // Compute covariance and correlation matrices
  //****************************************************
  if ((statisticalOptions.covMatrixCompute ()) ||
      (statisticalOptions.corrMatrixCompute())) {
    this->computeCovCorrMatrices(statisticalOptions,
                                 passedOfs);
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "All statistics took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing statistics for chain " << m_name
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeMeanVars(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                           passedOfs,
  V*                                       subMeanPtr,
  V*                                       subSampleVarPtr,
  V*                                       subPopulVarPtr)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing mean, sample variance and population variance"
                            << std::endl;
  }

  V subChainMean(m_vectorSpace.zeroVector());
  this->subMean(0,
                this->subSequenceSize(),
                subChainMean);

  V subChainSampleVariance(m_vectorSpace.zeroVector());
  this->subSampleVariance(0,
                          this->subSequenceSize(),
                          subChainMean,
                         subChainSampleVariance);

  if ((m_env.displayVerbosity() >= 5) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeMeanVars()"
                            << ": subChainMean.sizeLocal() = "           << subChainMean.sizeLocal()
                            << ", subChainMean = "                  << subChainMean
                            << ", subChainSampleVariance.sizeLocal() = " << subChainSampleVariance.sizeLocal()
                            << ", subChainSampleVariance = "        << subChainSampleVariance
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nEstimated variance of sample mean for the whole chain " << m_name
                            << ", under independence assumption:"
                            << std::endl;
  }
  V estimatedVarianceOfSampleMean(subChainSampleVariance);
  estimatedVarianceOfSampleMean /= (double) this->subSequenceSize();
  bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
  estimatedVarianceOfSampleMean.setPrintHorizontally(false);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << estimatedVarianceOfSampleMean
                            << std::endl;
  }
  estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);

  V subChainPopulationVariance(m_vectorSpace.zeroVector());
  this->subPopulationVariance(0,
                              this->subSequenceSize(),
                              subChainMean,
                              subChainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Sub Mean and variances took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nSub mean, sample std, population std"
                            << std::endl;
    char line[512];
    sprintf(line,"%s%4s%s%9s%s%9s%s",
	    "Parameter",
            " ",
            "Mean",
            " ",
            "SampleStd",
            " ",
            "Popul.Std");
    *m_env.subDisplayFile() << line;

    for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e",
              m_vectorSpace.localComponentName(i).c_str(), /*.*/
              " ",
	      subChainMean[i],
              " ",
              std::sqrt(subChainSampleVariance[i]),
              " ",
              std::sqrt(subChainPopulationVariance[i]));
      *m_env.subDisplayFile() << line;
    }
    *m_env.subDisplayFile() << std::endl;
  }

  if (subMeanPtr     ) *subMeanPtr      = subChainMean;
  if (subSampleVarPtr) *subSampleVarPtr = subChainSampleVariance;
  if (subPopulVarPtr ) *subPopulVarPtr  = subChainPopulationVariance;

  if (m_env.numSubEnvironments() > 1) {
    // Write unified min-max
    if (m_env.subDisplayFile()) {
      if (m_vectorSpace.numOfProcsForStorage() == 1) {
        V unifiedChainMean(m_vectorSpace.zeroVector());
        this->unifiedMean(0,
                          this->subSequenceSize(),
                          unifiedChainMean);

        V unifiedChainSampleVariance(m_vectorSpace.zeroVector());
        this->unifiedSampleVariance(0,
                                    this->subSequenceSize(),
                                    unifiedChainMean,
                                    unifiedChainSampleVariance);

        V unifiedChainPopulationVariance(m_vectorSpace.zeroVector());
        this->unifiedPopulationVariance(0,
                                        this->subSequenceSize(),
                                        unifiedChainMean,
                                        unifiedChainPopulationVariance);

        if (m_env.inter0Rank() == 0) {
          *m_env.subDisplayFile() << "\nUnif mean, sample std, population std"
                                  << std::endl;
          char line[512];
          sprintf(line,"%s%4s%s%9s%s%9s%s",
	          "Parameter",
                  " ",
                  "Mean",
                  " ",
                  "SampleStd",
                  " ",
                  "Popul.Std");
          *m_env.subDisplayFile() << line;

          for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
            sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e",
                    m_vectorSpace.localComponentName(i).c_str(), /*.*/
                    " ",
	            unifiedChainMean[i],
                    " ",
                    std::sqrt(unifiedChainSampleVariance[i]),
                    " ",
                    std::sqrt(unifiedChainPopulationVariance[i]));
            *m_env.subDisplayFile() << line;
          }
          *m_env.subDisplayFile() << std::endl;
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.fullRank(),
                            "uqBaseVectorSequenceClass<V,M>::computeMeanVars()",
                            "unified min-max writing, parallel vectors not supported yet");
      }
    } // if subDisplayFile
  } // if numSubEnvs > 1

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeBMM(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing variance of sample mean through BMM"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeBMM(): lengths for batchs in BMM =";
    for (unsigned int i = 0; i < statisticalOptions.bmmLengths().size(); ++i) {
      *m_env.subDisplayFile() << " " << statisticalOptions.bmmLengths()[i];
    }
    *m_env.subDisplayFile() << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfBMM(initialPosForStatistics.size(),statisticalOptions.bmmLengths().size());
  for (unsigned int i = 0; i < _2dArrayOfBMM.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfBMM.numCols(); ++j) {
      _2dArrayOfBMM.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  V bmmVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
      unsigned int batchLength = statisticalOptions.bmmLengths()[batchLengthId];
      this->bmm(initialPos,
                batchLength,
                bmmVec);
      _2dArrayOfBMM(initialPosId,batchLengthId) = bmmVec;
    }
  }

  if (m_env.subDisplayFile()) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      *m_env.subDisplayFile() << "\nEstimated variance of sample mean, through batch means method, for subchain beggining at position " << initialPosForStatistics[initialPosId]
                              << " (each column corresponds to a batch length)"
                              << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subDisplayFile() << line;
      for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.bmmLengths()[batchLengthId]);
        *m_env.subDisplayFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.localComponentName(i).c_str() /*.*/);
        *m_env.subDisplayFile() << line;
        for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfBMM(initialPosId,batchLengthId)[i]);
          *m_env.subDisplayFile() << line;
        }
      }
      *m_env.subDisplayFile() << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain BMM took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeFFT(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing FFT of chain on parameter of id = " << statisticalOptions.fftParamId()
                            << std::endl;
  }

  std::vector<std::complex<double> > forwardResult(0,std::complex<double>(0.,0.));
  std::vector<std::complex<double> > inverseResult(0,std::complex<double>(0.,0.));
  uqFftClass<std::complex<double> > fftObj(m_env);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->fftForward(initialPosition,
                     statisticalOptions.fftSize(),
                     statisticalOptions.fftParamId(),
                     forwardResult);

    if (statisticalOptions.fftWrite() && passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << "_fftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << 1
             << ","                                                                                                              << forwardResult.size()
             << ");"
             << std::endl;
      for (unsigned int j = 0; j < forwardResult.size(); ++j) {
        ofsvar << m_name << "_fftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << 1
               << ","                                                                                               << j+1
               << ") = "                                                                                            << forwardResult[j].real()
               << " + i*"                                                                                           << forwardResult[j].imag()
               << ";"
               << std::endl;
      }
    } // if write

    if (statisticalOptions.fftTestInversion()) {
      fftObj.inverse(forwardResult,
                     statisticalOptions.fftSize(),
                     inverseResult);
      if (statisticalOptions.fftWrite() && passedOfs) {
        std::ofstream& ofsvar = *passedOfs;
        ofsvar << m_name << "_iftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << 1
               << ","                                                                                                              << inverseResult.size()
               << ");"
               << std::endl;
        for (unsigned int j = 0; j < inverseResult.size(); ++j) {
          ofsvar << m_name << "_iftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << 1
                 << ","                                                                                               << j+1
                 << ") = "                                                                                            << inverseResult[j].real()
                 << " + i*"                                                                                           << inverseResult[j].imag()
                 << ";"
                 << std::endl;
        }
      } // if write
    }
  } // for initialPosId

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain FFT took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computePSD(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing PSD of chain on parameter of id = " << statisticalOptions.psdParamId()
                            << std::endl;
  }

  std::vector<double> psdResult(0,0.);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->psd(initialPosition,
              statisticalOptions.psdNumBlocks(),
              statisticalOptions.psdHopSizeRatio(),
              statisticalOptions.psdParamId(),
              psdResult);

    if (statisticalOptions.psdWrite() && passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << "_psdInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << 1
             << ","                                                                                                              << psdResult.size()
             << ");"
             << std::endl;
      for (unsigned int j = 0; j < psdResult.size(); ++j) {
        ofsvar << m_name << "_psdInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << 1
               << ","                                                                                                      << j+1
               << ") = "                                                                                                   << psdResult[j]
               << ";"
               << std::endl;
      }
    } // if write
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain PSD took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computePSDAtZero(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing PSD at frequency zero for all parameters"
                            << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfPSDAtZero(initialPosForStatistics.size(),statisticalOptions.psdAtZeroNumBlocks().size());
  for (unsigned int i = 0; i < _2dArrayOfPSDAtZero.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfPSDAtZero.numCols(); ++j) {
      _2dArrayOfPSDAtZero.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  V psdVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      unsigned int numBlocks = statisticalOptions.psdAtZeroNumBlocks()[numBlocksId];
      this->psdAtZero(initialPosition,
                      numBlocks,
                      statisticalOptions.psdAtZeroHopSizeRatio(),
                      psdVec);
      _2dArrayOfPSDAtZero(initialPosId,numBlocksId) = psdVec;
    }
  }

  // Display PSD at frequency zero
  if ((statisticalOptions.psdAtZeroDisplay()) && (m_env.subDisplayFile())) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      *m_env.subDisplayFile() << "\nComputed PSD at frequency zero for subchain beggining at position " << initialPos
                              << ", so effective data size = " << this->subSequenceSize() - initialPos
                              << " (each column corresponds to a number of blocks)"
                              << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subDisplayFile() << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        *m_env.subDisplayFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.localComponentName(i).c_str() /*.*/);
        *m_env.subDisplayFile() << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]);
          *m_env.subDisplayFile() << line;
        }
      }
      *m_env.subDisplayFile() << std::endl;
    }
  }

  // Display estimated variance of sample mean through PSD
  if (/*(statisticalOptions.psdAtZeroDisplay()) &&*/ (m_env.subDisplayFile())) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      *m_env.subDisplayFile() << "\nEstimated variance of sample mean, through psd, for subchain beggining at position " << initialPos
                              << ", so effective data size = " << this->subSequenceSize() - initialPos
                              << " (each column corresponds to a number of blocks)"
                              << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subDisplayFile() << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        *m_env.subDisplayFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.localComponentName(i).c_str() /*.*/);
        *m_env.subDisplayFile() << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  2.*M_PI*_2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]/(double) (this->subSequenceSize() - initialPos));
          *m_env.subDisplayFile() << line;
        }
      }
      *m_env.subDisplayFile() << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain PSD at frequency zero took " << tmpRunTime 
                            << " seconds"
                            << std::endl;
  }

  // Write PSD at frequency zero
  if (statisticalOptions.psdAtZeroWrite() && passedOfs) {
    std::ofstream& ofsvar = *passedOfs;
    ofsvar << m_name << "_psdAtZeroNumBlocks_sub" << m_env.subIdString() << " = zeros(" << 1
           << ","                                                                       << statisticalOptions.psdAtZeroNumBlocks().size()
           << ");"
           << std::endl;
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      ofsvar << m_name << "_psdAtZeroNumBlocks_sub" << m_env.subIdString() << "(" << 1
             << ","                                                               << numBlocksId+1
             << ") = "                                                            << statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]
             << ";"
             << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofsvar << m_name << "_psdAtZeroInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                                                                                                    << statisticalOptions.psdAtZeroNumBlocks().size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          ofsvar << m_name << "_psdAtZeroInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << i+1
                 << ","                                                                                                            << numBlocksId+1
                 << ") = "                                                                                                         << _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]
                 << ";"
                 << std::endl;
        }
      }
    }
  } 

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeGeweke(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing Geweke coefficients"
                            << std::endl;
  }

  std::vector<V*> vectorOfGeweke(initialPosForStatistics.size(),NULL);
  V gewVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->geweke(initialPosition,
                 statisticalOptions.gewekeNaRatio(),
                 statisticalOptions.gewekeNbRatio(),
                 gewVec);
    vectorOfGeweke[initialPosId] = new V(gewVec);
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nComputed Geweke coefficients with 10% and 50% percentages"
                            << " (each column corresponds to a different initial position on the full chain)"
                            << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    *m_env.subDisplayFile() << line;
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      sprintf(line,"%10s%3d",
              " ",
              initialPosForStatistics[initialPosId]);
      *m_env.subDisplayFile() << line;
    }

    for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
      sprintf(line,"\n%9.9s",
              m_vectorSpace.localComponentName(i).c_str() /*.*/);
      *m_env.subDisplayFile() << line;
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        sprintf(line,"%2s%11.4e",
                " ",
                (*(vectorOfGeweke[initialPosId]))[i]);
        *m_env.subDisplayFile() << line;
      }
    }
    *m_env.subDisplayFile() << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain Geweke took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeMeanStacc(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing mean statistical accuracy"
                            << std::endl;
  }

  std::vector<V*> vectorOfMeanStacc(initialPosForStatistics.size(),NULL);
  V meanStaccVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->meanStacc(initialPosition,
                    meanStaccVec);
    vectorOfMeanStacc[initialPosId] = new V(meanStaccVec);
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nComputed mean statistical accuracy"
                            << " (each column corresponds to a different initial position on the full chain)"
                            << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    *m_env.subDisplayFile() << line;
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      sprintf(line,"%10s%3d",
              " ",
              initialPosForStatistics[initialPosId]);
      *m_env.subDisplayFile() << line;
    }

    for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
      sprintf(line,"\n%9.9s",
              m_vectorSpace.localComponentName(i).c_str() /*.*/);
      *m_env.subDisplayFile() << line;
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        sprintf(line,"%2s%11.4e",
                " ",
                (*(vectorOfMeanStacc[initialPosId]))[i]);
        *m_env.subDisplayFile() << line;
      }
    }
    *m_env.subDisplayFile() << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain mean statistical accuracy took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaDef(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing autocorrelation coefficients (via def)"
                            << std::endl;
  }

  if (statisticalOptions.autoCorrDisplay() && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaDef(): lags for autocorrelation (via def) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      *m_env.subDisplayFile() << " " << lagsForCorrs[i];
    }
    *m_env.subDisplayFile() << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
  for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
      _2dArrayOfAutoCorrs.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  //V corrVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      unsigned int lag = lagsForCorrs[lagId];
      this->autoCorrViaDef(initialPos,
                           this->subSequenceSize()-initialPos,
                           lag,
                           _2dArrayOfAutoCorrs(initialPosId,lagId));
      //_2dArrayOfAutoCorrs(initialPosId,lagId) = corrVec;
    }
  }

  // It is not practical to compute the variance of sample mean by computing the autocorrelations via definition for each lag
  // The code computes the variance of sample mean by computing the autocorrelations via fft, below, in another routine

  if ((statisticalOptions.autoCorrDisplay()) && (m_env.subDisplayFile())) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      *m_env.subDisplayFile() << "\nComputed autocorrelation coefficients (via def), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                              << " (each column corresponds to a different lag)"
                              << std::endl;
      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subDisplayFile() << line;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        sprintf(line,"%10s%3d",
                " ",
                lagsForCorrs[lagId]);
        *m_env.subDisplayFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.localComponentName(i).c_str() /*.*/);
        *m_env.subDisplayFile() << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
          *m_env.subDisplayFile() << line;
        }
      }
      *m_env.subDisplayFile() << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain autocorrelation (via def) took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.autoCorrWrite() && passedOfs) {
    std::ofstream& ofsvar = *passedOfs;
    ofsvar << m_name << "_corrViaDefLags_sub" << m_env.subIdString() << " = zeros(" << 1
           << ","                                                                   << lagsForCorrs.size()
           << ");"
           << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofsvar << m_name << "_corrViaDefLags_sub" << m_env.subIdString() << "(" << 1
             << ","                                                           << lagId+1
             << ") = "                                                        << lagsForCorrs[lagId]
             << ";"
             << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofsvar << m_name << "_corrViaDefInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                                                                                                     << lagsForCorrs.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofsvar << m_name << "_corrViaDefInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << i+1
                 << ","                                                                                                             << lagId+1
                 << ") = "                                                                                                          << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
                 << ";"
                 << std::endl;
        }
      }
    }
  } 

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaFFT(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing autocorrelation coefficients (via fft)"
                            << std::endl;
  }

  if (statisticalOptions.autoCorrDisplay() && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaFFT(): lags for autocorrelation (via fft) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      *m_env.subDisplayFile() << " " << lagsForCorrs[i];
     }
     *m_env.subDisplayFile() << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
  for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
      _2dArrayOfAutoCorrs.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  std::vector<V*> corrVecs(lagsForCorrs.size(),NULL);
  std::vector<V*> corrSumVecs(initialPosForStatistics.size(),NULL);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    corrSumVecs[initialPosId] = new V(m_vectorSpace.zeroVector()) /*.*/;
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      corrVecs[lagId] = new V(m_vectorSpace.zeroVector()) /*.*/;
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaFFT()"
                              << ": about to call chain.autoCorrViaFft()"
                              << " with initialPos = "      << initialPos
                              << ", numPos = "              << this->subSequenceSize()-initialPos
                              << ", lagsForCorrs.size() = " << lagsForCorrs.size()
                              << ", corrVecs.size() = "     << corrVecs.size()
                              << std::endl;
    }
    this->autoCorrViaFft(initialPos,
                         this->subSequenceSize()-initialPos, // Use all possible data positions
                         lagsForCorrs,
                         corrVecs);
    this->autoCorrViaFft(initialPos,
                         this->subSequenceSize()-initialPos, // Use all possible data positions
                         (unsigned int) (1.0 * (double) (this->subSequenceSize()-initialPos)), // CHECK
                         *corrSumVecs[initialPosId]); // Sum of all possibly computable autocorrelations, not only the asked ones in lagsForCorrs
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      _2dArrayOfAutoCorrs(initialPosId,lagId) = *(corrVecs[lagId]);
    }
  }
  for (unsigned int j = 0; j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }

  if (statisticalOptions.autoCorrDisplay()) {
    V subChainMean                    (m_vectorSpace.zeroVector());
    V subChainSampleVariance          (m_vectorSpace.zeroVector());
    V estimatedVarianceOfSampleMean(m_vectorSpace.zeroVector());
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];

      this->subMean(initialPos,
                    this->subSequenceSize()-initialPos,
                    subChainMean);

      this->subSampleVariance(initialPos,
                              this->subSequenceSize()-initialPos,
                              subChainMean,
                              subChainSampleVariance);

      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "\nEstimated variance of sample mean, through autocorrelation (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                                << std::endl;
      }
      estimatedVarianceOfSampleMean.cwSet(-1.); // Yes, '-1' because the autocorrelation at lag 0, which values '+1', is already counted in the sum
      estimatedVarianceOfSampleMean += 2.* (*corrSumVecs[initialPosId]);
      estimatedVarianceOfSampleMean *= subChainSampleVariance;
      estimatedVarianceOfSampleMean /= (double) (this->subSequenceSize() - initialPos);
      bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
      estimatedVarianceOfSampleMean.setPrintHorizontally(false);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << estimatedVarianceOfSampleMean
                                << std::endl;
      }
      estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);

      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "\nComputed autocorrelation coefficients (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                                << " (each column corresponds to a different lag)"
                                << std::endl;

        char line[512];
        sprintf(line,"%s",
                "Parameter");
        *m_env.subDisplayFile() << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%10s%3d",
                  " ",
                  lagsForCorrs[lagId]);
          *m_env.subDisplayFile() << line;
        }

        for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
          sprintf(line,"\n%9.9s",
                  m_vectorSpace.localComponentName(i).c_str() /*.*/);
          *m_env.subDisplayFile() << line;
          for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
                    _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
            *m_env.subDisplayFile() << line;
          }
        }
        *m_env.subDisplayFile() << std::endl;
      }
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain autocorrelation (via fft) took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.autoCorrWrite() && passedOfs) {
    std::ofstream& ofsvar = *passedOfs;
    ofsvar << m_name << "_corrViaFftLags_sub" << m_env.subIdString() << " = zeros(" << 1
           << ","                                                                   << lagsForCorrs.size()
           << ");"
           << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofsvar << m_name << "_corrViaFftLags_sub" << m_env.subIdString() << "(" << 1
             << ","                                                           << lagId+1
             << ") = "                                                        << lagsForCorrs[lagId]
             << ";"
             << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofsvar << m_name << "_corrViaFftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                                                                                                     << lagsForCorrs.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofsvar << m_name << "_corrViaFftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << i+1
                 << ","                                                                                                             << lagId+1
                 << ") = "                                                                                                          << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
                 << ";"
                 << std::endl;
        }
      }
    }
  } 

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeFilterParams(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs,
  unsigned int&                         initialPos,
  unsigned int&                         spacing)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing filter parameters for chain " << m_name << " ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  bool okSituation = ((passedOfs == NULL                            ) ||
                      ((passedOfs != NULL) && (m_env.subRank() >= 0)));
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.fullRank(),
                      "uqBaseVectorSequenceClass<V,M>::computeFilterParams()",
                      "unexpected combination of file pointer and subRank");

  //initialPos = 0;
  spacing    = 1;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing filter parameters for chain " << m_name
                            << ": initialPos = " << initialPos
                            << ", spacing = "    << spacing
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde( // Use the whole chain
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing histogram and/or cdf stacc and/or KDE for chain " << m_name << " ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;

  //****************************************************
  // Compute MIN and MAX: for histograms and KDE
  //****************************************************
  double tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing min and max for histograms and KDE"
                            << std::endl;
  }

  V statsMinPositions(m_vectorSpace.zeroVector());
  V statsMaxPositions(m_vectorSpace.zeroVector());
  this->subMinMax(0, // Use the whole chain
                  statsMinPositions,
                  statsMaxPositions);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nComputed min values and max values for chain " << m_name
                            << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    *m_env.subDisplayFile() << line;

    sprintf(line,"%9s%s%9s%s",
            " ",
            "min",
            " ",
            "max");
    *m_env.subDisplayFile() << line;

    for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
      sprintf(line,"\n%8.8s",
              m_vectorSpace.localComponentName(i).c_str() /*.*/);
      *m_env.subDisplayFile() << line;

      sprintf(line,"%2s%11.4e%2s%11.4e",
              " ",
              statsMinPositions[i],
              " ",
              statsMaxPositions[i]);
      *m_env.subDisplayFile() << line;
    }
    *m_env.subDisplayFile() << std::endl;
  }

  V unifiedStatsMinPositions(statsMinPositions);
  V unifiedStatsMaxPositions(statsMaxPositions);
  if (m_env.numSubEnvironments() > 1) {
    // Compute unified min-max
    this->unifiedMinMax(0, // Use the whole chain
                        unifiedStatsMinPositions,
                        unifiedStatsMaxPositions);

    // Write unified min-max
    if (m_env.subDisplayFile()) {
      if (m_vectorSpace.numOfProcsForStorage() == 1) {
        if (m_env.inter0Rank() == 0) {
          *m_env.subDisplayFile() << "\nComputed unified min values and max values for chain " << m_name
                                  << std::endl;

          char line[512];
          sprintf(line,"%s",
                  "Parameter");
          *m_env.subDisplayFile() << line;

          sprintf(line,"%9s%s%9s%s",
                  " ",
                  "min",
                  " ",
                  "max");
          *m_env.subDisplayFile() << line;

          for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
            sprintf(line,"\n%8.8s",
                    m_vectorSpace.localComponentName(i).c_str() /*.*/);
            *m_env.subDisplayFile() << line;

            sprintf(line,"%2s%11.4e%2s%11.4e",
                    " ",
                    unifiedStatsMinPositions[i],
                    " ",
                    unifiedStatsMaxPositions[i]);
            *m_env.subDisplayFile() << line;
          }
          *m_env.subDisplayFile() << std::endl;
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.fullRank(),
                            "uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()",
                            "unified min-max writing, parallel vectors not supported yet");
      }
    } // if subDisplayFile
  } // if numSubEnvs > 1

  //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...
  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Chain min and max took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  //****************************************************
  // Compute histograms
  //****************************************************
  if ((statisticalOptions.histCompute()            ) &&
      (statisticalOptions.histNumInternalBins() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                              << "\nComputing histograms"
                              << std::endl;
    }

    std::string subCoreName_HistCenters((std::string)(    "_HistCenters_sub")+m_env.subIdString());
    std::string uniCoreName_HistCenters((std::string)("_unifHistCenters_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_HistCenters = uniCoreName_HistCenters;

    std::string subCoreName_HistQuantts((std::string)(    "_HistQuantts_sub")+m_env.subIdString());
    std::string uniCoreName_HistQuantts((std::string)("_unifHistQuantts_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_HistQuantts = uniCoreName_HistQuantts;

    for (unsigned int i = 0; i < statsMaxPositions.sizeLocal(); ++i) {
      statsMaxPositions[i] *= (1. + 1.e-15);
    }

    // Compute histograms
    std::vector<V*> histCentersForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    std::vector<V*> histQuanttsForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    this->subHistogram(0, // Use the whole chain
                       statsMinPositions,
                       statsMaxPositions,
                       histCentersForAllBins,
                       histQuanttsForAllBins);

    // Write histograms
    if (passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << subCoreName_HistCenters << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                              << histCentersForAllBins.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int j = 0; j < histCentersForAllBins.size(); ++j) {
           ofsvar << m_name << subCoreName_HistCenters << "(" << i+1
                  << ","                                      << j+1
                  << ") = "                                   << (*(histCentersForAllBins[j]))[i]
                  << ";"
                  << std::endl;
        }
      }

      ofsvar << m_name << subCoreName_HistQuantts << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                              << histQuanttsForAllBins.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int j = 0; j < histQuanttsForAllBins.size(); ++j) {
           ofsvar << m_name << subCoreName_HistQuantts << "(" << i+1
                  << ","                                      << j+1
                  << ") = "                                   << (*(histQuanttsForAllBins[j]))[i]
                  << ";"
                  << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < histQuanttsForAllBins.size(); ++i) {
      if (histQuanttsForAllBins[i] != NULL) delete histQuanttsForAllBins[i];
    }
    for (unsigned int i = 0; i < histCentersForAllBins.size(); ++i) {
      if (histCentersForAllBins[i] != NULL) delete histCentersForAllBins[i];
    }

    std::vector<V*> unifiedHistCentersForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    std::vector<V*> unifiedHistQuanttsForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    if (m_env.numSubEnvironments() > 1) {
      // Compute unified histogram
      this->unifiedHistogram(0, // Use the whole chain
                             unifiedStatsMinPositions,
                             unifiedStatsMaxPositions,
                             unifiedHistCentersForAllBins,
                             unifiedHistQuanttsForAllBins);

      // Write unified histogram
      if (passedOfs) {
        if (m_vectorSpace.numOfProcsForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            std::ofstream& ofsvar = *passedOfs;
            ofsvar << m_name << uniCoreName_HistCenters << " = zeros(" << this->vectorSizeLocal() /*.*/
                   << ","                                              << unifiedHistCentersForAllBins.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedHistCentersForAllBins.size(); ++j) {
                 ofsvar << m_name << uniCoreName_HistCenters << "(" << i+1
                        << ","                                      << j+1
                        << ") = "                                   << (*(unifiedHistCentersForAllBins[j]))[i]
                        << ";"
                        << std::endl;
              }
            }

            ofsvar << m_name << uniCoreName_HistQuantts << " = zeros(" << this->vectorSizeLocal() /*.*/
                   << ","                                              << unifiedHistQuanttsForAllBins.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedHistQuanttsForAllBins.size(); ++j) {
                 ofsvar << m_name << uniCoreName_HistQuantts << "(" << i+1
                        << ","                                      << j+1
                        << ") = "                                   << (*(unifiedHistQuanttsForAllBins[j]))[i]
                        << ";"
                        << std::endl;
              }
            }
          }
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()",
                              "unified histogram writing, parallel vectors not supported yet");
        }
      } // if passedOfs

      for (unsigned int i = 0; i < unifiedHistQuanttsForAllBins.size(); ++i) {
        if (unifiedHistQuanttsForAllBins[i] != NULL) delete unifiedHistQuanttsForAllBins[i];
      }
      for (unsigned int i = 0; i < unifiedHistCentersForAllBins.size(); ++i) {
        if (unifiedHistCentersForAllBins[i] != NULL) delete unifiedHistCentersForAllBins[i];
      }
    } // if numSubEnvs > 1

    //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...
    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "Chain histograms took " << tmpRunTime
                              << " seconds"
                              << std::endl;
    }
  }

  //****************************************************
  // Compute cdf statistical accuracy
  //****************************************************
  if ((statisticalOptions.cdfStaccCompute()             ) &&
      (statisticalOptions.cdfStaccNumEvalPositions() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                              << "\nComputing cdf statistical accuracy"
                              << std::endl;
    }

    std::vector<V*> cdfStaccEvalPositions(statisticalOptions.cdfStaccNumEvalPositions(),NULL);
    uqMiscComputePositionsBetweenMinMax(statsMinPositions,
                                        statsMaxPositions,
                                        cdfStaccEvalPositions);

    std::vector<V*> cdfStaccValues(statisticalOptions.cdfStaccNumEvalPositions(),NULL);
    this->subCdfStacc(0, // Use the whole chain
                      cdfStaccEvalPositions,
                      cdfStaccValues);

    //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...
    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "Chain cdf statistical accuracy took " << tmpRunTime
                              << " seconds"
                              << std::endl;
    }
  }

  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((statisticalOptions.kdeCompute()             ) &&
      (statisticalOptions.kdeNumEvalPositions() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                              << "\nComputing KDE"
                              << std::endl;
    }

    std::string subCoreName_GaussianKdePositions((std::string)(    "_GkdePosits_sub")+m_env.subIdString());
    std::string uniCoreName_GaussianKdePositions((std::string)("_unifGkdePosits_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_GaussianKdePositions = uniCoreName_GaussianKdePositions; // avoid temporarily (see '< -1' below)

    std::string subCoreName_GaussianKdeScaleVec ((std::string)(    "_GkdeScalev_sub")+m_env.subIdString());
    std::string uniCoreName_GaussianKdeScaleVec ((std::string)("_unifGkdeScalev_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_GaussianKdeScaleVec = uniCoreName_GaussianKdeScaleVec; // avoid temporarily (see '< -1' below)

    std::string subCoreName_GaussianKdeValues   ((std::string)(    "_GkdeValues_sub")+m_env.subIdString());
    std::string uniCoreName_GaussianKdeValues   ((std::string)("_unifGkdeValues_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_GaussianKdeValues = uniCoreName_GaussianKdeValues; // avoid temporarily (see '< -1' below)

    V iqrVec(m_vectorSpace.zeroVector());
    this->subInterQuantileRange(0, // Use the whole chain
                                iqrVec);

    V gaussianKdeScaleVec(m_vectorSpace.zeroVector());
    this->subScalesForKDE(0, // Use the whole chain
                          iqrVec,
                          1,
                          gaussianKdeScaleVec);

    std::vector<V*> kdeEvalPositions(statisticalOptions.kdeNumEvalPositions(),NULL);
    uqMiscComputePositionsBetweenMinMax(statsMinPositions,
                                        statsMaxPositions,
                                        kdeEvalPositions);

    std::vector<V*> gaussianKdeDensities(statisticalOptions.kdeNumEvalPositions(),NULL);
    this->subGaussianKDE(0, // Use the whole chain
                         gaussianKdeScaleVec,
                         kdeEvalPositions,
                         gaussianKdeDensities);

    // Write iqr
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "\nComputed inter quantile ranges for chain " << m_name
                              << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subDisplayFile() << line;

      sprintf(line,"%9s%s",
              " ",
              "iqr");
      *m_env.subDisplayFile() << line;

      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        sprintf(line,"\n%8.8s",
                m_vectorSpace.localComponentName(i).c_str() /*.*/);
        *m_env.subDisplayFile() << line;

        sprintf(line,"%2s%11.4e",
                " ",
                iqrVec[i]);
        *m_env.subDisplayFile() << line;
      }
      *m_env.subDisplayFile() << std::endl;
    }

    // Write KDE
    if (passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << subCoreName_GaussianKdePositions << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                                       << kdeEvalPositions.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int j = 0; j < kdeEvalPositions.size(); ++j) {
          ofsvar << m_name << subCoreName_GaussianKdePositions << "(" << i+1
                 << ","                                               << j+1
                 << ") = "                                            << (*(kdeEvalPositions[j]))[i]
                 << ";"
                 << std::endl;
        }
      }

      ofsvar << m_name << subCoreName_GaussianKdeScaleVec << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        ofsvar << m_name << subCoreName_GaussianKdeScaleVec << "(" << i+1
               << ") = "                                           << gaussianKdeScaleVec[i]
               << ";"
               << std::endl;
      }

      ofsvar << m_name << subCoreName_GaussianKdeValues << " = zeros(" << this->vectorSizeLocal() /*.*/
             << ","                                                    << gaussianKdeDensities.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
        for (unsigned int j = 0; j < gaussianKdeDensities.size(); ++j) {
          ofsvar << m_name << subCoreName_GaussianKdeValues << "(" << i+1
                 << ","                                            << j+1
                 << ") = "                                         << (*(gaussianKdeDensities[j]))[i]
                 << ";"
                 << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < gaussianKdeDensities.size(); ++i) {
      if (gaussianKdeDensities[i] != NULL) delete gaussianKdeDensities[i];
    }
    for (unsigned int i = 0; i < kdeEvalPositions.size(); ++i) {
      if (kdeEvalPositions[i] != NULL) delete kdeEvalPositions[i];
    }

    if ((int) m_env.numSubEnvironments() > 1) { // avoid code temporarily
      // Compute unified KDE
      V unifiedIqrVec(m_vectorSpace.zeroVector());
      this->unifiedInterQuantileRange(0, // Use the whole chain
                                      unifiedIqrVec);
      //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

      V unifiedGaussianKdeScaleVec(m_vectorSpace.zeroVector());
      this->unifiedScalesForKDE(0, // Use the whole chain
                                unifiedIqrVec,
                                1,
                                unifiedGaussianKdeScaleVec);
      //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

      std::vector<V*> unifiedKdeEvalPositions(statisticalOptions.kdeNumEvalPositions(),NULL);
      uqMiscComputePositionsBetweenMinMax(unifiedStatsMinPositions,
                                          unifiedStatsMaxPositions,
                                          unifiedKdeEvalPositions);

      std::vector<V*> unifiedGaussianKdeDensities(statisticalOptions.kdeNumEvalPositions(),NULL);
      this->unifiedGaussianKDE(0, // Use the whole chain
                               unifiedGaussianKdeScaleVec,
                               unifiedKdeEvalPositions,
                               unifiedGaussianKdeDensities);
      //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

      // Write unified iqr
      if (m_env.subDisplayFile()) {
        if (m_vectorSpace.numOfProcsForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            *m_env.subDisplayFile() << "\nComputed unified inter quantile ranges for chain " << m_name
                                    << std::endl;

            char line[512];
            sprintf(line,"%s",
                    "Parameter");
            *m_env.subDisplayFile() << line;

            sprintf(line,"%9s%s",
                    " ",
                    "iqr");
            *m_env.subDisplayFile() << line;

            for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
              sprintf(line,"\n%8.8s",
                      m_vectorSpace.localComponentName(i).c_str() /*.*/);
              *m_env.subDisplayFile() << line;

              sprintf(line,"%2s%11.4e",
                      " ",
                      unifiedIqrVec[i]);
              *m_env.subDisplayFile() << line;
            }
            *m_env.subDisplayFile() << std::endl;
          }
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()",
                              "unified iqr writing, parallel vectors not supported yet");
        }
      }

      // Write unified KDE
      if (passedOfs) {
        if (m_vectorSpace.numOfProcsForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            std::ofstream& ofsvar = *passedOfs;
            ofsvar << m_name << uniCoreName_GaussianKdePositions << " = zeros(" << this->vectorSizeLocal() /*.*/
                   << ","                                                       << unifiedKdeEvalPositions.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedKdeEvalPositions.size(); ++j) {
                ofsvar << m_name << uniCoreName_GaussianKdePositions << "(" << i+1
                       << ","                                               << j+1
                       << ") = "                                            << (*(unifiedKdeEvalPositions[j]))[i]
                       << ";"
                       << std::endl;
              }
            }

            ofsvar << m_name << uniCoreName_GaussianKdeScaleVec << " = zeros(" << this->vectorSizeLocal() /*.*/
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
              ofsvar << m_name << uniCoreName_GaussianKdeScaleVec << "(" << i+1
                     << ") = "                                           << unifiedGaussianKdeScaleVec[i]
                     << ";"
                     << std::endl;
            }

            ofsvar << m_name << uniCoreName_GaussianKdeValues << " = zeros(" << this->vectorSizeLocal() /*.*/
                   << ","                                                    << unifiedGaussianKdeDensities.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedGaussianKdeDensities.size(); ++j) {
                ofsvar << m_name << uniCoreName_GaussianKdeValues << "(" << i+1
                       << ","                                            << j+1
                       << ") = "                                         << (*(unifiedGaussianKdeDensities[j]))[i]
                       << ";"
                       << std::endl;
              }
            }
          }
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()",
                              "unified KDE writing, parallel vectors not supported yet");
        }
      }

      for (unsigned int i = 0; i < unifiedGaussianKdeDensities.size(); ++i) {
        if (unifiedGaussianKdeDensities[i] != NULL) delete unifiedGaussianKdeDensities[i];
      }
      for (unsigned int i = 0; i < unifiedKdeEvalPositions.size(); ++i) {
        if (unifiedKdeEvalPositions[i] != NULL) delete unifiedKdeEvalPositions[i];
      }
    }

    //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...
    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "Chain KDE took " << tmpRunTime
                              << " seconds"
                              << std::endl;
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing histogram and/or cdf stacc and/or KDE for chain " << m_name
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices( // Use the whole chain
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing covariance and correlation matrices for chain " << m_name << " ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  //int iRC = UQ_OK_RC;
  //struct timeval timevalTmp;
  M* covarianceMatrix = new M(m_env,
                              m_vectorSpace.map(),        // number of rows
                              m_vectorSpace.dimGlobal()); // number of cols
  M* correlationMatrix = new M(m_env,
                               m_vectorSpace.map(),        // number of rows
                               m_vectorSpace.dimGlobal()); // number of cols

  uqComputeCovCorrMatricesBetweenVectorSequences(*this,
                                                 *this,
                                                 this->subSequenceSize(),
                                                 *covarianceMatrix,
                                                 *correlationMatrix);

  if (m_env.subDisplayFile()) {
    if (m_vectorSpace.numOfProcsForStorage() == 1) {
      // Only unified covariance matrix is written. And only one processor writes it.
      if (m_env.inter0Rank() == 0) {
        *m_env.subDisplayFile() << "\nuqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices"
                                << ", chain " << m_name
                                << ": contents of covariance matrix are\n" << *covarianceMatrix
                                << std::endl;

        *m_env.subDisplayFile() << "\nuqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices"
                                << ", chain " << m_name
                                << ": contents of correlation matrix are\n" << *correlationMatrix
                                << std::endl;
      }
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.fullRank(),
                          "uqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices()",
                          "parallel vectors not supported yet");
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing covariance and correlation matrices for chain " << m_name
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}

template <class P_V, class P_M, class Q_V, class Q_M>
void
uqComputeCovCorrMatricesBetweenVectorSequences(
  const uqBaseVectorSequenceClass<P_V,P_M>& subPSeq,
  const uqBaseVectorSequenceClass<Q_V,Q_M>& subQSeq,
        unsigned int                        subNumSamples,
        P_M&                                pqCovMatrix,
        P_M&                                pqCorrMatrix)
{
  // Check input data consistency
  const uqBaseEnvironmentClass& env = subPSeq.vectorSpace().zeroVector().env();

  bool useOnlyInter0Comm = (subPSeq.vectorSpace().numOfProcsForStorage() == 1) &&
                           (subQSeq.vectorSpace().numOfProcsForStorage() == 1);

  UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                      env.fullRank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "parallel vectors not supported yet");

  unsigned int numRowsLocal = subPSeq.vectorSpace().dimLocal();
  unsigned int numCols = subQSeq.vectorSpace().dimGlobal();

  UQ_FATAL_TEST_MACRO((numRowsLocal != pqCovMatrix.numRowsLocal()) || (numCols != pqCovMatrix.numCols()),
                      env.fullRank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "inconsistent dimensions for covariance matrix");

  UQ_FATAL_TEST_MACRO((numRowsLocal != pqCorrMatrix.numRowsLocal()) || (numCols != pqCorrMatrix.numCols()),
                      env.fullRank(),
                      "uqComputeCorrelationBetweenVectorSequences()",
                      "inconsistent dimensions for correlation matrix");

  UQ_FATAL_TEST_MACRO((subNumSamples > subPSeq.subSequenceSize()) || (subNumSamples > subQSeq.subSequenceSize()),
                      env.fullRank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "subNumSamples is too large");

  // For both P and Q vector sequences: fill them
  P_V tmpP(subPSeq.vectorSpace().zeroVector());
  Q_V tmpQ(subQSeq.vectorSpace().zeroVector());

  // For both P and Q vector sequences: compute the unified mean
  P_V unifiedMeanP(subPSeq.vectorSpace().zeroVector());
  subPSeq.unifiedMean(0,subNumSamples,unifiedMeanP);

  Q_V unifiedMeanQ(subQSeq.vectorSpace().zeroVector());
  subQSeq.unifiedMean(0,subNumSamples,unifiedMeanQ);

  // Compute "sub" covariance matrix
  for (unsigned i = 0; i < numRowsLocal; ++i) {
    for (unsigned j = 0; j < numCols; ++j) {
      pqCovMatrix(i,j) = 0.;
    }
  }
  for (unsigned k = 0; k < subNumSamples; ++k) {
    // For both P and Q vector sequences: get the difference (wrt the unified mean) in them
    subPSeq.getPositionValues(k,tmpP);
    tmpP -= unifiedMeanP;

    subQSeq.getPositionValues(k,tmpQ);
    tmpQ -= unifiedMeanQ;

    for (unsigned i = 0; i < numRowsLocal; ++i) {
      for (unsigned j = 0; j < numCols; ++j) {
        pqCovMatrix(i,j) += tmpP[i]*tmpQ[j];
      }
    }
  }

  // For both P and Q vector sequences: compute the unified variance
  P_V unifiedSampleVarianceP(subPSeq.vectorSpace().zeroVector());
  subPSeq.unifiedSampleVariance(0,
                                subNumSamples,
                                unifiedMeanP,
                                unifiedSampleVarianceP);

  Q_V unifiedSampleVarianceQ(subQSeq.vectorSpace().zeroVector());
  subQSeq.unifiedSampleVariance(0,
                                subNumSamples,
                                unifiedMeanQ,
                                unifiedSampleVarianceQ);

  // Compute unified covariance matrix
  if (useOnlyInter0Comm) {
    if (env.inter0Rank() >= 0) {
      unsigned int unifiedNumSamples = 0;
      int mpiRC = MPI_Allreduce((void *) &subNumSamples, (void *) &unifiedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          env.fullRank(),
                          "uqComputeCovCorrMatricesBetweenVectorSequences()",
                          "failed MPI_Allreduce() for subNumSamples");

      for (unsigned i = 0; i < numRowsLocal; ++i) {
        for (unsigned j = 0; j < numCols; ++j) {
          double aux = 0.;
          int mpiRC = MPI_Allreduce((void *) &pqCovMatrix(i,j), (void *) &aux, (int) 1, MPI_DOUBLE, MPI_SUM, env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              env.fullRank(),
                              "uqComputeCovCorrMatricesBetweenVectorSequences()",
                              "failed MPI_Allreduce() for a matrix position");
          pqCovMatrix(i,j) = aux/((double) (unifiedNumSamples-1)); // Yes, '-1' in order to compensate for the 'N-1' denominator factor in the calculations of sample variances above (whose square roots will be used below)
        }
      }

      for (unsigned i = 0; i < numRowsLocal; ++i) {
        for (unsigned j = 0; j < numCols; ++j) {
          pqCorrMatrix(i,j) = pqCovMatrix(i,j)/std::sqrt(unifiedSampleVarianceP[i])/std::sqrt(unifiedSampleVarianceQ[j]);
          if (((pqCorrMatrix(i,j) + 1.) < -1.e-8) ||
              ((pqCorrMatrix(i,j) - 1.) >  1.e-8)) {
            if (env.inter0Rank() == 0) {
              std::cerr << "In uqComputeCovCorrMatricesBetweenVectorSequences()"
                        << ": fullRank = "            << env.fullRank()
                        << ", i = "                   << i
                        << ", j = "                   << j
                        << ", pqCorrMatrix(i,j)+1 = " << pqCorrMatrix(i,j)+1.
                        << ", pqCorrMatrix(i,j)-1 = " << pqCorrMatrix(i,j)-1.
                        << std::endl;
            }
            env.inter0Comm().Barrier();
          }
          UQ_FATAL_TEST_MACRO(((pqCorrMatrix(i,j) + 1.) < -1.e-8) ||
                              ((pqCorrMatrix(i,j) - 1.) >  1.e-8),
                               env.fullRank(),
                               "uqComputeCovCorrMatricesBetweenVectorSequences()",
                               "computed correlation is out of range");
        }
      }
    }
    else {
      // Node not in the 'inter0' communicator: do nothing extra
    }
  }
  else {
    UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                        env.fullRank(),
                        "uqComputeCovCorrMatricesBetweenVectorSequences()",
                        "parallel vectors not supported yet (2)");
  }

  return;
}
#endif // __UQ_VECTOR_SEQUENCE_H__
