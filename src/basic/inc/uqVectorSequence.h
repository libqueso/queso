//-----------------------------------------------------------------------Bl-
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

#ifndef __UQ_VECTOR_SEQUENCE_H__
#define __UQ_VECTOR_SEQUENCE_H__

#undef UQ_CODE_HAS_MONITORS

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

  virtual  unsigned int             subSequenceSize             () const = 0;
           unsigned int             unifiedSequenceSize         () const;
           unsigned int             vectorSizeLocal             () const;
           unsigned int             vectorSizeGlobal            () const;
  const    uqVectorSpaceClass<V,M>& vectorSpace                 () const;
  const    std::string&             name                        () const;
           void                     setName                     (const std::string& newName);
           void                     clear                       ();
  const    V&                       subMinPlain                 () const;
  const    V&                       unifiedMinPlain             () const;
  const    V&                       subMaxPlain                 () const;
  const    V&                       unifiedMaxPlain             () const;
  const    V&                       subMeanPlain                () const;
  const    V&                       unifiedMeanPlain            () const;
  const    V&                       subMedianPlain              () const;
  const    V&                       unifiedMedianPlain          () const;
  const    V&                       subSampleVariancePlain      () const;
  const    V&                       unifiedSampleVariancePlain  () const;
  const    uqBoxSubsetClass<V,M>&   subBoxPlain                 () const;
  const    uqBoxSubsetClass<V,M>&   unifiedBoxPlain             () const;
           void                     deleteStoredVectors         ();

           void                     append                      (const uqBaseVectorSequenceClass<V,M>& src,
                                                                       unsigned int                    initialPos,
                                                                       unsigned int                    numPos);             /* This routine deletes all stored computed vectors */

           double                   subPositionsOfMaximum       (const uqScalarSequenceClass<double>& subCorrespondingScalarValues,
                                                                 uqBaseVectorSequenceClass<V,M>&      subPositionsOfMaximum);
           double                   unifiedPositionsOfMaximum   (const uqScalarSequenceClass<double>& subCorrespondingScalarValues,
                                                                 uqBaseVectorSequenceClass<V,M>&      unifiedPositionsOfMaximum);

  virtual  void                     resizeSequence              (unsigned int newSubSequenceSize) = 0;                      /* This routine deletes all stored computed vectors */
  virtual  void                     resetValues                 (unsigned int initialPos, unsigned int numPos) = 0;         /* This routine deletes all stored computed vectors */
  virtual  void                     erasePositions              (unsigned int initialPos, unsigned int numPos) = 0;         /* This routine deletes all stored computed vectors */
  virtual  void                     getPositionValues           (unsigned int posId,       V& vec) const = 0;
  virtual  void                     setPositionValues           (unsigned int posId, const V& vec) = 0;                     /* This routine deletes all stored computed vectors */
           void                     setGaussian                 (const V& meanVec, const V& stdDevVec); /* This routine deletes all stored computed vectors */
           void                     setUniform                  (const V& aVec,    const V& bVec     ); /* This routine deletes all stored computed vectors */
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  virtual  void                     subUniformlySampledMdf      (const V&                       numEvaluationPointsVec,
                                                                 uqArrayOfOneDGridsClass <V,M>& mdfGrids,
                                                                 uqArrayOfOneDTablesClass<V,M>& mdfValues) const = 0;
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  virtual  void                     subUniformlySampledCdf      (const V&                       numEvaluationPointsVec,
                                                                 uqArrayOfOneDGridsClass <V,M>& cdfGrids,
                                                                 uqArrayOfOneDTablesClass<V,M>& cdfValues) const = 0;
  virtual  void                     unifiedUniformlySampledCdf  (const V&                       numEvaluationPointsVec,
                                                                 uqArrayOfOneDGridsClass <V,M>& unifiedCdfGrids,
                                                                 uqArrayOfOneDTablesClass<V,M>& unifieddfValues) const = 0;
#endif
  virtual  void                     subMeanExtra                (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 V&                                       meanVec) const = 0;
  virtual  void                     unifiedMeanExtra            (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 V&                                       unifiedMeanVec) const = 0;
  virtual  void                     subMedianExtra              (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 V&                                       medianVec) const = 0;
  virtual  void                     unifiedMedianExtra          (unsigned int                             initialPos,
                                                                 unsigned int                             localNumPos,
                                                                 V&                                       unifiedMedianVec) const = 0;
  virtual  void                     subSampleVarianceExtra      (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const V&                                 meanVec,
                                                                 V&                                       samVec) const = 0;
  virtual  void                     unifiedSampleVarianceExtra  (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const V&                                 unifiedMeanVec,
                                                                 V&                                       unifiedSamVec) const = 0;
  virtual  void                     subPopulationVariance       (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const V&                                 meanVec,
                                                                 V&                                       popVec) const = 0;
  virtual  void                     unifiedPopulationVariance   (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const V&                                 unifiedMeanVec,
                                                                 V&                                       unifiedPopVec) const = 0;

  virtual  void                     autoCovariance              (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const V&                                 meanVec,
                                                                 unsigned int                             lag,
                                                                 V&                                       covVec) const = 0;
  virtual  void                     autoCorrViaDef              (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 unsigned int                             lag,
                                                                 V&                                       corrVec) const = 0;
  virtual  void                     autoCorrViaFft              (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const std::vector<unsigned int>&         lags,
                                                                 std::vector<V*>&                         corrVecs) const = 0;
  virtual  void                     autoCorrViaFft              (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 unsigned int                             numSum,
                                                                 V&                                       autoCorrsSumVec) const = 0;
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  virtual  void                     bmm                         (unsigned int                             initialPos,
                                                                 unsigned int                             batchLength,
                                                                 V&                                       bmmVec) const = 0;
  virtual  void                     fftForward                  (unsigned int                             initialPos,
                                                                 unsigned int                             fftSize,
                                                                 unsigned int                             paramId,
                                                                 std::vector<std::complex<double> >&      fftResult) const = 0;
//virtual  void                     fftInverse                  (unsigned int                             fftSize) = 0;
  virtual  void                     psd                         (unsigned int                             initialPos,
                                                                 unsigned int                             numBlocks,
                                                                 double                                   hopSizeRatio,
                                                                 unsigned int                             paramId,
                                                                 std::vector<double>&                     psdResult) const = 0;
  virtual  void                     psdAtZero                   (unsigned int                             initialPos,
                                                                 unsigned int                             numBlocks,
                                                                 double                                   hopSizeRatio,
                                                                 V&                                       psdVec) const = 0;
  virtual  void                     geweke                      (unsigned int                             initialPos,
                                                                 double                                   ratioNa,
                                                                 double                                   ratioNb,
                                                                 V&                                       gewVec) const = 0;
  virtual  void                     meanStacc                   (unsigned int                             initialPos,
                                                                 V&                                       meanStaccVec) const = 0;
  virtual  void                     subCdfStacc                 (unsigned int                             initialPos,
                                                                 const std::vector<V*>&                   evalPositionsVecs,
                                                                 std::vector<V*>&                         cdfStaccVecs) const = 0;
#endif
  virtual  void                     subMinMaxExtra              (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 V&                                       minVec,
                                                                 V&                                       maxVec) const = 0;
  virtual  void                     unifiedMinMaxExtra          (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 V&                                       unifiedMinVec,
                                                                 V&                                       unifiedMaxVec) const = 0;
  virtual  void                     subHistogram                (unsigned int                             initialPos,
                                                                 const V&                                 minVec,
                                                                 const V&                                 maxVec,
                                                                 std::vector<V*>&                         centersForAllBins,
                                                                 std::vector<V*>&                         quanttsForAllBins) const = 0;
  virtual  void                     unifiedHistogram            (unsigned int                             initialPos,
                                                                 const V&                                 unifiedMinVec,
                                                                 const V&                                 unifiedMaxVec,
                                                                 std::vector<V*>&                         unifiedCentersForAllBins,
                                                                 std::vector<V*>&                         unifiedQuanttsForAllBins) const = 0;
  virtual  void                     subInterQuantileRange       (unsigned int                             initialPos,
                                                                 V&                                       iqrVec) const = 0;
  virtual  void                     unifiedInterQuantileRange   (unsigned int                             initialPos,
                                                                V&                                       unifiedIqrVec) const = 0;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
           void                     computeHistCdfstaccKde      (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 std::ofstream*                           passedOfs);
#endif
  virtual  void                     subScalesForKde             (unsigned int                             initialPos,
                                                                 const V&                                 iqrVec,
                                                                 unsigned int                             kdeDimension,
                                                                 V&                                       scaleVec) const = 0;
  virtual  void                     unifiedScalesForKde         (unsigned int                             initialPos,
                                                                 const V&                                 unifiedIqrVec,
                                                                 unsigned int                             kdeDimension,
                                                                 V&                                       unifiedScaleVec) const = 0;
  virtual  void                     subGaussian1dKde            (unsigned int                             initialPos,
                                                                 const V&                                 scaleVec,
                                                                 const std::vector<V*>&                   evaluationParamVecs,
                                                                 std::vector<V*>&                         densityVecs) const = 0;
  virtual  void                     unifiedGaussian1dKde        (unsigned int                             initialPos,
                                                                 const V&                                 unifiedScaleVec,
                                                                 const std::vector<V*>&                   unifiedEvaluationParamVecs,
                                                                 std::vector<V*>&                         unifiedDensityVecs) const = 0;
  virtual  void                     subWriteContents            (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 const std::string&                       fileName,
                                                                 const std::string&                       fileType,
                                                                 const std::set<unsigned int>&            allowedSubEnvIds) const = 0;
  virtual  void                     subWriteContents            (unsigned int                             initialPos,
                                                                 unsigned int                             numPos,
                                                                 std::ofstream&                           ofsvar,
                                                                 const std::string&                       fileType) const = 0;
  virtual  void                     unifiedWriteContents        (const std::string&                       fileName,
                                                                 const std::string&                       fileType) const = 0;
  virtual  void                     unifiedReadContents         (const std::string&                       fileName,
                                                                 const std::string&                       fileType,
                                                                 const unsigned int                       subSequenceSize) = 0;
  virtual  void                     select                      (const std::vector<unsigned int>&         idsOfUniquePositions) = 0;
  virtual  void                     filter                      (unsigned int                             initialPos,
                                                                 unsigned int                             spacing) = 0;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
           void                     computeStatistics           (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 std::ofstream*                           passedOfs);
#endif
           void                     computeFilterParams         (std::ofstream*                           passedOfs,
                                                                 unsigned int   &                         initialPos,
                                                                 unsigned int   &                         spacing);

  virtual  double                   estimateConvBrooksGelman    (unsigned int                             initialPos,
                                                                 unsigned int                             numPos) const = 0;

  virtual  void                     extractScalarSeq            (unsigned int                             initialPos,
                                                                 unsigned int                             spacing,
                                                                 unsigned int                             numPos,
                                                                 unsigned int                             paramId,
                                                                 uqScalarSequenceClass<double>&           scalarSeq) const = 0;
protected:
           void                     copy                        (const uqBaseVectorSequenceClass<V,M>&    src); /* This routine deletes all stored computed vectors */
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
           void                     computeMeanVars             (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 std::ofstream*                           passedOfs,
                                                                 V*                                       subMeanPtr,
                                                                 V*                                       subMedianPtr,
                                                                 V*                                       subSampleVarPtr,
                                                                 V*                                       subPopulVarPtr);
#endif
#ifdef UQ_CODE_HAS_MONITORS
  virtual  void                     subMeanMonitorAlloc         (unsigned int numberOfMonitorPositions) = 0;
  virtual  void                     subMeanInter0MonitorAlloc   (unsigned int numberOfMonitorPositions) = 0;
  virtual  void                     unifiedMeanMonitorAlloc     (unsigned int numberOfMonitorPositions) = 0;

  virtual  void                     subMeanMonitorRun           (unsigned int monitorPosition,
                                                                 V&           subMeanVec,
                                                                 V&           subMeanCltStd) = 0;
  virtual  void                     subMeanInter0MonitorRun     (unsigned int monitorPosition,
                                                                 V&           subMeanInter0Mean,
                                                                 V&           subMeanInter0Clt95,
                                                                 V&           subMeanInter0Empirical90,
                                                                 V&           subMeanInter0Min,
                                                                 V&           subMeanInter0Max) = 0;
  virtual  void                     unifiedMeanMonitorRun       (unsigned int monitorPosition,
                                                                 V&           unifiedMeanVec,
                                                                 V&           unifiedMeanCltStd) = 0;

  virtual  void                     subMeanMonitorStore         (unsigned int i,
                                                                 unsigned int monitorPosition,
                                                                 const V&     subMeanVec,
                                                                 const V&     subMeanCltStd) = 0;
  virtual  void                     subMeanInter0MonitorStore   (unsigned int i,
                                                                 unsigned int monitorPosition,
                                                                 const V&     subMeanInter0Mean,
                                                                 const V&     subMeanInter0Clt95,
                                                                 const V&     subMeanInter0Empirical90,
                                                                 const V&     subMeanInter0Min,
                                                                 const V&     subMeanInter0Max) = 0;
  virtual  void                     unifiedMeanMonitorStore     (unsigned int i,
                                                                 unsigned int monitorPosition,
                                                                 V&           unifiedMeanVec,
                                                                 V&           unifiedMeanCltStd) = 0;

  virtual  void                     subMeanMonitorWrite         (std::ofstream& ofs) = 0;
  virtual  void                     subMeanInter0MonitorWrite   (std::ofstream& ofs) = 0;
  virtual  void                     unifiedMeanMonitorWrite     (std::ofstream& ofs) = 0;

  virtual  void                     subMeanMonitorFree          () = 0;
  virtual  void                     subMeanInter0MonitorFree    () = 0;
  virtual  void                     unifiedMeanMonitorFree      () = 0;

           void                     computeMeanEvolution        (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 std::ofstream*                           passedOfs);
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
           void                     computeBMM                  (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 std::ofstream*                           passedOfs);
           void                     computeFFT                  (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 std::ofstream*                           passedOfs);
           void                     computePSD                  (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 std::ofstream*                           passedOfs);
           void                     computePSDAtZero            (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 std::ofstream*                           passedOfs);
           void                     computeGeweke               (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 std::ofstream*                           passedOfs);
           void                     computeMeanStacc            (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 std::ofstream*                           passedOfs);
#endif
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
           void                     computeAutoCorrViaDef       (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 const std::vector<unsigned int>&         lagsForCorrs,
                                                                 std::ofstream*                           passedOfs);
           void                     computeAutoCorrViaFFT       (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 const std::vector<unsigned int>&         initialPosForStatistics,
                                                                 const std::vector<unsigned int>&         lagsForCorrs,
                                                                 std::ofstream*                           passedOfs);
           void                     computeCovCorrMatrices      (const uqSequenceStatisticalOptionsClass& statisticalOptions,
                                                                 std::ofstream*                           passedOfs);
#endif
  virtual  void                     extractRawData              (unsigned int                             initialPos,
                                                                 unsigned int                             spacing,
                                                                 unsigned int                             numPos,
                                                                 unsigned int                             paramId,
                                                                 std::vector<double>&                     rawData) const = 0;

  const uqBaseEnvironmentClass&  m_env;
  const uqVectorSpaceClass<V,M>& m_vectorSpace;
  std::string                    m_name;

  mutable uqFftClass<double>*    m_fftObj;
  mutable V*                     m_subMinPlain;
  mutable V*                     m_unifiedMinPlain;
  mutable V*                     m_subMaxPlain;
  mutable V*                     m_unifiedMaxPlain;
  mutable V*                     m_subMeanPlain;
  mutable V*                     m_unifiedMeanPlain;
  mutable V*                     m_subMedianPlain;
  mutable V*                     m_unifiedMedianPlain;
  mutable V*                     m_subSampleVariancePlain;
  mutable V*                     m_unifiedSampleVariancePlain;
  mutable uqBoxSubsetClass<V,M>* m_subBoxPlain;
  mutable uqBoxSubsetClass<V,M>* m_unifiedBoxPlain;
};

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::uqBaseVectorSequenceClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   subSequenceSize,
  const std::string&             name)
  :
  m_env                       (vectorSpace.env()),
  m_vectorSpace               (vectorSpace),
  m_name                      (name),
  m_fftObj                    (new uqFftClass<double>(m_env)),
  m_subMinPlain               (NULL),
  m_unifiedMinPlain           (NULL),
  m_subMaxPlain               (NULL),
  m_unifiedMaxPlain           (NULL),
  m_subMeanPlain              (NULL),
  m_unifiedMeanPlain          (NULL),
  m_subMedianPlain            (NULL),
  m_unifiedMedianPlain        (NULL),
  m_subSampleVariancePlain    (NULL),
  m_unifiedSampleVariancePlain(NULL),
  m_subBoxPlain               (NULL),
  m_unifiedBoxPlain           (NULL)
{
  if (subSequenceSize) {}; // just to avoid compiler warning
}

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::~uqBaseVectorSequenceClass()
{
  //clear();
  this->deleteStoredVectors();
  if (m_fftObj != NULL) delete m_fftObj;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::copy(const uqBaseVectorSequenceClass<V,M>& src)
{
  // FIX ME: should check environments as well ???

  UQ_FATAL_TEST_MACRO(m_vectorSpace.dimLocal() != src.m_vectorSpace.dimLocal(),
                      m_env.worldRank(),
                      "uqSequenceOfVectorsClass<V,M>::copy()",
                      "incompatible vector space dimensions");

  m_name = src.m_name;
  this->deleteStoredVectors();

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
      m_env.inter0Comm().Allreduce((void *) &subNumSamples, (void *) &unifiedNumSamples, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                   "uqBaseVectorSequenceClass<V,M>::unifiedSequenceSize()",
                                   "failed MPI.Allreduce() for unifiedSequenceSize()");
    }
    else {
      // Node not in the 'inter0' communicator
      unifiedNumSamples = this->subSequenceSize();
    }
  }
  else {
    UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                        m_env.worldRank(),
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
uqBaseVectorSequenceClass<V,M>::subMinPlain() const
{
  if (m_subMinPlain == NULL) {
    m_subMinPlain = m_vectorSpace.newVector();
    if (m_subMaxPlain == NULL) m_subMaxPlain = m_vectorSpace.newVector();
    subMinMaxExtra(0,this->subSequenceSize(),*m_subMinPlain,*m_subMaxPlain);
  }

  return *m_subMinPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMinPlain() const
{
  if (m_unifiedMinPlain == NULL) {
    m_unifiedMinPlain = m_vectorSpace.newVector();
    if (m_unifiedMaxPlain == NULL) m_unifiedMaxPlain = m_vectorSpace.newVector();
    unifiedMinMaxExtra(0,this->subSequenceSize(),*m_unifiedMinPlain,*m_unifiedMaxPlain);
  }

  return *m_unifiedMinPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subMaxPlain() const
{
  if (m_subMaxPlain == NULL) {
    if (m_subMinPlain == NULL) m_subMinPlain = m_vectorSpace.newVector();
    m_subMaxPlain = m_vectorSpace.newVector();
    subMinMaxExtra(0,this->subSequenceSize(),*m_subMinPlain,*m_subMaxPlain);
  }

  return *m_subMaxPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMaxPlain() const
{
  if (m_unifiedMaxPlain == NULL) {
    if (m_unifiedMinPlain == NULL) m_unifiedMinPlain = m_vectorSpace.newVector();
    m_unifiedMaxPlain = m_vectorSpace.newVector();
    unifiedMinMaxExtra(0,this->subSequenceSize(),*m_unifiedMinPlain,*m_unifiedMaxPlain);
  }

  return *m_unifiedMaxPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subMeanPlain() const
{
  if (m_subMeanPlain == NULL) {
    m_subMeanPlain = m_vectorSpace.newVector();
    subMeanExtra(0,subSequenceSize(),*m_subMeanPlain);
  }

  return *m_subMeanPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMeanPlain() const
{
  if (m_unifiedMeanPlain == NULL) {
    m_unifiedMeanPlain = m_vectorSpace.newVector();
    unifiedMeanExtra(0,subSequenceSize(),*m_unifiedMeanPlain);
  }

  return *m_unifiedMeanPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subMedianPlain() const
{
  if (m_subMedianPlain == NULL) {
    m_subMedianPlain = m_vectorSpace.newVector();
    *m_subMedianPlain = subMedianExtra(0,subSequenceSize());
  }

  return *m_subMedianPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedMedianPlain() const
{
  if (m_unifiedMedianPlain == NULL) {
    m_unifiedMedianPlain = m_vectorSpace.newVector();
    *m_unifiedMedianPlain = unifiedMedianExtra(0,subSequenceSize());
  }

  return *m_unifiedMedianPlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::subSampleVariancePlain() const
{
  if (m_subSampleVariancePlain == NULL) {
    m_subSampleVariancePlain = m_vectorSpace.newVector();
    subSampleVarianceExtra(0,subSequenceSize(),subMeanPlain(),*m_subSampleVariancePlain);
  }

  return *m_subSampleVariancePlain;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::unifiedSampleVariancePlain() const
{
  if (m_unifiedSampleVariancePlain == NULL) {
    m_unifiedSampleVariancePlain = m_vectorSpace.newVector();
    unifiedSampleVarianceExtra(0,subSequenceSize(),unifiedMeanPlain(),*m_unifiedSampleVariancePlain);
  }

  return *m_unifiedSampleVariancePlain;
}

template <class V, class M>
const uqBoxSubsetClass<V,M>&
uqBaseVectorSequenceClass<V,M>::subBoxPlain() const
{
  if (m_subBoxPlain == NULL) {
    m_subBoxPlain = new uqBoxSubsetClass<V,M>(m_name.c_str(),
                                              m_vectorSpace,
                                              this->subMinPlain(),
                                              this->subMaxPlain());
  }

  return *m_subBoxPlain;
}

template <class V, class M>
const uqBoxSubsetClass<V,M>&
uqBaseVectorSequenceClass<V,M>::unifiedBoxPlain() const
{
  if (m_unifiedBoxPlain == NULL) {
    m_unifiedBoxPlain = new uqBoxSubsetClass<V,M>(m_name.c_str(),
                                                  m_vectorSpace,
                                                  this->unifiedMinPlain(),
                                                  this->unifiedMaxPlain());
  }

  return *m_unifiedBoxPlain;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::deleteStoredVectors()
{
  if (m_subMinPlain) {
    delete m_subMinPlain;
    m_subMinPlain                = NULL;
  }
  if (m_unifiedMinPlain) {
    delete m_unifiedMinPlain;
    m_unifiedMinPlain            = NULL;
  }
  if (m_subMaxPlain) {
    delete m_subMaxPlain;
    m_subMaxPlain                = NULL;
  }
  if (m_unifiedMaxPlain) {
    delete m_unifiedMaxPlain;
    m_unifiedMaxPlain            = NULL;
  }
  if (m_subMeanPlain) {
    delete m_subMeanPlain;
    m_subMeanPlain               = NULL;
  }
  if (m_unifiedMeanPlain) {
    delete m_unifiedMeanPlain;
    m_unifiedMeanPlain           = NULL;
  }
  if (m_subMedianPlain) {
    delete m_subMedianPlain;
    m_subMedianPlain             = NULL;
  }
  if (m_unifiedMedianPlain) {
    delete m_unifiedMedianPlain;
    m_unifiedMedianPlain         = NULL;
  }
  if (m_subSampleVariancePlain) {
    delete m_subSampleVariancePlain;
    m_subSampleVariancePlain     = NULL;
  }
  if (m_unifiedSampleVariancePlain) {
    delete m_unifiedSampleVariancePlain;
    m_unifiedSampleVariancePlain = NULL;
  }
  if (m_subBoxPlain) {
    delete m_subBoxPlain;
    m_subBoxPlain     = NULL;
  }
  if (m_unifiedBoxPlain) {
    delete m_unifiedBoxPlain;
    m_unifiedBoxPlain = NULL;
  }

  return;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::append(
  const uqBaseVectorSequenceClass<V,M>& src,
  unsigned int                          initialPos,
  unsigned int                          numPos)
{
  UQ_FATAL_TEST_MACRO((src.subSequenceSize() < (initialPos+1)),
                      m_env.worldRank(),
                      "uqBaseVectorSequence<V,M>::append()",
                      "initialPos is too big");

  UQ_FATAL_TEST_MACRO((src.subSequenceSize() < (initialPos+numPos)),
                      m_env.worldRank(),
                      "uqBaseVectorSequence<V,M>::append()",
                      "numPos is too big");

  this->deleteStoredVectors();
  unsigned int currentSize = this->subSequenceSize();
  this->resizeSequence(currentSize+numPos);
  V tmpVec(src.vectorSpace().zeroVector());
  for (unsigned int i = 0; i < numPos; ++i) {
    src.getPositionValues(initialPos+i,tmpVec);
    this->setPositionValues(currentSize+i,tmpVec);
  }

  return;
}

template <class V, class M>
double
uqBaseVectorSequenceClass<V,M>::subPositionsOfMaximum(
  const uqScalarSequenceClass<double>& subCorrespondingScalarValues,
  uqBaseVectorSequenceClass<V,M>&      subPositionsOfMaximum)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqBaseVectorSequenceClass<V,M>::subPositionsOfMaximum()"
                            << ": subCorrespondingScalarValues,subSequenceSize() = " << subCorrespondingScalarValues.subSequenceSize()
                            << ", this->subSequenceSize = " << this->subSequenceSize()
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(subCorrespondingScalarValues.subSequenceSize() != this->subSequenceSize(),
                      m_env.worldRank(),
                      "uqBaseVectorSequenceClass<V,M>::subPositionsOfMaximum()",
                      "invalid input");

  double subMaxValue = subCorrespondingScalarValues.subMaxPlain();
  unsigned int iMax = subCorrespondingScalarValues.subSequenceSize();

  unsigned int subNumPos = 0;
  for (unsigned int i = 0; i < iMax; ++i) {
    if (subCorrespondingScalarValues[i] == subMaxValue) {
      subNumPos++;
    }
  }

  V tmpVec(this->vectorSpace().zeroVector());
  subPositionsOfMaximum.resizeSequence(subNumPos);
  unsigned int j = 0;
  for (unsigned int i = 0; i < iMax; ++i) {
    if (subCorrespondingScalarValues[i] == subMaxValue) {
      this->getPositionValues                (i,tmpVec);
      subPositionsOfMaximum.setPositionValues(j,tmpVec);
      j++;
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqBaseVectorSequenceClass<V,M>::subPositionsOfMaximum()"
                            << std::endl;
  }

  return subMaxValue;
}

template <class V, class M>
double
uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum( // rr0
  const uqScalarSequenceClass<double>& subCorrespondingScalarValues,
  uqBaseVectorSequenceClass<V,M>&      unifiedPositionsOfMaximum)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()"
                            << ": subCorrespondingScalarValues,subSequenceSize() = " << subCorrespondingScalarValues.subSequenceSize()
                            << ", this->subSequenceSize = " << this->subSequenceSize()
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(subCorrespondingScalarValues.subSequenceSize() != this->subSequenceSize(),
                      m_env.worldRank(),
                      "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                      "invalid input");

  double subMaxValue = subCorrespondingScalarValues.subMaxPlain();

  //******************************************************************
  // Get overall max
  //******************************************************************
  double unifiedMaxValue;
  std::vector<double> sendbufPos(1,0.);
  for (unsigned int i = 0; i < sendbufPos.size(); ++i) {
    sendbufPos[i] = subMaxValue;
  }
  m_env.inter0Comm().Allreduce((void *) &sendbufPos[0], (void *) &unifiedMaxValue, (int) sendbufPos.size(), uqRawValue_MPI_DOUBLE, uqRawValue_MPI_MAX,
                               "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                               "failed MPI.Allreduce() for max");

  unsigned int iMax = subCorrespondingScalarValues.subSequenceSize();
  int subNumPos = 0; // Yes, 'int', due to MPI to be used soon
  for (unsigned int i = 0; i < iMax; ++i) {
    if (subCorrespondingScalarValues[i] == unifiedMaxValue) {
      subNumPos++;
    }
  }

  V tmpVec(this->vectorSpace().zeroVector());
  unifiedPositionsOfMaximum.resizeSequence(subNumPos);
  unsigned int j = 0;
  for (unsigned int i = 0; i < iMax; ++i) {
    if (subCorrespondingScalarValues[i] == unifiedMaxValue) {
      this->getPositionValues                    (i,tmpVec);
      unifiedPositionsOfMaximum.setPositionValues(j,tmpVec);
      j++;
    }
  }

  std::vector<int> auxBuf(1,0);
  int unifiedNumPos = 0; // Yes, 'int', due to MPI to be used soon
  auxBuf[0] = subNumPos;
  m_env.inter0Comm().Allreduce((void *) &auxBuf[0], (void *) &unifiedNumPos, (int) auxBuf.size(), uqRawValue_MPI_INT, uqRawValue_MPI_SUM,
                               "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                               "failed MPI.Allreduce() for sum");
  unifiedPositionsOfMaximum.resizeSequence(unifiedNumPos);

  //******************************************************************
  // Use MPI_Gatherv for number of positions
  //******************************************************************
  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();

  std::vector<int> recvcntsPos(Np,0); // '0' is NOT the correct value for recvcntsPos[0]
  m_env.inter0Comm().Gather((void *) &subNumPos, 1, uqRawValue_MPI_INT, (void *) &recvcntsPos[0], (int) 1, uqRawValue_MPI_INT, 0,
                            "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                            "failed MPI.Gatherv()");
  if (m_env.inter0Rank() == 0) {
    UQ_FATAL_TEST_MACRO(recvcntsPos[0] != (int) subNumPos,
                        m_env.worldRank(),
                        "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                        "failed MPI.Gather() result at proc 0 (recvcntsPos[0])");
  }

  std::vector<int> displsPos(Np,0);
  for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, from '1' on
    displsPos[nodeId] = displsPos[nodeId-1] + recvcntsPos[nodeId-1];
  }
  if (m_env.inter0Rank() == 0) {
    UQ_FATAL_TEST_MACRO(unifiedNumPos != (displsPos[Np - 1] + recvcntsPos[Np - 1]),
                        m_env.worldRank(),
                        "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                        "failed MPI.Gather() result at proc 0 (unifiedNumPos)");
  }

  //******************************************************************
  // Use MPI_Gatherv for number of doubles
  //******************************************************************
  unsigned int dimSize = m_vectorSpace.dimLocal();
  int subNumDbs = subNumPos * dimSize; // Yes, 'int', due to MPI to be used soon
  std::vector<int> recvcntsDbs(Np,0); // '0' is NOT the correct value for recvcntsDbs[0]
  m_env.inter0Comm().Gather((void *) &subNumDbs, 1, uqRawValue_MPI_INT, (void *) &recvcntsDbs[0], (int) 1, uqRawValue_MPI_INT, 0,
                            "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                            "failed MPI.Gatherv()");
  if (m_env.inter0Rank() == 0) {
    UQ_FATAL_TEST_MACRO(recvcntsDbs[0] != (int) subNumDbs,
                        m_env.worldRank(),
                        "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                        "failed MPI.Gather() result at proc 0 (recvcntsDbs[0])");
  }

  std::vector<int> displsDbs(Np,0);
  for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, from '1' on
    displsDbs[nodeId] = displsDbs[nodeId-1] + recvcntsDbs[nodeId-1];
  }
  if (m_env.inter0Rank() == 0) {
    UQ_FATAL_TEST_MACRO(((int) (unifiedNumPos*dimSize)) != (displsDbs[Np - 1] + recvcntsDbs[Np - 1]),
                        m_env.worldRank(),
                        "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                        "failed MPI.Gather() result at proc 0 (unifiedNumPos*dimSize)");
  }

  //******************************************************************
  // Prepare counters and buffers for gatherv of maximum positions
  //******************************************************************
  std::vector<double> sendbufDbs(subNumDbs,0.);
  for (unsigned int i = 0; i < (unsigned int) subNumPos; ++i) {
    unifiedPositionsOfMaximum.getPositionValues(i,tmpVec);
    for (unsigned int j = 0; j < dimSize; ++j) {
      sendbufDbs[i*dimSize + j] = tmpVec[j];
    }
  }

  std::vector<double> recvbufDbs(unifiedNumPos * dimSize);

  m_env.inter0Comm().Gatherv((void *) &sendbufDbs[0], (int) subNumDbs, uqRawValue_MPI_DOUBLE, (void *) &recvbufDbs[0], (int *) &recvcntsDbs[0], (int *) &displsDbs[0], uqRawValue_MPI_DOUBLE, 0,
                             "uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()",
                             "failed MPI.Gatherv()");

  //******************************************************************
  // Transfer data from 'recvbuf' to 'unifiedPositionsOfMaximum'
  //******************************************************************
  if (m_env.inter0Rank() == (int) 0) {
    for (unsigned int i = 0; i < (unsigned int) unifiedNumPos; ++i) {
      for (unsigned int j = 0; j < dimSize; ++j) {
        tmpVec[j] = recvbufDbs[i*dimSize + j];
      }
      unifiedPositionsOfMaximum.setPositionValues(i,tmpVec);
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqBaseVectorSequenceClass<V,M>::unifiedPositionsOfMaximum()"
                            << std::endl;
  }

  return unifiedMaxValue;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::setGaussian(const V& meanVec, const V& stdDevVec)
{
  V gaussianVector(m_vectorSpace.zeroVector());
  for (unsigned int j = 0; j < this->subSequenceSize(); ++j) {
    gaussianVector.cwSetGaussian(meanVec,stdDevVec);
    this->setPositionValues(j,gaussianVector);
  }

  this->deleteStoredVectors();

  return;
}


template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::setUniform(const V& aVec, const V& bVec)
{
  V uniformVector(m_vectorSpace.zeroVector());
  for (unsigned int j = 0; j < this->subSequenceSize(); ++j) {
    uniformVector.cwSetUniform(aVec,bVec);
    this->setPositionValues(j,uniformVector);
  }

  this->deleteStoredVectors();

  return;
}

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeStatistics(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                           passedOfs)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing statistics for chain '" << m_name << "' ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  bool okSituation = ((passedOfs == NULL                            ) ||
                      ((passedOfs != NULL) && (m_env.subRank() >= 0)));
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.worldRank(),
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
  // Compute mean, median, sample std, population std
  //****************************************************
  this->computeMeanVars(statisticalOptions,
                        passedOfs,
                        NULL,
                        NULL,
                        NULL,
                        NULL);

#ifdef UQ_CODE_HAS_MONITORS
  if (statisticalOptions.meanMonitorPeriod() != 0) {
    this->computeMeanEvolution(statisticalOptions,
                               passedOfs);
  }
#endif

  //****************************************************
  // Compute variance of sample mean through the 'batch means method' (BMM)
  //****************************************************
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.bmmRun()               ) &&
      (initialPosForStatistics.size()         > 0) &&
      (statisticalOptions.bmmLengths().size() > 0)) { 
    this->computeBMM(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }
#endif
  //****************************************************
  // Compute FFT of chain, for one parameter only
  //****************************************************
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.fftCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    this->computeFFT(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }
#endif
  //****************************************************
  // Compute power spectral density (PSD) of chain, for one parameter only
  //****************************************************
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.psdCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    this->computePSD(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }
#endif
  //****************************************************
  // Compute power spectral density (PSD) of chain at zero frequency
  //****************************************************
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.psdAtZeroCompute()             ) &&
      (initialPosForStatistics.size()                 > 0) &&
      (statisticalOptions.psdAtZeroNumBlocks().size() > 0)) { 
    this->computePSDAtZero(statisticalOptions,
                           initialPosForStatistics,
                           passedOfs);
  }
#endif
  //****************************************************
  // Compute Geweke
  //****************************************************
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.gewekeCompute()) &&
      (initialPosForStatistics.size() > 0)) {
    this->computeGeweke(statisticalOptions,
                        initialPosForStatistics,
                        passedOfs);
  }
#endif
  //****************************************************
  // Compute mean statistical accuracy
  //****************************************************
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.meanStaccCompute()) &&
      (initialPosForStatistics.size() > 0   )) {
    this->computeMeanStacc(statisticalOptions,
                           initialPosForStatistics,
                           passedOfs);
  }
#endif
  // Set lags for the computation of chain autocorrelations
  std::vector<unsigned int> lagsForCorrs(statisticalOptions.autoCorrNumLags(),1);
  for (unsigned int i = 1; i < lagsForCorrs.size(); ++i) {
    lagsForCorrs[i] = statisticalOptions.autoCorrSecondLag() + (i-1)*statisticalOptions.autoCorrLagSpacing();
  }

  //****************************************************
  // Compute autocorrelation coefficients via definition
  //****************************************************
  if ((statisticalOptions.autoCorrComputeViaDef()) &&
      (initialPosForStatistics.size() > 0        ) &&
      (lagsForCorrs.size()            > 0        )) { 
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
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((statisticalOptions.histCompute    ()) ||
      (statisticalOptions.cdfStaccCompute()) ||
      (statisticalOptions.kdeCompute     ())) {
#else
    if (statisticalOptions.kdeCompute()) {
#endif
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
    *m_env.subDisplayFile() << "All statistics of chain '" << m_name << "'"
                            << " took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing statistics for chain '" << m_name << "'"
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}
#endif // #ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeMeanVars(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                           passedOfs,
  V*                                       subMeanPtr,
  V*                                       subMedianPtr,
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
  this->subMeanExtra(0,
                     this->subSequenceSize(),
                     subChainMean);

  V subChainMedian(m_vectorSpace.zeroVector());
  this->subMedianExtra(0,
                     this->subSequenceSize(),
                     subChainMedian);

  V subChainSampleVariance(m_vectorSpace.zeroVector());
  this->subSampleVarianceExtra(0,
                               this->subSequenceSize(),
                               subChainMean,
                               subChainSampleVariance);

  if ((m_env.displayVerbosity() >= 5) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeMeanVars()"
                            << ": subChainMean.sizeLocal() = "           << subChainMean.sizeLocal()
                            << ", subChainMean = "                       << subChainMean
                            << ", subChainMedian = "                     << subChainMedian
                            << ", subChainSampleVariance.sizeLocal() = " << subChainSampleVariance.sizeLocal()
                            << ", subChainSampleVariance = "             << subChainSampleVariance
                            << std::endl;
  }

  V estimatedVarianceOfSampleMean(subChainSampleVariance);
  estimatedVarianceOfSampleMean /= (double) this->subSequenceSize();
  bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
  estimatedVarianceOfSampleMean.setPrintHorizontally(false);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nEstimated variance of sample mean for the whole chain '" << m_name << "'"
                            << ", under independence assumption:"
                            << "\n"
                            << estimatedVarianceOfSampleMean
                            << std::endl;
  }
  estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);

  V estimatedStdOfSampleMean(estimatedVarianceOfSampleMean);
  estimatedStdOfSampleMean.cwSqrt();
  savedVectorPrintState = estimatedStdOfSampleMean.getPrintHorizontally();
  estimatedStdOfSampleMean.setPrintHorizontally(false);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nEstimated standard deviation of sample mean for the whole chain '" << m_name << "'"
                            << ", under independence assumption:"
                            << "\n"
                            << estimatedStdOfSampleMean
                            << std::endl;
  }
  estimatedStdOfSampleMean.setPrintHorizontally(savedVectorPrintState);

  V subChainPopulationVariance(m_vectorSpace.zeroVector());
  this->subPopulationVariance(0,
                              this->subSequenceSize(),
                              subChainMean,
                              subChainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Sub Mean, median, and variances took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nSub mean, median, sample std, population std"
                            << std::endl;
    char line[512];
    sprintf(line,"%s%4s%s%6s%s%9s%s%9s%s",
	    "Parameter",
            " ",
            "Mean",
            " ",
            "Median",
            " ",
            "SampleStd",
            " ",
            "Popul.Std");
    *m_env.subDisplayFile() << line;

    for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e%2s%11.4e",
              m_vectorSpace.localComponentName(i).c_str(), /*.*/
              " ",
	      subChainMean[i],
              " ",
	      subChainMedian[i],
              " ",
              std::sqrt(subChainSampleVariance[i]),
              " ",
              std::sqrt(subChainPopulationVariance[i]));
      *m_env.subDisplayFile() << line;
    }
    *m_env.subDisplayFile() << std::endl;
  }

  if (subMeanPtr     ) *subMeanPtr      = subChainMean;
  if (subMedianPtr   ) *subMedianPtr    = subChainMedian;
  if (subSampleVarPtr) *subSampleVarPtr = subChainSampleVariance;
  if (subPopulVarPtr ) *subPopulVarPtr  = subChainPopulationVariance;

  if (m_env.numSubEnvironments() > 1) {
    // Write unified min-max
    if (m_vectorSpace.numOfProcsForStorage() == 1) {
      V unifiedChainMean(m_vectorSpace.zeroVector());
      this->unifiedMeanExtra(0,
                             this->subSequenceSize(),
                             unifiedChainMean);

      V unifiedChainMedian(m_vectorSpace.zeroVector());
      this->unifiedMedianExtra(0,
                               this->subSequenceSize(),
                               unifiedChainMedian);

      V unifiedChainSampleVariance(m_vectorSpace.zeroVector());
      this->unifiedSampleVarianceExtra(0,
                                       this->subSequenceSize(),
                                       unifiedChainMean,
                                       unifiedChainSampleVariance);

      V unifiedChainPopulationVariance(m_vectorSpace.zeroVector());
      this->unifiedPopulationVariance(0,
                                      this->subSequenceSize(),
                                      unifiedChainMean,
                                      unifiedChainPopulationVariance);

      if (m_env.inter0Rank() == 0) {
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "\nUnif mean, median, sample std, population std"
                                  << std::endl;
          char line[512];
          sprintf(line,"%s%4s%s%6s%s%9s%s%9s%s",
	          "Parameter",
                  " ",
                  "Mean",
                  " ",
                  "Median",
                  " ",
                  "SampleStd",
                  " ",
                  "Popul.Std");
          *m_env.subDisplayFile() << line;

          for (unsigned int i = 0; i < this->vectorSizeLocal() /*.*/; ++i) {
            sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e%2s%11.4e",
                    m_vectorSpace.localComponentName(i).c_str(), /*.*/
                    " ",
	            unifiedChainMean[i],
                    " ",
	            unifiedChainMedian[i],
                    " ",
                    std::sqrt(unifiedChainSampleVariance[i]),
                    " ",
                    std::sqrt(unifiedChainPopulationVariance[i]));
            *m_env.subDisplayFile() << line;
          }
          *m_env.subDisplayFile() << std::endl;
        } // if subDisplayFile
      }
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqBaseVectorSequenceClass<V,M>::computeMeanVars()",
                          "unified min-max writing, parallel vectors not supported yet");
    }
  } // if numSubEnvs > 1

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nEnded computing mean, sample variance and population variance"
                            << "\n-----------------------------------------------------"
                            << std::endl;
  }

  return;
}
#endif // #ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#ifdef UQ_CODE_HAS_MONITORS
template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeMeanEvolution(
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                           passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing mean evolution"
                            << std::endl;
  }

  unsigned int monitorPeriod = statisticalOptions.meanMonitorPeriod();
  unsigned int iMin = 0;
  unsigned int iMax = (unsigned int) ( ((double) this->subSequenceSize())/((double) monitorPeriod) );

  if (monitorPeriod == 1) {
    iMin = 1;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "  Sub sequence size = "                << this->subSequenceSize()
                            << "\n  Monitor period = "                 << monitorPeriod
                            << "\n  Number of monitoring positions = " << iMax
                            << std::endl;
  }

  this->subMeanMonitorAlloc(iMax);
  if (m_env.numSubEnvironments() > 1) {
    this->subMeanInter0MonitorAlloc(iMax);
    this->unifiedMeanMonitorAlloc(iMax);
  }

  for (unsigned int i = iMin; i < iMax; ++i) {
    unsigned int currentMonitoredFinalPosition = (i+1)*monitorPeriod - 1; // Yes, '-1'
    V subMeanVec   (m_vectorSpace.zeroVector());
    V subMeanCltStd(m_vectorSpace.zeroVector());
    this->subMeanMonitorRun(currentMonitoredFinalPosition,
                            subMeanVec,
                            subMeanCltStd);
    this->subMeanMonitorStore(i,
                              currentMonitoredFinalPosition,
                              subMeanVec,
                              subMeanCltStd);

    if (m_env.numSubEnvironments() > 1) {
      V subMeanInter0Mean       (m_vectorSpace.zeroVector());
      V subMeanInter0Clt95      (m_vectorSpace.zeroVector());
      V subMeanInter0Empirical90(m_vectorSpace.zeroVector());
      V subMeanInter0Min        (m_vectorSpace.zeroVector());
      V subMeanInter0Max        (m_vectorSpace.zeroVector());
      this->subMeanInter0MonitorRun(currentMonitoredFinalPosition,
                                    subMeanInter0Mean,
                                    subMeanInter0Clt95,
                                    subMeanInter0Empirical90,
                                    subMeanInter0Min,
                                    subMeanInter0Max);
      this->subMeanInter0MonitorStore(i,
                                      currentMonitoredFinalPosition,
                                      subMeanInter0Mean,
                                      subMeanInter0Clt95,
                                      subMeanInter0Empirical90,
                                      subMeanInter0Min,
                                      subMeanInter0Max);

      V unifiedMeanVec   (m_vectorSpace.zeroVector());
      V unifiedMeanCltStd(m_vectorSpace.zeroVector());
      this->unifiedMeanMonitorRun(currentMonitoredFinalPosition,
                                  unifiedMeanVec,
                                  unifiedMeanCltStd);
      this->unifiedMeanMonitorStore(i,
                                    currentMonitoredFinalPosition,
                                    unifiedMeanVec,
                                    unifiedMeanCltStd);
    }
  }

  if (passedOfs) {
    this->subMeanMonitorWrite(*passedOfs);
    if (m_env.numSubEnvironments() > 1) {
      this->subMeanInter0MonitorWrite(*passedOfs);
      this->unifiedMeanMonitorWrite(*passedOfs); // Yes, call 'unified' even though not all nodes might have 'passedOfs != NULL' ????????????? ernesto
    }
  }

  this->subMeanMonitorFree();
  if (m_env.numSubEnvironments() > 1) {
    this->subMeanInter0MonitorFree();
    this->unifiedMeanMonitorFree();
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Mean evolution took " << tmpRunTime
                            << " seconds"
                            << std::endl;
  }

  return;
}
#endif

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
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
#endif // #ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
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
    V subChainMean                 (m_vectorSpace.zeroVector());
    V subChainSampleVariance       (m_vectorSpace.zeroVector());
    V estimatedVarianceOfSampleMean(m_vectorSpace.zeroVector());
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];

      this->subMeanExtra(initialPos,
                         this->subSequenceSize()-initialPos,
                         subChainMean);

      this->subSampleVarianceExtra(initialPos,
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
#endif // #ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeFilterParams(
  std::ofstream* passedOfs,
  unsigned int&  initialPos,
  unsigned int&  spacing)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing filter parameters for chain '" << m_name << "' ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  bool okSituation = ((passedOfs == NULL                            ) ||
                      ((passedOfs != NULL) && (m_env.subRank() >= 0)));
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.worldRank(),
                      "uqBaseVectorSequenceClass<V,M>::computeFilterParams()",
                      "unexpected combination of file pointer and subRank");

  //initialPos = 0;
  spacing    = 1;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing filter parameters for chain '" << m_name << "'"
                            << ": initialPos = " << initialPos
                            << ", spacing = "    << spacing
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde( // Use the whole chain
  const uqSequenceStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                           passedOfs)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing histogram and/or cdf stacc and/or Kde for chain '" << m_name << "' ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;

  //****************************************************
  // Compute MIN and MAX: for histograms and Kde
  //****************************************************
  double tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\nComputing min and max for histograms and Kde"
                            << std::endl;
  }

  V statsMinPositions(m_vectorSpace.zeroVector());
  V statsMaxPositions(m_vectorSpace.zeroVector());
  this->subMinMaxExtra(0, // Use the whole chain
                       this->subSequenceSize(),
                       statsMinPositions,
                       statsMaxPositions);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\nComputed min values and max values for chain '" << m_name << "'"
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
    this->unifiedMinMaxExtra(0, // Use the whole chain
                             this->subSequenceSize(),
                             unifiedStatsMinPositions,
                             unifiedStatsMaxPositions);

    // Write unified min-max
    if (m_env.subDisplayFile()) {
      if (m_vectorSpace.numOfProcsForStorage() == 1) {
        if (m_env.inter0Rank() == 0) {
          *m_env.subDisplayFile() << "\nComputed unified min values and max values for chain '" << m_name << "'"
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
                            m_env.worldRank(),
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
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
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
                              m_env.worldRank(),
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
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
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
#endif
  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((statisticalOptions.kdeCompute()             ) &&
      (statisticalOptions.kdeNumEvalPositions() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                              << "\nComputing Kde"
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
    this->subScalesForKde(0, // Use the whole chain
                          iqrVec,
                          1,
                          gaussianKdeScaleVec);

    std::vector<V*> kdeEvalPositions(statisticalOptions.kdeNumEvalPositions(),NULL);
    uqMiscComputePositionsBetweenMinMax(statsMinPositions,
                                        statsMaxPositions,
                                        kdeEvalPositions);

    std::vector<V*> gaussianKdeDensities(statisticalOptions.kdeNumEvalPositions(),NULL);
    this->subGaussian1dKde(0, // Use the whole chain
                           gaussianKdeScaleVec,
                           kdeEvalPositions,
                           gaussianKdeDensities);

    // Write iqr
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "\nComputed inter quantile ranges for chain '" << m_name << "'"
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

    // Write Kde
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() > 10)) { // output debug
      *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()"
                              << ", for chain '" << m_name << "'"
                              << ", about to write sub kde to ofstream with pointer = " << passedOfs
                              << std::endl;
    }
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
      // Compute unified Kde
      V unifiedIqrVec(m_vectorSpace.zeroVector());
      this->unifiedInterQuantileRange(0, // Use the whole chain
                                      unifiedIqrVec);
      //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

      V unifiedGaussianKdeScaleVec(m_vectorSpace.zeroVector());
      this->unifiedScalesForKde(0, // Use the whole chain
                                unifiedIqrVec,
                                1,
                                unifiedGaussianKdeScaleVec);
      //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

      std::vector<V*> unifiedKdeEvalPositions(statisticalOptions.kdeNumEvalPositions(),NULL);
      uqMiscComputePositionsBetweenMinMax(unifiedStatsMinPositions,
                                          unifiedStatsMaxPositions,
                                          unifiedKdeEvalPositions);

      std::vector<V*> unifiedGaussianKdeDensities(statisticalOptions.kdeNumEvalPositions(),NULL);
      this->unifiedGaussian1dKde(0, // Use the whole chain
                                 unifiedGaussianKdeScaleVec,
                                 unifiedKdeEvalPositions,
                                 unifiedGaussianKdeDensities);
      //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

      // Write unified iqr
      if (m_env.subDisplayFile()) {
        if (m_vectorSpace.numOfProcsForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            *m_env.subDisplayFile() << "\nComputed unified inter quantile ranges for chain '" << m_name << "'"
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
                              m_env.worldRank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()",
                              "unified iqr writing, parallel vectors not supported yet");
        }
      }

      // Write unified Kde
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() > 10)) { // output debug
        *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()"
                                << ", for chain '" << m_name << "'"
                                << ", about to write unified kde to ofstream with pointer = " << passedOfs
                                << std::endl;
      }
      if (passedOfs) {
        if (m_vectorSpace.numOfProcsForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            std::ofstream& ofsvar = *passedOfs;
            ofsvar << m_name << uniCoreName_GaussianKdePositions << " = zeros(" << this->vectorSizeLocal() /*.*/
                   << ","                                                       << unifiedKdeEvalPositions.size()
                   << ");"
                   << std::endl;
            if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() > 10)) { // output debug
              *m_env.subDisplayFile() << "In uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()"
                                      << ", for chain '" << m_name << "'"
                                      << ": just wrote '... = zeros(.,.);' line to output file, which has pointer = " << passedOfs
                                      << std::endl;
            }
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
                              m_env.worldRank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistCdfstaccKde()",
                              "unified Kde writing, parallel vectors not supported yet");
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
      *m_env.subDisplayFile() << "Chain Kde took " << tmpRunTime
                              << " seconds"
                              << std::endl;
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing histogram and/or cdf stacc and/or Kde for chain '" << m_name << "'"
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
  std::ofstream*                           passedOfs)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Computing covariance and correlation matrices for chain '" << m_name << "' ..."
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
                          m_env.worldRank(),
                          "uqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices()",
                          "parallel vectors not supported yet");
    }
  }

  delete correlationMatrix;
  delete covarianceMatrix;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished computing covariance and correlation matrices for chain '" << m_name << "'"
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return;
}
#endif // #ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS

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
                      env.worldRank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "parallel vectors not supported yet");

  unsigned int numRowsLocal = subPSeq.vectorSpace().dimLocal();
  unsigned int numCols = subQSeq.vectorSpace().dimGlobal();

  UQ_FATAL_TEST_MACRO((numRowsLocal != pqCovMatrix.numRowsLocal()) || (numCols != pqCovMatrix.numCols()),
                      env.worldRank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "inconsistent dimensions for covariance matrix");

  UQ_FATAL_TEST_MACRO((numRowsLocal != pqCorrMatrix.numRowsLocal()) || (numCols != pqCorrMatrix.numCols()),
                      env.worldRank(),
                      "uqComputeCorrelationBetweenVectorSequences()",
                      "inconsistent dimensions for correlation matrix");

  UQ_FATAL_TEST_MACRO((subNumSamples > subPSeq.subSequenceSize()) || (subNumSamples > subQSeq.subSequenceSize()),
                      env.worldRank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "subNumSamples is too large");

  // For both P and Q vector sequences: fill them
  P_V tmpP(subPSeq.vectorSpace().zeroVector());
  Q_V tmpQ(subQSeq.vectorSpace().zeroVector());

  // For both P and Q vector sequences: compute the unified mean
  P_V unifiedMeanP(subPSeq.vectorSpace().zeroVector());
  subPSeq.unifiedMeanExtra(0,subNumSamples,unifiedMeanP);

  Q_V unifiedMeanQ(subQSeq.vectorSpace().zeroVector());
  subQSeq.unifiedMeanExtra(0,subNumSamples,unifiedMeanQ);

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
  subPSeq.unifiedSampleVarianceExtra(0,
                                     subNumSamples,
                                     unifiedMeanP,
                                     unifiedSampleVarianceP);

  Q_V unifiedSampleVarianceQ(subQSeq.vectorSpace().zeroVector());
  subQSeq.unifiedSampleVarianceExtra(0,
                                     subNumSamples,
                                     unifiedMeanQ,
                                     unifiedSampleVarianceQ);

  // Compute unified covariance matrix
  if (useOnlyInter0Comm) {
    if (env.inter0Rank() >= 0) {
      unsigned int unifiedNumSamples = 0;
      env.inter0Comm().Allreduce((void *) &subNumSamples, (void *) &unifiedNumSamples, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                 "uqComputeCovCorrMatricesBetweenVectorSequences()",
                                 "failed MPI.Allreduce() for subNumSamples");

      for (unsigned i = 0; i < numRowsLocal; ++i) {
        for (unsigned j = 0; j < numCols; ++j) {
          double aux = 0.;
          env.inter0Comm().Allreduce((void *) &pqCovMatrix(i,j), (void *) &aux, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                                     "uqComputeCovCorrMatricesBetweenVectorSequences()",
                                     "failed MPI.Allreduce() for a matrix position");
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
                        << ": worldRank = "            << env.worldRank()
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
                               env.worldRank(),
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
                        env.worldRank(),
                        "uqComputeCovCorrMatricesBetweenVectorSequences()",
                        "parallel vectors not supported yet (2)");
  }

  return;
}
#endif // __UQ_VECTOR_SEQUENCE_H__
