/* uq/libs/basic/inc/uqVectorSequence.h
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

#ifndef __UQ_VECTOR_SEQUENCE_H__
#define __UQ_VECTOR_SEQUENCE_H__

#include <uqVectorSpace.h>
#include <uqScalarSequence.h>
#include <uqChainStatisticalOptions.h>
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
                                     unsigned int                   sequenceSize,
                                     const std::string&             name);
  virtual ~uqBaseVectorSequenceClass();

  virtual  unsigned int             sequenceSize        () const = 0;
           unsigned int             vectorSize          () const;
  const    uqVectorSpaceClass<V,M>& vectorSpace         () const;
  const    std::string&             name                () const;
           void                     setName             (const std::string& newName);
           void                     clear               ();
  const    V&                       minValues           () const;
  const    V&                       maxValues           () const;
  const    V&                       meanValues          () const;
  const    V&                       sampleVarianceValues() const;
           void                     deleteStoredVectors ();

  virtual  void                     resizeSequence      (unsigned int newSequenceSize) = 0;
  virtual  void                     resetValues         (unsigned int initialPos, unsigned int numPos) = 0;
  virtual  void                     erasePositions      (unsigned int initialPos, unsigned int numPos) = 0;
  virtual  void                     getPositionValues   (unsigned int posId,       V& vec) const = 0;
  virtual  void                     setPositionValues   (unsigned int posId, const V& vec) = 0;
  virtual  void                     setGaussian         (const gsl_rng* rng, const V& meanVec, const V& stdDevVec) = 0;
  virtual  void                     setUniform          (const gsl_rng* rng, const V& aVec,    const V& bVec     ) = 0;
  virtual  void                     uniformlySampledMdf (const V&                       numEvaluationPointsVec,
                                                         uqArrayOfOneDGridsClass <V,M>& mdfGrids,
                                                         uqArrayOfOneDTablesClass<V,M>& mdfValues) const = 0;
  virtual  void                     uniformlySampledCdf (const V&                       numEvaluationPointsVec,
                                                         uqArrayOfOneDGridsClass <V,M>& cdfGrids,
                                                         uqArrayOfOneDTablesClass<V,M>& cdfValues) const = 0;

  virtual  void                     mean                (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         V&                                    meanVec) const = 0;
  virtual  void                     sampleVariance      (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         const V&                              meanVec,
                                                         V&                                    samVec) const = 0;
  virtual  void                     populationVariance  (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         const V&                              meanVec,
                                                         V&                                    popVec) const = 0;
  virtual  void                     autoCovariance      (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         const V&                              meanVec,
                                                         unsigned int                          lag,
                                                         V&                                    covVec) const = 0;

  virtual  void                     autoCorrViaDef      (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         unsigned int                          lag,
                                                         V&                                    corrVec) const = 0;
  virtual  void                     autoCorrViaFft      (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         const std::vector<unsigned int>&      lags,
                                                         std::vector<V*>&                      corrVecs) const = 0;
  virtual  void                     autoCorrViaFft      (unsigned int                          initialPos,
                                                         unsigned int                          numPos,
                                                         unsigned int                          numSum,
                                                         V&                                    autoCorrsSumVec) const = 0;
  virtual  void                     bmm                 (unsigned int                          initialPos,
                                                         unsigned int                          batchLength,
                                                         V&                                    bmmVec) const = 0;
  virtual  void                     fftForward          (unsigned int                          initialPos,
                                                         unsigned int                          fftSize,
                                                         unsigned int                          paramId,
                                                         std::vector<std::complex<double> >&   fftResult) const = 0;
//virtual  void                     fftInverse          (unsigned int fftSize) = 0;
  virtual  void                     psd                 (unsigned int                          initialPos,
                                                         unsigned int                          numBlocks,
                                                         double                                hopSizeRatio,
                                                         unsigned int                          paramId,
                                                         std::vector<double>&                  psdResult) const = 0;
  virtual  void                     psdAtZero           (unsigned int                          initialPos,
                                                         unsigned int                          numBlocks,
                                                         double                                hopSizeRatio,
                                                         V&                                    psdVec) const = 0;
  virtual  void                     geweke              (unsigned int                          initialPos,
                                                         double                                ratioNa,
                                                         double                                ratioNb,
                                                         V&                                    gewVec) const = 0;
  virtual  void                     minMax              (unsigned int                          initialPos,
                                                         V&                                    minVec,
                                                         V&                                    maxVec) const = 0;
  virtual  void                     histogram           (unsigned int                          initialPos,
                                                         const V&                              minVec,
                                                         const V&                              maxVec,
                                                         std::vector<V*>&                      centersForAllBins,
                                                         std::vector<V*>&                      binsForAllParams) const = 0;
  virtual  void                     interQuantileRange  (unsigned int                          initialPos,
                                                         V&                                    iqrVec) const = 0;
  virtual  void                     scalesForKDE        (unsigned int                          initialPos,
                                                         const V&                              iqrVec,
                                                         V&                                    scaleVec) const = 0;
  virtual  void                     gaussianKDE         (const V&                              evaluationParamVec,
                                                         V&                                    densityVec) const = 0;
  virtual  void                     gaussianKDE         (unsigned int                          initialPos,
                                                         const V&                              scaleVec,
                                                         const std::vector<V*>&                evaluationParamVecs,
                                                         std::vector<V*>&                      densityVecs) const = 0;
  virtual  void                     printContents       (std::ofstream&                        ofs) const = 0;
  virtual  void                     select              (const std::vector<unsigned int>&      idsOfUniquePositions) = 0;
  virtual  void                     filter              (unsigned int                          initialPos,
                                                         unsigned int                          spacing) = 0;

           void                     computeStatistics   (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         std::ofstream*                        passedOfs);

           void                     computeFilterParams (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         std::ofstream*                        passedOfs,
                                                         unsigned int&                         initialPos,
                                                         unsigned int&                         spacing);
protected:
           void                     computeMeanVars     (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         std::ofstream*                        passedOfs,
                                                         V*                                    meanPtr,
                                                         V*                                    sampleVarPtr,
                                                         V*                                    populVarPtr);
           void                     computeBMM          (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         std::ofstream*                        passedOfs);
           void                     computeFFT          (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         std::ofstream*                        passedOfs);
           void                     computePSD          (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         std::ofstream*                        passedOfs);
           void                     computePSDAtZero    (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         std::ofstream*                        passedOfs);
           void                     computeGeweke       (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         std::ofstream*                        passedOfs);
           void                     computeCorrViaDef   (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         const std::vector<unsigned int>&      lagsForCorrs,
                                                         std::ofstream*                        passedOfs);
           void                     computeCorrViaFFT   (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         const std::vector<unsigned int>&      initialPosForStatistics,
                                                         const std::vector<unsigned int>&      lagsForCorrs,
                                                         std::ofstream*                        passedOfs);
           void                     computeHistKde      (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                         std::ofstream*                        passedOfs);

  virtual  void                     extractScalarSeq    (unsigned int                          initialPos,
                                                         unsigned int                          spacing,
                                                         unsigned int                          numPos,
                                                         unsigned int                          paramId,
                                                         uqScalarSequenceClass<double>&        scalarSeq) const = 0;
  virtual  void                     extractRawData      (unsigned int                          initialPos,
                                                         unsigned int                          spacing,
                                                         unsigned int                          numPos,
                                                         unsigned int                          paramId,
                                                         std::vector<double>&                  rawData) const = 0;

  const uqBaseEnvironmentClass&      m_env;
  const uqVectorSpaceClass<V,M>& m_vectorSpace;
  std::string                    m_name;

  mutable uqFftClass<double>*    m_fftObj;
  mutable V*                     m_minValues;
  mutable V*                     m_maxValues;
  mutable V*                     m_meanValues;
  mutable V*                     m_sampleVarianceValues;
};

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::uqBaseVectorSequenceClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   sequenceSize,
  const std::string&             name)
  :
  m_env                 (vectorSpace.env()),
  m_vectorSpace         (vectorSpace),
  m_name                (name),
  m_fftObj              (new uqFftClass<double>(m_env)),
  m_minValues           (NULL),
  m_maxValues           (NULL),
  m_meanValues          (NULL),
  m_sampleVarianceValues(NULL)
{
}

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::~uqBaseVectorSequenceClass()
{
  //clear();
  if (m_sampleVarianceValues) delete m_sampleVarianceValues;
  if (m_meanValues          ) delete m_meanValues;
  if (m_maxValues           ) delete m_maxValues;
  if (m_minValues           ) delete m_minValues;
  if (m_fftObj != NULL      ) delete m_fftObj;
}

template <class V, class M>
unsigned int
uqBaseVectorSequenceClass<V,M>::vectorSize() const
{
  return m_vectorSpace.dim();
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
  unsigned int numPos = this->sequenceSize();
  if (numPos) {
    this->resetValues(0,numPos);
    this->resizeSequence(0);
  }

 return;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::minValues() const
{
  if (m_minValues == NULL) {
    m_minValues = m_vectorSpace.newVector();
    if (m_maxValues == NULL) m_maxValues = m_vectorSpace.newVector();
    minMax(0,*m_minValues,*m_maxValues);
  }

  return *m_minValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::maxValues() const
{
  if (m_maxValues == NULL) {
    if (m_minValues == NULL) m_minValues = m_vectorSpace.newVector();
    m_maxValues = m_vectorSpace.newVector();
    minMax(0,*m_minValues,*m_maxValues);
  }

  return *m_maxValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::meanValues() const
{
  if (m_meanValues == NULL) {
    m_meanValues = m_vectorSpace.newVector();
    mean(0,sequenceSize(),*m_meanValues);
  }

  return *m_meanValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::sampleVarianceValues() const
{
  if (m_sampleVarianceValues == NULL) {
    m_sampleVarianceValues = m_vectorSpace.newVector();
    sampleVariance(0,sequenceSize(),meanValues(),*m_sampleVarianceValues);
  }

  return *m_sampleVarianceValues;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::deleteStoredVectors()
{
  delete m_minValues;
  delete m_maxValues;
  delete m_meanValues;
  delete m_sampleVarianceValues;

  m_minValues            = NULL;
  m_maxValues            = NULL;
  m_meanValues           = NULL;
  m_sampleVarianceValues = NULL;

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeStatistics(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Computing statistics for chain " << m_name << " ..."
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  // Set initial positions for the computation of chain statistics
  std::vector<unsigned int> initialPosForStatistics(statisticalOptions.initialDiscardedPortions().size(),0);
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    initialPosForStatistics[i] = (unsigned int) (statisticalOptions.initialDiscardedPortions()[i] * (double) this->sequenceSize());
  }
  std::cout << "In uqBaseVectorSequenceClass<V,M>::computeStatistics(): initial positions for statistics =";
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    std::cout << " " << initialPosForStatistics[i];
  }
  std::cout << std::endl;

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

  // Set lags for the computation of chain autocorrelations
  std::vector<unsigned int> lagsForCorrs(statisticalOptions.corrNumLags(),1);
  for (unsigned int i = 1; i < lagsForCorrs.size(); ++i) {
    lagsForCorrs[i] = statisticalOptions.corrSecondLag() + (i-1)*statisticalOptions.corrLagSpacing();
  }

  //****************************************************
  // Compute autocorrelation coefficients via definition
  //****************************************************
  if ((statisticalOptions.corrComputeViaDef()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    this->computeCorrViaDef(statisticalOptions,
                            initialPosForStatistics,
                            lagsForCorrs,
                            passedOfs);
  }

  //****************************************************
  // Compute autocorrelation coefficients via FFT
  //****************************************************
  if ((statisticalOptions.corrComputeViaFft()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    this->computeCorrViaFFT(statisticalOptions,
                            initialPosForStatistics,
                            lagsForCorrs,
                            passedOfs);
  }

  //****************************************************
  // Compute histogram and/or Kde
  //****************************************************
  if ((statisticalOptions.histCompute()) ||
      (statisticalOptions.kdeCompute() )) {
    this->computeHistKde(statisticalOptions,
                         passedOfs);
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "All statistics took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
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
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs,
  V*                                    meanPtr,
  V*                                    sampleVarPtr,
  V*                                    populVarPtr)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing mean, sample variance and population variance"
              << std::endl;
  }

  V chainMean(m_vectorSpace.zeroVector());
  this->mean(0,
             this->sequenceSize(),
             chainMean);

  V chainSampleVariance(m_vectorSpace.zeroVector());
  this->sampleVariance(0,
                       this->sequenceSize(),
                       chainMean,
                       chainSampleVariance);

  if (m_env.rank() == 0) {
    std::cout << "\nEstimated variance of sample mean for the whole chain " << m_name
              << ", under independence assumption:"
              << std::endl;
  }
  V estimatedVarianceOfSampleMean(chainSampleVariance);
  estimatedVarianceOfSampleMean /= (double) this->sequenceSize();
  bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
  estimatedVarianceOfSampleMean.setPrintHorizontally(false);
  std::cout << estimatedVarianceOfSampleMean;
  estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);
  if (m_env.rank() == 0) {
    std::cout << std::endl;
  }

  V chainPopulationVariance(m_vectorSpace.zeroVector());
  this->populationVariance(0,
                           this->sequenceSize(),
                           chainMean,
                           chainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Mean and variances took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  if (m_env.rank() == 0) {
    std::cout << "\nMean, sample std, population std"
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
    std::cout << line;

    for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e",
              m_vectorSpace.componentName(i).c_str(), /*.*/
              " ",
	      chainMean[i],
              " ",
              sqrt(chainSampleVariance[i]),
              " ",
              sqrt(chainPopulationVariance[i]));
      std::cout << line;
    }
    std::cout << std::endl;
  }

  if (meanPtr     ) *meanPtr      = chainMean;
  if (sampleVarPtr) *sampleVarPtr = chainSampleVariance;
  if (populVarPtr ) *populVarPtr  = chainPopulationVariance;

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeBMM(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing variance of sample mean through BMM"
              << std::endl;
  }

  std::cout << "In uqBaseVectorSequenceClass<V,M>::computeBMM(): lengths for batchs in BMM =";
  for (unsigned int i = 0; i < statisticalOptions.bmmLengths().size(); ++i) {
    std::cout << " " << statisticalOptions.bmmLengths()[i];
  }
  std::cout << std::endl;

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

  if (m_env.rank() == 0) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      std::cout << "\nEstimated variance of sample mean, through batch means method, for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << " (each column corresponds to a batch length)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.bmmLengths()[batchLengthId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        std::cout << line;
        for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfBMM(initialPosId,batchLengthId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain BMM took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeFFT(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
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
      std::ofstream& ofs = *passedOfs;
      ofs << m_name << "_fft_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << 1
          << ","                                                                                 << forwardResult.size()
          << ");"
          << std::endl;
      for (unsigned int j = 0; j < forwardResult.size(); ++j) {
        ofs << m_name << "_fft_initPos" << initialPosForStatistics[initialPosId] << "(" << 1
            << ","                                                                         << j+1
            << ") = "                                                                      << forwardResult[j].real()
            << " + i*"                                                                     << forwardResult[j].imag()
            << ";"
            << std::endl;
      }
    } // if write

    if (statisticalOptions.fftTestInversion()) {
      fftObj.inverse(forwardResult,
                     statisticalOptions.fftSize(),
                     inverseResult);
      if (statisticalOptions.fftWrite() && passedOfs) {
        std::ofstream& ofs = *passedOfs;
        ofs << m_name << "_inv_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << 1
            << ","                                                                                 << inverseResult.size()
            << ");"
            << std::endl;
        for (unsigned int j = 0; j < inverseResult.size(); ++j) {
          ofs << m_name << "_inv_initPos" << initialPosForStatistics[initialPosId] << "(" << 1
              << ","                                                                         << j+1
              << ") = "                                                                      << inverseResult[j].real()
              << " + i*"                                                                     << inverseResult[j].imag()
              << ";"
              << std::endl;
        }
      } // if write
    }
  } // for initialPosId

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain FFT took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computePSD(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
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
      std::ofstream& ofs = *passedOfs;
      ofs << m_name << "_psd_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << 1
          << ","                                                                                 << psdResult.size()
          << ");"
          << std::endl;
      for (unsigned int j = 0; j < psdResult.size(); ++j) {
        ofs << m_name << "_psd_initPos" << initialPosForStatistics[initialPosId] << "(" << 1
            << ","                                                                         << j+1
            << ") = "                                                                      << psdResult[j]
            << ";"
            << std::endl;
      }
    } // if write
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain PSD took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computePSDAtZero(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
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
  if ((statisticalOptions.psdAtZeroDisplay()) && (m_env.rank() == 0)) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      std::cout << "\nComputed PSD at frequency zero for subchain beggining at position " << initialPos
                << ", so effective data size = " << this->sequenceSize() - initialPos
                << " (each column corresponds to a number of blocks)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        std::cout << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  // Display estimated variance of sample mean through PSD
  if (/*(statisticalOptions.psdAtZeroDisplay()) &&*/ (m_env.rank() == 0)) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      std::cout << "\nEstimated variance of sample mean, through psd, for subchain beggining at position " << initialPos
                << ", so effective data size = " << this->sequenceSize() - initialPos
                << " (each column corresponds to a number of blocks)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        std::cout << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  2.*M_PI*_2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]/(double) (this->sequenceSize() - initialPos));
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain PSD at frequency zero took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  // Write PSD at frequency zero
  if (statisticalOptions.psdAtZeroWrite() && passedOfs) {
    std::ofstream& ofs = *passedOfs;
    ofs << m_name << "_psdAtZero_numBlocks = zeros(" << 1
        << ","                                          << statisticalOptions.psdAtZeroNumBlocks().size()
        << ");"
        << std::endl;
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      ofs << m_name << "_psdAtZero_numBlocks(" << 1
          << ","                                  << numBlocksId+1
          << ") = "                               << statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]
          << ";"
          << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofs << m_name << "_psdAtZero_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << this->vectorSize() /*.*/
          << ","                                                                                       << statisticalOptions.psdAtZeroNumBlocks().size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          ofs << m_name << "_psdAtZero_initPos" << initialPosForStatistics[initialPosId] << "(" << i+1
              << ","                                                                               << numBlocksId+1
              << ") = "                                                                            << _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]
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
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
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

  if (m_env.rank() == 0) {
    std::cout << "\nComputed Geweke coefficients with 10% and 50% percentages"
              << " (each column corresponds to a different initial position on the full chain)"
              << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    std::cout << line;
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      sprintf(line,"%10s%3d",
              " ",
              initialPosForStatistics[initialPosId]);
      std::cout << line;
    }

    for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
      sprintf(line,"\n%9.9s",
              m_vectorSpace.componentName(i).c_str() /*.*/);
      std::cout << line;
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        sprintf(line,"%2s%11.4e",
                " ",
                (*(vectorOfGeweke[initialPosId]))[i]);
        std::cout << line;
      }
    }
    std::cout << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain Geweke took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeCorrViaDef(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing autocorrelation coefficients (via def)"
              << std::endl;
  }

  if (statisticalOptions.corrDisplay() && (m_env.rank() == 0)) {
    std::cout << "In uqBaseVectorSequenceClass<V,M>::computeCorrViaDef(): lags for autocorrelation (via def) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      std::cout << " " << lagsForCorrs[i];
    }
    std::cout << std::endl;
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
                           this->sequenceSize()-initialPos,
                           lag,
                           _2dArrayOfAutoCorrs(initialPosId,lagId));
      //_2dArrayOfAutoCorrs(initialPosId,lagId) = corrVec;
    }
  }

  // It is not practical to compute the variance of sample mean by computing the autocorrelations via definition for each lag
  // The code computes the variance of sample mean by computing the autocorrelations via fft, below, in another routine

  if ((statisticalOptions.corrDisplay()) && (m_env.rank() == 0)) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      std::cout << "\nComputed autocorrelation coefficients (via def), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << " (each column corresponds to a different lag)"
                << std::endl;
      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        sprintf(line,"%10s%3d",
                " ",
                lagsForCorrs[lagId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        std::cout << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain autocorrelation (via def) took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.corrWrite() && passedOfs) {
    std::ofstream& ofs = *passedOfs;
    ofs << m_name << "_corrViaDef_lags = zeros(" << 1
        << ","                                      << lagsForCorrs.size()
        << ");"
        << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofs << m_name << "_corrViaDef_lags(" << 1
          << ","                              << lagId+1
          << ") = "                           << lagsForCorrs[lagId]
          << ";"
          << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofs << m_name << "_corrViaDef_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << this->vectorSize() /*.*/
          << ","                                                                                        << lagsForCorrs.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofs << m_name << "_corrViaDef_initPos" << initialPosForStatistics[initialPosId] << "(" << i+1
              << ","                                                                                << lagId+1
              << ") = "                                                                             << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
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
uqBaseVectorSequenceClass<V,M>::computeCorrViaFFT(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing autocorrelation coefficients (via fft)"
              << std::endl;
  }

  if (statisticalOptions.corrDisplay() && (m_env.rank() == 0)) {
    std::cout << "In uqBaseVectorSequenceClass<V,M>::computeCorrViaFFT(): lags for autocorrelation (via fft) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      std::cout << " " << lagsForCorrs[i];
     }
     std::cout << std::endl;
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
    if (m_env.rank() == 0) {
      std::cout << "In uqBaseVectorSequenceClass<V,M>::computeCorrViaFFT()"
                << ": about to call chain.autoCorrViaFft()"
                << " with initialPos = "      << initialPos
                << ", numPos = "              << this->sequenceSize()-initialPos
                << ", lagsForCorrs.size() = " << lagsForCorrs.size()
                << ", corrVecs.size() = "     << corrVecs.size()
                << std::endl;
    }
    this->autoCorrViaFft(initialPos,
                         this->sequenceSize()-initialPos, // Use all possible data positions
                         lagsForCorrs,
                         corrVecs);
    this->autoCorrViaFft(initialPos,
                         this->sequenceSize()-initialPos, // Use all possible data positions
                         (unsigned int) (1.0 * (double) (this->sequenceSize()-initialPos)), // CHECK
                         *corrSumVecs[initialPosId]); // Sum of all possibly computable autocorrelations, not only the asked ones in lagsForCorrs
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      _2dArrayOfAutoCorrs(initialPosId,lagId) = *(corrVecs[lagId]);
    }
  }
  for (unsigned int j = 0; j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }

  if ((statisticalOptions.corrDisplay()) && (m_env.rank() == 0)) {
    V chainMean                    (m_vectorSpace.zeroVector());
    V chainSampleVariance          (m_vectorSpace.zeroVector());
    V estimatedVarianceOfSampleMean(m_vectorSpace.zeroVector());
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];

      this->mean(initialPos,
                 this->sequenceSize()-initialPos,
                 chainMean);

      this->sampleVariance(initialPos,
                           this->sequenceSize()-initialPos,
                           chainMean,
                           chainSampleVariance);

      std::cout << "\nEstimated variance of sample mean, through autocorrelation (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << std::endl;
      estimatedVarianceOfSampleMean.cwSet(-1.); // Yes, '-1' because the autocorrelation at lag 0, which values '+1', is already counted in the sum
      estimatedVarianceOfSampleMean += 2.* (*corrSumVecs[initialPosId]);
      estimatedVarianceOfSampleMean *= chainSampleVariance;
      estimatedVarianceOfSampleMean /= (double) (this->sequenceSize() - initialPos);
      bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
      estimatedVarianceOfSampleMean.setPrintHorizontally(false);
      std::cout << estimatedVarianceOfSampleMean;
      estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);
      if (m_env.rank() == 0) {
        std::cout << std::endl;
      }

      std::cout << "\nComputed autocorrelation coefficients (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << " (each column corresponds to a different lag)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        sprintf(line,"%10s%3d",
                " ",
                lagsForCorrs[lagId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        std::cout << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain autocorrelation (via fft) took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.corrWrite() && passedOfs) {
    std::ofstream& ofs = *passedOfs;
    ofs << m_name << "_corrViaFft_lags = zeros(" << 1
        << ","                                      << lagsForCorrs.size()
        << ");"
        << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofs << m_name << "_corrViaFft_lags(" << 1
          << ","                              << lagId+1
          << ") = "                           << lagsForCorrs[lagId]
          << ";"
          << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofs << m_name << "_corrViaFft_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << this->vectorSize() /*.*/
          << ","                                                                                        << lagsForCorrs.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofs << m_name << "_corrViaFft_initPos" << initialPosForStatistics[initialPosId] << "(" << i+1
              << ","                                                                                << lagId+1
              << ") = "                                                                             << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
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
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs,
  unsigned int&                         initialPos,
  unsigned int&                         spacing)
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Computing filter parameters for chain " << m_name << " ..."
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  initialPos = 0;
  spacing    = 1;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
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
uqBaseVectorSequenceClass<V,M>::computeHistKde( // Use the whole chain
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Computing histogram and/or KDE for chain " << m_name << " ..."
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
  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing min and max for histograms and KDE"
              << std::endl;
  }

  V statsMinPositions(m_vectorSpace.zeroVector());
  V statsMaxPositions(m_vectorSpace.zeroVector());
  this->minMax(0, // Use the whole chain
               statsMinPositions,
               statsMaxPositions);

  if (m_env.rank() == 0) {
    std::cout << "\nComputed min values and max values for chain " << m_name
              << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    std::cout << line;

    sprintf(line,"%9s%s%9s%s",
            " ",
            "min",
            " ",
            "max");
    std::cout << line;

    for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
      sprintf(line,"\n%8.8s",
              m_vectorSpace.componentName(i).c_str() /*.*/);
      std::cout << line;

      sprintf(line,"%2s%11.4e%2s%11.4e",
              " ",
              statsMinPositions[i],
              " ",
              statsMaxPositions[i]);
      std::cout << line;
    }
    std::cout << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain min and max took " << tmpRunTime
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
    if (m_env.rank() == 0) {
      std::cout << "\n-----------------------------------------------------"
                << "\nComputing histograms"
                << std::endl;
    }

    std::vector<V*> histCentersForAllBins(0);
    std::vector<V*> histBinsForAllParams(0);

    for (unsigned int i = 0; i < statsMaxPositions.size(); ++i) {
      statsMaxPositions[i] *= (1. + 1.e-15);
    }

    histCentersForAllBins.resize(statisticalOptions.histNumInternalBins()+2,NULL);
    histBinsForAllParams.resize (statisticalOptions.histNumInternalBins()+2,NULL);
    this->histogram(0, // Use the whole chain
                    statsMinPositions,
                    statsMaxPositions,
                    histCentersForAllBins,
                    histBinsForAllParams);

    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.rank() == 0) {
      std::cout << "Chain histograms took " << tmpRunTime
                << " seconds"
                << std::endl;
    }

    // Write histograms
    // plot(queso_centersOfHistBins(1,:)',queso_histBins(1,:)','r-');
    if (passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << m_name << "_centersOfHistBins = zeros(" << this->vectorSize() /*.*/
          << ","                                        << histCentersForAllBins.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < histCentersForAllBins.size(); ++j) {
           ofs << m_name << "_centersOfHistBins(" << i+1
               << ","                                << j+1
               << ") = "                             << (*(histCentersForAllBins[j]))[i]
               << ";"
               << std::endl;
        }
      }

      ofs << m_name << "_histBins = zeros(" << this->vectorSize() /*.*/
          << ","                               << histBinsForAllParams.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < histBinsForAllParams.size(); ++j) {
           ofs << m_name << "_histBins(" << i+1
               << ","                       << j+1
               << ") = "                    << (*(histBinsForAllParams[j]))[i]
               << ";"
               << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < histBinsForAllParams.size(); ++i) {
      if (histBinsForAllParams[i] != NULL) delete histBinsForAllParams[i];
    }
    for (unsigned int i = 0; i < histCentersForAllBins.size(); ++i) {
      if (histCentersForAllBins[i] != NULL) delete histCentersForAllBins[i];
    }
  }

  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((statisticalOptions.kdeCompute()             ) &&
      (statisticalOptions.kdeNumEvalPositions() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.rank() == 0) {
      std::cout << "\n-----------------------------------------------------"
                << "\nComputing KDE"
                << std::endl;
    }

    std::vector<V*> kdeEvalPositions(0);
    V               gaussianKdeScaleVec(m_vectorSpace.zeroVector());
    std::vector<V*> gaussianKdeDensities(0);

    kdeEvalPositions.resize(statisticalOptions.kdeNumEvalPositions(),NULL);
    uqMiscComputePositionsBetweenMinMax(statsMinPositions,
                                        statsMaxPositions,
                                        kdeEvalPositions);

    V iqrVec(m_vectorSpace.zeroVector());
    this->interQuantileRange(0, // Use the whole chain
                             iqrVec);

    if (m_env.rank() == 0) {
      std::cout << "\nComputed inter quantile ranges for chain " << m_name
                  << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;

      sprintf(line,"%9s%s",
              " ",
              "iqr");
      std::cout << line;

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%8.8s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        std::cout << line;

        sprintf(line,"%2s%11.4e",
                " ",
                iqrVec[i]);
        std::cout << line;
      }
      std::cout << std::endl;
    }

    this->scalesForKDE(0, // Use the whole chain
                       iqrVec,
                       gaussianKdeScaleVec);

    gaussianKdeDensities.resize(statisticalOptions.kdeNumEvalPositions(),NULL);
    this->gaussianKDE(0, // Use the whole chain
                      gaussianKdeScaleVec,
                      kdeEvalPositions,
                      gaussianKdeDensities);

    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.rank() == 0) {
      std::cout << "Chain KDE took " << tmpRunTime
                << " seconds"
                << std::endl;
    }

    // Write estimations of probability densities
    // hold
    // plot(queso_kdeEvalPositions(1,:)',7*queso_gaussianKdeDensities(1,:)','r-');
    if (passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << m_name << "_kdeEvalPositions = zeros(" << this->vectorSize() /*.*/
          << ","                                       << kdeEvalPositions.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < kdeEvalPositions.size(); ++j) {
          ofs << m_name << "_kdeEvalPositions(" << i+1
              << ","                               << j+1
              << ") = "                            << (*(kdeEvalPositions[j]))[i]
              << ";"
              << std::endl;
        }
      }

      ofs << m_name << "_gaussianKdeScaleVec = zeros(" << this->vectorSize() /*.*/
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        ofs << m_name << "_gaussianKdeScaleVec(" << i+1
            << ") = "                               << gaussianKdeScaleVec[i]
            << ";"
            << std::endl;
      }

      ofs << m_name << "_gaussianKdeDensities = zeros(" << this->vectorSize() /*.*/
          << ","                                           << gaussianKdeDensities.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < gaussianKdeDensities.size(); ++j) {
          ofs << m_name << "_gaussianKdeDensities(" << i+1
              << ","                                   << j+1
              << ") = "                                << (*(gaussianKdeDensities[j]))[i]
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
  }

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\n Finished computing histogram and/or KDE for chain " << m_name
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  return;
}

#endif // __UQ_VECTOR_SEQUENCE_H__
