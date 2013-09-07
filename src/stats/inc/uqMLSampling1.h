//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef __UQ_MULTI_LEVEL_SAMPLING1_H__
#define __UQ_MULTI_LEVEL_SAMPLING1_H__

#define ML_NEW_CODE_2009_12_29

#include <uqMLSamplingOptions.h>
#include <uqMetropolisHastingsSG1.h>
#include <uqFiniteDistribution.h>
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#ifdef QUESO_HAS_GLPK
#include <glpk.h>
#endif
#include <sys/time.h>
#include <fstream>

#define ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA 6

//---------------------------------------------------------

namespace QUESO {

// aqui 1
#ifdef QUESO_HAS_GLPK
struct BIP_routine_struct {
  const BaseEnvironment* env;
  unsigned int                  currLevel;
};

void BIP_routine(glp_tree *tree, void *info);
#endif

//---------------------------------------------------------

struct ExchangeInfoStruct
{
  int          originalNodeOfInitialPosition;
  unsigned int originalIndexOfInitialPosition;
  int          finalNodeOfInitialPosition;
  unsigned int numberOfPositions;
};

//---------------------------------------------------------

template <class P_V>
struct BalancedLinkedChainControlStruct
{
  P_V*         initialPosition; // KEY new: add logPrior and logLikelihood
  unsigned int numberOfPositions;
};

template <class P_V>
struct BalancedLinkedChainsPerNodeStruct
{
  std::vector<BalancedLinkedChainControlStruct<P_V> > balLinkedChains;
};

//---------------------------------------------------------

struct UnbalancedLinkedChainControlStruct
{
  unsigned int initialPositionIndexInPreviousChain; // KEY new: add logPrior and logLikelihood
  unsigned int numberOfPositions;
};

struct UnbalancedLinkedChainsPerNodeStruct
{
  std::vector<UnbalancedLinkedChainControlStruct> unbLinkedChains;
};

//---------------------------------------------------------

/*! A templated class that represents a Multi Level sampler.
 */
template <class P_V,class P_M>
class MLSampling
{
public:
  /*! Constructor: */
  MLSampling(/*! Prefix                  */ const char*                               prefix,
                    /*! The prior rv            */ const BaseVectorRV      <P_V,P_M>& priorRv,
                    /*! The likelihood function */ const BaseScalarFunction<P_V,P_M>& likelihoodFunction);
  /*! Destructor: */
 ~MLSampling();

  /*! Operation to generate the chain */
  void   generateSequence (BaseVectorSequence<P_V,P_M>& workingChain,
                           ScalarSequence<double>*      workingLogLikelihoodValues,
                           ScalarSequence<double>*      workingLogTargetValues);
  double logEvidence      () const;
  double meanLogLikelihood() const;
  double eig              () const;

  void   print           (std::ostream& os) const;

private:
  // Methods available at uqMLSampling2.h
  void   checkpointML                  (double                                          currExponent,                       // input
                                        double                                          currEta,                            // input
                                        const SequenceOfVectors<P_V,P_M>&        currChain,                          // input
                                        const ScalarSequence<double>&            currLogLikelihoodValues,            // input
                                        const ScalarSequence<double>&            currLogTargetValues);               // input

  void   restartML                     (double&                                         currExponent,                       // output
                                        double&                                         currEta,                            // output
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // output
                                        ScalarSequence<double>&                  currLogTargetValues);               // output

  void   generateSequence_Level0_all   (const MLSamplingLevelOptions&            currOptions,                        // input
                                        unsigned int&                                   unifiedRequestedNumSamples,         // output
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // output
                                        ScalarSequence<double>&                  currLogTargetValues);               // output

  void   generateSequence_Step01_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int&                                   unifiedRequestedNumSamples);        // output

  void   generateSequence_Step02_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // input/output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // input/output
                                        ScalarSequence<double>&                  currLogTargetValues,                // input/output
                                        SequenceOfVectors<P_V,P_M>&              prevChain,                          // output
                                        ScalarSequence<double>&                  prevLogLikelihoodValues,            // output
                                        ScalarSequence<double>&                  prevLogTargetValues,                // output
                                        unsigned int&                                   indexOfFirstWeight,                 // output
                                        unsigned int&                                   indexOfLastWeight);                 // output

  void   generateSequence_Step03_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                        double                                          prevExponent,                       // input
                                        double                                          failedExponent,                     // input // gpmsa
                                        double&                                         currExponent,                       // output
                                        ScalarSequence<double>&                  weightSequence);                    // output

  void   generateSequence_Step04_inter0(const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        const ScalarSequence<double>&            weightSequence,                     // input
                                        P_M&                                            unifiedCovMatrix);                  // output

  void   generateSequence_Step05_inter0(unsigned int                                    unifiedRequestedNumSamples,         // input
                                        const ScalarSequence<double>&            weightSequence,                     // input
                                        std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only,    // output
                                        std::vector<double>&                            unifiedWeightStdVectorAtProc0Only); // output

  void   generateSequence_Step06_all   (const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        bool&                                           useBalancedChains,                  // output
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // output

  void   generateSequence_Step07_inter0(bool                                            useBalancedChains,                  // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        UnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl,              // (possible) output
                                        const MLSamplingLevelOptions*            currOptions,                        // input
                                        const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec,                     // (possible) input/output
                                        BalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // (possible) output

  void   generateSequence_Step08_all   (BayesianJointPdf<P_V,P_M>&               currPdf,                            // input/output
                                        GenericVectorRV<P_V,P_M>&                currRv);                            // output

  void   generateSequence_Step09_all   (const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only,  // input
                                        const ScalarSequence<double>&            weightSequence,                     // input
                                        double                                          prevEta,                            // input
                                        const GenericVectorRV<P_V,P_M>&          currRv,                             // input
                                        MLSamplingLevelOptions*                  currOptions,                        // input (changed temporarily internally)
                                        P_M&                                            unifiedCovMatrix,                   // input/output
                                        double&                                         currEta);                           // output

  void   generateSequence_Step10_all   (MLSamplingLevelOptions&                  currOptions,                        // input (changed temporarily internally)
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const GenericVectorRV  <P_V,P_M>&        currRv,                             // input
                                        bool                                            useBalancedChains,                  // input
                                        const UnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,              // input // Round Rock
                                        unsigned int                                    indexOfFirstWeight,                 // input // Round Rock
                                        const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input // Round Rock
                                        const BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,                // input // Round Rock
                                        SequenceOfVectors      <P_V,P_M>&        currChain,                          // output
                                        double&                                         cumulativeRawChainRunTime,          // output
                                        unsigned int&                                   cumulativeRawChainRejections,       // output
                                        ScalarSequence         <double>*         currLogLikelihoodValues,            // output
                                        ScalarSequence         <double>*         currLogTargetValues);               // output

  void   generateSequence_Step11_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int                                    unifiedRequestedNumSamples,         // input
                                        unsigned int                                    cumulativeRawChainRejections,       // input
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // input/output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // input/output
                                        ScalarSequence<double>&                  currLogTargetValues,                // input/output
                                        unsigned int&                                   unifiedNumberOfRejections);         // output

  // Methods available at uqMLSampling3.h
  void   sampleIndexes_proc0           (unsigned int                                    unifiedRequestedNumSamples,         // input
                                        const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only,  // input
                                        std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only);   // output

  bool   decideOnBalancedChains_all    (const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // output

  void   prepareBalLinkedChains_inter0 (const MLSamplingLevelOptions*            currOptions,                        // input
                                        const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec,                     // input/output
                                        BalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // output

  void   prepareUnbLinkedChains_inter0 (unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        UnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl);             // output

  void   generateBalLinkedChains_all   (MLSamplingLevelOptions&                  inputOptions,                       // input, only m_rawChainSize changes
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const GenericVectorRV  <P_V,P_M>&        rv,                                 // input
                                        const BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,                // input // Round Rock
                                        SequenceOfVectors      <P_V,P_M>&        workingChain,                       // output
                                        double&                                         cumulativeRunTime,                  // output
                                        unsigned int&                                   cumulativeRejections,               // output
                                        ScalarSequence         <double>*         currLogLikelihoodValues,            // output
                                        ScalarSequence         <double>*         currLogTargetValues);               // output

  void   generateUnbLinkedChains_all   (MLSamplingLevelOptions&                  inputOptions,                       // input, only m_rawChainSize changes
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const GenericVectorRV  <P_V,P_M>&        rv,                                 // input
                                        const UnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,              // input // Round Rock
                                        unsigned int                                    indexOfFirstWeight,                 // input // Round Rock
                                        const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input // Round Rock
                                        SequenceOfVectors      <P_V,P_M>&        workingChain,                       // output
                                        double&                                         cumulativeRunTime,                  // output
                                        unsigned int&                                   cumulativeRejections,               // output
                                        ScalarSequence         <double>*         currLogLikelihoodValues,            // output
                                        ScalarSequence         <double>*         currLogTargetValues);               // output

#ifdef QUESO_HAS_GLPK
  void   solveBIP_proc0                (std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // input/output
#endif

  void   justBalance_proc0             (const MLSamplingLevelOptions*            currOptions,                        // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // input/output

  void   mpiExchangePositions_inter0   (const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        const std::vector<ExchangeInfoStruct>&        exchangeStdVec,                     // input
                                        const std::vector<unsigned int>&                finalNumChainsPerNode,              // input
                                        const std::vector<unsigned int>&                finalNumPositionsPerNode,           // input
                                        BalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // output

  // Private variables
  const BaseEnvironment&             m_env;
  const BaseVectorRV      <P_V,P_M>& m_priorRv;
  const BaseScalarFunction<P_V,P_M>& m_likelihoodFunction;
  const VectorSpace       <P_V,P_M>& m_vectorSpace;
        VectorSet         <P_V,P_M>* m_targetDomain;

        MLSamplingOptions            m_options;

        unsigned int                        m_currLevel;          // restart
        unsigned int                        m_currStep;
        double                              m_debugExponent;
	std::vector<double>                 m_logEvidenceFactors; // restart
        double                              m_logEvidence;
        double                              m_meanLogLikelihood;
        double                              m_eig;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const MLSampling<P_V,P_M>& obj);
}  // End namespace QUESO

#include <uqMLSampling3.h>  // We should do something about this
#include <uqMLSampling2.h>  // We should do something about this

namespace QUESO {
template<class P_V,class P_M>
MLSampling<P_V,P_M>::MLSampling(
  const char*                               prefix,
  const BaseVectorRV      <P_V,P_M>& priorRv,            
  const BaseScalarFunction<P_V,P_M>& likelihoodFunction)
  :
  m_env               (priorRv.env()),
  m_priorRv           (priorRv),
  m_likelihoodFunction(likelihoodFunction),
  m_vectorSpace       (m_priorRv.imageSet().vectorSpace()),
  m_targetDomain      (InstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet())),
  m_options           (m_env,prefix),
  m_currLevel         (0),
  m_currStep          (0),
  m_debugExponent     (0.),
  m_logEvidenceFactors(0),
  m_logEvidence       (0.),
  m_meanLogLikelihood (0.),
  m_eig               (0.)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::constructor()"
                            << std::endl;
  }

  m_options.scanOptionsValues();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::constructor()"
                            << std::endl;
  }
}

template<class P_V,class P_M>
MLSampling<P_V,P_M>::~MLSampling()
{
  if (m_targetDomain) delete m_targetDomain;
}

template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence(
  BaseVectorSequence<P_V,P_M>& workingChain,
  ScalarSequence<double>*      workingLogLikelihoodValues,
  ScalarSequence<double>*      workingLogTargetValues)
{
  struct timeval timevalRoutineBegin;
  int iRC = 0;
  iRC = gettimeofday(&timevalRoutineBegin, NULL);
  if (iRC) {}; // just to remove compiler warning

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::generateSequence()"
                            << ", at  "   << ctime(&timevalRoutineBegin.tv_sec)
                            << ", after " << timevalRoutineBegin.tv_sec - m_env.timevalBegin().tv_sec
                            << " seconds from queso environment instatiation..."
                            << std::endl;
  }

  //***********************************************************
  // Declaration of Variables
  //***********************************************************
  double                            currExponent                   = 0.;   // restate
  double                            currEta                        = 1.;   // restate
  unsigned int                      currUnifiedRequestedNumSamples = 0;
  SequenceOfVectors<P_V,P_M> currChain              (m_vectorSpace, // restate
                                                            0,
                                                            m_options.m_prefix+"curr_chain");
  ScalarSequence<double>     currLogLikelihoodValues(m_env,         // restate
                                                            0,
                                                            "");
  ScalarSequence<double>     currLogTargetValues    (m_env,         // restate
                                                            0,
                                                            "");

  bool stopAtEndOfLevel = false;
  char levelPrefix[256];

  //***********************************************************
  // Take care of first level (level '0')
  //***********************************************************
  MLSamplingLevelOptions defaultLevelOptions(m_env,(m_options.m_prefix + "default_").c_str());
  defaultLevelOptions.scanOptionsValues(NULL);

  MLSamplingLevelOptions lastLevelOptions(m_env,(m_options.m_prefix + "last_").c_str());
  lastLevelOptions.scanOptionsValues(&defaultLevelOptions);

  if (m_options.m_restartInput_baseNameForFiles != ".") {
    restartML(currExponent,            // output
              currEta,                 // output
              currChain,               // output
              currLogLikelihoodValues, // output
              currLogTargetValues);    // output

    if (currExponent == 1.) {
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      if (lastLevelOptions.m_rawChainComputeStats) {
        FilePtrSetStruct filePtrSet;
        m_env.openOutputFile(lastLevelOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                             lastLevelOptions.m_dataOutputAllowedSet,
                             false,
                             filePtrSet);

        currChain.computeStatistics(*lastLevelOptions.m_rawChainStatisticalOptionsObj,
                                    filePtrSet.ofsVar);

        m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
      }
#endif
      // Compute MLE and MAP
      // rr0

      if (lastLevelOptions.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
        currChain.unifiedWriteContents              (lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType);
        currLogLikelihoodValues.unifiedWriteContents(lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType);
        currLogTargetValues.unifiedWriteContents    (lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType);
      }

      if (lastLevelOptions.m_filteredChainGenerate) {
        FilePtrSetStruct filePtrSet;
        m_env.openOutputFile(lastLevelOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                             lastLevelOptions.m_dataOutputAllowedSet,
                             false,
                             filePtrSet);

        unsigned int filterInitialPos = (unsigned int) (lastLevelOptions.m_filteredChainDiscardedPortion * (double) currChain.subSequenceSize());
        unsigned int filterSpacing    = lastLevelOptions.m_filteredChainLag;
        if (filterSpacing == 0) {
          currChain.computeFilterParams(filePtrSet.ofsVar,
                                        filterInitialPos,
                                        filterSpacing);
        }

        // Filter positions from the converged portion of the chain
        currChain.filter(filterInitialPos,
                         filterSpacing);
        currChain.setName(lastLevelOptions.m_prefix + "filtChain");

        currLogLikelihoodValues.filter(filterInitialPos,
                                       filterSpacing);
        currLogLikelihoodValues.setName(lastLevelOptions.m_prefix + "filtLogLikelihood");

        currLogTargetValues.filter(filterInitialPos,
                                   filterSpacing);
        currLogTargetValues.setName(lastLevelOptions.m_prefix + "filtLogTarget");

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        if (lastLevelOptions.m_filteredChainComputeStats) {
          currChain.computeStatistics(*lastLevelOptions.m_filteredChainStatisticalOptionsObj,
                                      filePtrSet.ofsVar);
        }
#endif
        // Compute MLE and MAP
        // rr0
        m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);

        if (lastLevelOptions.m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
          currChain.unifiedWriteContents              (lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
          currLogLikelihoodValues.unifiedWriteContents(lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
          currLogTargetValues.unifiedWriteContents    (lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
        }
      } // if (lastLevelOptions.m_filteredChainGenerate)
    }
  }
  else {
    sprintf(levelPrefix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    MLSamplingLevelOptions currOptions(m_env,(m_options.m_prefix + levelPrefix).c_str());
    currOptions.scanOptionsValues(&defaultLevelOptions);

    generateSequence_Level0_all(currOptions,                    // input
                                currUnifiedRequestedNumSamples, // output
                                currChain,                      // output
                                currLogLikelihoodValues,        // output
                                currLogTargetValues);           // output

    stopAtEndOfLevel = currOptions.m_stopAtEnd;
    bool performCheckpoint = stopAtEndOfLevel;
    if (m_options.m_restartOutput_levelPeriod > 0) {
      performCheckpoint = performCheckpoint || ( ((m_currLevel + 1) % m_options.m_restartOutput_levelPeriod) == 0 );
    }
    if (performCheckpoint) {
      checkpointML(currExponent,            // input
                   currEta,                 // input
                   currChain,               // input
                   currLogLikelihoodValues, // input
                   currLogTargetValues);    // input
    }
  }
  //std::cout << "In QUESO: end of level 0. Exiting on purpose" << std::endl;
  //exit(1);

  double minLogLike = -INFINITY;
  double maxLogLike =  INFINITY;
  currLogLikelihoodValues.subMinMaxExtra(0,
                                         currLogLikelihoodValues.subSequenceSize(),
                                         minLogLike,
                                         maxLogLike);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                            << ": at end of level "  << m_currLevel+LEVEL_REF_ID
                            << ", sub minLogLike = " << minLogLike
                            << ", sub maxLogLike = " << maxLogLike
                            << std::endl;
  }

  m_env.fullComm().Barrier();

  minLogLike = -INFINITY;
  maxLogLike =  INFINITY;
  currLogLikelihoodValues.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                             0,
                                             currLogLikelihoodValues.subSequenceSize(),
                                             minLogLike,
                                             maxLogLike);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                            << ": at end of level "      << m_currLevel+LEVEL_REF_ID
                            << ", unified minLogLike = " << minLogLike
                            << ", unified maxLogLike = " << maxLogLike
                            << std::endl;
  }

  //***********************************************************
  // Take care of next levels
  //***********************************************************
  while ((currExponent     <  1.   ) && // begin level while
         (stopAtEndOfLevel == false)) {
    m_currLevel++; // restate

    struct timeval timevalLevelBegin;
    iRC = 0;
    iRC = gettimeofday(&timevalLevelBegin, NULL);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level " << m_currLevel+LEVEL_REF_ID
                              << ", at  "   << ctime(&timevalLevelBegin.tv_sec)
                              << ", after " << timevalLevelBegin.tv_sec - timevalRoutineBegin.tv_sec
                              << " seconds from entering the routine"
                              << ", after " << timevalLevelBegin.tv_sec - m_env.timevalBegin().tv_sec
                              << " seconds from queso environment instatiation"
                              << std::endl;
    }

    iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);
    double       cumulativeRawChainRunTime    = 0.;
    unsigned int cumulativeRawChainRejections = 0;

    bool   tryExponentEta = true; // gpmsa
    double failedExponent = 0.;   // gpmsa
    double failedEta      = 0.;   // gpmsa

    MLSamplingLevelOptions*            currOptions           = NULL;  // step 1
    SequenceOfVectors<P_V,P_M>*        prevChain             = NULL;  // step 2
    unsigned int                              indexOfFirstWeight    = 0;     // step 2
    unsigned int                              indexOfLastWeight     = 0;     // step 2
    P_M*                                      unifiedCovMatrix      = NULL;  // step 4
    bool                                      useBalancedChains     = false; // step 6
    BalancedLinkedChainsPerNodeStruct<P_V>* balancedLinkControl   = NULL;  // step 7
    UnbalancedLinkedChainsPerNodeStruct*    unbalancedLinkControl = NULL;  // step 7
    BayesianJointPdf<P_V,P_M>*         currPdf               = NULL;  // step 8
    GenericVectorRV<P_V,P_M>*          currRv                = NULL;  // step 8

    unsigned int exponentEtaTriedAmount = 0;
    while (tryExponentEta) { // gpmsa
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In IMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", beginning 'do-while(tryExponentEta)"
                                << ": failedExponent = " << failedExponent
                                << ", failedEta = "      << failedEta
                                << std::endl;
      }

    //***********************************************************
    // Step 1 of 11: read options
    //***********************************************************
    m_currStep = 1;
    sprintf(levelPrefix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    currOptions = new MLSamplingLevelOptions(m_env,(m_options.m_prefix + levelPrefix).c_str());
    currOptions->scanOptionsValues(&defaultLevelOptions);

    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step01_inter0(currOptions,                     // input
                                     currUnifiedRequestedNumSamples); // output
    }

    //***********************************************************
    // Step 2 of 11: save [chain and corresponding target pdf values] from previous level
    //***********************************************************
    m_currStep = 2;
    double       prevExponent                   = currExponent;
    double       prevEta                        = currEta;
    unsigned int prevUnifiedRequestedNumSamples = currUnifiedRequestedNumSamples;
    prevChain = new SequenceOfVectors<P_V,P_M>(m_vectorSpace,
                                                      0,
                                                      m_options.m_prefix+"prev_chain");
    ScalarSequence<double> prevLogLikelihoodValues(m_env,0,"");
    ScalarSequence<double> prevLogTargetValues    (m_env,0,"");

    indexOfFirstWeight = 0;
    indexOfLastWeight  = 0;

    //std::cout << "m_env.inter0Rank() = " << m_env.inter0Rank() << std::endl;
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step02_inter0(currOptions,             // input
                                     currChain,               // input/output // restate
                                     currLogLikelihoodValues, // input/output // restate
                                     currLogTargetValues,     // input/output // restate
                                     *prevChain,              // output
                                     prevLogLikelihoodValues, // output
                                     prevLogTargetValues,     // output
                                     indexOfFirstWeight,      // output
                                     indexOfLastWeight);      // output
    }

    //***********************************************************
    // Step 3 of 11: compute [currExponent and sequence of weights] for current level
    //               update 'm_logEvidenceFactors'
    //***********************************************************
    m_currStep = 3;
    ScalarSequence<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step03_inter0(currOptions,             // input
                                     prevLogLikelihoodValues, // input
                                     prevExponent,            // input
                                     failedExponent,          // input // gpmsa
                                     currExponent,            // output
                                     weightSequence);         // output
    }

    // All nodes in 'subComm' should have the same 'currExponent'
    m_env.subComm().Bcast((void *) &currExponent, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm', important
                          "MLSampling<P_V,P_M>::generateSequence()",
                          "failed MPI.Bcast() for currExponent");
    m_debugExponent = currExponent;

    if (currExponent == 1.) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": copying 'last' level options to current options" // In all nodes of 'subComm', important
                                << std::endl;
      }
      delete currOptions;
      currOptions = &lastLevelOptions;

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": after copying 'last' level options to current options, the current options are"
                                << "\n" << *currOptions
                                << std::endl;
      }

      if (m_env.inter0Rank() >= 0) {
        // It is necessary to recompute 'currUnifiedRequestedNumSamples' because
        // 'currOptions' has just been replaced by 'lastLevelOptions'
        unsigned int tmpSize = currOptions->m_rawChainSize;
        m_env.inter0Comm().Allreduce((void *) &tmpSize, (void *) &currUnifiedRequestedNumSamples, (int) 1, RawValue_MPI_UNSIGNED, RawValue_MPI_SUM,
                                     "MLSampling<P_V,P_M>::generateSequence()",
                                     "failed MPI.Allreduce() for requested num samples in step 3");
      }
    }

    //***********************************************************
    // Step 4 of 11: create covariance matrix for current level
    //***********************************************************
    m_currStep = 4;
    P_V oneVec(m_vectorSpace.zeroVector());
    oneVec.cwSet(1.);

    unifiedCovMatrix = NULL;
    if (m_env.inter0Rank() >= 0) {
      unifiedCovMatrix = m_vectorSpace.newMatrix();
    }
    else {
      unifiedCovMatrix = new P_M(oneVec);
    }

    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step04_inter0(*prevChain,         // input
                                     weightSequence,     // input
                                     *unifiedCovMatrix); // output
    }

    //***********************************************************
    // Step 5 of 11: create *unified* finite distribution for current level
    //***********************************************************
    m_currStep = 5;
    std::vector<unsigned int> unifiedIndexCountersAtProc0Only(0);
    std::vector<double>       unifiedWeightStdVectorAtProc0Only(0); // KAUST, to check
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step05_inter0(currUnifiedRequestedNumSamples,     // input
                                     weightSequence,                     // input
                                     unifiedIndexCountersAtProc0Only,    // output
                                     unifiedWeightStdVectorAtProc0Only); // output
    }

    //***********************************************************
    // Step 6 of 11: decide on using balanced chains or not
    //***********************************************************
    m_currStep = 6;
    useBalancedChains = false;
    std::vector<ExchangeInfoStruct> exchangeStdVec(0);
    // All processors should call this routine in order to have the same decision value
    generateSequence_Step06_all(currOptions,                     // input
                                indexOfFirstWeight,              // input
                                indexOfLastWeight,               // input
                                unifiedIndexCountersAtProc0Only, // input
                                useBalancedChains,               // output
                                exchangeStdVec);                 // output

    //***********************************************************
    // Step 7 of 11: plan for number of linked chains for each node so that all
    //               nodes generate the closest possible to the same number of positions
    //***********************************************************
    m_currStep = 7;
    balancedLinkControl   = new BalancedLinkedChainsPerNodeStruct<P_V>();
    unbalancedLinkControl = new UnbalancedLinkedChainsPerNodeStruct   ();
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step07_inter0(useBalancedChains,               // input
                                     indexOfFirstWeight,              // input
                                     indexOfLastWeight,               // input
                                     unifiedIndexCountersAtProc0Only, // input
                                     *unbalancedLinkControl,          // (possible) output
                                     currOptions,                     // input
                                     *prevChain,                      // input
                                     exchangeStdVec,                  // (possible) input/output
                                     *balancedLinkControl);           // (possible) output
    }

    // aqui: prevChain might not be needed anymore. Delete it to save memory

    //***********************************************************
    // Step 8 of 11: create vector RV for current level
    //***********************************************************
    m_currStep = 8;
    currPdf = new BayesianJointPdf<P_V,P_M> (m_options.m_prefix.c_str(),
                                                    m_priorRv.pdf(),
                                                    m_likelihoodFunction,
                                                    currExponent,
                                                    *m_targetDomain);

    currRv = new GenericVectorRV<P_V,P_M> (m_options.m_prefix.c_str(),
                                                  *m_targetDomain);

    // All nodes should set 'currRv'
    generateSequence_Step08_all(*currPdf,
                                *currRv);

    //***********************************************************
    // Step 9 of 11: scale unified covariance matrix until min <= rejection rate <= max
    //***********************************************************
    m_currStep = 9;
    generateSequence_Step09_all(*prevChain,                        // input
                                indexOfFirstWeight,                // input
                                indexOfLastWeight,                 // input
                                unifiedWeightStdVectorAtProc0Only, // input
                                weightSequence,                    // input
                                prevEta,                           // input
                                *currRv,                           // input
                                currOptions,                       // input (changed temporarily internally)
                                *unifiedCovMatrix,                 // input/output
                                currEta);                          // output

    tryExponentEta = false; // gpmsa
    if ((currOptions->m_minAcceptableEta > 0.     ) && // gpmsa
        (currEta < currOptions->m_minAcceptableEta)) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", preparing to retry ExponentEta"
                                << ": currExponent = " << currExponent
                                << ", currEta = "      << currEta
                                << std::endl;
      }
      exponentEtaTriedAmount++;
      tryExponentEta = true;
      failedExponent = currExponent;
      failedEta      = currEta;

      // "Return" to previous level
      delete currRv;                 // Step 8
      currRv = NULL;                 // Step 8
      delete currPdf;                // Step 8
      currPdf = NULL;                // Step 8

      delete balancedLinkControl;    // Step 7
      balancedLinkControl   = NULL;  // Step 7
      delete unbalancedLinkControl;  // Step 7
      unbalancedLinkControl = NULL;  // Step 7

      delete unifiedCovMatrix;       // Step 4
      unifiedCovMatrix = NULL;       // Step 4

      currExponent                   = prevExponent;
      currEta                        = 1.; // prevEta;
      currUnifiedRequestedNumSamples = prevUnifiedRequestedNumSamples;

      currChain.clear();             // Step 2
      currChain = (*prevChain);      // Step 2
      delete prevChain;              // Step 2
      prevChain = NULL;              // Step 2

      currLogLikelihoodValues        = prevLogLikelihoodValues;
      currLogTargetValues            = prevLogTargetValues;
    }
    } // while (tryExponentEta) // gpmsa

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) { // gpmsa
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", exited 'do-while(tryExponentEta)"
                              << ", failedExponent = " << failedExponent
                              << ", failedEta = "      << failedEta
                              << std::endl;
    }

    //***********************************************************
    // Step 10 of 11: sample vector RV of current level
    //***********************************************************
    m_currStep = 10;
    // All nodes should call here
    generateSequence_Step10_all(*currOptions,                 // input (changed temporarily internally)
                                *unifiedCovMatrix,            // input
                                *currRv,                      // input
                                useBalancedChains,            // input
                                *unbalancedLinkControl,       // input // Round Rock
                                indexOfFirstWeight,           // input // Round Rock
                                *prevChain,                   // input // Round Rock
                                *balancedLinkControl,         // input // Round Rock
                                currChain,                    // output
                                cumulativeRawChainRunTime,    // output
                                cumulativeRawChainRejections, // output
                                &currLogLikelihoodValues,     // output // likelihood is important
                                &currLogTargetValues);        // output);

    //***********************************************************
    // Perform checkpoint if necessary
    //***********************************************************
    stopAtEndOfLevel = currOptions->m_stopAtEnd;
    bool performCheckpoint = stopAtEndOfLevel;
    if (m_options.m_restartOutput_levelPeriod > 0) {
      performCheckpoint = performCheckpoint || ( ((m_currLevel + 1) % m_options.m_restartOutput_levelPeriod) == 0 );
      if (currExponent == 1.) {
        performCheckpoint = true;
      }
    }
    if (performCheckpoint) {
      checkpointML(currExponent,            // input
                   currEta,                 // input
                   currChain,               // input
                   currLogLikelihoodValues, // input
                   currLogTargetValues);    // input
    }

    //***********************************************************
    // Just free some memory
    //***********************************************************
    {
      delete unifiedCovMatrix;

      for (unsigned int i = 0; i < balancedLinkControl->balLinkedChains.size(); ++i) {
        UQ_FATAL_TEST_MACRO(balancedLinkControl->balLinkedChains[i].initialPosition == NULL,
                            m_env.worldRank(),
                            "MLSampling<P_V,P_M>::generateSequence()",
                            "Initial position pointer in step 9 should not be NULL");
        delete balancedLinkControl->balLinkedChains[i].initialPosition;
        balancedLinkControl->balLinkedChains[i].initialPosition = NULL;
      }
      balancedLinkControl->balLinkedChains.clear();
    }

    //***********************************************************
    // Step 11 of 11: filter chain if requested
    //***********************************************************
    m_currStep = 11;
    unsigned int unifiedNumberOfRejections = 0;
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step11_inter0(currOptions,                      // input
                                     currUnifiedRequestedNumSamples,   // input
                                     cumulativeRawChainRejections,     // input
                                     currChain,                        // input/output
                                     currLogLikelihoodValues,          // input/output
                                     currLogTargetValues,              // input/output
                                     unifiedNumberOfRejections);       // output
    }

    minLogLike = -INFINITY;
    maxLogLike =  INFINITY;
    currLogLikelihoodValues.subMinMaxExtra(0,
                                           currLogLikelihoodValues.subSequenceSize(),
                                           minLogLike,
                                           maxLogLike);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": at end of level "  << m_currLevel+LEVEL_REF_ID
                              << ", sub minLogLike = " << minLogLike
                              << ", sub maxLogLike = " << maxLogLike
                              << std::endl;
    }

    m_env.fullComm().Barrier();

    minLogLike = -INFINITY;
    maxLogLike =  INFINITY;
    currLogLikelihoodValues.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                               0,
                                               currLogLikelihoodValues.subSequenceSize(),
                                               minLogLike,
                                               maxLogLike);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": at end of level "      << m_currLevel+LEVEL_REF_ID
                              << ", unified minLogLike = " << minLogLike
                              << ", unified maxLogLike = " << maxLogLike
                              << std::endl;
    }

    //***********************************************************
    // Prepare to end current level
    //***********************************************************
    double levelRunTime = MiscGetEllapsedSeconds(&timevalLevel);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level "                   << m_currLevel+LEVEL_REF_ID
                              << ", having generated "               << currChain.subSequenceSize()
                              << " chain positions"
                              << ", cumulativeRawChainRunTime = "    << cumulativeRawChainRunTime << " seconds"
                              << ", total level time = "             << levelRunTime              << " seconds"
                              << ", cumulativeRawChainRejections = " << cumulativeRawChainRejections
                              << " (" << 100.*((double) cumulativeRawChainRejections)/((double) currOptions->m_rawChainSize)
                              << "% at this processor)"
                              << " (" << 100.*((double) unifiedNumberOfRejections)/((double) currUnifiedRequestedNumSamples)
                              << "% over all processors)"
                              << ", stopAtEndOfLevel = " << stopAtEndOfLevel
                              << std::endl;
    }

    if (m_env.inter0Rank() >= 0) {
      double minCumulativeRawChainRunTime = 0.;
      m_env.inter0Comm().Allreduce((void *) &cumulativeRawChainRunTime, (void *) &minCumulativeRawChainRunTime, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_MIN,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for min cumulative raw chain run time");

      double maxCumulativeRawChainRunTime = 0.;
      m_env.inter0Comm().Allreduce((void *) &cumulativeRawChainRunTime, (void *) &maxCumulativeRawChainRunTime, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_MAX,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for max cumulative raw chain run time");

      double avgCumulativeRawChainRunTime = 0.;
      m_env.inter0Comm().Allreduce((void *) &cumulativeRawChainRunTime, (void *) &avgCumulativeRawChainRunTime, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for sum cumulative raw chain run time");
      avgCumulativeRawChainRunTime /= ((double) m_env.inter0Comm().NumProc());

      double minLevelRunTime = 0.;
      m_env.inter0Comm().Allreduce((void *) &levelRunTime, (void *) &minLevelRunTime, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_MIN,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for min level run time");

      double maxLevelRunTime = 0.;
      m_env.inter0Comm().Allreduce((void *) &levelRunTime, (void *) &maxLevelRunTime, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_MAX,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for max level run time");

      double avgLevelRunTime = 0.;
      m_env.inter0Comm().Allreduce((void *) &levelRunTime, (void *) &avgLevelRunTime, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for sum level run time");
      avgLevelRunTime /= ((double) m_env.inter0Comm().NumProc());

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level "               << m_currLevel+LEVEL_REF_ID
                                << ": min cumul seconds = " << minCumulativeRawChainRunTime
                                << ", avg cumul seconds = " << avgCumulativeRawChainRunTime
                                << ", max cumul seconds = " << maxCumulativeRawChainRunTime
                                << ", min level seconds = " << minLevelRunTime
                                << ", avg level seconds = " << avgLevelRunTime
                                << ", max level seconds = " << maxLevelRunTime
                                << std::endl;
      }
    }

    if (currExponent != 1.) delete currOptions;

    struct timeval timevalLevelEnd;
    iRC = 0;
    iRC = gettimeofday(&timevalLevelEnd, NULL);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "Getting at the end of level " << m_currLevel+LEVEL_REF_ID
                              << ", as part of a 'while' on levels"
                              << ", at  "   << ctime(&timevalLevelEnd.tv_sec)
                              << ", after " << timevalLevelEnd.tv_sec - timevalRoutineBegin.tv_sec
                              << " seconds from entering the routine"
                              << ", after " << timevalLevelEnd.tv_sec - m_env.timevalBegin().tv_sec
                              << " seconds from queso environment instatiation"
                              << std::endl;
    }
  } // end of level while

  //UQ_FATAL_TEST_MACRO((currExponent < 1.),
  //                    m_env.worldRank(),
  //                    "MLSampling<P_V,P_M>::generateSequence()",
  //                    "exponent has not achieved value '1' even after exiting level loop");

  //***********************************************************
  // Compute information gain
  // ln( \pi(D|M) ) = E[ln( \pi(D|\theta,M) )] - E[ln( \pi(\theta|D,M) / \pi(\theta|M) )]
  //***********************************************************
  if (m_env.inter0Rank() >= 0) { // KAUST
    UQ_FATAL_TEST_MACRO((m_currLevel != m_logEvidenceFactors.size()),
                        m_env.worldRank(),
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "invalid m_currLevel at the exit of the level loop");
    m_logEvidence = 0.;
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      m_logEvidence += m_logEvidenceFactors[i];
    }

#if 1 // prudenci-2012-07-06
    m_meanLogLikelihood = currLogLikelihoodValues.unifiedMeanPlain(m_vectorSpace.numOfProcsForStorage() == 1);
#else
    m_meanLogLikelihood = currLogLikelihoodValues.unifiedMeanExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                                                   0,
                                                                   currLogLikelihoodValues.subSequenceSize());
#endif

    m_eig = m_meanLogLikelihood - m_logEvidence;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ", log(evidence) = "     << m_logEvidence
                              << ", evidence = "          << exp(m_logEvidence)
                              << ", meanLogLikelihood = " << m_meanLogLikelihood
                              << ", eig = "               << m_eig
                              << std::endl;
    }
  }

  m_env.subComm().Bcast((void *) &m_logEvidence, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm'
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "failed MPI.Bcast() for m_logEvidence");

  m_env.subComm().Bcast((void *) &m_meanLogLikelihood, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm'
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "failed MPI.Bcast() for m_meanLogLikelihood");

  m_env.subComm().Bcast((void *) &m_eig, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm'
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "failed MPI.Bcast() for m_eig");

  //***********************************************************
  // Prepare to return
  //***********************************************************
  workingChain.clear();
  workingChain.resizeSequence(currChain.subSequenceSize());
  P_V auxVec(m_vectorSpace.zeroVector());
  for (unsigned int i = 0; i < workingChain.subSequenceSize(); ++i) {
    if (m_env.inter0Rank() >= 0) {
      currChain.getPositionValues(i,auxVec);
    }
    workingChain.setPositionValues(i,auxVec);
  }

  if (workingLogLikelihoodValues) *workingLogLikelihoodValues = currLogLikelihoodValues;
  if (workingLogTargetValues    ) *workingLogTargetValues     = currLogTargetValues;

  struct timeval timevalRoutineEnd;
  iRC = 0;
  iRC = gettimeofday(&timevalRoutineEnd, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence()"
                            << ", at  "   << ctime(&timevalRoutineEnd.tv_sec)
                            << ", after " << timevalRoutineEnd.tv_sec - timevalRoutineBegin.tv_sec
                            << " seconds from entering the routine"
                            << ", after " << timevalRoutineEnd.tv_sec - m_env.timevalBegin().tv_sec
                            << " seconds from queso environment instatiation"
                            << std::endl;
  }

  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const MLSampling<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}

template <class P_V,class P_M>
double MLSampling<P_V,P_M>::logEvidence() const
{
  return m_logEvidence;
}

template <class P_V,class P_M>
double MLSampling<P_V,P_M>::meanLogLikelihood() const
{
  return m_meanLogLikelihood;
}

template <class P_V,class P_M>
double MLSampling<P_V,P_M>::eig() const
{
  return m_eig;
}

}  // End namespace QUESO

#endif // __UQ_MULTI_LEVEL_SAMPLING1_H__
