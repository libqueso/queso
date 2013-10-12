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

#ifndef UQ_MULTI_LEVEL_SAMPLING_H
#define UQ_MULTI_LEVEL_SAMPLING_H

#define ML_NEW_CODE_2009_12_29

#include <queso/MLSamplingOptions.h>
#include <queso/MetropolisHastingsSG1.h>
#include <queso/FiniteDistribution.h>
#include <queso/VectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/MarkovChainPositionData.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/SequenceOfVectors.h>
#include <queso/ArrayOfSequences.h>
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

}  // End namespace QUESO

#endif // UQ_MULTI_LEVEL_SAMPLING_H
