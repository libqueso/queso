//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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

//---------------------------------------------------------

// aqui 1
#ifdef QUESO_HAS_GLPK
struct BIP_routine_struct {
  const uqBaseEnvironmentClass* env;
  unsigned int                  currLevel;
};

void BIP_routine(glp_tree *tree, void *info);
#endif

//---------------------------------------------------------

template <class P_V, class P_M>
struct uqMLStateStruct
{
  uqMLStateStruct(const uqVectorSpaceClass<P_V,P_M>& vectorSpace);
 ~uqMLStateStruct();

  const uqVectorSpaceClass<P_V,P_M>& m_vectorSpace;
  const uqBaseEnvironmentClass&      m_env;
  unsigned int                       m_vectorSpaceDimension;
  int                                m_level; // Yes, 'int', in order to store value '-1'
  std::vector<double>                m_logEvidenceFactors;
  double                             m_exponent;
  double                             m_eta;
  unsigned int                       m_unifiedRequestedNumSamples;
  uqSequenceOfVectorsClass<P_V,P_M>  m_chain;
  uqScalarSequenceClass<double>      m_logLikelihoodValues;
  uqScalarSequenceClass<double>      m_logTargetValues;
};

template<class P_V,class P_M>
uqMLStateStruct<P_V,P_M>::uqMLStateStruct(const uqVectorSpaceClass<P_V,P_M>& vectorSpace)
  :
  m_vectorSpace               (vectorSpace),
  m_env                       (vectorSpace.env()),
  m_vectorSpaceDimension      (0),
  m_level                     (-1), // Yes, '-1'
  m_logEvidenceFactors        (0),
  m_exponent                  (0.),
  m_eta                       (1.),
  m_unifiedRequestedNumSamples(0),
  m_chain                     (m_vectorSpace,0,""),
  m_logLikelihoodValues       (m_env,0,""),
  m_logTargetValues           (m_env,0,"")
{
}

template<class P_V,class P_M>
uqMLStateStruct<P_V,P_M>::~uqMLStateStruct()
{
}

//---------------------------------------------------------

struct uqExchangeInfoStruct
{
  int          originalNodeOfInitialPosition;
  unsigned int originalIndexOfInitialPosition;
  int          finalNodeOfInitialPosition;
  unsigned int numberOfPositions;
};

//---------------------------------------------------------

template <class P_V>
struct uqBalancedLinkedChainControlStruct
{
  P_V*         initialPosition;
  unsigned int numberOfPositions;
};

template <class P_V>
struct uqBalancedLinkedChainsPerNodeStruct
{
  std::vector<uqBalancedLinkedChainControlStruct<P_V> > balLinkedChains;
};

//---------------------------------------------------------

struct uqUnbalancedLinkedChainControlStruct
{
  unsigned int initialPositionIndexInPreviousChain;
  unsigned int numberOfPositions;
};

struct uqUnbalancedLinkedChainsPerNodeStruct
{
  std::vector<uqUnbalancedLinkedChainControlStruct> unbLinkedChains;
};

//---------------------------------------------------------

/*! A templated class that represents a Multi Level sampler.
 */
template <class P_V,class P_M>
class uqMLSamplingClass
{
public:
  /*! Constructor: */
  uqMLSamplingClass(/*! Prefix                  */ const char*                               prefix,
                    /*! The prior rv            */ const uqBaseVectorRVClass      <P_V,P_M>& priorRv,
                    /*! The likelihood function */ const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction);
  /*! Destructor: */
 ~uqMLSamplingClass();

  /*! Operation to generate the chain */
  void   generateSequence(uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
                          uqScalarSequenceClass<double>*      workingLogLikelihoodValues,
                          uqScalarSequenceClass<double>*      workingLogTargetValues);

  void   print           (std::ostream& os) const;

private:
  // Methods available at uqMLSampling2.h
  void   restartML                     (uqMLStateStruct<P_V,P_M>&                       restartPrevState,                   // output
                                        uqMLStateStruct<P_V,P_M>&                       restartCurrState);                  // output

  void   checkpointML                  (double                                          currExponent,                       // input
                                        double                                          currEta,                            // input
                                        unsigned int                                    currUnifiedRequestedNumSamples,     // input
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        currChain,                          // input
                                        const uqScalarSequenceClass<double>&            currLogLikelihoodValues,            // input
                                        const uqScalarSequenceClass<double>&            currLogTargetValues);               // input

  void   checkpointML                  (double                                          prevExponent,                       // input
                                        double                                          prevEta,                            // input
                                        unsigned int                                    prevUnifiedRequestedNumSamples,     // input
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input
                                        const uqScalarSequenceClass<double>&            prevLogLikelihoodValues,            // input
                                        const uqScalarSequenceClass<double>&            prevLogTargetValues,                // input
                                        double                                          currExponent,                       // input
                                        double                                          currEta,                            // input
                                        unsigned int                                    currUnifiedRequestedNumSamples,     // input
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        currChain,                          // input
                                        const uqScalarSequenceClass<double>&            currLogLikelihoodValues,            // input
                                        const uqScalarSequenceClass<double>&            currLogTargetValues);               // input

  void   generateSequence_Level0_all   (const uqMLSamplingLevelOptionsClass&            defaultLevelOptions,                // input
                                        unsigned int&                                   unifiedRequestedNumSamples,         // output
                                        uqSequenceOfVectorsClass<P_V,P_M>&              currChain,                          // output
                                        uqScalarSequenceClass<double>&                  currLogLikelihoodValues,            // output
                                        uqScalarSequenceClass<double>&                  currLogTargetValues);               // output

  void   generateSequence_Step01_inter0(const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        unsigned int&                                   unifiedRequestedNumSamples);        // output

  void   generateSequence_Step02_inter0(const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        uqSequenceOfVectorsClass<P_V,P_M>&              currChain,                          // input/output
                                        uqScalarSequenceClass<double>&                  currLogLikelihoodValues,            // input/output
                                        uqScalarSequenceClass<double>&                  currLogTargetValues,                // input/output
                                        uqSequenceOfVectorsClass<P_V,P_M>&              prevChain,                          // output
                                        uqScalarSequenceClass<double>&                  prevLogLikelihoodValues,            // output
                                        uqScalarSequenceClass<double>&                  prevLogTargetValues,                // output
                                        unsigned int&                                   indexOfFirstWeight,                 // output
                                        unsigned int&                                   indexOfLastWeight);                 // output

  void   generateSequence_Step03_inter0(const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        const uqScalarSequenceClass<double>&            prevLogLikelihoodValues,            // input
                                        double                                          prevExponent,                       // input
                                        double&                                         currExponent,                       // output
                                        uqScalarSequenceClass<double>&                  weightSequence);                    // output

  void   generateSequence_Step04_inter0(const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input
                                        const uqScalarSequenceClass<double>&            weightSequence,                     // input
                                        P_M&                                            unifiedCovMatrix);                  // output

  void   generateSequence_Step05_inter0(unsigned int                                    unifiedRequestedNumSamples,         // input
                                        const uqScalarSequenceClass<double>&            weightSequence,                     // input
                                        std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only,    // output
                                        std::vector<double>&                            unifiedWeightStdVectorAtProc0Only); // output

  void   generateSequence_Step06_all   (const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        bool&                                           useBalancedChains,                  // output
                                        std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                    // output

  void   generateSequence_Step07_inter0(bool                                            useBalancedChains,                  // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        uqUnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl,              // (possible) output
                                        const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input
                                        std::vector<uqExchangeInfoStruct>&              exchangeStdVec,                     // (possible) input/output
                                        uqBalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // (possible) output

  void   generateSequence_Step08_all   (uqBayesianJointPdfClass<P_V,P_M>&               currPdf,                            // input/output
                                        uqGenericVectorRVClass<P_V,P_M>&                currRv);                            // output

  void   generateSequence_Step09_all   (const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only,  // input
                                        const uqScalarSequenceClass<double>&            weightSequence,                     // input
                                        double                                          prevEta,                            // input
                                        const uqGenericVectorRVClass<P_V,P_M>&          currRv,                             // input
                                        uqMLSamplingLevelOptionsClass*                  currOptions,                        // input (changed temporarily internally)
                                        P_M&                                            unifiedCovMatrix,                   // input/output
                                        double&                                         currEta);                           // output

  void   generateSequence_Step10_all   (uqMLSamplingLevelOptionsClass&                  currOptions,                        // input (changed temporarily internally)
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const uqGenericVectorRVClass  <P_V,P_M>&        currRv,                             // input
                                        bool                                            useBalancedChains,                  // input
                                        const uqUnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,              // input // Round Rock
                                        unsigned int                                    indexOfFirstWeight,                 // input // Round Rock
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input // Round Rock
                                        const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,                // input // Round Rock
                                        uqSequenceOfVectorsClass      <P_V,P_M>&        currChain,                          // output
                                        double&                                         cumulativeRawChainRunTime,          // output
                                        unsigned int&                                   cumulativeRawChainRejections,       // output
                                        uqScalarSequenceClass         <double>*         currLogLikelihoodValues,            // output
                                        uqScalarSequenceClass         <double>*         currLogTargetValues);               // output

  void   generateSequence_Step11_inter0(const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        unsigned int                                    unifiedRequestedNumSamples,         // input
                                        unsigned int                                    cumulativeRawChainRejections,       // input
                                        uqSequenceOfVectorsClass<P_V,P_M>&              currChain,                          // input/output
                                        uqScalarSequenceClass<double>&                  currLogLikelihoodValues,            // input/output
                                        uqScalarSequenceClass<double>&                  currLogTargetValues,                // input/output
                                        unsigned int&                                   unifiedNumberOfRejections);         // output

  // Methods available at uqMLSampling3.h
  void   sampleIndexes_proc0           (unsigned int                                    unifiedRequestedNumSamples,         // input
                                        const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only,  // input
                                        std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only);   // output

  bool   decideOnBalancedChains_all    (const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                    // output

  void   prepareBalLinkedChains_inter0 (const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input
                                        std::vector<uqExchangeInfoStruct>&              exchangeStdVec,                     // input/output
                                        uqBalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // output

  void   prepareUnbLinkedChains_inter0 (unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        uqUnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl);             // output

  void   generateBalLinkedChains_all   (uqMLSamplingLevelOptionsClass&                  inputOptions,                       // input, only m_rawChainSize changes
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const uqGenericVectorRVClass  <P_V,P_M>&        rv,                                 // input
                                        const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,                // input // Round Rock
                                        uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,                       // output
                                        double&                                         cumulativeRunTime,                  // output
                                        unsigned int&                                   cumulativeRejections,               // output
                                        uqScalarSequenceClass         <double>*         currLogLikelihoodValues,            // output
                                        uqScalarSequenceClass         <double>*         currLogTargetValues);               // output

  void   generateUnbLinkedChains_all   (uqMLSamplingLevelOptionsClass&                  inputOptions,                       // input, only m_rawChainSize changes
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const uqGenericVectorRVClass  <P_V,P_M>&        rv,                                 // input
                                        const uqUnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,              // input // Round Rock
                                        unsigned int                                    indexOfFirstWeight,                 // input // Round Rock
                                        const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input // Round Rock
                                        uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,                       // output
                                        double&                                         cumulativeRunTime,                  // output
                                        unsigned int&                                   cumulativeRejections,               // output
                                        uqScalarSequenceClass         <double>*         currLogLikelihoodValues,            // output
                                        uqScalarSequenceClass         <double>*         currLogTargetValues);               // output

#ifdef QUESO_HAS_GLPK
  void   solveBIP_proc0                (std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                    // input/output
#endif

  void   justBalance_proc0             (const uqMLSamplingLevelOptionsClass*            currOptions,                        // input
                                        std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                    // input/output

  void   mpiExchangePositions_inter0   (const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                          // input
                                        const std::vector<uqExchangeInfoStruct>&        exchangeStdVec,                     // input
                                        const std::vector<unsigned int>&                finalNumChainsPerNode,              // input
                                        const std::vector<unsigned int>&                finalNumPositionsPerNode,           // input
                                        uqBalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // output

  // Private variables
  const uqBaseEnvironmentClass&             m_env;
  const uqBaseVectorRVClass      <P_V,P_M>& m_priorRv;
  const uqBaseScalarFunctionClass<P_V,P_M>& m_likelihoodFunction;
  const uqVectorSpaceClass       <P_V,P_M>& m_vectorSpace;
        uqVectorSetClass         <P_V,P_M>* m_targetDomain;

        uqMLSamplingOptionsClass            m_options;

        unsigned int                        m_currLevel;
        unsigned int                        m_currStep;
        double                              m_debugExponent;
	std::vector<double>                 m_logEvidenceFactors;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj);

#include <uqMLSampling3.h>
#include <uqMLSampling2.h>

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::uqMLSamplingClass(
  const char*                               prefix,
  const uqBaseVectorRVClass      <P_V,P_M>& priorRv,            
  const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction)
  :
  m_env               (priorRv.env()),
  m_priorRv           (priorRv),
  m_likelihoodFunction(likelihoodFunction),
  m_vectorSpace       (m_priorRv.imageSet().vectorSpace()),
  m_targetDomain      (uqInstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet())),
  m_options           (m_env,prefix),
  m_currLevel         (0),
  m_currStep          (0),
  m_debugExponent     (0.),
  m_logEvidenceFactors(0)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::constructor()"
                            << std::endl;
  }

  m_options.scanOptionsValues();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::constructor()"
                            << std::endl;
  }
}

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::~uqMLSamplingClass()
{
  if (m_targetDomain) delete m_targetDomain;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  uqScalarSequenceClass<double>*      workingLogLikelihoodValues,
  uqScalarSequenceClass<double>*      workingLogTargetValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateSequence()..."
                            << std::endl;
  }

  // ernesto
  if (m_options.m_restartInputFileName != ".") {
    workingChain.unifiedReadContents(m_options.m_restartInputFileName,
                                     m_options.m_restartInputFileType,
                                     m_options.m_restartChainSize);

    uqMLSamplingLevelOptionsClass defaultLevelOptions(m_env,(m_options.m_prefix + "default_").c_str());
    defaultLevelOptions.scanOptionsValues(NULL);

    uqMLSamplingLevelOptionsClass lastLevelOptions(m_env,(m_options.m_prefix + "last_").c_str());
    lastLevelOptions.scanOptionsValues(&defaultLevelOptions);

        if (lastLevelOptions.m_rawChainComputeStats) {
          uqFilePtrSetStruct filePtrSet;
          m_env.openOutputFile(lastLevelOptions.m_dataOutputFileName,
                               UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                               lastLevelOptions.m_dataOutputAllowedSet,
                               false,
                               filePtrSet);

          workingChain.computeStatistics(*lastLevelOptions.m_rawChainStatisticalOptionsObj,
                                         filePtrSet.ofsVar);

          m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
        }

        if (lastLevelOptions.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
          workingChain.unifiedWriteContents              (lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType); // KAUST5
          //currLogLikelihoodValues.unifiedWriteContents(lastLevelOptions.m_rawChainDataOutputFileName);
          //currLogTargetValues.unifiedWriteContents    (lastLevelOptions.m_rawChainDataOutputFileName);
        }

        if (lastLevelOptions.m_filteredChainGenerate) {
          uqFilePtrSetStruct filePtrSet;
          m_env.openOutputFile(lastLevelOptions.m_dataOutputFileName,
                               UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                               lastLevelOptions.m_dataOutputAllowedSet,
                               false,
                               filePtrSet);

          unsigned int filterInitialPos = (unsigned int) (lastLevelOptions.m_filteredChainDiscardedPortion * (double) workingChain.subSequenceSize());
          unsigned int filterSpacing    = lastLevelOptions.m_filteredChainLag;
          if (filterSpacing == 0) {
            workingChain.computeFilterParams(*lastLevelOptions.m_filteredChainStatisticalOptionsObj,
                                             filePtrSet.ofsVar,
                                             filterInitialPos,
                                             filterSpacing);
          }

          // Filter positions from the converged portion of the chain
          workingChain.filter(filterInitialPos,
                           filterSpacing);
          workingChain.setName(lastLevelOptions.m_prefix + "filtChain");

          //currLogLikelihoodValues.filter(filterInitialPos,
          //                               filterSpacing);
          //currLogLikelihoodValues.setName(lastLevelOptions.m_prefix + "filtLogLikelihood");

          //currLogTargetValues.filter(filterInitialPos,
          //                           filterSpacing);
          //currLogTargetValues.setName(lastLevelOptions.m_prefix + "filtLogTarget");

          if (lastLevelOptions.m_filteredChainComputeStats) {
            workingChain.computeStatistics(*lastLevelOptions.m_filteredChainStatisticalOptionsObj,
                                           filePtrSet.ofsVar);
          }

          m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);

          if (lastLevelOptions.m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
            workingChain.unifiedWriteContents              (lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
            //currLogLikelihoodValues.unifiedWriteContents(lastLevelOptions.m_filteredChainDataOutputFileName);
            //currLogTargetValues.unifiedWriteContents    (lastLevelOptions.m_filteredChainDataOutputFileName);
          }
        } // if (lastLevelOptions.m_filteredChainGenerate)
  }
  else {
  //***********************************************************
  // Declaration of Variables
  //***********************************************************

  double                            currExponent = 0.;
  double                            currEta      = 1.;
  unsigned int                      currUnifiedRequestedNumSamples = 0;
  uqSequenceOfVectorsClass<P_V,P_M> currChain(m_vectorSpace,
                                              0,
                                              m_options.m_prefix+"curr_chain");
  uqScalarSequenceClass<double>     currLogLikelihoodValues(m_env,0,"");
  uqScalarSequenceClass<double>     currLogTargetValues    (m_env,0,"");

  //***********************************************************
  // Take care of first level (level '0')
  //***********************************************************
  uqMLSamplingLevelOptionsClass defaultLevelOptions(m_env,(m_options.m_prefix + "default_").c_str());
  defaultLevelOptions.scanOptionsValues(NULL);

  uqMLSamplingLevelOptionsClass lastLevelOptions(m_env,(m_options.m_prefix + "last_").c_str());
  lastLevelOptions.scanOptionsValues(&defaultLevelOptions);

  uqMLStateStruct<P_V,P_M> restartPrevState(m_vectorSpace);
  uqMLStateStruct<P_V,P_M> restartCurrState(m_vectorSpace);

  if (m_options.m_restartInputFileName != ".") {
    restartML(restartPrevState,  // output
              restartCurrState); // output

    UQ_FATAL_TEST_MACRO(restartCurrState.m_level == -1,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "restartCurrState.m_level should not be '-1'");

    if (restartCurrState.m_level > 0) {
      UQ_FATAL_TEST_MACRO(restartPrevState.m_level == -1,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "restartPrevState.m_level should not be '-1'");
    }

    if (restartCurrState.m_exponent < 1.) {
      // Ok
      //m_vectorSpaceDimension
      m_currLevel                    = 0;
      //m_logEvidenceFactors           = 0;
      currExponent                   = 0;
      currEta                        = 0;
      currUnifiedRequestedNumSamples = 0;
      //currChain                      = 0;
      //currLogLikelihoodValues        = 0;
      //currLogTargetValues            = 0;
    }
    else if (restartCurrState.m_exponent == 1.) {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "restart exponent == 1 is not supported yet");
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "restart exponent > 1 should not happen");
    }
  }
  else {
    generateSequence_Level0_all(defaultLevelOptions,            // input
                                currUnifiedRequestedNumSamples, // output
                                currChain,                      // output
                                currLogLikelihoodValues,        // output
                                currLogTargetValues);           // output
  }
  //std::cout << "In QUESO: end of level 0. Exiting on purpose" << std::endl;
  //exit(1);

  //***********************************************************
  // Take care of next levels
  //***********************************************************
  bool stopAtEndOfLevel = false;
  char tmpSufix[256];
  do { // begin level while
    m_currLevel++;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level " << m_currLevel+LEVEL_REF_ID
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);
    double       cumulativeRawChainRunTime    = 0.;
    unsigned int cumulativeRawChainRejections = 0;

    //***********************************************************
    // Step 1 of 11: read options
    //***********************************************************
    m_currStep = 1;
    sprintf(tmpSufix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass* currOptions = new uqMLSamplingLevelOptionsClass(m_env,(m_options.m_prefix + tmpSufix).c_str());
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
    uqSequenceOfVectorsClass<P_V,P_M> prevChain(m_vectorSpace,
                                                0,
                                                m_options.m_prefix+"prev_chain");
    uqScalarSequenceClass<double> prevLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> prevLogTargetValues    (m_env,0,"");

    unsigned int indexOfFirstWeight = 0;
    unsigned int indexOfLastWeight  = 0;

    //std::cout << "m_env.inter0Rank() = " << m_env.inter0Rank() << std::endl;
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step02_inter0(currOptions,             // input
                                     currChain,               // input/output
                                     currLogLikelihoodValues, // input/output
                                     currLogTargetValues,     // input/output
                                     prevChain,               // output
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
    uqScalarSequenceClass<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step03_inter0(currOptions,             // input
                                     prevLogLikelihoodValues, // input
                                     prevExponent,            // input
                                     currExponent,            // output
                                     weightSequence);         // output
    }

    // All nodes in 'subComm' should have the same 'currExponent'
    int mpiRC = MPI_Bcast((void *) &currExponent, (int) 1, MPI_DOUBLE, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "failed MPI_Bcast() for currExponent");
    m_debugExponent = currExponent;

    if (currExponent == 1.) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": copying 'last' level options to current options" // In all nodes of 'subComm', important
                                << std::endl;
      }
      delete currOptions;
      currOptions = &lastLevelOptions;

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
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
        mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &currUnifiedRequestedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for requested num samples in step 3");
      }
    }

    //***********************************************************
    // Step 4 of 11: create covariance matrix for current level
    //***********************************************************
    m_currStep = 4;
    P_V oneVec(m_vectorSpace.zeroVector());
    oneVec.cwSet(1.);

    P_M* unifiedCovMatrix = NULL;
    if (m_env.inter0Rank() >= 0) {
      unifiedCovMatrix = m_vectorSpace.newMatrix();
    }
    else {
      unifiedCovMatrix = new P_M(oneVec);
    }

    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step04_inter0(prevChain,          // input
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
    bool useBalancedChains = false;
    std::vector<uqExchangeInfoStruct> exchangeStdVec(0);
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
    uqBalancedLinkedChainsPerNodeStruct<P_V> balancedLinkControl;
    uqUnbalancedLinkedChainsPerNodeStruct    unbalancedLinkControl;
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step07_inter0(useBalancedChains,               // input
                                     indexOfFirstWeight,              // input
                                     indexOfLastWeight,               // input
                                     unifiedIndexCountersAtProc0Only, // input
                                     unbalancedLinkControl,           // (possible) output
                                     currOptions,                     // input
                                     prevChain,                       // input
                                     exchangeStdVec,                  // (possible) input/output
                                     balancedLinkControl);            // (possible) output
    }

    // aqui: prevChain might not be needed anymore. Delete it to save memory

    //***********************************************************
    // Step 8 of 11: create vector RV for current level
    //***********************************************************
    m_currStep = 8;
    uqBayesianJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                             m_priorRv.pdf(),
                                             m_likelihoodFunction,
                                             currExponent,
                                             *m_targetDomain);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           *m_targetDomain);

    // All nodes should set 'currRv'
    generateSequence_Step08_all(currPdf,
                                currRv);

    //***********************************************************
    // Step 9 of 11: scale unified covariance matrix until min <= rejection rate <= max
    //***********************************************************
    m_currStep = 9;
    generateSequence_Step09_all(prevChain,                         // input
                                indexOfFirstWeight,                // input
                                indexOfLastWeight,                 // input
                                unifiedWeightStdVectorAtProc0Only, // input
                                weightSequence,                    // input
                                prevEta,                           // input
                                currRv,                            // input
                                currOptions,                       // input (changed temporarily internally)
                                *unifiedCovMatrix,                 // input/output
                                currEta);                          // output

    //***********************************************************
    // Step 10 of 11: sample vector RV of current level
    //***********************************************************
    m_currStep = 10;
    // All nodes should call here
    generateSequence_Step10_all(*currOptions,                 // input (changed temporarily internally)
                                *unifiedCovMatrix,            // input
                                currRv,                       // input
                                useBalancedChains,            // input
                                unbalancedLinkControl,        // input // Round Rock
                                indexOfFirstWeight,           // input // Round Rock
                                prevChain,                    // input // Round Rock
                                balancedLinkControl,          // input // Round Rock
                                currChain,                    // output
                                cumulativeRawChainRunTime,    // output
                                cumulativeRawChainRejections, // output
                                &currLogLikelihoodValues,     // output // likelihood is important
                                &currLogTargetValues);        // output);

    if (currOptions->m_checkpointOutputFileName != ".") {
      checkpointML(prevExponent,                   // input
                   prevEta,                        // input
                   prevUnifiedRequestedNumSamples, // input
                   prevChain,                      // input
                   prevLogLikelihoodValues,        // input
                   prevLogTargetValues,            // input
                   currExponent,                   // input
                   currEta,                        // input
                   currUnifiedRequestedNumSamples, // input
                   currChain,                      // input
                   currLogLikelihoodValues,        // input
                   currLogTargetValues);           // input
    }

    //***********************************************************
    // Just free some memory
    //***********************************************************
    {
      delete unifiedCovMatrix;

      for (unsigned int i = 0; i < balancedLinkControl.balLinkedChains.size(); ++i) {
        UQ_FATAL_TEST_MACRO(balancedLinkControl.balLinkedChains[i].initialPosition == NULL,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "Initial position pointer in step 9 should not be NULL");
        delete balancedLinkControl.balLinkedChains[i].initialPosition;
        balancedLinkControl.balLinkedChains[i].initialPosition = NULL;
      }
      balancedLinkControl.balLinkedChains.clear();
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

    //***********************************************************
    // Prepare to end current level
    //***********************************************************
    stopAtEndOfLevel = currOptions->m_stopAtEnd;
    double levelRunTime = uqMiscGetEllapsedSeconds(&timevalLevel);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
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
      mpiRC = MPI_Allreduce((void *) &cumulativeRawChainRunTime, (void *) &minCumulativeRawChainRunTime, (int) 1, MPI_DOUBLE, MPI_MIN, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for min cumulative raw chain run time");

      double maxCumulativeRawChainRunTime = 0.;
      mpiRC = MPI_Allreduce((void *) &cumulativeRawChainRunTime, (void *) &maxCumulativeRawChainRunTime, (int) 1, MPI_DOUBLE, MPI_MAX, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for max cumulative raw chain run time");

      double avgCumulativeRawChainRunTime = 0.;
      mpiRC = MPI_Allreduce((void *) &cumulativeRawChainRunTime, (void *) &avgCumulativeRawChainRunTime, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for sum cumulative raw chain run time");
      avgCumulativeRawChainRunTime /= ((double) m_env.inter0Comm().NumProc());

      double minLevelRunTime = 0.;
      mpiRC = MPI_Allreduce((void *) &levelRunTime, (void *) &minLevelRunTime, (int) 1, MPI_DOUBLE, MPI_MIN, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for min level run time");

      double maxLevelRunTime = 0.;
      mpiRC = MPI_Allreduce((void *) &levelRunTime, (void *) &maxLevelRunTime, (int) 1, MPI_DOUBLE, MPI_MAX, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for max level run time");

      double avgLevelRunTime = 0.;
      mpiRC = MPI_Allreduce((void *) &levelRunTime, (void *) &avgLevelRunTime, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for sum level run time");
      avgLevelRunTime /= ((double) m_env.inter0Comm().NumProc());

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
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
  } while ((currExponent     <  1.   ) &&
           (stopAtEndOfLevel == false)); // end of level while

  UQ_FATAL_TEST_MACRO((currExponent < 1.),
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                      "exponent has not achieved value '1' even after exiting level loop");

  if (m_env.inter0Rank() >= 0) { // KAUST
    UQ_FATAL_TEST_MACRO((m_currLevel != m_logEvidenceFactors.size()),
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "invalid m_currLevel at the exit of the level loop");
    double logEvidence = 0.;
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      logEvidence += m_logEvidenceFactors[i];
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ", log(evidence) = " << logEvidence
                              << ", evidence = "      << exp(logEvidence)
                              << std::endl;
    }
  }

  //***********************************************************
  // Prepare to return
  //***********************************************************
  workingChain.clear();
  workingChain.resizeSequence(currChain.subSequenceSize());
  P_V auxVec(m_vectorSpace.zeroVector());
  for (unsigned int i = 0; i < workingChain.subSequenceSize(); ++i) {
    currChain.getPositionValues(i,auxVec);
    workingChain.setPositionValues(i,auxVec);
  }

  if (workingLogLikelihoodValues) *workingLogLikelihoodValues = currLogLikelihoodValues;
  if (workingLogTargetValues    ) *workingLogTargetValues     = currLogTargetValues;
  } // ernesto

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateSequence()"
                            << std::endl;
  }

  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING1_H__
