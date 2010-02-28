/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyrightg (C) 2008 The PECOS Development Team
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

#ifndef __UQ_MULTI_LEVEL_SAMPLING_H__
#define __UQ_MULTI_LEVEL_SAMPLING_H__

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
#include <glpk.h>
#include <sys/time.h>
#include <fstream>

//---------------------------------------------------------

// aqui 1
struct BIP_routine_struct {
  const uqBaseEnvironmentClass* env;
  unsigned int                  currLevel;
};

void BIP_routine(glp_tree *tree, void *info);

struct uqExchangeInfoStruct
{
  int          originalNodeOfInitialPosition;
  unsigned int originalIndexOfInitialPosition;
  int          finalNodeOfInitialPosition;
  unsigned int numberOfPositions;
};

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
  void   sampleIndexesAtProc0   (unsigned int                                    unifiedRequestedNumSamples,        // input
                                 const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only, // input
                                 std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only);  // output

  bool   decideOnBalancedChains (uqMLSamplingLevelOptionsClass*                  currOptions,                       // input
                                 unsigned int                                    indexOfFirstWeight,                // input
                                 unsigned int                                    indexOfLastWeight,                 // input
                                 const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,   // input
                                 std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                   // output

  void   prepareBalLinkedChains (uqMLSamplingLevelOptionsClass*                  currOptions,                       // input
                                 const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                         // input
                                 std::vector<uqExchangeInfoStruct>&              exchangeStdVec,                    // input/output
                                 uqBalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);              // output

  void   prepareUnbLinkedChains (unsigned int                                    indexOfFirstWeight,                // input
                                 unsigned int                                    indexOfLastWeight,                 // input
                                 const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,   // input
                                 uqUnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl);            // output

  void   generateBalLinkedChains(uqMLSamplingLevelOptionsClass&                  inputOptions,                      // input, only m_rawChainSize changes
                                 const P_M&                                      unifiedCovMatrix,                  // input
                                 const uqGenericVectorRVClass  <P_V,P_M>&        rv,                                // input
                                 const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,               // input // Round Rock
                                 uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,                      // output
                                 double&                                         cumulativeRunTime,                 // output
                                 unsigned int&                                   cumulativeRejections,              // output
                                 uqScalarSequenceClass         <double>*         currLogLikelihoodValues,           // output
                                 uqScalarSequenceClass         <double>*         currLogTargetValues);              // output

  void   generateUnbLinkedChains(uqMLSamplingLevelOptionsClass&                  inputOptions,                      // input, only m_rawChainSize changes
                                 const P_M&                                      unifiedCovMatrix,                  // input
                                 const uqGenericVectorRVClass  <P_V,P_M>&        rv,                                // input
                                 const uqUnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,             // input // Round Rock
                                 unsigned int                                    indexOfFirstWeight,                // input // Round Rock
                                 const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                         // input // Round Rock
                                 uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,                      // output
                                 double&                                         cumulativeRunTime,                 // output
                                 unsigned int&                                   cumulativeRejections,              // output
                                 uqScalarSequenceClass         <double>*         currLogLikelihoodValues,           // output
                                 uqScalarSequenceClass         <double>*         currLogTargetValues);              // output

  void   solveBIPAtProc0        (std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                   // input/output

  void   justBalanceAtProc0     (uqMLSamplingLevelOptionsClass*                  currOptions,                       // input
                                 std::vector<uqExchangeInfoStruct>&              exchangeStdVec);                   // input/output

  void   mpiExchangePositions   (const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                         // input
                                 const std::vector<uqExchangeInfoStruct>&        exchangeStdVec,                    // input
                                 const std::vector<unsigned int>&                finalNumChainsPerNode,             // input
                                 const std::vector<unsigned int>&                finalNumPositionsPerNode,          // input
                                 uqBalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);              // output

  void   generateSequence_Step8 (const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                         // input
                                 unsigned int                                    indexOfFirstWeight,                // input
                                 unsigned int                                    indexOfLastWeight,                 // input
                                 const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only, // input
                                 const uqScalarSequenceClass<double>&            weightSequence,                    // input
                                 double                                          prevEta,                           // input
                                 const uqGenericVectorRVClass<P_V,P_M>&          currRv,                            // input
                                 uqMLSamplingLevelOptionsClass*                  currOptions,                       // input (changed temporarily internally)
                                 P_M&                                            unifiedCovMatrix,                  // input/output
                                 double&                                         currEta);                          // output

  const uqBaseEnvironmentClass&             m_env;
  const uqBaseVectorRVClass      <P_V,P_M>& m_priorRv;
  const uqBaseScalarFunctionClass<P_V,P_M>& m_likelihoodFunction;
  const uqVectorSpaceClass       <P_V,P_M>& m_vectorSpace;
        uqVectorSetClass         <P_V,P_M>* m_targetDomain;

        uqMLSamplingOptionsClass            m_options;

        unsigned int                        m_currLevel;
        unsigned int                        m_currStep;
	std::vector<double>                 m_logEvidenceFactors;
        double                              m_logEvidence;
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
  m_logEvidenceFactors(0),
  m_logEvidence       (0.)
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

  //***********************************************************
  // Declaration of Variables
  //***********************************************************
  unsigned int unifiedRequestedNumSamples = 0;

  double                            currExponent = 0.;
  double                            currEta      = 1.;
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

  char tmpSufix[256];

  {
    sprintf(tmpSufix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass currOptions(m_env,(m_options.m_prefix + tmpSufix).c_str());
    currOptions.scanOptionsValues(&defaultLevelOptions);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level "              << m_currLevel+LEVEL_REF_ID
                              << ", currOptions.m_rawChainSize = " << currOptions.m_rawChainSize // Ok to use rawChainSize
                              << ", currExponent = "               << currExponent
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);

    if (m_env.inter0Rank() >= 0) {
      unsigned int tmpSize = currOptions.m_rawChainSize;
      int mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &unifiedRequestedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for requested num samples in level 0");
    }
    else {
      unifiedRequestedNumSamples = currOptions.m_rawChainSize;
    }

    currChain.setName              (currOptions.m_prefix + "rawChain");
    currLogLikelihoodValues.setName(currOptions.m_prefix + "rawLogLikelihood");
    currLogTargetValues.setName    (currOptions.m_prefix + "rawLogTarget");

    currChain.resizeSequence              (currOptions.m_rawChainSize); // Ok to use rawChainSize
    currLogLikelihoodValues.resizeSequence(currOptions.m_rawChainSize); // Ok to use rawChainSize
    currLogTargetValues.resizeSequence    (currOptions.m_rawChainSize); // Ok to use rawChainSize

    P_V auxVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 0; i < currChain.subSequenceSize(); ++i) {
      //std::cout << "In QUESO: before prior realizer with i = " << i << std::endl;
      m_priorRv.realizer().realization(auxVec);
      currChain.setPositionValues(i,auxVec);
      // KAUST: all nodes should call here
      currLogLikelihoodValues[i] = m_likelihoodFunction.lnValue(auxVec,NULL,NULL,NULL,NULL);  // likelihood is important
      currLogTargetValues[i]     = m_priorRv.pdf().lnValue(auxVec,NULL,NULL,NULL,NULL) + currLogLikelihoodValues[i];
      //std::cout << "In QUESO: currLogTargetValues[" << i << "] = " << currLogTargetValues[i] << std::endl;
    }

    if (m_env.inter0Rank() >= 0) { // KAUST
      if (currOptions.m_rawChainComputeStats) {
        std::ofstream* genericOfsVar = NULL;
        m_env.openOutputFile(currOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                             currOptions.m_dataOutputAllowedSet,
                             false,
                             genericOfsVar);

        currChain.computeStatistics(*currOptions.m_rawChainStatisticalOptions,
                                    genericOfsVar);

        //genericOfsVar->close();
        delete genericOfsVar;
      }

      if (currOptions.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
        currChain.unifiedWriteContents              (currOptions.m_rawChainDataOutputFileName);
        currLogLikelihoodValues.unifiedWriteContents(currOptions.m_rawChainDataOutputFileName);
        currLogTargetValues.unifiedWriteContents    (currOptions.m_rawChainDataOutputFileName);
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;

        //unsigned int numZeros = 0;
        //for (unsigned int i = 0; i < currTargetValues.subSequenceSize(); ++i) {
        //  *m_env.subDisplayFile() << "currTargetValues[" << i
        //                          << "] = " << currTargetValues[i]
        //                          << std::endl;
        //  if (currTargetValues[i] == 0.) numZeros++;
        //}
        //*m_env.subDisplayFile() << "Number of zeros in currTargetValues = " << numZeros
        //                        << std::endl;
      }
    } // KAUST

    UQ_FATAL_TEST_MACRO((currChain.subSequenceSize() != currOptions.m_rawChainSize), // Ok to use rawChainSize
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "currChain (first one) has been generated with invalid size");
    double levelRunTime = uqMiscGetEllapsedSeconds(&timevalLevel);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << m_currLevel+LEVEL_REF_ID
                              << " after " << levelRunTime << " seconds"
                              << std::endl;
    }
  }

  //std::cout << "In QUESO: end of level 0. Exiting on purpose" << std::endl;
  //exit(1);

  //***********************************************************
  // Take care of next levels
  //***********************************************************
  while (currExponent < 1.) {
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
    // Step 1 of 9: read options
    //***********************************************************
    m_currStep = 1;
    sprintf(tmpSufix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass* currOptions = new uqMLSamplingLevelOptionsClass(m_env,(m_options.m_prefix + tmpSufix).c_str());
    currOptions->scanOptionsValues(&defaultLevelOptions);
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 1 of 9"
                                << std::endl;
      }

      unsigned int tmpSize = currOptions->m_rawChainSize;
      // This computed 'unifiedRequestedNumSamples' needs to be recomputed only at the last
      // level, when 'currOptions' is replaced by 'lastLevelOptions' (see step 3 of 9)
      int mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &unifiedRequestedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Allreduce() for requested num samples in step 1");

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", currOptions->m_rawChainSize = " << currOptions->m_rawChainSize // Ok to use rawChainSize
                                << std::endl;
      }
    }

    //***********************************************************
    // Step 2 of 9: save [chain and corresponding target pdf values] from previous level
    //***********************************************************
    m_currStep++;
    double prevExponent = currExponent;
    double prevEta      = currEta;
    uqSequenceOfVectorsClass<P_V,P_M> prevChain(m_vectorSpace,
                                                0,
                                                m_options.m_prefix+"prev_chain");
    uqScalarSequenceClass<double> prevLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> prevLogTargetValues    (m_env,0,"");

    unsigned int indexOfFirstWeight = 0;
    unsigned int indexOfLastWeight  = 0;

    //std::cout << "m_env.inter0Rank() = " << m_env.inter0Rank() << std::endl;
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 2 of 9"
                                << std::endl;
      }

      prevChain = currChain;
      currChain.clear();
      currChain.setName(currOptions->m_prefix + "rawChain");

      prevLogLikelihoodValues = currLogLikelihoodValues; // likelihood is important
      prevLogTargetValues     = currLogTargetValues;

      currLogLikelihoodValues.clear();
      currLogLikelihoodValues.setName(currOptions->m_prefix + "rawLogLikelihood");

      currLogTargetValues.clear();
      currLogTargetValues.setName(currOptions->m_prefix + "rawLogTarget");

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        P_V prevPosition(m_vectorSpace.zeroVector());
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ":"
                                << std::endl;
        for (unsigned int i = 0; i < prevChain.subSequenceSize(); ++i) {
          prevChain.getPositionValues(i,prevPosition);
          *m_env.subDisplayFile() << "  prevChain[" << i
                                  << "] = " << prevPosition
                                  << ", prevLogLikelihoodValues[" << i
                                  << "] = " << prevLogLikelihoodValues[i]
                                  << ", prevLogTargetValues[" << i
                                  << "] = " << prevLogTargetValues[i]
                                  << std::endl;
        }
      }
#endif

      unsigned int quantity1 = prevChain.unifiedSequenceSize();
      unsigned int quantity2 = currChain.unifiedSequenceSize();
      unsigned int quantity3 = prevLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity4 = currLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity5 = prevLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity6 = currLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": prevChain.unifiedSequenceSize() = " << quantity1
                                << ", currChain.unifiedSequenceSize() = " << quantity2
                                << ", prevLogLikelihoodValues.unifiedSequenceSize() = " << quantity3
                                << ", currLogLikelihoodValues.unifiedSequenceSize() = " << quantity4
                                << ", prevLogTargetValues.unifiedSequenceSize() = " << quantity5
                                << ", currLogTargetValues.unifiedSequenceSize() = " << quantity6
                                << std::endl;
      }

      UQ_FATAL_TEST_MACRO((prevChain.subSequenceSize() != prevLogLikelihoodValues.subSequenceSize()),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "different sizes between previous chain and previous sequence of likelihood values");

      UQ_FATAL_TEST_MACRO((prevChain.subSequenceSize() != prevLogTargetValues.subSequenceSize()),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "different sizes between previous chain and previous sequence of target values");

      // Set 'indexOfFirstWeight' and 'indexOfLastWeight' // KAUST
      indexOfFirstWeight = 0;
      indexOfLastWeight  = indexOfFirstWeight + prevChain.subSequenceSize()-1;
      {
        //std::cout << "m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc() << std::endl;
        int r = m_env.inter0Rank();
        //std::cout << "r = " << r << std::endl;
        m_env.inter0Comm().Barrier();
        unsigned int auxUint = 0;
        if (r > 0) {
          MPI_Status status;
	  //std::cout << "Rank " << r << " is entering MPI_Recv()" << std::endl;
          int mpiRC = MPI_Recv((void*) &auxUint, 1, MPI_UNSIGNED, r-1, r-1, m_env.inter0Comm().Comm(), &status);
	  //std::cout << "Rank " << r << " received auxUint = " << auxUint << std::endl;
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Recv()");
          indexOfFirstWeight = auxUint;
          indexOfLastWeight = indexOfFirstWeight + prevChain.subSequenceSize()-1;
        }
        if (r < (m_env.inter0Comm().NumProc()-1)) {
          auxUint = indexOfLastWeight + 1;
	  //std::cout << "Rank " << r << " is sending auxUint = " << auxUint << std::endl;
          int mpiRC = MPI_Send((void*) &auxUint, 1, MPI_UNSIGNED, r+1, r, m_env.inter0Comm().Comm());
	  //std::cout << "Rank " << r << " sent auxUint = " << auxUint << std::endl;
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Send()");
        }
        m_env.inter0Comm().Barrier();
      }
    } // end of step 2

    //***********************************************************
    // Step 3 of 9: compute [currExponent and sequence of weights] for current level
    //***********************************************************
    m_currStep++;
    uqScalarSequenceClass<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 3 of 9"
                                << std::endl;
      }

      std::vector<double> exponents(2,0.);
      exponents[0] = prevExponent;
      exponents[1] = 1.;

      double nowExponent = 1.; // Try '1.' right away
      double nowEffectiveSizeRatio = 0.; // To be computed

      unsigned int nowAttempt = 0;
      bool testResult = false;
      double meanEffectiveSizeRatio = .5*(currOptions->m_minEffectiveSizeRatio + currOptions->m_maxEffectiveSizeRatio);
      uqScalarSequenceClass<double> omegaLnDiffSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");

      double nowUnifiedEvidenceLnFactor = 0.;
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ": entering loop for computing next exponent"
                                  << ", with nowAttempt = " << nowAttempt
                                  << std::endl;
        }

        if (nowAttempt > 0) {
          if (nowEffectiveSizeRatio > meanEffectiveSizeRatio) {
            exponents[0] = nowExponent;
          }
          else {
            exponents[1] = nowExponent;
          }
          nowExponent = .5*(exponents[0] + exponents[1]);
        }
        double auxExponent = nowExponent;
        if (prevExponent != 0.) {
          auxExponent /= prevExponent;
          auxExponent -= 1.;
        }
        double subWeightRatioSum     = 0.;
        double unifiedWeightRatioSum = 0.;

        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] = prevLogLikelihoodValues[i]*auxExponent; // likelihood is important
        }

        double unifiedOmegaLnMax = 0.;
        double unifiedOmegaLnMin = 0.;
        omegaLnDiffSequence.unifiedMinMax(m_vectorSpace.numOfProcsForStorage() == 1, // KAUST3
                                          0,
                                          omegaLnDiffSequence.subSequenceSize(),
                                          unifiedOmegaLnMin,
                                          unifiedOmegaLnMax);
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] -= unifiedOmegaLnMax;
          weightSequence[i] = exp(omegaLnDiffSequence[i]);
          subWeightRatioSum += weightSequence[i];
        }
        int mpiRC = MPI_Allreduce((void *) &subWeightRatioSum, (void *) &unifiedWeightRatioSum, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for weight ratio sum");

        nowUnifiedEvidenceLnFactor = log(unifiedWeightRatioSum) + unifiedOmegaLnMax - log(weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1));

        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] /= unifiedWeightRatioSum;
          effectiveSampleSize += weightSequence[i]*weightSequence[i];
          //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          //  *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
          //                          << ", level "                 << m_currLevel+LEVEL_REF_ID
          //                          << ", step "                  << m_currStep
          //                          << ": i = "                   << i
          //                          << ", effectiveSampleSize = " << effectiveSampleSize
          //                          << std::endl;
          //}
        }

#if 0 // For debug only
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ":"
                                  << std::endl;
        }
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "  weightSequence[" << i
                                    << "] = "              << weightSequence[i]
                                    << std::endl;
          }
        }
#endif

        double subQuantity = effectiveSampleSize;
        effectiveSampleSize = 0.;
        mpiRC = MPI_Allreduce((void *) &subQuantity, (void *) &effectiveSampleSize, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for effective sample size");

        effectiveSampleSize = 1./effectiveSampleSize;
        nowEffectiveSizeRatio = effectiveSampleSize/((double) weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1));
        UQ_FATAL_TEST_MACRO((nowEffectiveSizeRatio > (1.+1.e-8)),
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "effective sample size ratio cannot be > 1");
        //UQ_FATAL_TEST_MACRO((nowEffectiveSizeRatio < (1.-1.e-8)),
        //                    m_env.fullRank(),
        //                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
        //                    "effective sample size ratio cannot be < 1");

        //bool aux1 = (nowEffectiveSizeRatio == meanEffectiveSizeRatio);
        bool aux2 = (nowExponent == 1.                             )
                    &&
                    (nowEffectiveSizeRatio > meanEffectiveSizeRatio);
        bool aux3 = (nowEffectiveSizeRatio >= currOptions->m_minEffectiveSizeRatio)
                    &&
                    (nowEffectiveSizeRatio <= currOptions->m_maxEffectiveSizeRatio);
        testResult = aux2 || aux3;

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "                   << m_currLevel+LEVEL_REF_ID
                                  << ", step "                    << m_currStep
                                  << ": nowAttempt = "            << nowAttempt
                                  << ", prevExponent = "          << prevExponent
                                  << ", exponents[0] = "          << exponents[0]
                                  << ", nowExponent = "           << nowExponent
                                  << ", exponents[1] = "          << exponents[1]
                                  << ", effectiveSampleSize = "   << effectiveSampleSize
                                  << ", weightSequenceSize = "    << weightSequence.subSequenceSize()
                                  << ", minEffectiveSizeRatio = " << currOptions->m_minEffectiveSizeRatio
                                  << ", nowEffectiveSizeRatio = " << nowEffectiveSizeRatio
                                  << ", maxEffectiveSizeRatio = " << currOptions->m_maxEffectiveSizeRatio
                                  << std::endl;
        }
        nowAttempt++;

        // Make sure all nodes in 'inter0Comm' have the same value of 'nowExponent'
        uqMiscCheckForSameValueInAllNodes(nowExponent,
                                          0.,
                                          m_env.inter0Comm(),
                                          "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, testResult");

        // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
        uqMiscCheckForSameValueInAllNodes(testResult,
                                          0.,
                                          m_env.inter0Comm(),
                                          "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, testResult");
      } while (testResult == false);
      currExponent = nowExponent;
      m_logEvidenceFactors.push_back(nowUnifiedEvidenceLnFactor);

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                << ", step "                                   << m_currStep
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", currExponent = "                         << currExponent
                                << ", effective ratio = "                      << nowEffectiveSizeRatio
                                << ", log(evidence factor) = "                 << m_logEvidenceFactors[m_logEvidenceFactors.size()-1]
                                << ", evidence factor = "                      << exp(m_logEvidenceFactors[m_logEvidenceFactors.size()-1])
                                << std::endl;

        //unsigned int numZeros = 0;
        //for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        //  *m_env.subDisplayFile() << "weightSequence[" << i
        //                          << "] = " << weightSequence[i]
        //                         << std::endl;
        //  if (weightSequence[i] == 0.) numZeros++;
        //}
        //*m_env.subDisplayFile() << "Number of zeros in weightSequence = " << numZeros
        //                        << std::endl;
      }

      // Make sure all nodes in 'inter0Comm' have the same value of 'logEvidenceFactor'
      uqMiscCheckForSameValueInAllNodes(m_logEvidenceFactors[m_logEvidenceFactors.size()-1],
                                        1.e-16,
                                        m_env.inter0Comm(),
                                        "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, logEvidenceFactor");
    } // end of step 3 // if (m_env.inter0Rank() >= 0)

    // KAUST: all nodes in 'subComm' should have the same 'currExponent'
    int mpiRC = MPI_Bcast((void *) &currExponent, (int) 1, MPI_DOUBLE, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "failed MPI_Bcast() for currExponent");

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

      if (m_env.inter0Rank() >= 0) {
        // It is necessary to recompute 'unifiedRequestedNumSamples' because
        // 'currOptions' has just been replaced by 'lastLevelOptions'
        unsigned int tmpSize = currOptions->m_rawChainSize;
        int mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &unifiedRequestedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for requested num samples in step 3");
      }
    }

    //***********************************************************
    // Step 4 of 9: create covariance matrix for current level
    //***********************************************************
    m_currStep++;
    P_V oneVec(m_vectorSpace.zeroVector());
    oneVec.cwSet(1.);

    P_M* unifiedCovMatrix = NULL;
    if (m_env.inter0Rank() >= 0) {
      unifiedCovMatrix = m_vectorSpace.newMatrix();
    }
    else {
      unifiedCovMatrix = new P_M(oneVec);
    }
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 4 of 9"
                                << std::endl;
      }

      P_V auxVec(m_vectorSpace.zeroVector());
      P_V weightedMeanVec(m_vectorSpace.zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        weightedMeanVec += weightSequence[i]*auxVec;
      }

      P_V diffVec(m_vectorSpace.zeroVector());
      P_M subCovMatrix(m_vectorSpace.zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        diffVec = auxVec - weightedMeanVec;
        subCovMatrix += weightSequence[i]*matrixProduct(diffVec,diffVec);
      }

      for (unsigned int i = 0; i < unifiedCovMatrix->numRowsLocal(); ++i) { // KAUST5
        for (unsigned int j = 0; j < unifiedCovMatrix->numCols(); ++j) {
          double localValue = subCovMatrix(i,j);
          double sumValue = 0.;
          if (m_env.inter0Rank() >= 0) {
            int mpiRC = MPI_Allreduce((void *) &localValue, (void *) &sumValue, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
            UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                                m_env.fullRank(),
                                "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                "failed MPI_Allreduce() for cov matrix");
          }
          else {
            sumValue = localValue;
          }
          (*unifiedCovMatrix)(i,j) = sumValue;
        }
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "              << m_currLevel+LEVEL_REF_ID
                                << ", step "               << m_currStep
                                << ": unifiedCovMatrix = " << *unifiedCovMatrix
                                << std::endl;
      }

    } // end of step 4

    //***********************************************************
    // Step 5 of 9: create *unified* finite distribution for current level
    //***********************************************************
    m_currStep++;
    std::vector<unsigned int> unifiedIndexCountersAtProc0Only(0);
    std::vector<double>       unifiedWeightStdVectorAtProc0Only(0); // KAUST, to check
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 5 of 9"
                                << std::endl;
      }

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", before weightSequence.getUnifiedContentsAtProc0Only()"
                                << ":"
                                << std::endl;
      }
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << ", weightSequence[" << i
                                  << "] = "              << weightSequence[i]
                                  << std::endl;
        }
      }
#endif

      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", after weightSequence.getUnifiedContentsAtProc0Only()"
                                << ":"
                                << std::endl;
      }
      for (unsigned int i = 0; i < unifiedWeightStdVectorAtProc0Only.size(); ++i) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "  unifiedWeightStdVectorAtProc0Only[" << i
                                  << "] = "                                 << unifiedWeightStdVectorAtProc0Only[i]
                                  << std::endl;
        }
      }
#endif
      sampleIndexesAtProc0(unifiedRequestedNumSamples,        // input
                           unifiedWeightStdVectorAtProc0Only, // input
                           unifiedIndexCountersAtProc0Only);  // output

      unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if (m_env.inter0Rank() == 0) {
        UQ_FATAL_TEST_MACRO(unifiedIndexCountersAtProc0Only.size() != auxUnifiedSize,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "wrong output from sampleIndexesAtProc0() in step 5");
      }
    } // end of step 5

    //***********************************************************
    // Step 6 of 9: plan for number of linked chains for each node so that all
    //              nodes generate the closest possible to the same number of positions
    //***********************************************************
    m_currStep++;
    uqBalancedLinkedChainsPerNodeStruct<P_V> balancedLinkControl;
    uqUnbalancedLinkedChainsPerNodeStruct    unbalancedLinkControl; // KAUST

    unsigned int Np = 0;
    if (m_env.inter0Rank() == 0) { // Yes, '== 0'
      Np = (unsigned int) m_env.inter0Comm().NumProc();
    }
    std::vector<uqExchangeInfoStruct> exchangeStdVec(0);

    // All processors should call this routine in order to have the same decision value
    bool useBalancedChains = decideOnBalancedChains(currOptions,                     // input
                                                    indexOfFirstWeight,              // input
                                                    indexOfLastWeight,               // input
                                                    unifiedIndexCountersAtProc0Only, // input
                                                    exchangeStdVec);                 // output

    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 6 of 9"
                                << std::endl;
      }

      if (useBalancedChains) {
        prepareBalLinkedChains(currOptions,                     // input
                               prevChain,                       // input
                               exchangeStdVec,                  // input/output
                               balancedLinkControl);            // output
        // aqui: prevChain might not be needed anymore. Delete it to save memory
      }
      else {
        prepareUnbLinkedChains(indexOfFirstWeight,              // input
                               indexOfLastWeight,               // input
                               unifiedIndexCountersAtProc0Only, // input
                               unbalancedLinkControl);          // output
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": balancedLinkControl.balLinkedChains.size() = "   << balancedLinkControl.balLinkedChains.size()
                                << ", unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                                << std::endl;
      }
    } // end of step 6

    //***********************************************************
    // Step 7 of 9: create vector RV for current level
    //***********************************************************
    m_currStep++;
    uqBayesianJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                             m_priorRv.pdf(),
                                             m_likelihoodFunction,
                                             currExponent,
                                             *m_targetDomain);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           *m_targetDomain);

    //if (m_env.inter0Rank() >= 0) // KAUST
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 7 of 9"
                                << std::endl;
      }

      currRv.setPdf(currPdf);
    } // end of step 7

    //***********************************************************
    // Step 8 of 9: scale unified covariance matrix until min <= rejection rate <= max
    //***********************************************************
    m_currStep++;
    if (currOptions->m_scaleCovMatrix == false) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": skipping step 8 of 9"
                                << std::endl;
      }
    }
    else {
      generateSequence_Step8(prevChain,                         // input
                             indexOfFirstWeight,                // input
                             indexOfLastWeight,                 // input
                             unifiedWeightStdVectorAtProc0Only, // input
                             weightSequence,                    // input
                             prevEta,                           // input
                             currRv,                            // input
                             currOptions,                       // input (changed temporarily internally)
                             *unifiedCovMatrix,                 // input/output
                             currEta);                          // output
    } // end of step 8

    //***********************************************************
    // Step 9 of 9: sample vector RV of current level
    //***********************************************************
    m_currStep++;
    unsigned int unifiedNumberOfRejections = 0;
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 9 of 9"
                                << std::endl;
      }

      // KAUST: all nodes should call here
      bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
      unsigned int savedRawChainSize          = currOptions->m_rawChainSize; // Ok to use rawChainSize
      bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
      bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;

      currOptions->m_totallyMute           = true;
      currOptions->m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
      currOptions->m_rawChainComputeStats  = false;
      currOptions->m_filteredChainGenerate = false;

      // KAUST: all nodes should call here
      if (useBalancedChains) {
        generateBalLinkedChains(*currOptions,                 // input, only m_rawChainSize changes
                                *unifiedCovMatrix,            // input
                                currRv,                       // input
                                balancedLinkControl,          // input // Round Rock
                                currChain,                    // output
                                cumulativeRawChainRunTime,    // output
                                cumulativeRawChainRejections, // output
                                &currLogLikelihoodValues,     // output // likelihood is important
                                &currLogTargetValues);        // output
      }
      else {
        generateUnbLinkedChains(*currOptions,                 // input, only m_rawChainSize changes
                                *unifiedCovMatrix,            // input
                                currRv,                       // input
                                unbalancedLinkControl,        // input // Round Rock
                                indexOfFirstWeight,           // input // Round Rock
                                prevChain,                    // input // Round Rock
                                currChain,                    // output
                                cumulativeRawChainRunTime,    // output
                                cumulativeRawChainRejections, // output
                                &currLogLikelihoodValues,     // output // likelihood is important
                                &currLogTargetValues);        // output
      }

      // KAUST: all nodes should call here
      currOptions->m_totallyMute           = savedTotallyMute;
      currOptions->m_rawChainSize          = savedRawChainSize;
      currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
      currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

      for (unsigned int i = 0; i < balancedLinkControl.balLinkedChains.size(); ++i) {
        UQ_FATAL_TEST_MACRO(balancedLinkControl.balLinkedChains[i].initialPosition == NULL,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "Initial position pointer in step 9 should not be NULL");
        delete balancedLinkControl.balLinkedChains[i].initialPosition;
        balancedLinkControl.balLinkedChains[i].initialPosition = NULL;
      }
      balancedLinkControl.balLinkedChains.clear();

      if (m_env.inter0Rank() >= 0) { // KAUST
        if (currOptions->m_rawChainComputeStats) {
          std::ofstream* genericOfsVar = NULL;
          m_env.openOutputFile(currOptions->m_dataOutputFileName,
                               UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                               currOptions->m_dataOutputAllowedSet,
                               false,
                               genericOfsVar);

          currChain.computeStatistics(*currOptions->m_rawChainStatisticalOptions,
                                      genericOfsVar);

          //genericOfsVar->close();
          delete genericOfsVar;
        }

        if (currOptions->m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
          currChain.unifiedWriteContents              (currOptions->m_rawChainDataOutputFileName); // KAUST5
          currLogLikelihoodValues.unifiedWriteContents(currOptions->m_rawChainDataOutputFileName);
          currLogTargetValues.unifiedWriteContents    (currOptions->m_rawChainDataOutputFileName);
        }

        if (currOptions->m_filteredChainGenerate) {
          std::ofstream* genericOfsVar = NULL;
          m_env.openOutputFile(currOptions->m_dataOutputFileName,
                               UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                               currOptions->m_dataOutputAllowedSet,
                               false,
                               genericOfsVar);

          unsigned int filterInitialPos = (unsigned int) (currOptions->m_filteredChainDiscardedPortion * (double) currChain.subSequenceSize());
          unsigned int filterSpacing    = currOptions->m_filteredChainLag;
          if (filterSpacing == 0) {
            currChain.computeFilterParams(*currOptions->m_filteredChainStatisticalOptions,
                                          genericOfsVar,
                                          filterInitialPos,
                                          filterSpacing);
          }

          // Filter positions from the converged portion of the chain
          currChain.filter(filterInitialPos,
                           filterSpacing);
          currChain.setName(currOptions->m_prefix + "filtChain");

          currLogLikelihoodValues.filter(filterInitialPos,
                                         filterSpacing);
          currLogLikelihoodValues.setName(currOptions->m_prefix + "filtLogLikelihood");

          currLogTargetValues.filter(filterInitialPos,
                                     filterSpacing);
          currLogTargetValues.setName(currOptions->m_prefix + "filtLogTarget");

          if (currOptions->m_filteredChainComputeStats) {
            currChain.computeStatistics(*currOptions->m_filteredChainStatisticalOptions,
                                      genericOfsVar);
          }

          //genericOfsVar->close();
          delete genericOfsVar;

          if (currOptions->m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
            currChain.unifiedWriteContents              (currOptions->m_filteredChainDataOutputFileName);
            currLogLikelihoodValues.unifiedWriteContents(currOptions->m_filteredChainDataOutputFileName);
            currLogTargetValues.unifiedWriteContents    (currOptions->m_filteredChainDataOutputFileName);
          }
        } // if (currOptions->m_filteredChainGenerate)

        if (currOptions->m_filteredChainGenerate) {
          // Do not check
        }
        else {
          // Check if unified size of generated chain matches the unified requested size // KAUST
          unsigned int tmpSize = currChain.subSequenceSize();
          unsigned int unifiedGeneratedNumSamples = 0;
          unsigned int mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &unifiedGeneratedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Allreduce() for generated num samples in step 9");
          //std::cout << "unifiedGeneratedNumSamples = "   << unifiedGeneratedNumSamples
          //          << ", unifiedRequestedNumSamples = " << unifiedRequestedNumSamples
          //          << std::endl;
          UQ_FATAL_TEST_MACRO(unifiedGeneratedNumSamples != unifiedRequestedNumSamples,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "currChain (linked one) has been generated with invalid size");
        }

        // Compute unified number of rejections
        mpiRC = MPI_Allreduce((void *) &cumulativeRawChainRejections, (void *) &unifiedNumberOfRejections, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for number of rejections");
      } // if (m_env.inter0Rank() >= 0) // KAUST

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }
    } // end of step 9

    //***********************************************************
    // Prepare to end current level
    //***********************************************************
    delete unifiedCovMatrix;

    double levelRunTime = uqMiscGetEllapsedSeconds(&timevalLevel);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level "                   << m_currLevel+LEVEL_REF_ID
                              << " after "                           << levelRunTime << " seconds"
                              << ", cumulativeRawChainRunTime = "    << cumulativeRawChainRunTime << " seconds"
                              << ", cumulativeRawChainRejections = " << cumulativeRawChainRejections
                              << " (" << 100.*((double) cumulativeRawChainRejections)/((double) currOptions->m_rawChainSize)
                              << "% at this processor)"
                              << " (" << 100.*((double) unifiedNumberOfRejections)/((double) unifiedRequestedNumSamples)
                              << "% over all processors)"
                              << std::endl;
    }
    if (currExponent != 1.) delete currOptions;
  } // end of level while

  UQ_FATAL_TEST_MACRO((currExponent < 1),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                      "exponent has not achieved value '1' even after exiting level loop");

  if (m_env.inter0Rank() >= 0) { // KAUST
    UQ_FATAL_TEST_MACRO((m_currLevel != m_logEvidenceFactors.size()),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "invalid m_currLevel at the exit of the level loop");
    m_logEvidence = 0.;
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      m_logEvidence += m_logEvidenceFactors[i];
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ", log(evidence) = " << m_logEvidence
                              << ", evidence = "      << exp(m_logEvidence)
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
#endif // __UQ_MULTI_LEVEL_SAMPLING_H__
