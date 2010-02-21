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
  void   generateSequence(uqBaseVectorSequenceClass<P_V,P_M>&             workingChain,
                          uqScalarSequenceClass<double>*                  workingLogLikelihoodValues,
                          uqScalarSequenceClass<double>*                  workingLogTargetValues);

  void   print           (std::ostream& os) const;

private:
  void   sampleIndexesAtProc0   (unsigned int                                    currLevel,                         // input
                                 unsigned int                                    unifiedRequestedNumSamples,        // input
                                 const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only, // input
                                 std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only);  // output

  void   prepareBalLinkedChains (unsigned int                                    currLevel,                       // input
                                 const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                       // input
                                 unsigned int                                    indexOfFirstWeight,              // input
                                 unsigned int                                    indexOfLastWeight,               // input
                                 const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only, // input
                                 uqBalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);            // output

  void   prepareUnbLinkedChains (unsigned int                                    currLevel,                       // input
                                 unsigned int                                    indexOfFirstWeight,              // input
                                 unsigned int                                    indexOfLastWeight,               // input
                                 const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only, // input
                                 uqUnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl);          // output

  void   generateBalLinkedChains(unsigned int                                    currLevel,               // input
                                 uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
                                 const P_M&                                      unifiedCovMatrix,        // input
                                 const uqGenericVectorRVClass  <P_V,P_M>&        rv,                      // input
                                 const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,     // input // Round Rock
                                 uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,            // output
                                 double&                                         cumulativeRunTime,       // output
                                 unsigned int&                                   cumulativeRejections,    // output
                                 uqScalarSequenceClass         <double>*         currLogLikelihoodValues, // output
                                 uqScalarSequenceClass         <double>*         currLogTargetValues);    // output

  void   generateUnbLinkedChains(unsigned int                                    currLevel,               // input
                                 uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
                                 const P_M&                                      unifiedCovMatrix,        // input
                                 const uqGenericVectorRVClass  <P_V,P_M>&        rv,                      // input
                                 const uqUnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,   // input // Round Rock
                                 unsigned int                                    indexOfFirstWeight,      // input // Round Rock
                                 const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,               // input // Round Rock
                                 uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,            // output
                                 double&                                         cumulativeRunTime,       // output
                                 unsigned int&                                   cumulativeRejections,    // output
                                 uqScalarSequenceClass         <double>*         currLogLikelihoodValues, // output
                                 uqScalarSequenceClass         <double>*         currLogTargetValues);    // output

  void   solveBIPAtProc0       (unsigned int                       currLevel,
                                const std::vector<unsigned int>&   origNumChainsPerNode,
                                const std::vector<unsigned int>&   origNumPositionsPerNode,
                                std::vector<unsigned int>&         finalNumChainsPerNode,
                                std::vector<unsigned int>&         finalNumPositionsPerNode,
                                std::vector<uqExchangeInfoStruct>& exchangeStdVec);

  void   mpiExchangePositions  (unsigned int                              currLevel,
                                const uqSequenceOfVectorsClass<P_V,P_M>&  prevChain,
                                const std::vector<uqExchangeInfoStruct>&  exchangeStdVec,
                                const std::vector<unsigned int>&          finalNumChainsPerNode,
                                const std::vector<unsigned int>&          finalNumPositionsPerNode,
                                uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl);

  const uqBaseEnvironmentClass&             m_env;
  const uqBaseVectorRVClass      <P_V,P_M>& m_priorRv;
  const uqBaseScalarFunctionClass<P_V,P_M>& m_likelihoodFunction;
  const uqVectorSpaceClass       <P_V,P_M>& m_vectorSpace;
        uqVectorSetClass         <P_V,P_M>* m_targetDomain;

        uqMLSamplingOptionsClass            m_options;

	std::vector<double>                 m_logEvidenceFactors;
        double                              m_logEvidence;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj);

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

  unsigned int currLevel = 0;
  {
    sprintf(tmpSufix,"%d_",currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass currOptions(m_env,(m_options.m_prefix + tmpSufix).c_str());
    currOptions.scanOptionsValues(&defaultLevelOptions);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level "              << currLevel+LEVEL_REF_ID
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
                                << ", level " << currLevel+LEVEL_REF_ID
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
                              << ": ending level " << currLevel+LEVEL_REF_ID
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
    currLevel++;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level " << currLevel+LEVEL_REF_ID
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
    sprintf(tmpSufix,"%d_",currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass* currOptions = new uqMLSamplingLevelOptionsClass(m_env,(m_options.m_prefix + tmpSufix).c_str());
    currOptions->scanOptionsValues(&defaultLevelOptions);
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
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
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ", currOptions->m_rawChainSize = " << currOptions->m_rawChainSize // Ok to use rawChainSize
                                << std::endl;
      }
    }

    //***********************************************************
    // Step 2 of 9: save [chain and corresponding target pdf values] from previous level
    //***********************************************************
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
                                << ", level " << currLevel+LEVEL_REF_ID
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
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ", step 2"
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
                                << ", level " << currLevel+LEVEL_REF_ID
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
    uqScalarSequenceClass<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
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
                                  << ", level " << currLevel+LEVEL_REF_ID
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
          //                          << ", level "                 << currLevel+LEVEL_REF_ID
          //                          << ": i = "                   << i
          //                          << ", effectiveSampleSize = " << effectiveSampleSize
          //                          << std::endl;
          //}
        }

#if 0 // For debug only
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "         << currLevel+LEVEL_REF_ID
                                  << ", step 3"
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
                                  << ", level "                   << currLevel+LEVEL_REF_ID
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
                                << ", level "                                  << currLevel+LEVEL_REF_ID
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
                                << ", level " << currLevel+LEVEL_REF_ID
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
                                << ", level " << currLevel+LEVEL_REF_ID
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
                                << ", level "              << currLevel+LEVEL_REF_ID
                                << ": unifiedCovMatrix = " << *unifiedCovMatrix
                                << std::endl;
      }

    } // end of step 4

    //***********************************************************
    // Step 5 of 9: create *unified* finite distribution for current level
    //***********************************************************
    std::vector<unsigned int> unifiedIndexCountersAtProc0Only(0);
    std::vector<double>       unifiedWeightStdVectorAtProc0Only(0); // KAUST, to check
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 5 of 9"
                                << std::endl;
      }

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "         << currLevel+LEVEL_REF_ID
                                << ", step 5, before weightSequence.getUnifiedContentsAtProc0Only()"
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
                                << ", level "         << currLevel+LEVEL_REF_ID
                                << ", step 5, after weightSequence.getUnifiedContentsAtProc0Only()"
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
      sampleIndexesAtProc0(currLevel,                         // input
                           unifiedRequestedNumSamples,        // input
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
    //               nodes generate the closest possible to the same number of positions
    //***********************************************************
    uqBalancedLinkedChainsPerNodeStruct<P_V> balancedLinkControl;
    uqUnbalancedLinkedChainsPerNodeStruct    unbalancedLinkControl; // KAUST
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 6 of 9"
                                << std::endl;
      }

      if ((currOptions->m_loadBalance      ) &&
	  (m_env.inter0Comm().NumProc() > 1)) {
        prepareBalLinkedChains(currLevel,                       // input
                               prevChain,                       // input
                               indexOfFirstWeight,              // input
                               indexOfLastWeight,               // input
                               unifiedIndexCountersAtProc0Only, // input
                               balancedLinkControl);            // output
        // aqui: precChain might not be needed anymore. Delete it to save memory
      }
      else {
        prepareUnbLinkedChains(currLevel,                       // input
                               indexOfFirstWeight,              // input
                               indexOfLastWeight,               // input
                               unifiedIndexCountersAtProc0Only, // input
                               unbalancedLinkControl);          // output
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": balancedLinkControl.balLinkedChains.size() = "   << balancedLinkControl.balLinkedChains.size()
                                << ", unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                                << std::endl;
      }
    } // end of step 6

    //***********************************************************
    // Step 7 of 9: create vector RV for current level
    //***********************************************************
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
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 7 of 9"
                                << std::endl;
      }

      currRv.setPdf(currPdf);
    } // end of step 7

    //***********************************************************
    // Step 8 of 9: scale unified covariance matrix until min <= rejection rate <= max
    //***********************************************************
    if (currOptions->m_scaleCovMatrix == false) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": skipping step 8 of 9"
                                << std::endl;
      }
    }
    else {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 8 of 9"
                                << std::endl;
      }

      double beforeEta           = prevEta;
      double beforeRejectionRate = 0.;               // To be updated
      bool   beforeRejectionRateIsBelowRange = true; // To be updated

      double nowEta           = prevEta;
      double nowRejectionRate = 0.;               // To be computed
      bool   nowRejectionRateIsBelowRange = true; // To be computed

      std::vector<double> etas(2,0.);
      etas[0] = beforeEta;
      etas[1] = 1.;

      std::vector<double> rejs(2,0.);
      rejs[0] = 0.; // To be computed
      rejs[1] = 0.; // To be computed

      unsigned int nowAttempt = 0;
      bool testResult = false;
      double meanRejectionRate = .5*(currOptions->m_minRejectionRate + currOptions->m_maxRejectionRate);
      bool useMiddlePointLogicForEta = false;
      P_M nowCovMatrix(*unifiedCovMatrix);
#if 0 // KAUST, to check
      std::vector<double> unifiedWeightStdVectorAtProc0Only(0);
      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);
#endif
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+LEVEL_REF_ID
                                  << ": entering loop for assessing rejection rate"
                                  << ", with nowAttempt = "  << nowAttempt
                                  << ", nowRejectionRate = " << nowRejectionRate
                                  << std::endl;
        }
        nowCovMatrix = *unifiedCovMatrix;

        if (nowRejectionRate < currOptions->m_minRejectionRate) {
          nowRejectionRateIsBelowRange = true;
        }
        else if (nowRejectionRate > currOptions->m_maxRejectionRate) {
          nowRejectionRateIsBelowRange = false;
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "nowRejectionRate should be out of the requested range at this point of the logic");
        }

        if (m_env.inter0Rank() >= 0) { // KAUST
          if (nowAttempt > 0) {
            if (useMiddlePointLogicForEta == false) {
              if (nowAttempt == 1) {
                // Ok, keep useMiddlePointLogicForEta = false
              }
              else if ((beforeRejectionRateIsBelowRange == true) &&
                       (nowRejectionRateIsBelowRange    == true)) {
                // Ok
              }
              else if ((beforeRejectionRateIsBelowRange == false) &&
                       (nowRejectionRateIsBelowRange    == false)) {
                // Ok
              }
              else if ((beforeRejectionRateIsBelowRange == true ) &&
                       (nowRejectionRateIsBelowRange    == false)) {
                useMiddlePointLogicForEta = true;

                // This is the first time the middle point logic will be used below
                etas[0] = std::min(beforeEta,nowEta);
                etas[1] = std::max(beforeEta,nowEta);

                if (etas[0] == beforeEta) {
                  rejs[0] = beforeRejectionRate;
                  rejs[1] = nowRejectionRate;
                }
                else {
                  rejs[0] = nowRejectionRate;
                  rejs[1] = beforeRejectionRate;
                }
              }
              else if ((beforeRejectionRateIsBelowRange == false) &&
                       (nowRejectionRateIsBelowRange    == true )) {
                useMiddlePointLogicForEta = true;

                // This is the first time the middle point logic will be used below
                etas[0] = std::min(beforeEta,nowEta);
                etas[1] = std::max(beforeEta,nowEta);
              }
              else {
                UQ_FATAL_TEST_MACRO(true,
                                    m_env.fullRank(),
                                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                    "before and now range flags are inconsistent");
              }
            } // if (useMiddlePointLogicForEta == false)

            beforeEta                       = nowEta;
            beforeRejectionRate             = nowRejectionRate;
            beforeRejectionRateIsBelowRange = nowRejectionRateIsBelowRange;
            if (useMiddlePointLogicForEta == false) {
              if (beforeRejectionRateIsBelowRange) nowEta *= 4.;
              else                                 nowEta /= 4.;
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                        << ", level " << currLevel+LEVEL_REF_ID
                                        << ": in loop for assessing rejection rate"
                                        << ", with nowAttempt = "  << nowAttempt
                                        << ", useMiddlePointLogicForEta = false"
                                        << ", nowEta just updated to value (to be tested) " << nowEta
                                        << std::endl;
              }
            }
            else {
              if (nowRejectionRate > meanRejectionRate) {
                if (rejs[0] > meanRejectionRate) {
                  etas[0] = nowEta;
                  etas[1] = etas[1];
                }
                else {
                  etas[0] = etas[0];
                  etas[1] = nowEta;
                }
              }
              else {
                if (rejs[0] < meanRejectionRate) {
                  etas[0] = nowEta;
                  etas[1] = etas[1];
                }
                else {
                  etas[0] = etas[0];
                  etas[1] = nowEta;
                }
              }
              nowEta = .5*(etas[0] + etas[1]);
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                        << ", level " << currLevel+LEVEL_REF_ID
                                        << ": in loop for assessing rejection rate"
                                        << ", with nowAttempt = " << nowAttempt
                                        << ", useMiddlePointLogicForEta = true"
                                        << ", nowEta just updated to value (to be tested) " << nowEta
                                        << ", etas[0] = " << etas[0]
                                        << ", etas[1] = " << etas[1]
                                        << std::endl;
              }
            }
          } // if (nowAttempt > 0)
        } // if (m_env.inter0Rank() >= 0) // KAUST

        nowCovMatrix *= nowEta;

        unsigned int originalSubNumSamples = 1 + (unsigned int) ( (1.-meanRejectionRate)/meanRejectionRate/currOptions->m_covRejectionRate/currOptions->m_covRejectionRate );

        if (m_env.inter0Rank() >= 0) { // KAUST
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level " << currLevel+LEVEL_REF_ID
                                    << ": in loop for assessing rejection rate"
                                    << ", about to sample "     << originalSubNumSamples << " indexes"
                                    << ", meanRejectionRate = " << meanRejectionRate
                                    << ", covRejectionRate = "  << currOptions->m_covRejectionRate
                                    << std::endl;
          }
        } // KAUST

        std::vector<unsigned int> nowUnifiedIndexCountersAtProc0Only(0); // It will be resized by 'sampleIndexesAtProc0()' below
        if (m_env.inter0Rank() >= 0) { // KAUST
          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          sampleIndexesAtProc0(currLevel,                           // input
                               tmpUnifiedNumSamples,                // input
                               unifiedWeightStdVectorAtProc0Only,   // input
                               nowUnifiedIndexCountersAtProc0Only); // output

          unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
          if (m_env.inter0Rank() == 0) {
            UQ_FATAL_TEST_MACRO(nowUnifiedIndexCountersAtProc0Only.size() != auxUnifiedSize,
                                m_env.fullRank(),
                                "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                "wrong output from sampleIndexesAtProc0() in step 8");
          }

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level " << currLevel+LEVEL_REF_ID
                                    << ": in loop for assessing rejection rate"
                                    << ", about to distribute sampled assessment indexes"
                                    << std::endl;
          }
        } // KAUST

        uqBalancedLinkedChainsPerNodeStruct<P_V> nowBalLinkControl;
        uqUnbalancedLinkedChainsPerNodeStruct    nowUnbLinkControl; // KAUST

        if (m_env.inter0Rank() >= 0) { // KAUST
          if ((currOptions->m_loadBalance      ) &&
              (m_env.inter0Comm().NumProc() > 1)) {
            prepareBalLinkedChains(currLevel,                          // input
                                   prevChain,                          // input
                                   indexOfFirstWeight,                 // input
                                   indexOfLastWeight,                  // input
                                   nowUnifiedIndexCountersAtProc0Only, // input
                                   nowBalLinkControl);                 // output
          }
          else {
            prepareUnbLinkedChains(currLevel,                          // input
                                   indexOfFirstWeight,                 // input
                                   indexOfLastWeight,                  // input
                                   nowUnifiedIndexCountersAtProc0Only, // input
                                   nowUnbLinkControl);                 // output
          }
        } // KAUST

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+LEVEL_REF_ID
                                  << ": in loop for assessing rejection rate"
                                  << ", about to generate assessment chain"
                                  << std::endl;
        }

        uqSequenceOfVectorsClass<P_V,P_M> nowChain(m_vectorSpace,
                                                   0,
                                                   m_options.m_prefix+"now_chain");
        double       nowRunTime    = 0.;
        unsigned int nowRejections = 0;

        // KAUST: all nodes should call here
        bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
        unsigned int savedRawChainSize          = currOptions->m_rawChainSize; // Ok to use rawChainSize
        bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
        bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;
        unsigned int savedDrMaxNumExtraStages   = currOptions->m_drMaxNumExtraStages;
        unsigned int savedAmAdaptInterval       = currOptions->m_amAdaptInterval;

        currOptions->m_totallyMute           = true;
        currOptions->m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
        currOptions->m_rawChainComputeStats  = false;
        currOptions->m_filteredChainGenerate = false;
        currOptions->m_drMaxNumExtraStages   = 0;
        currOptions->m_amAdaptInterval       = 0;

        // KAUST: all nodes in 'subComm' should call here, important
        if ((currOptions->m_loadBalance    ) &&
            (m_env.numSubEnvironments() > 1)) { // Cannot use 'm_env.inter0Comm().NumProc()' because not all nodes at this point of the code belong to 'inter0Comm'
          generateBalLinkedChains(currLevel,          // input
                                  *currOptions,       // input, only m_rawChainSize changes
                                  nowCovMatrix,       // input
                                  currRv,             // input
                                  nowBalLinkControl,  // input // Round Rock
                                  nowChain,           // output 
                                  nowRunTime,         // output
                                  nowRejections,      // output
                                  NULL,               // output
                                  NULL);              // output
        }
        else {
          generateUnbLinkedChains(currLevel,          // input
                                  *currOptions,       // input, only m_rawChainSize changes
                                  nowCovMatrix,       // input
                                  currRv,             // input
                                  nowUnbLinkControl,  // input // Round Rock
                                  indexOfFirstWeight, // input // Round Rock
                                  prevChain,          // input // Round Rock
                                  nowChain,           // output 
                                  nowRunTime,         // output
                                  nowRejections,      // output
                                  NULL,               // output
                                  NULL);              // output
        }

        // KAUST: all nodes should call here
        currOptions->m_totallyMute           = savedTotallyMute;
        currOptions->m_rawChainSize          = savedRawChainSize;
        currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
        currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME
        currOptions->m_drMaxNumExtraStages   = savedDrMaxNumExtraStages;
        currOptions->m_amAdaptInterval       = savedAmAdaptInterval;

        for (unsigned int i = 0; i < nowBalLinkControl.balLinkedChains.size(); ++i) {
          UQ_FATAL_TEST_MACRO(nowBalLinkControl.balLinkedChains[i].initialPosition == NULL,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "Initial position pointer in step 8 should not be NULL");
          delete nowBalLinkControl.balLinkedChains[i].initialPosition;
          nowBalLinkControl.balLinkedChains[i].initialPosition = NULL;
        }
        nowBalLinkControl.balLinkedChains.clear();

        if (m_env.inter0Rank() >= 0) { // KAUST
          // If only one cov matrix is used, then the rejection should be assessed among all inter0Comm nodes // KAUST3
          unsigned int nowUnifiedRejections = 0;
          int mpiRC = MPI_Allreduce((void *) &nowRejections, (void *) &nowUnifiedRejections, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Allreduce() for now rejections");

#if 0 // Round Rock 2009 12 29
          unsigned int tmpUnifiedNumSamples = 0;
          mpiRC = MPI_Allreduce((void *) &tmpSubNumSamples, (void *) &tmpUnifiedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Allreduce() for num samples in step 8");
#endif

          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          nowRejectionRate = ((double) nowUnifiedRejections) / ((double) tmpUnifiedNumSamples);

          //bool aux1 = (nowRejectionRate == meanRejectionRate);
          bool aux2 = (nowRejectionRate >= currOptions->m_minRejectionRate)
                      &&
                      (nowRejectionRate <= currOptions->m_maxRejectionRate);
          testResult = aux2;

          // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
          uqMiscCheckForSameValueInAllNodes(testResult,
                                            0.,
                                            m_env.inter0Comm(),
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 8, testResult");
        } // if (m_env.inter0Rank() >= 0) { // KAUST

        // KAUST: all nodes in 'subComm' should have the same 'testResult'
        unsigned int tmpUint = (unsigned int) testResult;
        int mpiRC = MPI_Bcast((void *) &tmpUint, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Bcast() for testResult");
        testResult = (bool) tmpUint;

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "              << currLevel+LEVEL_REF_ID
                                  << ": in loop for assessing rejection rate"
                                  << ", nowAttempt = "       << nowAttempt
                                  << ", beforeEta = "        << beforeEta
                                  << ", etas[0] = "          << etas[0]
                                  << ", nowEta = "           << nowEta
                                  << ", etas[1] = "          << etas[1]
                                  << ", minRejectionRate = " << currOptions->m_minRejectionRate
                                  << ", nowRejectionRate = " << nowRejectionRate
                                  << ", maxRejectionRate = " << currOptions->m_maxRejectionRate
                                  << std::endl;
        }
        nowAttempt++;

        if (m_env.inter0Rank() >= 0) { // KAUST
          // Make sure all nodes in 'inter0Comm' have the same value of 'nowEta'
          uqMiscCheckForSameValueInAllNodes(nowEta,
                                            1.e-16,
                                            m_env.inter0Comm(),
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 8, testResult");
        }
      } while (testResult == false);
      currEta = nowEta;
      if (currEta != 1.) {
        *unifiedCovMatrix *= currEta;
      }

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << currLevel+LEVEL_REF_ID
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", currEta = "                              << currEta
                                << ", assessed rejection rate = "              << nowRejectionRate
                                << std::endl;
      }
    } // end of step 8

    //***********************************************************
    // Step 9 of 9: sample vector RV of current level
    //***********************************************************
    unsigned int unifiedNumberOfRejections = 0;
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
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
      if ((currOptions->m_loadBalance    ) &&
          (m_env.numSubEnvironments() > 1)) { // Cannot use 'm_env.inter0Comm().NumProc()' because not all nodes at this point of the code belong to 'inter0Comm'
        generateBalLinkedChains(currLevel,                    // input
                                *currOptions,                 // input, only m_rawChainSize changes
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
        generateUnbLinkedChains(currLevel,                    // input
                                *currOptions,                 // input, only m_rawChainSize changes
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
                                << ", level " << currLevel+LEVEL_REF_ID
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
                              << ": ending level "                   << currLevel+LEVEL_REF_ID
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
    UQ_FATAL_TEST_MACRO((currLevel != m_logEvidenceFactors.size()),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "invalid currLevel at the exit of the level loop");
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

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::sampleIndexesAtProc0(
  unsigned int               currLevel,                         // input
  unsigned int               unifiedRequestedNumSamples,        // input
  const std::vector<double>& unifiedWeightStdVectorAtProc0Only, // input
  std::vector<unsigned int>& unifiedIndexCountersAtProc0Only)   // output
{
  if (m_env.inter0Rank() != 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::sampleIndexesAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": unifiedRequestedNumSamples = "               << unifiedRequestedNumSamples
                            << ", unifiedWeightStdVectorAtProc0Only.size() = " << unifiedWeightStdVectorAtProc0Only.size()
                            << std::endl;
  }

#if 0 // For debug only
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::sampleIndexesAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ":"
                            << std::endl;
    unsigned int numZeros = 0;
    for (unsigned int i = 0; i < unifiedWeightStdVectorAtProc0Only.size(); ++i) {
      *m_env.subDisplayFile() << "  unifiedWeightStdVectorAtProc0Only[" << i
                              << "] = " << unifiedWeightStdVectorAtProc0Only[i]
                              << std::endl;
      if (unifiedWeightStdVectorAtProc0Only[i] == 0.) numZeros++;
    }
    *m_env.subDisplayFile() << "Number of zeros in unifiedWeightStdVectorAtProc0Only = " << numZeros
                            << std::endl;
  }
#endif
 
  if (m_env.inter0Rank() == 0) {
    unsigned int resizeSize = unifiedWeightStdVectorAtProc0Only.size();
    unifiedIndexCountersAtProc0Only.resize(resizeSize,0);

    // Generate 'unifiedRequestedNumSamples' samples from 'tmpFD'
    uqFiniteDistributionClass tmpFd(m_env,
                                    "",
                                    unifiedWeightStdVectorAtProc0Only);
    for (unsigned int i = 0; i < unifiedRequestedNumSamples; ++i) {
      unsigned int index = tmpFd.sample();
      unifiedIndexCountersAtProc0Only[index] += 1;
    }
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains( // EXTRA FOR LOAD BALANCE
  unsigned int                              currLevel,                       // input
  const uqSequenceOfVectorsClass<P_V,P_M>&  prevChain,                       // input
  unsigned int                              indexOfFirstWeight,              // input
  unsigned int                              indexOfLastWeight,               // input
  const std::vector<unsigned int>&          unifiedIndexCountersAtProc0Only, // input
  uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl)             // output
{
  if (m_env.inter0Rank() < 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": indexOfFirstWeight = " << indexOfFirstWeight
                            << ", indexOfLastWeight = "  << indexOfLastWeight
                            << std::endl;
  }

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();

  //////////////////////////////////////////////////////////////////////////
  // Gather information at proc 0, so that it can solve BIP (binary integer problem)
  //////////////////////////////////////////////////////////////////////////
  //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
  //                void *recvbuf, int recvcount, MPI_Datatype recvtype, 
  //                int root, MPI_Comm comm )
  unsigned int auxUInt = indexOfFirstWeight;
  std::vector<unsigned int> allFirstIndexes(Np,0); // '0' is already the correct value for recvcnts[0]
  int mpiRC = MPI_Gather((void *) &auxUInt, 1, MPI_INT, (void *) &allFirstIndexes[0], (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Gather() for first indexes");
  if (m_env.inter0Rank() == 0) {
    UQ_FATAL_TEST_MACRO(allFirstIndexes[0] != indexOfFirstWeight,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                        "failed MPI_Gather() result for first indexes, at proc 0");
  }

  auxUInt = indexOfLastWeight;
  std::vector<unsigned int> allLastIndexes(Np,0); // '0' is NOT the correct value for recvcnts[0]
  mpiRC = MPI_Gather((void *) &auxUInt, 1, MPI_INT, (void *) &allLastIndexes[0], (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Gather() for last indexes");
  if (m_env.inter0Rank() == 0) {
    //allLastIndexes[0] = indexOfLastWeight; // FIX ME: really necessary????
    UQ_FATAL_TEST_MACRO(allLastIndexes[0] != indexOfLastWeight,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                        "failed MPI_Gather() result for last indexes, at proc 0");
  }

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 assembles and solves BIP
  //////////////////////////////////////////////////////////////////////////
  std::vector<uqExchangeInfoStruct> exchangeStdVec(0);
  std::vector<unsigned int> origNumChainsPerNode    (Np,0);
  std::vector<unsigned int> origNumPositionsPerNode (Np,0);
  std::vector<unsigned int> finalNumChainsPerNode   (Np,0);
  std::vector<unsigned int> finalNumPositionsPerNode(Np,0);
  if (m_env.inter0Rank() == 0) {
    *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::prepareBalLinkdeChains()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": original distribution of unified indexes in 'inter0Comm' is as follows"
                            << std::endl;
    for (unsigned int r = 0; r < Np; ++r) {
      *m_env.subDisplayFile() << "  allFirstIndexes[" << r << "] = " << allFirstIndexes[r]
                              << "  allLastIndexes["  << r << "] = " << allLastIndexes[r]
                              << std::endl;
    }
    for (unsigned int r = 0; r < (Np-1); ++r) { // Yes, '-1'
      UQ_FATAL_TEST_MACRO(allFirstIndexes[r+1] != (allLastIndexes[r]+1),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                          "wrong indexes");
    }

    for (unsigned int r = 0; r < (Np-1); ++r) { // Yes, '-1'
      UQ_FATAL_TEST_MACRO(allFirstIndexes[r+1] != (allLastIndexes[r]+1),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                          "wrong indexes");
    }

    int r = 0;
    for (unsigned int i = 0; i < unifiedIndexCountersAtProc0Only.size(); ++i) {
      if ((allFirstIndexes[r] <= i) && // FIX ME: not a robust logic
          (i <= allLastIndexes[r] )) {
        // Ok
      }
      else {
        r++;
        if ((r < (int) Np           ) &&
            (allFirstIndexes[r] <= i) && 
            (i <= allLastIndexes[r] )) {
          // Ok
        }
        else {
	  std::cerr << "In uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()"
                    << ": i = " << i
                    << ", r = " << r
                    << ", allFirstIndexes[r] = " << allFirstIndexes[r]
                    << ", allLastIndexes[r] = "  << allLastIndexes[r]
                    << std::endl;
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                              "wrong indexes or 'r' got too large");
        }
      }
      if (unifiedIndexCountersAtProc0Only[i] != 0) {
        origNumChainsPerNode[r]++;
        origNumPositionsPerNode[r] += unifiedIndexCountersAtProc0Only[i];

        uqExchangeInfoStruct auxInfo;
        auxInfo.originalNodeOfInitialPosition  = r;
        auxInfo.originalIndexOfInitialPosition = i - allFirstIndexes[r];
        auxInfo.finalNodeOfInitialPosition     = -1; // Yes, '-1' for now, important
        auxInfo.numberOfPositions              = unifiedIndexCountersAtProc0Only[i];
        exchangeStdVec.push_back(auxInfo);
      }
      // FIX ME: swap trick to save memory
    }

    // Add extra logic: check if number of procs is too large + check if ratio max/min justifies optimization

    // Get final node responsible for a linked chain by solving BIP at node zero only
    solveBIPAtProc0(currLevel,
                    origNumChainsPerNode,
                    origNumPositionsPerNode,
                    finalNumChainsPerNode,
                    finalNumPositionsPerNode,
                    exchangeStdVec);
  } // if (m_env.inter0Rank() == 0)

  m_env.inter0Comm().Barrier();

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 now broadcasts the information on 'exchangeStdVec'
  //////////////////////////////////////////////////////////////////////////
  unsigned int exchangeStdVecSize = exchangeStdVec.size();
  mpiRC = MPI_Bcast((void *) &exchangeStdVecSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Bcast() for exchangeStdVec size");
  if (m_env.inter0Rank() > 0) exchangeStdVec.resize(exchangeStdVecSize);

  mpiRC = MPI_Bcast((void *) &exchangeStdVec[0], (int) (exchangeStdVecSize*sizeof(uqExchangeInfoStruct)), MPI_CHAR, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Bcast() for exchangeStdVec data");

  //////////////////////////////////////////////////////////////////////////
  // Sanity check
  //////////////////////////////////////////////////////////////////////////
  unsigned int Nc = exchangeStdVec.size();
  if (m_env.inter0Rank() > 0) {
    for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
      unsigned int nodeId = exchangeStdVec[chainId].finalNodeOfInitialPosition;
      finalNumPositionsPerNode[nodeId] += exchangeStdVec[chainId].numberOfPositions;
    }
  }
  unsigned int finalMinPosPerNode = *std::min_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  unsigned int finalMaxPosPerNode = *std::max_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  double finalRatioOfPosPerNode = ((double) finalMaxPosPerNode) / ((double)finalMinPosPerNode);
  //std::cout << m_env.fullRank() << ", finalRatioOfPosPerNode = " << finalRatioOfPosPerNode << std::endl;

  std::vector<double> auxBuf(1,0.);
  double minRatio = 0.;
  auxBuf[0] = finalRatioOfPosPerNode;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minRatio, (int) auxBuf.size(), MPI_DOUBLE, MPI_MIN, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Allreduce() for min");
  //std::cout << m_env.fullRank() << ", minRatio = " << minRatio << std::endl;
  UQ_FATAL_TEST_MACRO(minRatio != finalRatioOfPosPerNode,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed minRatio sanity check");

  double maxRatio = 0.;
  auxBuf[0] = finalRatioOfPosPerNode;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxRatio, (int) auxBuf.size(), MPI_DOUBLE, MPI_MAX, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Allreduce() for max");
  //std::cout << m_env.fullRank() << ", maxRatio = " << maxRatio << std::endl;
  UQ_FATAL_TEST_MACRO(maxRatio != finalRatioOfPosPerNode,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed maxRatio sanity check");

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 now broadcasts the information on 'finalNumChainsPerNode'
  //////////////////////////////////////////////////////////////////////////
  unsigned int finalNumChainsPerNodeSize = finalNumChainsPerNode.size();
  mpiRC = MPI_Bcast((void *) &finalNumChainsPerNodeSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Bcast() for finalNumChainsPerNode size");
  if (m_env.inter0Rank() > 0) finalNumChainsPerNode.resize(finalNumChainsPerNodeSize);

  mpiRC = MPI_Bcast((void *) &finalNumChainsPerNode[0], (int) finalNumChainsPerNodeSize, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains()",
                      "failed MPI_Bcast() for finalNumChainsPerNode data");

  //////////////////////////////////////////////////////////////////////////
  // Mpi exchange information between nodes and properly populate
  // balancedLinkControl.linkedChains at each node
  //////////////////////////////////////////////////////////////////////////
  mpiExchangePositions(currLevel,
                       prevChain,
                       exchangeStdVec,
                       finalNumChainsPerNode,
                       finalNumPositionsPerNode, // It is already valid at all nodes (not only at node 0) because of the sanity check above
                       balancedLinkControl);

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains(
  unsigned int                           currLevel,                       // input
  unsigned int                           indexOfFirstWeight,              // input
  unsigned int                           indexOfLastWeight,               // input
  const std::vector<unsigned int>&       unifiedIndexCountersAtProc0Only, // input
  uqUnbalancedLinkedChainsPerNodeStruct& unbalancedLinkControl)           // output
{
  if (m_env.inter0Rank() < 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": indexOfFirstWeight = " << indexOfFirstWeight
                            << ", indexOfLastWeight = "  << indexOfLastWeight
                            << std::endl;
  }

  unsigned int              subNumSamples = 0;
  std::vector<unsigned int> unifiedIndexCountersAtAllProcs(0);

  // All nodes in 'inter0Comm' should resize to the same size // KAUST3
  unsigned int resizeSize = unifiedIndexCountersAtProc0Only.size();
  int mpiRC = MPI_Bcast((void *) &resizeSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "failed MPI_Bcast() for resizeSize");
  unifiedIndexCountersAtAllProcs.resize(resizeSize,0);

  if (m_env.inter0Rank() == 0) unifiedIndexCountersAtAllProcs = unifiedIndexCountersAtProc0Only;

  // Broadcast index counters to all nodes
  mpiRC = MPI_Bcast((void *) &unifiedIndexCountersAtAllProcs[0], (int) unifiedIndexCountersAtAllProcs.size(), MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "failed MPI_Bcast() for unified index counters");
#if 0 // Use allgatherv ??? for subNumSamples instead
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::prepareUnbLinkedChains()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ":"
                            << std::endl;
    for (int r = 0; r < m_env.inter0Comm().NumProc(); ++r) {
      *m_env.subDisplayFile() << "  unifiedIndexCountersAtAllProcs[" << r << "] = " << unifiedIndexCountersAtAllProcs[r]
                              << std::endl;
    }
  }
#endif
  //for (unsigned int i = 0; i < unifiedIndexCountersAtAllProcs.size(); ++i) {
  //  *m_env.subDisplayFile() << "unifiedIndexCountersAtAllProcs[" << i
  //                          << "] = " << unifiedIndexCountersAtAllProcs[i]
  //                          << std::endl;
  //}

  // Use 'indexOfFirstWeight' and 'indexOfLastWeight' in order to update 'subNumSamples'
  UQ_FATAL_TEST_MACRO(indexOfFirstWeight >= unifiedIndexCountersAtAllProcs.size(),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "invalid indexOfFirstWeight");
  UQ_FATAL_TEST_MACRO(indexOfLastWeight >= unifiedIndexCountersAtAllProcs.size(),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "invalid indexOfLastWeight");
  subNumSamples = 0;
  for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
    subNumSamples += unifiedIndexCountersAtAllProcs[i];
  }

  std::vector<unsigned int> auxBuf(1,0);

  unsigned int minModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "failed MPI_Allreduce() for min");

  unsigned int maxModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "failed MPI_Allreduce() for max");

  unsigned int sumModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "failed MPI_Allreduce() for sum");

  //UQ_FATAL_TEST_MACRO(unifiedRequestedNumSamples != sumModifiedSubNumSamples,
  //                    m_env.fullRank(),
  //                    "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
  //                    "invalid state");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()"
                            << ", level "                                   << currLevel+LEVEL_REF_ID
                            << ": subNumSamples = "                         << subNumSamples
                            << ", unifiedIndexCountersAtAllProcs.size() = " << unifiedIndexCountersAtAllProcs.size()
                            << std::endl;
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()"
                            << ", level "                      << currLevel+LEVEL_REF_ID
                            << ": minModifiedSubNumSamples = " << minModifiedSubNumSamples
                            << ", avgModifiedSubNumSamples = " << ((double) sumModifiedSubNumSamples)/((double) m_env.inter0Comm().NumProc())
                            << ", maxModifiedSubNumSamples = " << maxModifiedSubNumSamples
                            << std::endl;
  }

  unsigned int numberOfPositionsToGuaranteeForNode = subNumSamples;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::prepareUnbLinkedChains()"
                            << ", level "                                 << currLevel+LEVEL_REF_ID
                            << ": numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                            << std::endl;
  }
  for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
//for (unsigned int i = 0; i < unifiedIndexCountersAtAllProcs.size(); ++i) { // KAUST4: important
    while (unifiedIndexCountersAtAllProcs[i] != 0) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 30)) {
        *m_env.subDisplayFile() << ", numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                                << ", unifiedIndexCountersAtAllProcs["        << i
                                << "] = "                                     << unifiedIndexCountersAtAllProcs[i]
                                << std::endl;
      }
      if (unifiedIndexCountersAtAllProcs[i] < numberOfPositionsToGuaranteeForNode) {
        uqUnbalancedLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = unifiedIndexCountersAtAllProcs[i];
        unbalancedLinkControl.unbLinkedChains.push_back(auxControl);

        numberOfPositionsToGuaranteeForNode -= unifiedIndexCountersAtAllProcs[i];
        unifiedIndexCountersAtAllProcs[i] = 0;
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] > 0                                   )) {
      //else { // KAUST4
        uqUnbalancedLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = numberOfPositionsToGuaranteeForNode;
        unbalancedLinkControl.unbLinkedChains.push_back(auxControl);

        unifiedIndexCountersAtAllProcs[i] -= numberOfPositionsToGuaranteeForNode;
        numberOfPositionsToGuaranteeForNode = 0;
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] == 0                                  )) {
        // Ok
      }
      else {
        UQ_FATAL_TEST_MACRO(true, // KAUST4
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                            "should never get here");
      }
    }
  }
  UQ_FATAL_TEST_MACRO(numberOfPositionsToGuaranteeForNode != 0, // subNumSamples, // KAUST4
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()",
                      "numberOfPositionsToGuaranteeForNode exited loop with wrong value");
  // FIX ME: swap trick to save memory

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains()"
                            << ", level "                          << currLevel+LEVEL_REF_ID
                            << ": unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains( // EXTRA FOR LOAD BALANCE
  unsigned int                                    currLevel,               // input
  uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
  const P_M&                                      unifiedCovMatrix,        // input
  const uqGenericVectorRVClass  <P_V,P_M>&        rv,                      // input
  const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,     // input // Round Rock
  uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,            // output
  double&                                         cumulativeRunTime,       // output
  unsigned int&                                   cumulativeRejections,    // output
  uqScalarSequenceClass         <double>*         currLogLikelihoodValues, // output
  uqScalarSequenceClass         <double>*         currLogTargetValues)     // output
{
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()"
                            << ": balancedLinkControl.balLinkedChains.size() = " << balancedLinkControl.balLinkedChains.size()
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = balancedLinkControl.balLinkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  int mpiRC = MPI_Bcast((void *) &chainIdMax, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()",
                      "failed MPI_Bcast() for chainIdMax");

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfPositions = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfPositions += balancedLinkControl.balLinkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm()); // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()",
                        "failed MPI_Allreduce() for min");

    unsigned int maxNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm()); // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()",
                        "failed MPI_Allreduce() for max");

    unsigned int sumNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm()); // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()",
                        "failed MPI_Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()"
                              << ", level "               << currLevel+LEVEL_REF_ID
                              << ": chainIdMax = "        << chainIdMax
                              << ", numberOfPositions = " << numberOfPositions
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()"
                              << ", level "                  << currLevel+LEVEL_REF_ID
                              << ": minNumberOfPositions = " << minNumberOfPositions
                              << ", avgNumberOfPositions = " << ((double) sumNumberOfPositions)/((double) m_env.inter0Comm().NumProc())
                              << ", maxNumberOfPositions = " << maxNumberOfPositions
                              << std::endl;
    }
  }
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      // aqui 4
      auxInitialPosition = *(balancedLinkControl.balLinkedChains[chainId].initialPosition); // Round Rock
      tmpChainSize = balancedLinkControl.balLinkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()"
                                << ", level "          << currLevel+LEVEL_REF_ID
                                << ": chainId = "      << chainId
                                << ", tmpChainSize = " << tmpChainSize
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm().Comm()); // Yes, 'subComm', important // KAUST
#if 0 // For debug only
    for (int r = 0; r < m_env.subComm().NumProc(); ++r) {
      if (r == m_env.subComm().MyPID()) {
	std::cout << "Vector 'auxInitialPosition at rank " << r
                  << " has contents "                      << auxInitialPosition
                  << std::endl;
      }
      m_env.subComm().Barrier();
    }
    sleep(1);
#endif

    // KAUST: all nodes in 'subComm' should have the same 'tmpChainSize'
    mpiRC = MPI_Bcast((void *) &tmpChainSize, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()",
                        "failed MPI_Bcast() for tmpChainSize");

    inputOptions.m_rawChainSize = tmpChainSize;
    uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    uqScalarSequenceClass<double> tmpLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> tmpLogTargetValues    (m_env,0,"");

    // KAUST: all nodes should call here
    uqMetropolisHastingsSGClass<P_V,P_M> mcSeqGenerator(inputOptions,
                                                        rv,
                                                        auxInitialPosition,
                                                        &unifiedCovMatrix);

    // KAUST: all nodes should call here
    mcSeqGenerator.generateSequence(tmpChain,
                                    &tmpLogLikelihoodValues, // likelihood is IMPORTANT
                                    &tmpLogTargetValues);
    cumulativeRunTime    += mcSeqGenerator.rawChainRunTime();
    cumulativeRejections += mcSeqGenerator.numRejections();

    if (m_env.inter0Rank() >= 0) {
      //for (unsigned int i = 0; i < tmpLogLikelihoodValues.subSequenceSize(); ++i) {
      //  std::cout << "tmpLogLikelihoodValues[" << i << "] = " << tmpLogLikelihoodValues[i]
      //            << ", tmpLogTargetValues["   << i << "] = " << tmpLogTargetValues[i]
      //            << std::endl;
      //}
        
      if ((m_env.subDisplayFile()             ) &&
          (m_env.displayVerbosity()   >= 0    ) &&
          (inputOptions.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateBalLinkedChains()"
                                << ", level "               << currLevel+LEVEL_REF_ID
                                << ", chainId = "           << chainId
                                << ": finished generating " << tmpChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      // KAUST5: what if workingChain ends up with different size in different nodes? Important
      workingChain.append              (tmpChain,              1,tmpChain.subSequenceSize()-1              ); // IMPORTANT: '1' in order to discard initial position
      if (currLogLikelihoodValues) {
        currLogLikelihoodValues->append(tmpLogLikelihoodValues,1,tmpLogLikelihoodValues.subSequenceSize()-1); // IMPORTANT: '1' in order to discard initial position
      }
      if (currLogTargetValues) {
        currLogTargetValues->append    (tmpLogTargetValues,    1,tmpLogTargetValues.subSequenceSize()-1    ); // IMPORTANT: '1' in order to discard initial position
      }
    }
  } // for

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains()"
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains(
  unsigned int                                 currLevel,               // input
  uqMLSamplingLevelOptionsClass&               inputOptions,            // input, only m_rawChainSize changes
  const P_M&                                   unifiedCovMatrix,        // input
  const uqGenericVectorRVClass  <P_V,P_M>&     rv,                      // input
  const uqUnbalancedLinkedChainsPerNodeStruct& unbalancedLinkControl,   // input // Round Rock
  unsigned int                                 indexOfFirstWeight,      // input // Round Rock
  const uqSequenceOfVectorsClass<P_V,P_M>&     prevChain,               // input // Round Rock
  uqSequenceOfVectorsClass      <P_V,P_M>&     workingChain,            // output
  double&                                      cumulativeRunTime,       // output
  unsigned int&                                cumulativeRejections,    // output
  uqScalarSequenceClass         <double>*      currLogLikelihoodValues, // output
  uqScalarSequenceClass         <double>*      currLogTargetValues)     // output
{
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()"
                            << ": unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                            << ", indexOfFirstWeight = "                           << indexOfFirstWeight
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = unbalancedLinkControl.unbLinkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  int mpiRC = MPI_Bcast((void *) &chainIdMax, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()",
                      "failed MPI_Bcast() for chainIdMax");

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfPositions = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfPositions += unbalancedLinkControl.unbLinkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()",
                        "failed MPI_Allreduce() for min");

    unsigned int maxNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()",
                        "failed MPI_Allreduce() for max");

    unsigned int sumNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()",
                        "failed MPI_Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()"
                              << ", level "               << currLevel+LEVEL_REF_ID
                              << ": chainIdMax = "        << chainIdMax
                              << ", numberOfPositions = " << numberOfPositions
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()"
                              << ", level "                  << currLevel+LEVEL_REF_ID
                              << ": minNumberOfPositions = " << minNumberOfPositions
                              << ", avgNumberOfPositions = " << ((double) sumNumberOfPositions)/((double) m_env.inter0Comm().NumProc())
                              << ", maxNumberOfPositions = " << maxNumberOfPositions
                              << std::endl;
    }
  }
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      unsigned int auxIndex = unbalancedLinkControl.unbLinkedChains[chainId].initialPositionIndexInPreviousChain - indexOfFirstWeight; // KAUST4 // Round Rock
      prevChain.getPositionValues(auxIndex,auxInitialPosition); // Round Rock
      tmpChainSize = unbalancedLinkControl.unbLinkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()"
                                << ", level "          << currLevel+LEVEL_REF_ID
                                << ": chainId = "      << chainId
                                << ", tmpChainSize = " << tmpChainSize
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm().Comm()); // Yes, 'subComm', important // KAUST
#if 0 // For debug only
    for (int r = 0; r < m_env.subComm().NumProc(); ++r) {
      if (r == m_env.subComm().MyPID()) {
	std::cout << "Vector 'auxInitialPosition at rank " << r
                  << " has contents "                      << auxInitialPosition
                  << std::endl;
      }
      m_env.subComm().Barrier();
    }
    sleep(1);
#endif

    // KAUST: all nodes in 'subComm' should have the same 'tmpChainSize'
    mpiRC = MPI_Bcast((void *) &tmpChainSize, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()",
                        "failed MPI_Bcast() for tmpChainSize");

    inputOptions.m_rawChainSize = tmpChainSize;
    uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    uqScalarSequenceClass<double> tmpLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> tmpLogTargetValues    (m_env,0,"");

    // KAUST: all nodes should call here
    uqMetropolisHastingsSGClass<P_V,P_M> mcSeqGenerator(inputOptions,
                                                        rv,
                                                        auxInitialPosition,
                                                        &unifiedCovMatrix);

    // KAUST: all nodes should call here
    mcSeqGenerator.generateSequence(tmpChain,
                                    &tmpLogLikelihoodValues, // likelihood is IMPORTANT
                                    &tmpLogTargetValues);
    cumulativeRunTime    += mcSeqGenerator.rawChainRunTime();
    cumulativeRejections += mcSeqGenerator.numRejections();

    if (m_env.inter0Rank() >= 0) {
      //for (unsigned int i = 0; i < tmpLogLikelihoodValues.subSequenceSize(); ++i) {
      //  std::cout << "tmpLogLikelihoodValues[" << i << "] = " << tmpLogLikelihoodValues[i]
      //            << ", tmpLogTargetValues["   << i << "] = " << tmpLogTargetValues[i]
      //            << std::endl;
      //}
        
      if ((m_env.subDisplayFile()             ) &&
          (m_env.displayVerbosity()   >= 0    ) &&
          (inputOptions.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateUnbLinkedChains()"
                                << ", level "               << currLevel+LEVEL_REF_ID
                                << ", chainId = "           << chainId
                                << ": finished generating " << tmpChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      // KAUST5: what if workingChain ends up with different size in different nodes? Important
      workingChain.append              (tmpChain,              1,tmpChain.subSequenceSize()-1              ); // IMPORTANT: '1' in order to discard initial position
      if (currLogLikelihoodValues) {
        currLogLikelihoodValues->append(tmpLogLikelihoodValues,1,tmpLogLikelihoodValues.subSequenceSize()-1); // IMPORTANT: '1' in order to discard initial position
      }
      if (currLogTargetValues) {
        currLogTargetValues->append    (tmpLogTargetValues,    1,tmpLogTargetValues.subSequenceSize()-1    ); // IMPORTANT: '1' in order to discard initial position
      }
    }
  } // for

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains()"
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0( // EXTRA FOR LOAD BALANCE 
  unsigned int                       currLevel,
  const std::vector<unsigned int>&   origNumChainsPerNode,
  const std::vector<unsigned int>&   origNumPositionsPerNode,
  std::vector<unsigned int>&         finalNumChainsPerNode,
  std::vector<unsigned int>&         finalNumPositionsPerNode,
  std::vector<uqExchangeInfoStruct>& exchangeStdVec)
{
  if (m_env.inter0Rank() != 0) return;

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  unsigned int Nc = exchangeStdVec.size();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": Np = " << Np
                            << ", Nc = " << Nc
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Instantiate BIP
  //////////////////////////////////////////////////////////////////////////
  glp_prob *lp; 
  lp = glp_create_prob(); 
  glp_set_prob_name(lp, "sample"); 

  //////////////////////////////////////////////////////////////////////////
  // Set rows and colums of BIP constraint matrix
  //////////////////////////////////////////////////////////////////////////
  unsigned int m = Nc+Np-1;
  unsigned int n = Nc*Np;
  unsigned int ne = Nc*Np + 2*Nc*(Np -1);

  glp_add_rows(lp, m); // Not 'm+1'
  for (int i = 1; i <= (int) Nc; ++i) {
    glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0); 
    glp_set_row_name(lp, i, ""); 
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0); 
    glp_set_row_name(lp, i, ""); 
  }
 
  glp_add_cols(lp, n); // Not 'n+1'
  for (int j = 1; j <= (int) n; ++j) {
    //glp_set_col_kind(lp, j, GLP_BV);
    glp_set_col_kind(lp, j, GLP_IV);
    glp_set_col_bnds(lp, j, GLP_DB, 0.0, 1.0); 
    glp_set_col_name(lp, j, ""); 
  }

  glp_set_obj_dir(lp, GLP_MIN); 
  for (int chainId = 0; chainId <= (int) (Nc-1); ++chainId) {
    glp_set_obj_coef(lp, (chainId*Np)+1, exchangeStdVec[chainId].numberOfPositions); 
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": finished setting BIP rows and cols"
                            << ", m = "  << m
                            << ", n = "  << n
                            << ", ne = " << ne
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Load constraint matrix
  //////////////////////////////////////////////////////////////////////////
  std::vector<int   > iVec(ne+1,0);
  std::vector<int   > jVec(ne+1,0);
  std::vector<double> aVec(ne+1,0.);
  int coefId = 1; // Yes, '1'
  for (int i = 1; i <= (int) Nc; ++i) {
    for (int j = 1; j <= (int) Np; ++j) {
      iVec[coefId] = i;
      jVec[coefId] = (i-1)*Np + j;
      aVec[coefId] = 1.;
      coefId++;
    }
  }
  for (int i = 1; i <= (int) (Np-1); ++i) {
    for (int j = 1; j <= (int) Nc; ++j) {
      iVec[coefId] = Nc+i;
      jVec[coefId] = (j-1)*Np + i;
      aVec[coefId] = -((double) exchangeStdVec[j-1].numberOfPositions);
      coefId++;

      iVec[coefId] = Nc+i;
      jVec[coefId] = (j-1)*Np + i + 1;
      aVec[coefId] = exchangeStdVec[j-1].numberOfPositions;
      coefId++;
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": finished setting BIP constraint matrix"
                            << ", ne = "     << ne
                            << ", coefId = " << coefId
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(coefId != (int) (ne+1),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "invalid final coefId");

  glp_load_matrix(lp, ne, &iVec[0], &jVec[0], &aVec[0]);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": finished loading BIP constraint matrix"
                            << ", glp_get_num_rows(lp) = " << glp_get_num_rows(lp)
                            << ", glp_get_num_cols(lp) = " << glp_get_num_cols(lp)
                            << ", glp_get_num_nz(lp) = "   << glp_get_num_nz(lp)
                            << ", glp_get_num_int(lp) = "  << glp_get_num_int(lp)
                            << ", glp_get_num_bin(lp) = "  << glp_get_num_bin(lp)
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Check BIP before solving it
  //////////////////////////////////////////////////////////////////////////
  UQ_FATAL_TEST_MACRO(glp_get_num_rows(lp) != (int) m, // Not 'm+1'
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "invalid number of rows");

  UQ_FATAL_TEST_MACRO(glp_get_num_cols(lp) != (int) n, // Not 'n+1'
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "invalid number of columnss");

  UQ_FATAL_TEST_MACRO(glp_get_num_nz(lp) != (int) ne,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "invalid number of nonzero constraint coefficients");

  UQ_FATAL_TEST_MACRO(glp_get_num_int(lp) != (int) n, // ????
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "invalid number of integer structural variables");

  UQ_FATAL_TEST_MACRO(glp_get_num_bin(lp) != (int) n,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "invalid number of binary structural variables");

  //////////////////////////////////////////////////////////////////////////
  // Set initial state
  //////////////////////////////////////////////////////////////////////////
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      int j = chainId*Np + nodeId + 1;
      if (nodeId == 0) {
        glp_set_col_stat(lp, j, GLP_BS);
      }
      else {
        glp_set_col_stat(lp, j, GLP_BS);
      }
    }
  }
#if 0
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      int j = chainId*Np + nodeId + 1;
      int initialState = glp_mip_col_val(lp, j);
      if (nodeId == 0) {
        UQ_FATAL_TEST_MACRO(initialState != 1,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                            "for nodeId = 0, initial state should be '1'");
      }
      else {
        UQ_FATAL_TEST_MACRO(initialState != 0,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                            "for nodeId > 0, initial state should be '0'");
      }
    }
  }
#endif
  for (int i = 1; i <= (int) Nc; ++i) {
    glp_set_row_stat(lp, i, GLP_NS); 
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    glp_set_row_stat(lp, i, GLP_BS); 
  }

  //glp_write_mps(lp, GLP_MPS_DECK, NULL, "nada.fixed_mps");
  //glp_write_mps(lp, GLP_MPS_FILE, NULL, "nada.free_mps" );
  //glp_write_lp (lp, NULL, "nada.cplex");

  //////////////////////////////////////////////////////////////////////////
  // Solve BIP
  //////////////////////////////////////////////////////////////////////////
  BIP_routine_struct BIP_routine_info;
  BIP_routine_info.env = &m_env;
  BIP_routine_info.currLevel = currLevel;

  glp_iocp BIP_params;
  glp_init_iocp(&BIP_params);
  BIP_params.presolve = GLP_ON;
  // aqui 2
  //BIP_params.binarize = GLP_ON;
  //BIP_params.cb_func = BIP_routine; 
  //BIP_params.cb_info = (void *) (&BIP_routine_info);
  int BIP_rc = glp_intopt(lp, &BIP_params);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": finished solving BIP"
                            << ", BIP_rc = " << BIP_rc
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(BIP_rc != 0,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "BIP returned rc != 0");

  //////////////////////////////////////////////////////////////////////////
  // Check BIP status after solution
  //////////////////////////////////////////////////////////////////////////
  int BIP_Status = glp_mip_status(lp);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": BIP_Status = " << BIP_Status
                            << std::endl;
  }

  switch (BIP_Status) {
    case GLP_OPT:
      // Ok 
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": BIP solution is optimal"
                                << std::endl;
      }
    break;

    case GLP_FEAS:
      // Ok 
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": BIP solution is guaranteed to be 'only' feasible"
                                << std::endl;
      }
    break;

    default:
      UQ_FATAL_TEST_MACRO(true,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                          "BIP has an undefined solution or has no solution");
    break;
  }

  for (int i = 1; i <= (int) Nc; ++i) {
    UQ_FATAL_TEST_MACRO(glp_mip_row_val(lp, i) != 1,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                        "row should have value 1 at solution");
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    UQ_FATAL_TEST_MACRO(glp_mip_row_val(lp, i) > 0,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                        "row should have value 0 or should be negative at solution");
  }

  //////////////////////////////////////////////////////////////////////////
  // Prepare output information, needed to for MPI distribution afterwards
  //////////////////////////////////////////////////////////////////////////
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      int j = chainId*Np + nodeId + 1;
      if (glp_mip_col_val(lp, j) == 0) {
        // Do nothing
      }
      else if (glp_mip_col_val(lp, j) == 1) {
        UQ_FATAL_TEST_MACRO(exchangeStdVec[chainId].finalNodeOfInitialPosition != -1, // important
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                            "chain has already been taken care of");
        exchangeStdVec[chainId].finalNodeOfInitialPosition = nodeId;
        finalNumChainsPerNode[nodeId]++;
        finalNumPositionsPerNode[nodeId] += exchangeStdVec[chainId].numberOfPositions;
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                            "control variable should be either '0' or '1'");
      }
    }
  }

  unsigned int origMinPosPerNode  = *std::min_element(origNumPositionsPerNode.begin(), origNumPositionsPerNode.end());
  unsigned int origMaxPosPerNode  = *std::max_element(origNumPositionsPerNode.begin(), origNumPositionsPerNode.end());
  unsigned int finalMinPosPerNode = *std::min_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  unsigned int finalMaxPosPerNode = *std::max_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": finished preparing output information"
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Printout solution information
  //////////////////////////////////////////////////////////////////////////
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": BIP solution gives the following redistribution"
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY"
                              << ", origNumChainsPerNode["     << nodeId << "] = " << origNumChainsPerNode[nodeId]
                              << ", origNumPositionsPerNode["  << nodeId << "] = " << origNumPositionsPerNode[nodeId]
                              << ", finalNumChainsPerNode["    << nodeId << "] = " << finalNumChainsPerNode[nodeId]
                              << ", finalNumPositionsPerNode[" << nodeId << "] = " << finalNumPositionsPerNode[nodeId]
                              << std::endl;
    }
    *m_env.subDisplayFile() << "  KEY"
                            << ", origRatioOfPosPerNode = "  << ((double) origMaxPosPerNode ) / ((double) origMinPosPerNode)
                            << ", finalRatioOfPosPerNode = " << ((double) finalMaxPosPerNode) / ((double)finalMinPosPerNode)
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Make sanity checks
  //////////////////////////////////////////////////////////////////////////
  UQ_FATAL_TEST_MACRO(glp_mip_obj_val(lp) != (double) finalNumPositionsPerNode[0],
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                      "Invalid objective value");

  for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, '1'
    UQ_FATAL_TEST_MACRO(finalNumPositionsPerNode[nodeId-1] < finalNumPositionsPerNode[nodeId],
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                        "Next node should have a number of positions equal or less than the current node");
  }

  for (int i = (int) (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    unsigned int nodeId = i - Nc;
    int diff = ((int) finalNumPositionsPerNode[nodeId]) - ((int) finalNumPositionsPerNode[nodeId-1]);
    UQ_FATAL_TEST_MACRO(glp_mip_row_val(lp, i) != diff,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIPAtProc0()",
                        "wrong state");
  }
    
  //////////////////////////////////////////////////////////////////////////
  // Free memory and return
  //////////////////////////////////////////////////////////////////////////
  glp_delete_prob(lp); 

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::solveBIPAtProc0()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::mpiExchangePositions( // EXTRA FOR LOAD BALANCE
  unsigned int                              currLevel,
  const uqSequenceOfVectorsClass<P_V,P_M>&  prevChain,
  const std::vector<uqExchangeInfoStruct>&  exchangeStdVec,
  const std::vector<unsigned int>&          finalNumChainsPerNode,
  const std::vector<unsigned int>&          finalNumPositionsPerNode,
  uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl)
{
  if (m_env.inter0Rank() < 0) return;

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  unsigned int Nc = exchangeStdVec.size();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSampling<P_V,P_M>::mpiExchangePositions()"
                            << ", level " << currLevel+LEVEL_REF_ID
                            << ": Np = " << Np
                            << ", Nc = " << Nc
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Each node performs:
  // --> a 'gatherv' for collecting all necessary initial positions from other nodes
  // --> a 'gatherv' for collecting all necessary chain lenghts from other nodes
  //////////////////////////////////////////////////////////////////////////
  for (unsigned int r = 0; r < Np; ++r) {
    //////////////////////////////////////////////////////////////////////////
    // Prepare some counters
    //////////////////////////////////////////////////////////////////////////
    unsigned int              numberOfInitialPositionsNodeRAlreadyHas = 0;
    std::vector<unsigned int> numberOfInitialPositionsNodeRHasToReceiveFromNode(Np,0);
    std::vector<unsigned int> indexesOfInitialPositionsNodeRHasToReceiveFromMe(0);

    unsigned int              sumOfChainLenghtsNodeRAlreadyHas = 0;
    std::vector<unsigned int> chainLenghtsNodeRHasToInherit(0);

    for (unsigned int i = 0; i < Nc; ++i) {
      if (exchangeStdVec[i].finalNodeOfInitialPosition == (int) r) {
        if (exchangeStdVec[i].originalNodeOfInitialPosition == (int) r) {
          numberOfInitialPositionsNodeRAlreadyHas++;
          sumOfChainLenghtsNodeRAlreadyHas += exchangeStdVec[i].numberOfPositions;
        }
        else {
          numberOfInitialPositionsNodeRHasToReceiveFromNode[exchangeStdVec[i].originalNodeOfInitialPosition]++;
          chainLenghtsNodeRHasToInherit.push_back(exchangeStdVec[i].numberOfPositions);
          if (m_env.inter0Rank() == exchangeStdVec[i].originalNodeOfInitialPosition) {
            indexesOfInitialPositionsNodeRHasToReceiveFromMe.push_back(exchangeStdVec[i].originalIndexOfInitialPosition);
          }
        }
      }
    }

    unsigned int totalNumberOfInitialPositionsNodeRHasToReceive = 0;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      totalNumberOfInitialPositionsNodeRHasToReceive += numberOfInitialPositionsNodeRHasToReceiveFromNode[nodeId];
    }

    unsigned int totalNumberOfChainLenghtsNodeRHasToInherit = chainLenghtsNodeRHasToInherit.size();
    unsigned int totalSumOfChainLenghtsNodeRHasToInherit = 0;
    for (unsigned int i = 0; i < totalNumberOfChainLenghtsNodeRHasToInherit; ++i) {
      totalSumOfChainLenghtsNodeRHasToInherit += chainLenghtsNodeRHasToInherit[i];
    }

    //////////////////////////////////////////////////////////////////////////
    // Printout important information
    //////////////////////////////////////////////////////////////////////////
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()"
                              << ", level " << currLevel+LEVEL_REF_ID
                              << ": r = "                                              << r
                              << ", finalNumChainsPerNode[r] = "                       << finalNumChainsPerNode[r]
                              << ", totalNumberOfInitialPositionsNodeRHasToReceive = " << totalNumberOfInitialPositionsNodeRHasToReceive
                              << ", numberOfInitialPositionsNodeRAlreadyHas = "        << numberOfInitialPositionsNodeRAlreadyHas
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()"
                              << ", level " << currLevel+LEVEL_REF_ID
                              << ": r = "                                       << r
                              << ", finalNumPositionsPerNode[r] = "             << finalNumPositionsPerNode[r]
                              << ", totalSumOfChainLenghtsNodeRHasToInherit = " << totalSumOfChainLenghtsNodeRHasToInherit
                              << ", sumOfChainLenghtsNodeRAlreadyHas = "        << sumOfChainLenghtsNodeRAlreadyHas
                              << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////
    // Make sanity checks
    //////////////////////////////////////////////////////////////////////////
    UQ_FATAL_TEST_MACRO(indexesOfInitialPositionsNodeRHasToReceiveFromMe.size() != numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()],
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()",
                        "inconsistent number of initial positions to send to node 'r'");

    UQ_FATAL_TEST_MACRO(finalNumChainsPerNode[r] != (totalNumberOfInitialPositionsNodeRHasToReceive + numberOfInitialPositionsNodeRAlreadyHas),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()",
                        "inconsistent number of chains in node 'r'");

    UQ_FATAL_TEST_MACRO(finalNumPositionsPerNode[r] != (totalSumOfChainLenghtsNodeRHasToInherit + sumOfChainLenghtsNodeRAlreadyHas),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()",
                        "inconsistent sum of chain lenghts in node 'r'");

    UQ_FATAL_TEST_MACRO(totalNumberOfInitialPositionsNodeRHasToReceive != totalNumberOfChainLenghtsNodeRHasToInherit,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()",
                        "inconsistent on total number of initial positions to receive in node 'r'");

    // Optimize use of memory (FIX ME: don't need to use swap here ????)
    indexesOfInitialPositionsNodeRHasToReceiveFromMe.resize(numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()]);
    chainLenghtsNodeRHasToInherit.resize                   (totalSumOfChainLenghtsNodeRHasToInherit);

    //////////////////////////////////////////////////////////////////////////
    // Prepare counters and buffers for gatherv of initial positions
    //////////////////////////////////////////////////////////////////////////
    unsigned int dimSize = m_vectorSpace.dimLocal();
    P_V auxInitialPosition(m_vectorSpace.zeroVector());
    std::vector<double> sendbuf(0);
    unsigned int sendcnt = 0;
    if (m_env.inter0Rank() != (int) r) {
      sendcnt = numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()] * dimSize;
      sendbuf.resize(sendcnt);
      for (unsigned int i = 0; i < numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()]; ++i) {
        unsigned int auxIndex = indexesOfInitialPositionsNodeRHasToReceiveFromMe[i];
        prevChain.getPositionValues(auxIndex,auxInitialPosition);
        for (unsigned int j = 0; j < dimSize; ++j) {
          sendbuf[i*dimSize + j] = auxInitialPosition[j];
        }
      }
    }

    std::vector<double> recvbuf(0);
    std::vector<int> recvcnts(Np,0); // '0' is already the correct value for recvcnts[r]
    if (m_env.inter0Rank() == (int) r) {
      recvbuf.resize(totalNumberOfInitialPositionsNodeRHasToReceive * dimSize);
      for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) { // Yes, from '0' on (for 'r', numberOf...ToReceiveFromNode[r] = 0 anyway)
        recvcnts[nodeId] = numberOfInitialPositionsNodeRHasToReceiveFromNode[nodeId/*m_env.inter0Rank()*/] * dimSize;
      }
    }

    std::vector<int> displs(Np,0);
    for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, from '1' on
      displs[nodeId] = displs[nodeId-1] + recvcnts[nodeId-1];
    }

    //int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
    //                void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, 
    //                int root, MPI_Comm comm )
    int mpiRC = 0;
#if 0
    if (m_env.inter0Rank() == r) {
      mpiRC = MPI_Gatherv(MPI_IN_PLACE, (int) sendcnt, MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], MPI_DOUBLE, r, m_env.inter0Comm().Comm()); // LOAD BALANCE
    }
    else {
      mpiRC = MPI_Gatherv((void *) &sendbuf[0], (int) sendcnt, MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], MPI_DOUBLE, r, m_env.inter0Comm().Comm()); // LOAD BALANCE
    }
#else
    mpiRC = MPI_Gatherv((void *) &sendbuf[0], (int) sendcnt, MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], MPI_DOUBLE, r, m_env.inter0Comm().Comm()); // LOAD BALANCE
#endif
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions()",
                        "failed MPI_Gatherv()");

    //////////////////////////////////////////////////////////////////////////
    // Make sanity checks
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // Transfer data from 'recvbuf' to 'balancedLinkControl'
    // Remember that finalNumChainsPerNode[r] = (totalNumberOfInitialPositionsNodeRHasToReceive + numberOfInitialPositionsNodeRAlreadyHas)
    // Remember that totalNumberOfInitialPositionsNodeRHasToReceive = totalNumberOfChainLenghtsNodeRHasToInherit
    //////////////////////////////////////////////////////////////////////////
    if (m_env.inter0Rank() == (int) r) {
      balancedLinkControl.balLinkedChains.resize(finalNumChainsPerNode[r]);
      unsigned int auxIndex = 0;

      for (unsigned int i = 0; i < Nc; ++i) {
        if ((exchangeStdVec[i].finalNodeOfInitialPosition    == (int) r) &&
            (exchangeStdVec[i].originalNodeOfInitialPosition == (int) r)) {
          prevChain.getPositionValues(exchangeStdVec[i].originalIndexOfInitialPosition,auxInitialPosition);
          balancedLinkControl.balLinkedChains[auxIndex].initialPosition = new P_V(auxInitialPosition);
          balancedLinkControl.balLinkedChains[auxIndex].numberOfPositions = exchangeStdVec[i].numberOfPositions;
          auxIndex++;
	}
      }

      for (unsigned int i = 0; i < totalNumberOfInitialPositionsNodeRHasToReceive; ++i) {
        for (unsigned int j = 0; j < dimSize; ++j) {
          auxInitialPosition[j] = recvbuf[i*dimSize + j];
        }
        balancedLinkControl.balLinkedChains[auxIndex].initialPosition = new P_V(auxInitialPosition);
        balancedLinkControl.balLinkedChains[auxIndex].numberOfPositions = chainLenghtsNodeRHasToInherit[i]; // aqui 3
        auxIndex++;
      }
    }

    m_env.inter0Comm().Barrier();
  } // for 'r'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::mpiExchangePositions()"
                            << ", level " << currLevel+LEVEL_REF_ID
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
