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

#ifndef __UQ_MULTI_LEVEL_SAMPLING_H__
#define __UQ_MULTI_LEVEL_SAMPLING_H__

#define ML_KAUST3

#include <uqMLSamplingOptions.h>
#include <uqFiniteDistribution.h>
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>

#define UQ_ML_SAMPLING_USES_OMEGA_LN

struct uqLinkedChainControlStruct
{
  unsigned int initialPositionIndexInPreviousChain;
  unsigned int numberOfPositions;
};

struct uqLinkedChainsPerNodeStruct
{
  std::vector<uqLinkedChainControlStruct> linkedChains;
};

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
  void   generateSequence   (uqBaseVectorSequenceClass<P_V,P_M>&             workingChain,
                             uqScalarSequenceClass<double>*                  workingLogLikelihoodValues,
                             uqScalarSequenceClass<double>*                  workingLogTargetValues);

  void   print              (std::ostream& os) const;

private:
  void   sampleIndexes      (unsigned int                                    currLevel,                         // input
                             unsigned int                                    unifiedRequestedNumSamples,        // input
                             const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only, // input
                             unsigned int                                    indexOfFirstWeight,                // input
                             unsigned int                                    indexOfLastWeight,                 // input
                             unsigned int&                                   modifiedSubNumSamples,             // output
                             std::vector<unsigned int>&                      unifiedIndexCountersAtAllProcs);   // output

  void   distribIndexSamples(unsigned int                                    currLevel,                      // input
                             unsigned int                                    subNumSamples,                  // input
                             unsigned int                                    indexOfFirstWeight,             // input
                             unsigned int                                    indexOfLastWeight,              // input
                             std::vector<unsigned int>&                      unifiedIndexCountersAtAllProcs, // input, modified
                             std::vector<uqLinkedChainsPerNodeStruct>&       linkControl);                   // output

  void   generateChain      (unsigned int                                    currLevel,               // input
                             uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
                             const std::vector<uqLinkedChainsPerNodeStruct>& linkControl,             // input
                             unsigned int                                    indexOfFirstWeight,      // input
                             const P_M&                                      unifiedCovMatrix,        // input
                             const uqGenericVectorRVClass  <P_V,P_M>&        rv,                      // input
                             const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,               // input
                             uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,            // output
                             double&                                         cumulativeRunTime,       // output
                             unsigned int&                                   cumulativeRejections,    // output
                             uqScalarSequenceClass         <double>*         currLogLikelihoodValues, // output
                             uqScalarSequenceClass         <double>*         currLogTargetValues);    // output

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
    // Step 1 of 10: read options
    //***********************************************************
    sprintf(tmpSufix,"%d_",currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass* currOptions = new uqMLSamplingLevelOptionsClass(m_env,(m_options.m_prefix + tmpSufix).c_str());
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 1 of 10"
                                << std::endl;
      }

      currOptions->scanOptionsValues(&defaultLevelOptions);
      unsigned int tmpSize = currOptions->m_rawChainSize;
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
    // Step 2 of 10: save [chain and corresponding target pdf values] from previous level
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
                                << ": beginning step 2 of 10"
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
    // Step 3 of 10: compute [currExponent and sequence of weights] for current level
    //***********************************************************
    uqScalarSequenceClass<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 3 of 10"
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
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN
      uqScalarSequenceClass<double> omegaLnDiffSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");

#ifdef ML_KAUST3
      double nowUnifiedEvidenceLnFactor = 0.;
#else
      double nowEvidenceLnFactor = 0.;
#endif

#else
      double nowEvidenceFactor = 0.;
#endif
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
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN

#ifdef ML_KAUST3
        double subWeightRatioSum     = 0.;
        double unifiedWeightRatioSum = 0.;
#else
        double weightRatioSum = 0.;
#endif

#else
        double weightSum = 0.;
#endif
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN
          omegaLnDiffSequence[i] = prevLogLikelihoodValues[i]*auxExponent; // likelihood is important
#else
          weightSequence[i] = exp(prevLogLikelihoodValues[i]*auxExponent);
          weightSum += weightSequence[i];
#endif
          //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          //  *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
          //                          << ", level "              << currLevel+LEVEL_REF_ID
          //                          << ": auxExponent = "      << auxExponent
          //                          << ", 'prev'TargetValues[" << i
          //                          << "] = "                  << prevTargetValues[i]
          //                          << ", weightSequence["     << i
          //                          << "] = "                  << weightSequence[i]
          //                          << ", weightSum = "        << weightSum
          //                          << std::endl;
          //}
        }
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN

#ifdef ML_KAUST3
        double unifiedOmegaLnMax = 0.;
        double unifiedOmegaLnMin = 0.;
        omegaLnDiffSequence.unifiedMinMax(m_vectorSpace.numOfProcsForStorage() == 1, // KAUST3
                                          0,
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
#else
        double omegaLnMax = 0.;
        double omegaLnMin = 0.;
        omegaLnDiffSequence.unifiedMinMax(m_vectorSpace.numOfProcsForStorage() == 1,
                                          0,
                                          omegaLnMin,
                                          omegaLnMax);
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] -= omegaLnMax;
          weightSequence[i] = exp(omegaLnDiffSequence[i]);
          weightRatioSum += weightSequence[i];
        }
        nowEvidenceLnFactor = log(weightRatioSum) + omegaLnMax - log(weightSequence.subSequenceSize());
#endif

#else
        nowEvidenceFactor = weightSum/prevChain.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
#endif
        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN

#ifdef ML_KAUST3
          weightSequence[i] /= unifiedWeightRatioSum;
#else
          weightSequence[i] /= weightRatioSum;
#endif

#else
          weightSequence[i] /= weightSum;
#endif
          effectiveSampleSize += weightSequence[i]*weightSequence[i];
          //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          //  *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
          //                          << ", level "                 << currLevel+LEVEL_REF_ID
          //                          << ": i = "                   << i
          //                          << ", effectiveSampleSize = " << effectiveSampleSize
          //                          << std::endl;
          //}
        }
#ifdef ML_KAUST3
        double subQuantity = effectiveSampleSize;
        effectiveSampleSize = 0.;
        mpiRC = MPI_Allreduce((void *) &subQuantity, (void *) &effectiveSampleSize, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for effective sample size");
#endif
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
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN

#ifdef ML_KAUST3
      m_logEvidenceFactors.push_back(nowUnifiedEvidenceLnFactor);
#else
      m_logEvidenceFactors.push_back(nowEvidenceLnFactor);
#endif

#else
      m_logEvidenceFactors.push_back(log(nowEvidenceFactor));
#endif

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

      if (currExponent == 1.) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+LEVEL_REF_ID
                                  << ": copying 'last' level options to current options"
                                  << std::endl;
        }
        delete currOptions;
        currOptions = &lastLevelOptions;
        unsigned int tmpSize = currOptions->m_rawChainSize;
        int mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &unifiedRequestedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for requested num samples in step 3");
      }

      // Make sure all nodes in 'inter0Comm' have the same value of 'logEvidenceFactor'
      uqMiscCheckForSameValueInAllNodes(m_logEvidenceFactors[m_logEvidenceFactors.size()-1],
                                        1.e-16,
                                        m_env.inter0Comm(),
                                        "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, logEvidenceFactor");
    } // end of step 3
    // KAUST: all nodes in 'subComm' should have the same 'currExponent'
    int mpiRC = MPI_Bcast((void *) &currExponent, (int) 1, MPI_DOUBLE, 0, m_env.subComm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "failed MPI_Bcast() for currExponent");

    //***********************************************************
    // Step 4 of 10: create covariance matrix for current level
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
                                << ": beginning step 4 of 10"
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
    // Step 5 of 10: create *unified* finite distribution for current level
    //***********************************************************
    std::vector<unsigned int> unifiedIndexCountersAtAllProcs(0); // It will be resized by 'sampleIndexes()' below
    unsigned int modifiedSubNumSamples = 0;
    std::vector<double> unifiedWeightStdVectorAtProc0Only(0); // KAUST, to check
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 5 of 10"
                                << std::endl;
      }

      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);
      sampleIndexes(currLevel,                         // input
                    unifiedRequestedNumSamples,        // input
                    unifiedWeightStdVectorAtProc0Only, // input
                    indexOfFirstWeight,                // input
                    indexOfLastWeight,                 // input
                    modifiedSubNumSamples,             // output
                    unifiedIndexCountersAtAllProcs);   // output

      UQ_FATAL_TEST_MACRO(unifiedIndexCountersAtAllProcs.size() != weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "wrong output from sampleIndexes() in step 5");
    } // end of step 5

    //***********************************************************
    // Step 6 of 10: plan for number of linked chains for each node so that
    //               each node generates the same number of positions
    //***********************************************************
    std::vector<uqLinkedChainsPerNodeStruct> linkControl(0); // KAUST
    if (m_env.inter0Rank() >= 0) { // KAUST
      linkControl.resize(m_env.inter0Comm().NumProc()); // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 6 of 10"
                                << std::endl;
      }

      distribIndexSamples(currLevel,                      // input
                          modifiedSubNumSamples,          // input
                          indexOfFirstWeight,             // input
                          indexOfLastWeight,              // input
                          unifiedIndexCountersAtAllProcs, // input, modified
                          linkControl);                   // output

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": linkControl[m_env.inter0Rank()].linkedChains.size() = " << linkControl[m_env.inter0Rank()].linkedChains.size()
                                << std::endl;
      }
    } // end of step 6

    //***********************************************************
    // Step 7 of 10: load balancing = subChainSize [indexes + corresponding positions] for each node
    //***********************************************************
    if (m_env.inter0Rank() >= 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 7 of 10"
                                << std::endl;
      }

      // KAUST5: important
      //if (m_env.inter0Comm().NumProc() > 1) {
      //  UQ_FATAL_TEST_MACRO(true,
      //                      m_env.fullRank(),
      //                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
      //                      "incomplete code for load balancing");
      //}
    } // end of step 7
  
    //***********************************************************
    // Step 8 of 10: create vector RV for current level
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
                                << ": beginning step 8 of 10"
                                << std::endl;
      }

      currRv.setPdf(currPdf);
    } // end of step 8

    //***********************************************************
    // Step 9 of 10: scale unified covariance matrix until min <= rejection rate <= max
    //***********************************************************
    if (currOptions->m_scaleCovMatrix == false) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": skipping step 9 of 10"
                                << std::endl;
      }
    }
    else {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 9 of 10"
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
        } // KAUST
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

        std::vector<unsigned int> nowUnifiedIndexCountersAtAllProcs(0); // It will be resized by 'sampleIndexes()' below

        unsigned int tmpSubNumSamples = 0;
        if (m_env.inter0Rank() >= 0) { // KAUST
          sampleIndexes(currLevel,                                          // input
                        originalSubNumSamples*m_env.inter0Comm().NumProc(), // input
                        unifiedWeightStdVectorAtProc0Only,                  // input
                        indexOfFirstWeight,                                 // input
                        indexOfLastWeight,                                  // input
                        tmpSubNumSamples,                                   // output
                        nowUnifiedIndexCountersAtAllProcs);                 // output

          UQ_FATAL_TEST_MACRO(nowUnifiedIndexCountersAtAllProcs.size() != weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1),
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "wrong output from sampleIndexes() in step 9");

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level " << currLevel+LEVEL_REF_ID
                                    << ": in loop for assessing rejection rate"
                                    << ", about to distribute sampled assessment indexes"
                                    << std::endl;
          }
        } // KAUST

        std::vector<uqLinkedChainsPerNodeStruct> nowLinkControl(0); // KAUST

        if (m_env.inter0Rank() >= 0) { // KAUST
          nowLinkControl.resize(m_env.inter0Comm().NumProc()); // KAUST
          distribIndexSamples(currLevel,                         // input
                              tmpSubNumSamples,                  // input
                              indexOfFirstWeight,                // input
                              indexOfLastWeight,                 // input
                              nowUnifiedIndexCountersAtAllProcs, // input, modified
                              nowLinkControl);                   // output
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
        currOptions->m_rawChainSize          = 0; // will be set inside generateChain()
        currOptions->m_rawChainComputeStats  = false;
        currOptions->m_filteredChainGenerate = false;
        currOptions->m_drMaxNumExtraStages   = 0;
        currOptions->m_amAdaptInterval       = 0;

        // KAUST: all nodes should call here
        generateChain(currLevel,          // input
                      *currOptions,       // input, only m_rawChainSize changes
                      nowLinkControl,     // input
                      indexOfFirstWeight, // input
                      nowCovMatrix,       // input
                      currRv,             // input
                      prevChain,          // input
                      nowChain,           // output 
                      nowRunTime,         // output
                      nowRejections,      // output
                      NULL,               // output
                      NULL);              // output

        // KAUST: all nodes should call here
        currOptions->m_totallyMute           = savedTotallyMute;
        currOptions->m_rawChainSize          = savedRawChainSize;
        currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
        currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME
        currOptions->m_drMaxNumExtraStages   = savedDrMaxNumExtraStages;
        currOptions->m_amAdaptInterval       = savedAmAdaptInterval;

        if (m_env.inter0Rank() >= 0) { // KAUST
          // If only one cov matrix is used, then the rejection should be assessed among all inter0Comm nodes // KAUST3
          unsigned int nowUnifiedRejections = 0;
          int mpiRC = MPI_Allreduce((void *) &nowRejections, (void *) &nowUnifiedRejections, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Allreduce() for now rejections");

          unsigned int unifiedNumSamples = 0;
          mpiRC = MPI_Allreduce((void *) &tmpSubNumSamples, (void *) &unifiedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Allreduce() for num samples in step 9");

          nowRejectionRate = ((double) nowUnifiedRejections) / ((double) unifiedNumSamples);

          //bool aux1 = (nowRejectionRate == meanRejectionRate);
          bool aux2 = (nowRejectionRate >= currOptions->m_minRejectionRate)
                      &&
                      (nowRejectionRate <= currOptions->m_maxRejectionRate);
          testResult = aux2;

          // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
          uqMiscCheckForSameValueInAllNodes(testResult,
                                            0.,
                                            m_env.inter0Comm(),
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 9, testResult");
        }

        // KAUST: all nodes in 'subComm' should have the same 'testResult'
        unsigned int tmpUint = (unsigned int) testResult;
        int mpiRC = MPI_Bcast((void *) &tmpUint, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm());
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
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 9, testResult");
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
    } // end of step 9

    //***********************************************************
    // Step 10 of 10: sample vector RV of current level
    //***********************************************************
    unsigned int unifiedNumberOfRejections = 0;
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 10 of 10"
                                << std::endl;
      }

      // KAUST: all nodes should call here
      bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
      unsigned int savedRawChainSize          = currOptions->m_rawChainSize; // Ok to use rawChainSize
      bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
      bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;

      currOptions->m_totallyMute           = true;
      currOptions->m_rawChainSize          = 0; // will be set inside generateChain()
      currOptions->m_rawChainComputeStats  = false;
      currOptions->m_filteredChainGenerate = false;

      // KAUST: all nodes should call here
      generateChain(currLevel,                    // input
                    *currOptions,                 // input, only m_rawChainSize changes
                    linkControl,                  // input
                    indexOfFirstWeight,           // input
                    *unifiedCovMatrix,            // input
                    currRv,                       // input
                    prevChain,                    // input
                    currChain,                    // output
                    cumulativeRawChainRunTime,    // output
                    cumulativeRawChainRejections, // output
                    &currLogLikelihoodValues,     // output // likelihood is important
                    &currLogTargetValues);        // output

      // KAUST: all nodes should call here
      currOptions->m_totallyMute           = savedTotallyMute;
      currOptions->m_rawChainSize          = savedRawChainSize;
      currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
      currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

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
        }

        // Check if unified size of generated chain matches the unified requested size // KAUST
        unsigned int tmpSize = currChain.subSequenceSize();
        unsigned int unifiedGeneratedNumSamples = 0;
        unsigned int mpiRC = MPI_Allreduce((void *) &tmpSize, (void *) &unifiedGeneratedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for generated num samples in step 10");
        UQ_FATAL_TEST_MACRO(unifiedGeneratedNumSamples != unifiedRequestedNumSamples,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "currChain (linked one) has been generated with invalid size");

        // Compute unified number of rejections
        mpiRC = MPI_Allreduce((void *) &cumulativeRawChainRejections, (void *) &unifiedNumberOfRejections, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "failed MPI_Allreduce() for number of rejections");
      } // KAUST

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }
    } // end of step 10

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
uqMLSamplingClass<P_V,P_M>::sampleIndexes(
  unsigned int               currLevel,                         // input
  unsigned int               unifiedRequestedNumSamples,        // input
  const std::vector<double>& unifiedWeightStdVectorAtProc0Only, // input
  unsigned int               indexOfFirstWeight,                // input
  unsigned int               indexOfLastWeight,                 // input
  unsigned int&              modifiedSubNumSamples,             // output
  std::vector<unsigned int>& unifiedIndexCountersAtAllProcs)    // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::sampleIndexes()"
                            << ": unifiedRequestedNumSamples = "               << unifiedRequestedNumSamples
                            << ", unifiedWeightStdVectorAtProc0Only.size() = " << unifiedWeightStdVectorAtProc0Only.size()
                            << ", indexOfFirstWeight = "                       << indexOfFirstWeight
                            << ", indexOfLastWeight = "                        << indexOfLastWeight
                            << std::endl;
  }

  // All nodes in 'inter0Comm' should resize to the same size // KAUST3
  unsigned int resizeSize = unifiedWeightStdVectorAtProc0Only.size();
  int mpiRC = MPI_Bcast((void *) &resizeSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                      "failed MPI_Bcast() for resizeSize");
  unifiedIndexCountersAtAllProcs.resize(resizeSize,0);

  if (m_env.inter0Rank() >= 0) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::sampleIndexes()"
                              << ", level " << currLevel+LEVEL_REF_ID
                              << ": unifiedWeightStdVectorAtProc0Only.size() = " << unifiedWeightStdVectorAtProc0Only.size()
                              << std::endl;

      //unsigned int numZeros = 0;
      //for (unsigned int i = 0; i < unifiedWeightStdVectorAtProc0Only.size(); ++i) {
      //  *m_env.subDisplayFile() << "unifiedWeightStdVectorAtProc0Only[" << i
      //                          << "] = " << unifiedWeightStdVectorAtProc0Only[i]
      //                          << std::endl;
      //  if (unifiedWeightStdVectorAtProc0Only[i] == 0.) numZeros++;
      //}
      //*m_env.subDisplayFile() << "Number of zeros in unifiedWeightStdVectorAtProc0Only = " << numZeros
      //                        << std::endl;
    }
 
    if (m_env.inter0Rank() == 0) {
      // Generate 'unifiedNumSamples' samples from 'tmpFD'
      uqFiniteDistributionClass tmpFd(m_env,
                                      "",
                                      unifiedWeightStdVectorAtProc0Only);
      for (unsigned int i = 0; i < unifiedRequestedNumSamples; ++i) {
        unsigned int index = tmpFd.sample();
        unifiedIndexCountersAtAllProcs[index] += 1;
      }
    }

    mpiRC = MPI_Bcast((void *) &unifiedIndexCountersAtAllProcs[0], (int) unifiedIndexCountersAtAllProcs.size(), MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                        "failed MPI_Bcast() for unified index counters");
#if 0 // Use allgatherv ??? for modifiedNumSamples instead
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::sampleIndexes()"
                              << ", level " << currLevel+LEVEL_REF_ID
                              << ":"
                              << std::endl;
      for (int r = 0; r < m_env.inter0Comm().NumProc(); ++r) {
        *m_env.subDisplayFile() << "  unifiedIndexCountersAtAllProcs[" << r << "] = " << unifiedIndexCountersAtAllProcs[r]
                                << std::endl;
      }
    }
#endif
    // Use 'indexOfFirstWeight' and 'indexOfLastWeight' in order to update 'modifiedSubNumSamples'
    UQ_FATAL_TEST_MACRO(indexOfFirstWeight >= unifiedIndexCountersAtAllProcs.size(),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                        "invalid indexOfFirstWeight");
    UQ_FATAL_TEST_MACRO(indexOfLastWeight >= unifiedIndexCountersAtAllProcs.size(),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                        "invalid indexOfLastWeight");
    modifiedSubNumSamples = 0;
    for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
      modifiedSubNumSamples += unifiedIndexCountersAtAllProcs[i];
    }
  }

  //for (unsigned int i = 0; i < unifiedIndexCountersAtAllProcs.size(); ++i) {
  //  *m_env.subDisplayFile() << "unifiedIndexCountersAtAllProcs[" << i
  //                          << "] = " << unifiedIndexCountersAtAllProcs[i]
  //                          << std::endl;
  //}

#if 0 // for debug only
  //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
  //                void *recvbuf, int recvcount, MPI_Datatype recvtype, 
  //                int root, MPI_Comm comm )
  std::vector<unsigned int> recvcnts(m_env.inter0Comm().NumProc(),0);
  mpiRC = MPI_Gather((void *) &modifiedSubNumSamples, 1, MPI_UNSIGNED, (void *) &recvcnts[0], (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                      "failed MPI_Gather()");
#endif

  std::vector<unsigned int> auxBuf(1,0);

  unsigned int minModifiedSubNumSamples = 0;
  auxBuf[0] = modifiedSubNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                      "failed MPI_Allreduce() for min");

  unsigned int maxModifiedSubNumSamples = 0;
  auxBuf[0] = modifiedSubNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                      "failed MPI_Allreduce() for max");

  unsigned int sumModifiedSubNumSamples = 0;
  auxBuf[0] = modifiedSubNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                      "failed MPI_Allreduce() for sum");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::sampleIndexes()"
                            << ", level "                                   << currLevel+LEVEL_REF_ID
                            << ": modifiedSubNumSamples = "                 << modifiedSubNumSamples
                            << ", unifiedIndexCountersAtAllProcs.size() = " << unifiedIndexCountersAtAllProcs.size()
                            << std::endl;
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::sampleIndexes()"
                            << ", level "                      << currLevel+LEVEL_REF_ID
                            << ": minModifiedSubNumSamples = " << minModifiedSubNumSamples
                            << ", avgModifiedSubNumSamples = " << ((double) sumModifiedSubNumSamples)/((double) m_env.inter0Comm().NumProc())
                            << ", maxModifiedSubNumSamples = " << maxModifiedSubNumSamples
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::distribIndexSamples(
  unsigned int                              currLevel,                      // input
  unsigned int                              subNumSamples,                  // input
  unsigned int                              indexOfFirstWeight,             // input
  unsigned int                              indexOfLastWeight,              // input
  std::vector<unsigned int>&                unifiedIndexCountersAtAllProcs, // input, modified
  std::vector<uqLinkedChainsPerNodeStruct>& linkControl)                    // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::distribIndexSamples()"
                            << ": subNumSamples = "                         << subNumSamples
                            << ", indexOfFirstWeight = "                    << indexOfFirstWeight
                            << ", indexOfLastWeight = "                     << indexOfLastWeight
                            << ", unifiedIndexCountersAtAllProcs.size() = " << unifiedIndexCountersAtAllProcs.size()
                            << std::endl;
  }

  unsigned int auxNode = (unsigned int) m_env.inter0Rank(); // 0 // KAUST4
  unsigned int numberOfPositionsToGuaranteeForNode = subNumSamples;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::distribIndexSamples()"
                            << ", level "                                 << currLevel+LEVEL_REF_ID
                            << ": numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                            << std::endl;
  }
  for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
//for (unsigned int i = 0; i < unifiedIndexCountersAtAllProcs.size(); ++i) { // KAUST4: important
    while (unifiedIndexCountersAtAllProcs[i] != 0) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 30)) {
        *m_env.subDisplayFile() << "auxNode = "                               << auxNode
                                << ", numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                                << ", unifiedIndexCountersAtAllProcs["        << i
                                << "] = "                                     << unifiedIndexCountersAtAllProcs[i]
                                << std::endl;
      }
      UQ_FATAL_TEST_MACRO(auxNode >= (unsigned int) m_env.inter0Comm().NumProc(),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                          "auxNode got too large");
      if (unifiedIndexCountersAtAllProcs[i] < numberOfPositionsToGuaranteeForNode) {
        uqLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = unifiedIndexCountersAtAllProcs[i];
        linkControl[auxNode].linkedChains.push_back(auxControl);

        numberOfPositionsToGuaranteeForNode -= unifiedIndexCountersAtAllProcs[i];
        unifiedIndexCountersAtAllProcs[i] = 0;
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] > 0                                   )) {
      //else { // KAUST4
        uqLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = numberOfPositionsToGuaranteeForNode;
        linkControl[auxNode].linkedChains.push_back(auxControl);

        unifiedIndexCountersAtAllProcs[i] -= numberOfPositionsToGuaranteeForNode;
        numberOfPositionsToGuaranteeForNode = 0;

        // Go to next node
        //auxNode++; // KAUST4
        //numberOfPositionsToGuaranteeForNode = subNumSamples; // KAUST4
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] == 0                                  )) {
        // Ok
      }
      else {
        UQ_FATAL_TEST_MACRO(true, // KAUST4
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                            "should never get here");
      }
    }
  }
  UQ_FATAL_TEST_MACRO(auxNode != (unsigned int) m_env.inter0Rank(), // m_env.inter0Comm().NumProc(), // KAUST4
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                      "auxNode exited loop with wrong value");
  UQ_FATAL_TEST_MACRO(numberOfPositionsToGuaranteeForNode != 0, // subNumSamples, // KAUST4
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                      "numberOfPositionsToGuaranteeForNode exited loop with wrong value");
  // FIX ME: swap trick to save memory

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::distribIndexSamples()"
                            << ", level "                << currLevel+LEVEL_REF_ID
                            << ": linkControl.size() = " << linkControl.size()
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateChain(
  unsigned int                                    currLevel,               // input
  uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
  const std::vector<uqLinkedChainsPerNodeStruct>& linkControl,             // input
  unsigned int                                    indexOfFirstWeight,      // input
  const P_M&                                      unifiedCovMatrix,        // input
  const uqGenericVectorRVClass  <P_V,P_M>&        rv,                      // input
  const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,               // input
  uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,            // output
  double&                                         cumulativeRunTime,       // output
  unsigned int&                                   cumulativeRejections,    // output
  uqScalarSequenceClass         <double>*         currLogLikelihoodValues, // output
  uqScalarSequenceClass         <double>*         currLogTargetValues)     // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateChain()"
                            << ": linkControl.size() = " << linkControl.size()
                            << ", indexOfFirstWeight = " << indexOfFirstWeight
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = linkControl[m_env.inter0Rank()].linkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  int mpiRC = MPI_Bcast((void *) &chainIdMax, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateChain()",
                      "failed MPI_Bcast() for chainIdMax");

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfUsefulSamples = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfUsefulSamples += linkControl[m_env.inter0Rank()].linkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfUsefulSamples = 0;
    auxBuf[0] = numberOfUsefulSamples;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minNumberOfUsefulSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateChain()",
                        "failed MPI_Allreduce() for min");

    unsigned int maxNumberOfUsefulSamples = 0;
    auxBuf[0] = numberOfUsefulSamples;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxNumberOfUsefulSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateChain()",
                        "failed MPI_Allreduce() for max");

    unsigned int sumNumberOfUsefulSamples = 0;
    auxBuf[0] = numberOfUsefulSamples;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumNumberOfUsefulSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateChain()",
                        "failed MPI_Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateChain()"
                              << ", level "                   << currLevel+LEVEL_REF_ID
                              << ": chainIdMax = "            << chainIdMax
                              << ", numberOfUsefulSamples = " << numberOfUsefulSamples
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateChain()"
                              << ", level "                      << currLevel+LEVEL_REF_ID
                              << ": minNumberOfUsefulSamples = " << minNumberOfUsefulSamples
                              << ", avgNumberOfUsefulSamples = " << ((double) sumNumberOfUsefulSamples)/((double) m_env.inter0Comm().NumProc())
                              << ", maxNumberOfUsefulSamples = " << maxNumberOfUsefulSamples
                              << std::endl;
    }
  }
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      unsigned int auxIndex = linkControl[m_env.inter0Rank()].linkedChains[chainId].initialPositionIndexInPreviousChain - indexOfFirstWeight; // KAUST4
      prevChain.getPositionValues(auxIndex,auxInitialPosition);
      tmpChainSize = linkControl[m_env.inter0Rank()].linkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::generateChain()"
                                << ": chainId = "      << chainId
                                << ", tmpChainSize = " << tmpChainSize
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm().Comm()); // KAUST
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
    mpiRC = MPI_Bcast((void *) &tmpChainSize, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateChain()",
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
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateChain()"
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
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateChain()"
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING_H__
