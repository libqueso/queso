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
  void   sampleIndexes      (unsigned int                                    subNumSamples,          // input
                             const std::vector<double>&                      unifiedWeightStdVector, // input
                             std::vector<unsigned int>&                      unifiedIndexCounters);  // output

  void   distribIndexSamples(unsigned int                                    subNumSamples,        // input
                             std::vector<unsigned int>&                      unifiedIndexCounters, // input, modified
                             std::vector<uqLinkedChainsPerNodeStruct>&       nodes);               // output

  void   generateChain      (uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
                             const std::vector<uqLinkedChainsPerNodeStruct>& nodes,                   // input
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
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level "              << currLevel+LEVEL_REF_ID
                              << ", currOptions.m_rawChainSize = " << currOptions.m_rawChainSize
                              << ", currExponent = "               << currExponent
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);

    currChain.setName              (currOptions.m_prefix + "rawChain");
    currLogLikelihoodValues.setName(currOptions.m_prefix + "rawLogLikelihood");
    currLogTargetValues.setName    (currOptions.m_prefix + "rawLogTarget");

    currChain.resizeSequence              (currOptions.m_rawChainSize);
    currLogLikelihoodValues.resizeSequence(currOptions.m_rawChainSize);
    currLogTargetValues.resizeSequence    (currOptions.m_rawChainSize);

    P_V auxVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 0; i < currChain.subSequenceSize(); ++i) {
      //std::cout << "In QUESO: before prior realizer with i = " << i << std::endl;
      m_priorRv.realizer().realization(auxVec);
      currChain.setPositionValues(i,auxVec);
      // KAUST: all processors should call here
      currLogLikelihoodValues[i] = m_likelihoodFunction.lnValue(auxVec,NULL,NULL,NULL,NULL);  // likelihood is important
      currLogTargetValues[i]     = m_priorRv.pdf().lnValue(auxVec,NULL,NULL,NULL,NULL) + currLogLikelihoodValues[i];
      //std::cout << "In QUESO: currLogTargetValues[" << i << "] = " << currLogTargetValues[i] << std::endl;
    }

    if (currOptions.m_rawChainComputeStats) {
      std::ofstream* genericOfsVar = NULL;
      m_env.openOutputFile(currOptions.m_dataOutputFileName,
                           UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                           currOptions.m_dataOutputAllowedSet,
                           false,
                           genericOfsVar);

      currChain.computeStatistics(*currOptions.m_rawChainStatisticalOptions,
                                  genericOfsVar);

      genericOfsVar->close();
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

    //UQ_FATAL_TEST_MACRO((currChain.subSequenceSize() != currOptions.m_rawChainSize),
    //                    m_env.fullRank(),
    //                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
    //                    "currChain (first one) has been generated with invalid size");
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
    if (m_env.inter0Rank() == 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 1 of 10"
                                << std::endl;
      }

      currOptions->scanOptionsValues(&defaultLevelOptions);

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ", currOptions->m_rawChainSize = " << currOptions->m_rawChainSize
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

    if (m_env.inter0Rank() == 0) { // KAUST
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

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": prevChain.unifiedSequenceSize() = " << prevChain.unifiedSequenceSize()
                                << ", currChain.unifiedSequenceSize() = " << currChain.unifiedSequenceSize()
                                << ", prevLogLikelihoodValues.unifiedSequenceSize() = " << prevLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1)
                                << ", currLogLikelihoodValues.unifiedSequenceSize() = " << currLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1)
                                << ", prevLogTargetValues.unifiedSequenceSize() = " << prevLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1)
                                << ", currLogTargetValues.unifiedSequenceSize() = " << currLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1)
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
    } // end of step 2

    //***********************************************************
    // Step 3 of 10: compute [currExponent and sequence of weights] for current level
    //***********************************************************
    uqScalarSequenceClass<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() == 0) { // KAUST
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
      double nowEvidenceLnFactor = 0.;
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
        double weightRatioSum = 0.;
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
        double omegaLnMax = omegaLnDiffSequence.subMax(0,omegaLnDiffSequence.subSequenceSize()); // FIX ME: unifiedMax()
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] -= omegaLnMax;
          weightSequence[i] = exp(omegaLnDiffSequence[i]);
          weightRatioSum += weightSequence[i];
        }
        nowEvidenceLnFactor = log(weightRatioSum) + omegaLnMax - log(weightSequence.subSequenceSize());
#else
        nowEvidenceFactor = weightSum/prevChain.unifiedSequenceSize(); // FIX ME: unified
#endif
        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN
          weightSequence[i] /= weightRatioSum;
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
        effectiveSampleSize = 1./effectiveSampleSize;
        nowEffectiveSizeRatio = effectiveSampleSize/((double) weightSequence.subSequenceSize()); // FIX ME: unified
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
      } while (testResult == false);
      currExponent = nowExponent;
#ifdef UQ_ML_SAMPLING_USES_OMEGA_LN
      m_logEvidenceFactors.push_back(nowEvidenceLnFactor);
#else
      m_logEvidenceFactors.push_back(log(nowEvidenceFactor));
#endif

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << currLevel+LEVEL_REF_ID
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1)
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
      }
    } // end of step 3

    //***********************************************************
    // Step 4 of 10: create covariance matrix for current level
    //***********************************************************
    P_M* unifiedCovMatrix = m_vectorSpace.newMatrix();
    if (m_env.inter0Rank() == 0) { // KAUST
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
      P_M* subCovMatrix = m_vectorSpace.newMatrix();
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        diffVec = auxVec - weightedMeanVec;
        *subCovMatrix += weightSequence[i]*matrixProduct(diffVec,diffVec);
      }

      for (unsigned int i = 0; i < unifiedCovMatrix->numRowsLocal(); ++i) {
        for (unsigned int j = 0; j < unifiedCovMatrix->numCols(); ++j) {
          double localValue = (*subCovMatrix)(i,j);
          double sumValue = 0.;
          if (m_env.inter0Rank() == 0) {
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
      delete subCovMatrix;

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
    std::vector<unsigned int> unifiedIndexCounters(0); // It will be resized by 'sampleIndexes()' below
    if (m_env.inter0Rank() == 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 5 of 10"
                                << std::endl;
      }

      std::vector<double> unifiedWeightStdVector(0);
      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVector);
      sampleIndexes(currOptions->m_rawChainSize, // input
                    unifiedWeightStdVector,      // input
                    unifiedIndexCounters);       // output

      UQ_FATAL_TEST_MACRO(unifiedIndexCounters.size() != weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "wrong output from sampleIndexes() in step 5");
    } // end of step 5

    //***********************************************************
    // Step 6 of 10: plan for number of linked chains for each node so that
    //               each node generates the same number of positions
    //***********************************************************
    std::vector<uqLinkedChainsPerNodeStruct> nodes(m_env.inter0Comm().NumProc());
    if (m_env.inter0Rank() == 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 6 of 10"
                                << std::endl;
      }

      distribIndexSamples(currOptions->m_rawChainSize, // input
                          unifiedIndexCounters,        // input, modified
                          nodes);                      // output

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": nodes[m_env.subId()].linkedChains.size() = " << nodes[m_env.subId()].linkedChains.size()
                                << std::endl;
      }
    } // end of step 6

    //***********************************************************
    // Step 7 of 10: load balancing = subChainSize [indexes + corresponding positions] for each node
    //***********************************************************
    if (m_env.inter0Rank() == 0) { // KAUST
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 7 of 10"
                                << std::endl;
      }

      if (m_env.inter0Comm().NumProc() > 1) {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "incomplete code for load balancing");
      }
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

    if (m_env.inter0Rank() == 0) { // KAUST
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

      std::vector<double> unifiedWeightStdVector(0);
      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVector);

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

        if (m_env.inter0Rank() == 0) { // KAUST
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

        unsigned int subNumSamples = 1 + (unsigned int) ( (1.-meanRejectionRate)/meanRejectionRate/currOptions->m_covRejectionRate/currOptions->m_covRejectionRate );

        if (m_env.inter0Rank() == 0) { // KAUST
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+LEVEL_REF_ID
                                  << ": in loop for assessing rejection rate"
                                  << ", about to sample " << subNumSamples << " indexes"
                                  << ", meanRejectionRate = " << meanRejectionRate
                                  << ", covRejectionRate = "  << currOptions->m_covRejectionRate
                                  << std::endl;
        }
        } // KAUST

        std::vector<unsigned int> nowUnifiedIndexCounters(0); // It will be resized by 'sampleIndexes()' below

        if (m_env.inter0Rank() == 0) { // KAUST
        sampleIndexes(subNumSamples,            // input
                      unifiedWeightStdVector,   // input
                      nowUnifiedIndexCounters); // output

        UQ_FATAL_TEST_MACRO(nowUnifiedIndexCounters.size() != weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1),
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

        std::vector<uqLinkedChainsPerNodeStruct> nowNodes(m_env.inter0Comm().NumProc());

        if (m_env.inter0Rank() == 0) { // KAUST
        distribIndexSamples(subNumSamples,           // input
                            nowUnifiedIndexCounters, // input, modified
                            nowNodes);               // output
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

        // KAUST: all processors should call here
        bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
        unsigned int savedRawChainSize          = currOptions->m_rawChainSize;
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

        // KAUST: all processors should call here
        generateChain(*currOptions,  // input, only m_rawChainSize changes
                      nowNodes,      // input
                      nowCovMatrix,  // input
                      currRv,        // input
                      prevChain,     // input
                      nowChain,      // output 
                      nowRunTime,    // output
                      nowRejections, // output
                      NULL,          // output
                      NULL);         // output

        // KAUST: all processors should call here
        currOptions->m_totallyMute           = savedTotallyMute;
        currOptions->m_rawChainSize          = savedRawChainSize;
        currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
        currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME
        currOptions->m_drMaxNumExtraStages   = savedDrMaxNumExtraStages;
        currOptions->m_amAdaptInterval       = savedAmAdaptInterval;

        if (m_env.inter0Rank() == 0) { // KAUST
        nowRejectionRate = ((double) nowRejections) / ((double) subNumSamples); // FIX ME: rejection among all inter0rank nodes ???
        //bool aux1 = (nowRejectionRate == meanRejectionRate);
        bool aux2 = (nowRejectionRate >= currOptions->m_minRejectionRate)
                    &&
                    (nowRejectionRate <= currOptions->m_maxRejectionRate);
        // KAUST2: all processors should have the same 'testResult'
        testResult = aux2;
        }

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
      } while (testResult == false);
      currEta = nowEta;
      if (currEta != 1.) {
        *unifiedCovMatrix *= currEta;
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << currLevel+LEVEL_REF_ID
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1)
                                << ", currEta = "                              << currEta
                                << ", assessed rejection rate = "              << nowRejectionRate
                                << std::endl;
      }
    } // end of step 9

    //***********************************************************
    // Step 10 of 10: sample vector RV of current level
    //***********************************************************
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 10 of 10"
                                << std::endl;
      }

      // KAUST: all processors should call here
      bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
      unsigned int savedRawChainSize          = currOptions->m_rawChainSize;
      bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
      bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;

      currOptions->m_totallyMute           = true;
      currOptions->m_rawChainSize          = 0; // will be set inside generateChain()
      currOptions->m_rawChainComputeStats  = false;
      currOptions->m_filteredChainGenerate = false;

      // KAUST: all processors should call here
      generateChain(*currOptions,                 // input, only m_rawChainSize changes
                    nodes,                        // input
                    *unifiedCovMatrix,            // input
                    currRv,                       // input
                    prevChain,                    // input
                    currChain,                    // output
                    cumulativeRawChainRunTime,    // output
                    cumulativeRawChainRejections, // output
                    &currLogLikelihoodValues,     // output // likelihood is important
                    &currLogTargetValues);        // output

      // KAUST: all processors should call here
      currOptions->m_totallyMute           = savedTotallyMute;
      currOptions->m_rawChainSize          = savedRawChainSize;
      currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
      currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

      if (m_env.inter0Rank() == 0) { // KAUST
      if (currOptions->m_rawChainComputeStats) {
        std::ofstream* genericOfsVar = NULL;
        m_env.openOutputFile(currOptions->m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                             currOptions->m_dataOutputAllowedSet,
                             false,
                             genericOfsVar);

        currChain.computeStatistics(*currOptions->m_rawChainStatisticalOptions,
                                    genericOfsVar);

        genericOfsVar->close();
      }

      if (currOptions->m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
        currChain.unifiedWriteContents              (currOptions->m_rawChainDataOutputFileName);
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

        genericOfsVar->close();

        if (currOptions->m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
          currChain.unifiedWriteContents              (currOptions->m_filteredChainDataOutputFileName);
          currLogLikelihoodValues.unifiedWriteContents(currOptions->m_filteredChainDataOutputFileName);
          currLogTargetValues.unifiedWriteContents    (currOptions->m_filteredChainDataOutputFileName);
        }
      }
      } // KAUST

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      //UQ_FATAL_TEST_MACRO((currChain.subSequenceSize() != currOptions->m_rawChainSize),
      //                    m_env.fullRank(),
      //                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
      //                    "currChain (a linked one) has been generated with invalid size");
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
                              << "%)" // FIX ME: unified
                              << std::endl;
    }

    // KAUST: all processors should have the same 'currExponent'
    if (currExponent != 1.) delete currOptions;
  } // end of level while

  UQ_FATAL_TEST_MACRO((currExponent < 1),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                      "exponent has not achieved value '1' even after maximum number of levels");

  UQ_FATAL_TEST_MACRO((currLevel != m_logEvidenceFactors.size()),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                      "exponent has not achieved value '1' even after maximum number of levels");

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
  unsigned int               subNumSamples,          // input
  const std::vector<double>& unifiedWeightStdVector, // input
  std::vector<unsigned int>& unifiedIndexCounters)   // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::sampleIndexes()..."
                            << std::endl;
  }

  unifiedIndexCounters.resize(unifiedWeightStdVector.size(),0);

  if (m_env.inter0Rank() == 0) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::sampleIndexes()"
                            //<< ", level " << currLevel+LEVEL_REF_ID
                              << ": unifiedWeightStdVector.size() = " << unifiedWeightStdVector.size()
                              << std::endl;

      //unsigned int numZeros = 0;
      //for (unsigned int i = 0; i < unifiedWeightStdVector.size(); ++i) {
      //  *m_env.subDisplayFile() << "unifiedWeightStdVector[" << i
      //                          << "] = " << unifiedWeightStdVector[i]
      //                          << std::endl;
      //  if (unifiedWeightStdVector[i] == 0.) numZeros++;
      //}
      //*m_env.subDisplayFile() << "Number of zeros in unifiedWeightStdVector = " << numZeros
      //                        << std::endl;
    }

    uqFiniteDistributionClass tmpFd(m_env,
                                    "",
                                    unifiedWeightStdVector);

    // Generate 'unifiedNumSamples' samples from 'tmpFD'
    unsigned int unifiedNumSamples = m_env.inter0Comm().NumProc() * subNumSamples;
    for (unsigned int i = 0; i < unifiedNumSamples; ++i) {
      unsigned int index = tmpFd.sample();
      unifiedIndexCounters[index] += 1;
    }

    int mpiRC = MPI_Bcast((void *) &unifiedIndexCounters[0], (int) unifiedIndexCounters.size(), MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::sampleIndexes()",
                        "failed MPI_Bcast() for unified index counters");
  }

  //for (unsigned int i = 0; i < unifiedIndexCounters.size(); ++i) {
  //  *m_env.subDisplayFile() << "unifiedIndexCounters[" << i
  //                          << "] = " << unifiedIndexCounters[i]
  //                          << std::endl;
  //}

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::sampleIndexes()"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::distribIndexSamples(
  unsigned int                              subNumSamples,        // input
  std::vector<unsigned int>&                unifiedIndexCounters, // input, modified
  std::vector<uqLinkedChainsPerNodeStruct>& nodes)                // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::distribIndexSamples()..."
                            << std::endl;
  }

  unsigned int auxNode = 0;
  unsigned int numberOfPositionsToGuaranteeForNode = subNumSamples;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::distribIndexSamples()"
                          //<< ", level " << currLevel+LEVEL_REF_ID
                            << ": numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                            << std::endl;
  }
  for (unsigned int i = 0; i < unifiedIndexCounters.size(); ++i) {
    while (unifiedIndexCounters[i] != 0) {
      //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 30)) {
      //  *m_env.subDisplayFile() << "auxNode = "                               << auxNode
      //                          << ", numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
      //                          << ", unifiedIndexCounters["                  << i
      //                          << "] = "                                     << unifiedIndexCounters[i]
      //                          << std::endl;
      //}
      UQ_FATAL_TEST_MACRO(auxNode >= (unsigned int) m_env.inter0Comm().NumProc(),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                          "auxNode got too large");
      if (unifiedIndexCounters[i] < numberOfPositionsToGuaranteeForNode) {
        uqLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = unifiedIndexCounters[i];
        nodes[auxNode].linkedChains.push_back(auxControl);

        numberOfPositionsToGuaranteeForNode -= unifiedIndexCounters[i];
        unifiedIndexCounters[i] = 0;
      }
      else {
        uqLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = numberOfPositionsToGuaranteeForNode;
        nodes[auxNode].linkedChains.push_back(auxControl);

        unifiedIndexCounters[i] -= numberOfPositionsToGuaranteeForNode;
        numberOfPositionsToGuaranteeForNode = 0;

        // Go to next node
        auxNode++;
        numberOfPositionsToGuaranteeForNode = subNumSamples;
      }
    }
  }
  UQ_FATAL_TEST_MACRO(auxNode != (unsigned int) m_env.inter0Comm().NumProc(),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                      "auxNode exited loop with wrong value");
  UQ_FATAL_TEST_MACRO(numberOfPositionsToGuaranteeForNode != subNumSamples,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::distribIndexSamples()",
                      "numberOfPositionsToGuaranteeForNode exited loop with wrong value");
  // FIX ME: swap trick to save memory

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::distribIndexSamples()"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateChain(
  uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
  const std::vector<uqLinkedChainsPerNodeStruct>& nodes,                   // input
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
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateChain()..."
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());

  // KAUST2: all processors should have the same 'chainIdMax'
  unsigned int chainIdMax = nodes[m_env.subId()].linkedChains.size();
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) { // FIX ME
      unsigned int auxIndex = nodes[m_env.subId()].linkedChains[chainId].initialPositionIndexInPreviousChain;
      prevChain.getPositionValues(auxIndex,auxInitialPosition); // FIX ME

      tmpChainSize = nodes[m_env.subId()].linkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
    }

    // KAUST2: all processors should have the same 'tmpChainSize'
    inputOptions.m_rawChainSize = tmpChainSize;
    uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    uqScalarSequenceClass<double> tmpLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> tmpLogTargetValues    (m_env,0,"");

    // KAUST: all processors should call here
    uqMetropolisHastingsSGClass<P_V,P_M> mcSeqGenerator(inputOptions,
                                                        rv,
                                                        auxInitialPosition,
                                                        &unifiedCovMatrix);

    // KAUST: all processors should call here
    mcSeqGenerator.generateSequence(tmpChain,
                                    &tmpLogLikelihoodValues, // likelihood is IMPORTANT
                                    &tmpLogTargetValues);
    cumulativeRunTime    += mcSeqGenerator.rawChainRunTime();
    cumulativeRejections += mcSeqGenerator.numRejections();

    if (m_env.inter0Rank() >= 0) { // FIX ME
      //for (unsigned int i = 0; i < tmpLogLikelihoodValues.subSequenceSize(); ++i) {
      //  std::cout << "tmpLogLikelihoodValues[" << i << "] = " << tmpLogLikelihoodValues[i]
      //            << ", tmpLogTargetValues["   << i << "] = " << tmpLogTargetValues[i]
      //            << std::endl;
      //}
        
      if ((m_env.subDisplayFile()             ) &&
          (m_env.displayVerbosity()   >= 0    ) &&
          (inputOptions.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateChain()"
                              //<< ", level "               << currLevel+LEVEL_REF_ID
                                << ", chainId = "           << chainId
                                << ": finished generating " << tmpChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      // FIX ME: unified
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

  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING_H__
