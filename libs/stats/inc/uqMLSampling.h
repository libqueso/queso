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
  void   generateSequence(uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
                          uqScalarSequenceClass<double>*      workingTargetValues,
                          uqScalarSequenceClass<double>*      workingLogTargetValues);

  void   print           (std::ostream& os) const;

private:
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
  uqScalarSequenceClass<double>*      workingTargetValues,
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
  uqSequenceOfVectorsClass<P_V,P_M> currChain(m_vectorSpace,
                                              0,
                                              m_options.m_prefix+"curr_chain");
  uqScalarSequenceClass<double>     currLogTargetValues(m_env,0);

  //***********************************************************
  // Take care of first level
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
    currChain.resizeSequence          (currOptions.m_rawChainSize);
    currLogTargetValues.resizeSequence(currOptions.m_rawChainSize);
    P_V auxVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 0; i < currChain.subSequenceSize(); ++i) {
      m_priorRv.realizer().realization(auxVec);
      currChain.setPositionValues(i,auxVec);
      currLogTargetValues[i] = -.5*m_likelihoodFunction.minus2LnValue(auxVec,NULL,NULL,NULL,NULL);
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
    double       cumulativeRawChainRunTime = 0.;
    unsigned int cumulativeNumRejections   = 0;

    //***********************************************************
    // Step 1 of 10: read options
    //***********************************************************
    sprintf(tmpSufix,"%d_",currLevel+LEVEL_REF_ID); // Yes, '+0'
    uqMLSamplingLevelOptionsClass* currOptions = new uqMLSamplingLevelOptionsClass(m_env,(m_options.m_prefix + tmpSufix).c_str());
    {
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
    uqSequenceOfVectorsClass<P_V,P_M> prevChain(m_vectorSpace,
                                                0,
                                                m_options.m_prefix+"prev_chain");
    uqScalarSequenceClass<double> prevLogTargetValues(m_env,0);

    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 2 of 10"
                                << std::endl;
      }

      prevChain = currChain;
      currChain.clear();
      currChain.setName(currOptions->m_prefix + "rawChain");

      prevLogTargetValues = currLogTargetValues;
      currLogTargetValues.clear();

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": prevChain.unifiedSequenceSize() = " << prevChain.unifiedSequenceSize()
                                << ", currChain.unifiedSequenceSize() = " << currChain.unifiedSequenceSize()
                                << ", prevLogTargetValues.unifiedSequenceSize() = " << prevLogTargetValues.unifiedSequenceSize()
                                << ", currLogTargetValues.unifiedSequenceSize() = " << currLogTargetValues.unifiedSequenceSize()
                                << std::endl;
      }

      UQ_FATAL_TEST_MACRO((prevChain.subSequenceSize() != prevLogTargetValues.subSequenceSize()),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "different sizes between previous chain and previous sequence of target values");
    } // end of step 2

    //***********************************************************
    // Step 3 of 10: compute [currExponent and sequence of weights] for current level
    //***********************************************************
    uqScalarSequenceClass<double> weightSequence(m_env,prevLogTargetValues.subSequenceSize());
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 3 of 10"
                                << std::endl;
      }

      std::vector<double> exponents(2,0.);
      exponents[0] = prevExponent;
      exponents[1] = 1.;

      currExponent = 1.; // Try '1.' right away
      double currEffectiveSizeRatio  = 0.; // To be computed

      unsigned int currAttempt = 0;
      bool testResult = false;
      double tmpEvidenceFactor = 0.;
      double meanEffectiveSizeRatio = .5*(currOptions->m_minEffectiveSizeRatio + currOptions->m_maxEffectiveSizeRatio);
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+LEVEL_REF_ID
                                  << ": entering loop for computing next exponent, with currAttempt = " << currAttempt
                                  << std::endl;
        }
        if (currAttempt > 0) {
          if (currEffectiveSizeRatio > meanEffectiveSizeRatio) {
            exponents[0] = currExponent;
          }
          else {
            exponents[1] = currExponent;
          }
          currExponent = .5*(exponents[0] + exponents[1]);
        }
        double auxExp = currExponent;
        if (prevExponent != 0.) {
          auxExp /= prevExponent;
          auxExp -= 1.;
        }
        double weightSum = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] = exp(prevLogTargetValues[i]*auxExp);
          weightSum += weightSequence[i];
          //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          //  *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
          //                          << ", level "              << currLevel+LEVEL_REF_ID
          //                          << ": auxExp = "           << auxExp
          //                          << ", 'prev'TargetValues[" << i
          //                          << "] = "                  << prevTargetValues[i]
          //                          << ", weightSequence["     << i
          //                          << "] = "                  << weightSequence[i]
          //                          << ", weightSum = "        << weightSum
          //                          << std::endl;
          //}
        }
        tmpEvidenceFactor = weightSum/prevChain.unifiedSequenceSize(); // FIX ME: unified
        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] /= weightSum;
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
        currEffectiveSizeRatio = effectiveSampleSize/((double) weightSequence.subSequenceSize()); // FIX ME: unified
        UQ_FATAL_TEST_MACRO((currEffectiveSizeRatio > (1.+1.e-8)),
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "effective sample size ratio cannot be > 1");

        bool aux1 = (currEffectiveSizeRatio == meanEffectiveSizeRatio);
        bool aux2 = (currExponent == 1.                             )
                    &&
                    (currEffectiveSizeRatio > meanEffectiveSizeRatio);
        bool aux3 = (currEffectiveSizeRatio >= currOptions->m_minEffectiveSizeRatio)
                    &&
                    (currEffectiveSizeRatio <= currOptions->m_maxEffectiveSizeRatio);
        testResult = aux1 || aux2 || aux3;

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "                    << currLevel+LEVEL_REF_ID
                                  << ": currAttempt = "            << currAttempt
                                  << ", previousExponent = "       << prevExponent
                                  << ", exponents[0] = "           << exponents[0]
                                  << ", attemptedExponent = "      << currExponent
                                  << ", exponents[1] = "           << exponents[1]
                                  << ", effectiveSampleSize = "    << effectiveSampleSize
                                  << ", weightSequenceSize = "     << weightSequence.subSequenceSize()
                                  << ", currEffectiveSizeRatio = " << currEffectiveSizeRatio
                                  << ", minEffectiveSizeRatio = "  << currOptions->m_minEffectiveSizeRatio
                                  << ", maxEffectiveSizeRatio = "  << currOptions->m_maxEffectiveSizeRatio
                                  << std::endl;
        }
        currAttempt++;
      } while (testResult == false);
      m_logEvidenceFactors.push_back(log(tmpEvidenceFactor));

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << currLevel+LEVEL_REF_ID
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << weightSequence.unifiedSequenceSize()
                                << ", currExponent = "                         << currExponent
                                << ", effective ratio = "                      << currEffectiveSizeRatio
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
    {
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
          int mpiRC = MPI_Allreduce((void *) &localValue, (void *) &sumValue, (int) 1, MPI_DOUBLE, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                              "failed MPI_Allreduce() for cov matrix");
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
    std::vector<unsigned int> unifiedIndexCounters(weightSequence.unifiedSequenceSize(),0);
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 5 of 10"
                                << std::endl;
      }

      std::vector<double> unifiedWeightStdVector(0);
      weightSequence.getUnifiedContentsAtProc0Only(unifiedWeightStdVector);

      if (m_env.inter0Rank() == 0) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+LEVEL_REF_ID
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

        unsigned int subChainSize = currOptions->m_rawChainSize;
        unsigned int unifiedChainSize = m_env.inter0Comm().NumProc() * subChainSize;

        // Generate 'unifiedChainSize' samples from 'tmpFD'
        for (unsigned int i = 0; i < unifiedChainSize; ++i) {
          unsigned int index = tmpFd.sample();
          unifiedIndexCounters[index] += 1;
        }
      }

      int mpiRC = MPI_Bcast((void *) &unifiedIndexCounters[0], (int) unifiedIndexCounters.size(), MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "failed MPI_Bcast() for unified index counters");

      //for (unsigned int i = 0; i < unifiedIndexCounters.size(); ++i) {
      //  *m_env.subDisplayFile() << "unifiedIndexCounters[" << i
      //                          << "] = " << unifiedIndexCounters[i]
      //                          << std::endl;
      //}
    } // end of step 5

    //***********************************************************
    // Step 6 of 10: plan for number of linked chains for each node so that
    //              each node generates the same number of positions
    //***********************************************************
    std::vector<uqLinkedChainsPerNodeStruct> nodes(m_env.inter0Comm().NumProc());
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 6 of 10"
                                << std::endl;
      }

      unsigned int auxNode = 0;
      unsigned int numberOfPositionsToGuaranteeForNode = currOptions->m_rawChainSize;
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
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
                              "uqMLSamplingClass<P_V,P_M>::generateSequence()",
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
            numberOfPositionsToGuaranteeForNode = currOptions->m_rawChainSize;
          }
        }
      }
      UQ_FATAL_TEST_MACRO(auxNode != (unsigned int) m_env.inter0Comm().NumProc(),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "auxNode exited loop with wrong value");
      UQ_FATAL_TEST_MACRO(numberOfPositionsToGuaranteeForNode != currOptions->m_rawChainSize,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "numberOfPositionsToGuaranteeForNode exited loop with wrong value");

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": nodes[m_env.subId()].linkedChains.size() = " << nodes[m_env.subId()].linkedChains.size()
                                << std::endl;
      }

      // FIX ME: swap trick to save memory
    } // end of step 6

    //***********************************************************
    // Step 7 of 10: load balancing = subChainSize [indexes + corresponding positions] for each node
    //***********************************************************
    {
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
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+LEVEL_REF_ID
                                << ": beginning step 10 of 10"
                                << std::endl;
      }

      double prevEta = 1.;
      std::vector<double> eta(2,0.);
      eta[0] = prevEta;
      eta[1] = 1.;

      double currEta = 1.; // Try '1.' right away
      double currRejectionRate  = 0.; // To be computed

      unsigned int currAttempt = 0;
      bool testResult = false;
      double tmpEvidenceFactor = 0.;
      double meanRejectionRate = .5*(currOptions->m_minRejectionRate + currOptions->m_maxRejectionRate);
      do {
        tmpEvidenceFactor = currRejectionRate - meanRejectionRate;
        currAttempt++;
        testResult = true;
      } while (testResult == false);

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << currLevel+LEVEL_REF_ID
                                << ": eta = " << 0.
                                << std::endl;
      }
      // HERE - ENHANCEMENT
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

      P_V auxInitialPosition(m_vectorSpace.zeroVector());
      bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
      unsigned int savedRawChainSize          = currOptions->m_rawChainSize;
      bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
      bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;
      if (m_env.inter0Rank() >= 0) { // FIX ME
        for (unsigned int chainId = 0; chainId < nodes[m_env.subId()].linkedChains.size(); ++chainId) {
          unsigned int auxIndex = nodes[m_env.subId()].linkedChains[chainId].initialPositionIndexInPreviousChain;
          prevChain.getPositionValues(auxIndex,auxInitialPosition); // FIX ME

          unsigned int auxNumPositions = nodes[m_env.subId()].linkedChains[chainId].numberOfPositions;
          currOptions->m_totallyMute           = true;
          currOptions->m_rawChainSize          = auxNumPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
          currOptions->m_rawChainComputeStats  = false;
          currOptions->m_filteredChainGenerate = false;

          uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_vectorSpace,
                                                     0,
                                                     m_options.m_prefix+"tmp_chain");
          uqScalarSequenceClass<double> tmpLogTargetValues(m_env,0);

          uqMetropolisHastingsSGClass<P_V,P_M> mcSeqGenerator(*currOptions,
                                                              currRv,
                                                              auxInitialPosition,
                                                              unifiedCovMatrix);

          mcSeqGenerator.generateSequence(tmpChain,
                                          NULL,
                                          &tmpLogTargetValues);
          cumulativeRawChainRunTime += mcSeqGenerator.rawChainRunTime();
          cumulativeNumRejections   += mcSeqGenerator.numRejections();
        
          if ((m_env.subDisplayFile()            ) &&
              (m_env.displayVerbosity() >= 0     ) &&
              (currOptions->m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level "               << currLevel+LEVEL_REF_ID
                                    << ", chainId = "           << chainId
                                    << ": finished generating " << tmpChain.subSequenceSize()
                                    << " chain positions"
                                    << std::endl;
          }

          // FIX ME: unified
          currChain.append          (tmpChain,          1,tmpChain.subSequenceSize()-1          ); // IMPORTANT: '1' in order to discard initial position
          currLogTargetValues.append(tmpLogTargetValues,1,tmpLogTargetValues.subSequenceSize()-1); // IMPORTANT: '1' in order to discard initial position
        }
      }
      currOptions->m_totallyMute           = savedTotallyMute;
      currOptions->m_rawChainSize          = savedRawChainSize;
      currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
      currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

      // HERE - ENHANCEMENT: write chain contents to the disk

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

      //std::cout << "In loop, pos 1" << std::endl;
      //sleep(1);

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

        currLogTargetValues.filter(filterInitialPos,
                                   filterSpacing);

        if (currOptions->m_filteredChainComputeStats) {
          currChain.computeStatistics(*currOptions->m_filteredChainStatisticalOptions,
                                      genericOfsVar);
        }

        genericOfsVar->close();
      }

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
                              << ": ending level " << currLevel+LEVEL_REF_ID
                              << " after " << levelRunTime << " seconds"
                              << ", cumulativeRawChainRunTime = " << cumulativeRawChainRunTime << " seconds"
                              << ", cumulativeNumRejections = "   << cumulativeNumRejections
                              << " (" << 100.*((double) cumulativeNumRejections)/((double) currOptions->m_rawChainSize)
                              << "%)" // FIX ME: unified
                              << std::endl;
    }
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

  if (workingTargetValues) {
    workingTargetValues->clear();
    workingTargetValues->resizeSequence(currLogTargetValues.subSequenceSize());
    for (unsigned int i = 0; i < workingTargetValues->subSequenceSize(); ++i) {
      (*workingTargetValues)[i] = exp(currLogTargetValues[i]);
    }
  }
  if (workingLogTargetValues) *workingLogTargetValues = currLogTargetValues;

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
