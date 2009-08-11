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
  uqMLSamplingClass(/*! Prefix                 */ const char*                         prefix,                  
                    /*! The source rv          */ const uqBaseVectorRVClass<P_V,P_M>& sourceRv,                
                    /*! Initial chain position */ const P_V&                          initialPosition,
                    /*! Proposal cov. matrix   */ const P_M*                          inputProposalCovMatrix);  
  /*! Destructor: */
 ~uqMLSamplingClass();

  /*! Operation to generate the chain */
 void   generateSequence(uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
                         uqScalarSequenceClass<double>*      workingTargetValues);

  void   print          (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass&       m_env;
  const uqBaseVectorRVClass<P_V,P_M>& m_sourceRv;
        P_V                           m_initialPosition;
  const P_M*                          m_initialProposalCovMatrix;
        bool                          m_nullInputProposalCovMatrix;

        uqMLSamplingOptionsClass      m_options;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj);

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::uqMLSamplingClass(
  const char*                         prefix,
  const uqBaseVectorRVClass<P_V,P_M>& sourceRv,
  const P_V&                          initialPosition,
  const P_M*                          inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_sourceRv                  (sourceRv),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (inputProposalCovMatrix),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_options                   (m_env,prefix)
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
  if (m_nullInputProposalCovMatrix) delete m_initialProposalCovMatrix;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  uqScalarSequenceClass<double>*      workingTargetValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateSequence()..."
                            << std::endl;
  }

  //***********************************************************
  // Declaration of Variables
  //***********************************************************
  double                            currExponent = m_options.m_initialExponent;
  uqSequenceOfVectorsClass<P_V,P_M> currChain(m_sourceRv.imageSet().vectorSpace(),
                                              0,
                                              m_options.m_prefix+"curr_chain");
  uqScalarSequenceClass<double>     currTargetValues(m_env,0);

  //***********************************************************
  // Take care of first level
  //***********************************************************
  unsigned int currLevel = 0;
  {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level " << currLevel+1
                              << ", m_options.m_levelOptions[currLevel]->m_rawChainSize = " << m_options.m_levelOptions[currLevel]->m_rawChainSize
                              << ", currExponent = "  << currExponent
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);

    uqPoweredJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                            m_sourceRv.pdf(),
                                            currExponent);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           m_sourceRv.pdf().domainSet());

    currRv.setPdf(currPdf);

    uqMarkovChainSGClass<P_V,P_M> mcSeqGenerator(*(m_options.m_levelOptions[currLevel]),
                                                 currRv,
                                                 m_initialPosition,
                                                 m_initialProposalCovMatrix);

    mcSeqGenerator.generateSequence(currChain,
                                    &currTargetValues);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ", level " << currLevel+1
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

    //UQ_FATAL_TEST_MACRO((currChain.subSequenceSize() != m_options.m_levelOptions[currLevel]->m_rawChainSize),
    //                    m_env.fullRank(),
    //                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
    //                    "currChain (first one) has been generated with invalid size");
    double levelRunTime = uqMiscGetEllapsedSeconds(&timevalLevel);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << currLevel+1
                              << " after " << levelRunTime << " seconds"
                              << std::endl;
    }
  }

  //***********************************************************
  // Take care of remaining levels
  //***********************************************************
  while ((currLevel    < (m_options.m_maxNumberOfLevels-1)) &&
         (currExponent < 1                                )) {
    currLevel++;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level " << currLevel+1
                              << ", m_options.m_levelOptions[currLevel]->m_rawChainSize = " << m_options.m_levelOptions[currLevel]->m_rawChainSize
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);

    //***********************************************************
    // Step 1 of 8: save [chain and corresponding target pdf values] from previous level
    //***********************************************************
    double prevExponent = currExponent;
    uqSequenceOfVectorsClass<P_V,P_M> prevChain(m_sourceRv.imageSet().vectorSpace(),
                                                0,
                                                m_options.m_prefix+"prev_chain");
    uqScalarSequenceClass<double> prevTargetValues(m_env,0);

    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 1"
                                << std::endl;
      }

      prevChain = currChain;
      currChain.clear();
      currChain.setName(m_options.m_levelOptions[currLevel]->m_prefix + "rawChain");

      prevTargetValues = currTargetValues;
      currTargetValues.clear();

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": prevChain.unifiedSequenceSize() = " << prevChain.unifiedSequenceSize()
                                << ", currChain.unifiedSequenceSize() = " << currChain.unifiedSequenceSize()
                                << ", prevTargetValues.unifiedSequenceSize() = " << prevTargetValues.unifiedSequenceSize()
                                << ", currTargetValues.unifiedSequenceSize() = " << currTargetValues.unifiedSequenceSize()
                                << std::endl;
      }

      UQ_FATAL_TEST_MACRO((prevChain.subSequenceSize() != prevTargetValues.subSequenceSize()),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "different sizes between previous chain and previous sequence of target values");
    } // end of step 1

    //***********************************************************
    // Step 2 of 8: create [currExponent and sequence of weights] for current level
    //***********************************************************
    uqScalarSequenceClass<double> weightSequence(m_env,prevTargetValues.subSequenceSize());
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 2"
                                << std::endl;
      }

      double exponentQuanta = std::min(1.,m_options.m_levelOptions[currLevel]->m_maxExponent) - prevExponent;
      exponentQuanta /= (double) m_options.m_levelOptions[currLevel]->m_maxNumberOfAttempts;

      unsigned int currAttempt = 0;
      bool testResult = false;
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+1
                                  << ": entering loop for computing next exponent, with currAttempt = " << currAttempt
                                  << std::endl;
        }
        currExponent = m_options.m_levelOptions[currLevel]->m_maxExponent - currAttempt*exponentQuanta;
        double auxExp = (currExponent/prevExponent) - 1.;
        double weightSum = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] = pow(prevTargetValues[i],auxExp);
          weightSum += weightSequence[i];
        }
        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] /= weightSum;
          effectiveSampleSize += 1/weightSequence[i]/weightSequence[i];
        }
        double auxRatio = effectiveSampleSize/((double) weightSequence.subSequenceSize());
        testResult = (auxRatio >= m_options.m_levelOptions[currLevel]->m_minEffectiveSizeRatio);
        currAttempt++;
      } while ((currAttempt < m_options.m_levelOptions[currLevel]->m_maxNumberOfAttempts) &&
               (testResult == false));

      UQ_FATAL_TEST_MACRO((testResult == false),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "test for next exponent failed even after maximum number of attempts");

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << weightSequence.unifiedSequenceSize()
                                << ", currExponent = "                         << currExponent
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
    } // end of step 2

    //***********************************************************
    // Step 3 of 8: create covariance matrix for current level
    //***********************************************************
    P_M* unifiedCovMatrix = m_sourceRv.imageSet().vectorSpace().newMatrix();
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 3"
                                << std::endl;
      }

      P_V auxVec(m_sourceRv.imageSet().vectorSpace().zeroVector());
      P_V weightedMeanVec(m_sourceRv.imageSet().vectorSpace().zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        weightedMeanVec += weightSequence[i]*auxVec;
      }

      P_V diffVec(m_sourceRv.imageSet().vectorSpace().zeroVector());
      P_M* subCovMatrix = m_sourceRv.imageSet().vectorSpace().newMatrix();
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
    } // end of step 3

    //***********************************************************
    // Step 4 of 8: create *unified* finite distribution for current level
    //***********************************************************
    std::vector<unsigned int> unifiedIndexCounters(weightSequence.unifiedSequenceSize(),0);
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 4"
                                << std::endl;
      }

      std::vector<double> unifiedWeightStdVector(0);
      weightSequence.getUnifiedContentsAtProc0Only(unifiedWeightStdVector);

      if (m_env.inter0Rank() == 0) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << currLevel+1
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

        unsigned int subChainSize = m_options.m_levelOptions[currLevel]->m_rawChainSize;
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
    } // end of step 4

    //***********************************************************
    // Step 5 of 8: plan for number of linked chains for each node so that
    //              each node generates the same number of positions
    //***********************************************************
    std::vector<uqLinkedChainsPerNodeStruct> nodes(m_env.inter0Comm().NumProc());
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 5"
                                << std::endl;
      }

      unsigned int auxNode = 0;
      unsigned int numberOfPositionsToGuaranteeForNode = m_options.m_levelOptions[currLevel]->m_rawChainSize;
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
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
            numberOfPositionsToGuaranteeForNode = m_options.m_levelOptions[currLevel]->m_rawChainSize;
          }
        }
      }
      UQ_FATAL_TEST_MACRO(auxNode != (unsigned int) m_env.inter0Comm().NumProc(),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "auxNode exited loop with wrong value");
      UQ_FATAL_TEST_MACRO(numberOfPositionsToGuaranteeForNode != m_options.m_levelOptions[currLevel]->m_rawChainSize,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "numberOfPositionsToGuaranteeForNode exited loop with wrong value");

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": nodes[m_env.subId()].linkedChains.size() = " << nodes[m_env.subId()].linkedChains.size()
                                << std::endl;
      }

      // FIX ME: swap trick to save memory
    } // end of step 5

    //***********************************************************
    // Step 6 of 8: load balancing = subChainSize [indexes + corresponding positions] for each node
    //***********************************************************
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 6"
                                << std::endl;
      }

      if (m_env.inter0Comm().NumProc() > 1) {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "incomplete code for load balancing");
      }
    } // end of step 6
  
    //***********************************************************
    // Step 7 of 8: create vector RV for current level
    //***********************************************************
    uqPoweredJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                            m_sourceRv.pdf(),
                                            currExponent);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           m_sourceRv.pdf().domainSet());

    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 7"
                                << std::endl;
      }

      currRv.setPdf(currPdf);
    } // end of step 7

    //***********************************************************
    // Step 8 of 8: sample vector RV of current level
    //***********************************************************
    {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": beginning step 8"
                                << std::endl;
      }

      P_V auxInitialPosition(m_sourceRv.imageSet().vectorSpace().zeroVector());
      unsigned int savedRawChainSize          = m_options.m_levelOptions[currLevel]->m_rawChainSize;
      bool         savedRawChainComputeStats  = m_options.m_levelOptions[currLevel]->m_rawChainComputeStats;
      bool         savedFilteredChainGenerate = m_options.m_levelOptions[currLevel]->m_filteredChainGenerate;
      if (m_env.inter0Rank() >= 0) { // FIX ME
        for (unsigned int chainId = 0; chainId < nodes[m_env.subId()].linkedChains.size(); ++chainId) {
          unsigned int auxIndex = nodes[m_env.subId()].linkedChains[chainId].initialPositionIndexInPreviousChain;
          prevChain.getPositionValues(auxIndex,auxInitialPosition); // FIX ME

          unsigned int auxNumPositions = nodes[m_env.subId()].linkedChains[chainId].numberOfPositions;
          m_options.m_levelOptions[currLevel]->m_rawChainSize          = auxNumPositions;
          m_options.m_levelOptions[currLevel]->m_rawChainComputeStats  = false;
          m_options.m_levelOptions[currLevel]->m_filteredChainGenerate = false;

          uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_sourceRv.imageSet().vectorSpace(),
                                                     0,
                                                     m_options.m_prefix+"tmp_chain");
          uqScalarSequenceClass<double> tmpTargetValues(m_env,0);

          uqMarkovChainSGClass<P_V,P_M> mcSeqGenerator(*(m_options.m_levelOptions[currLevel]),
                                                       currRv,
                                                       auxInitialPosition,
                                                       unifiedCovMatrix);

          mcSeqGenerator.generateSequence(tmpChain,
                                          &tmpTargetValues);
        
          if ((m_env.subDisplayFile()                                     ) &&
              (m_env.displayVerbosity() >= 0                              ) &&
              (m_options.m_levelOptions[currLevel]->m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level "               << currLevel+1
                                    << ", chainId = "           << chainId
                                    << ": finished generating " << tmpChain.subSequenceSize()
                                    << " chain positions"
                                    << std::endl;
          }

          currChain.append(tmpChain);
          currTargetValues.append(tmpTargetValues);
        }
      }
      m_options.m_levelOptions[currLevel]->m_rawChainSize          = savedRawChainSize;
      m_options.m_levelOptions[currLevel]->m_rawChainComputeStats  = savedRawChainComputeStats;
      m_options.m_levelOptions[currLevel]->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

      if (m_options.m_levelOptions[currLevel]->m_rawChainComputeStats) {
        std::ofstream* genericOfsVar = NULL;
        m_env.openOutputFile(m_options.m_levelOptions[currLevel]->m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                             m_options.m_levelOptions[currLevel]->m_dataOutputAllowedSet,
                             false,
                             genericOfsVar);

        currChain.computeStatistics(*m_options.m_levelOptions[currLevel]->m_rawChainStatisticalOptions,
                                    genericOfsVar);

        genericOfsVar->close();
      }

      if (m_options.m_levelOptions[currLevel]->m_filteredChainGenerate) {
        std::ofstream* genericOfsVar = NULL;
        m_env.openOutputFile(m_options.m_levelOptions[currLevel]->m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                             m_options.m_levelOptions[currLevel]->m_dataOutputAllowedSet,
                             false,
                             genericOfsVar);

        unsigned int filterInitialPos = (unsigned int) (m_options.m_levelOptions[currLevel]->m_filteredChainDiscardedPortion * (double) currChain.subSequenceSize());
        unsigned int filterSpacing    = m_options.m_levelOptions[currLevel]->m_filteredChainLag;
        if (filterSpacing == 0) {
          currChain.computeFilterParams(*m_options.m_levelOptions[currLevel]->m_filteredChainStatisticalOptions,
                                        genericOfsVar,
                                        filterInitialPos,
                                        filterSpacing);
        }

        // Filter positions from the converged portion of the chain
        currChain.filter(filterInitialPos,
                         filterSpacing);
        currChain.setName(m_options.m_levelOptions[currLevel]->m_prefix + "filtChain");

        currTargetValues.filter(filterInitialPos,
                                filterSpacing);

        if (m_options.m_levelOptions[currLevel]->m_filteredChainComputeStats) {
          currChain.computeStatistics(*m_options.m_levelOptions[currLevel]->m_filteredChainStatisticalOptions,
                                      genericOfsVar);
        }

        genericOfsVar->close();
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << currLevel+1
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      //UQ_FATAL_TEST_MACRO((currChain.subSequenceSize() != m_options.m_levelOptions[currLevel]->m_rawChainSize),
      //                    m_env.fullRank(),
      //                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
      //                    "currChain (a linked one) has been generated with invalid size");
    } // end of step 8

    //***********************************************************
    // Prepare to end current level
    //***********************************************************
    delete unifiedCovMatrix;

    double levelRunTime = uqMiscGetEllapsedSeconds(&timevalLevel);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << currLevel+1
                              << " after " << levelRunTime << " seconds"
                              << std::endl;
    }
  } // end of level loop

  UQ_FATAL_TEST_MACRO((currExponent < 1),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                      "exponent has not achieved value '1' even after maximum number of leves");

  //***********************************************************
  // Prepare to return
  //***********************************************************
  workingChain.clear();
  workingChain.resizeSequence(currChain.subSequenceSize());
  P_V auxVec(m_sourceRv.imageSet().vectorSpace().zeroVector());
  for (unsigned int i = 0; i < workingChain.subSequenceSize(); ++i) {
    currChain.getPositionValues(i,auxVec);
    workingChain.setPositionValues(i,auxVec);
  }

  if (workingTargetValues) *workingTargetValues = currTargetValues;

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
