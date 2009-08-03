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
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
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
                              << std::endl;
    }

    uqPoweredJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                            m_sourceRv.pdf(),
                                            currExponent);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           m_sourceRv.pdf().domainSet());

    currRv.setPdf(currPdf);

    uqMarkovChainSGClass<P_V,P_M> mcSeqGenerator(m_options.m_prefix.c_str(),
                                                 currRv,
                                                 m_initialPosition,
                                                 m_initialProposalCovMatrix);

    mcSeqGenerator.generateSequence(currChain,
                                    &currTargetValues);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << currLevel+1
                              << ", m_options.m_levelOptions[currLevel]->m_rawChainSize = " << m_options.m_levelOptions[currLevel]->m_rawChainSize
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

    // Step 1 of 8: save [chain and corresponding target pdf values] from previous level
    double prevExponent = currExponent;

    uqSequenceOfVectorsClass<P_V,P_M> prevChain(m_sourceRv.imageSet().vectorSpace(),
                                                0,
                                                m_options.m_prefix+"prev_chain");
    prevChain = currChain;
    currChain.clear();

    uqScalarSequenceClass<double> prevTargetValues(m_env,0);
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

    // Step 2 of 8: create [currExponent and sequence of weights] for current level
    uqScalarSequenceClass<double> weightSequence(m_env,prevTargetValues.subSequenceSize());
    double exponentQuanta = std::min(1.,m_options.m_levelOptions[currLevel]->m_maxExponent) - prevExponent;
    exponentQuanta /= (double) m_options.m_levelOptions[currLevel]->m_maxNumberOfAttempts;

    unsigned int currAttempt = 0;
    bool testResult = false;
    do {
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
                              << std::endl;
    }

    // Step 3 of 8: create covariance matrix for current level
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

    P_M* unifiedCovMatrix = m_sourceRv.imageSet().vectorSpace().newMatrix();
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

    // Step 4 of 8: create *unified* finite distribution for current level
    std::vector<double> unifiedSampleStdVector(weightSequence.unifiedSequenceSize(),0.);
    for (unsigned int i = 0; i < weightSequence.unifiedSequenceSize(); ++i) {
      unifiedSampleStdVector[i] = i;
    }

    std::vector<unsigned int> unifiedIndexCounters(weightSequence.unifiedSequenceSize(),0);
    std::vector<double> unifiedWeightStdVector(0);
    weightSequence.getUnifiedContentsAtProc0Only(unifiedWeightStdVector);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ", level " << currLevel+1
                              << ": unifiedWeightStdVector.size() = " << unifiedWeightStdVector.size()
                              << std::endl;
    }

    if (m_env.inter0Rank() == 0) {
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

    // Step 5 of 8: plan for number of linked chains for each node so that
    //              each node generates the same number of positions
    std::vector<uqLinkedChainsPerNodeStruct> nodes(m_env.inter0Comm().NumProc());
    unsigned int auxNode = 0;
    unsigned int numberOfPositionsToGuaranteeForNode = m_options.m_levelOptions[currLevel]->m_rawChainSize;
    for (unsigned int i = 0; i < unifiedIndexCounters.size(); ++i) {
      while (unifiedIndexCounters[i] != 0) {
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

          numberOfPositionsToGuaranteeForNode = 0;
          unifiedIndexCounters[i] -= numberOfPositionsToGuaranteeForNode;

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

    // FIX ME: swap trick to save memory

    // Step 6 of 8: load balancing = subChainSize [indexes + corresponding positions] for each node
    if (m_env.inter0Comm().NumProc() > 1) {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "incomplete code for load balancing");
    }
  
    // Step 7 of 8: create vector RV for current level
    uqPoweredJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                            m_sourceRv.pdf(),
                                            currExponent);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           m_sourceRv.pdf().domainSet());

    currRv.setPdf(currPdf);

    // Step 8 of 8: sample vector RV of current level
    P_V auxInitialPosition(m_sourceRv.imageSet().vectorSpace().zeroVector());
    if (m_env.inter0Rank() > 0) {
      unsigned int myInter0Rank = (unsigned int) m_env.inter0Rank();
      for (unsigned int chainId = 0; chainId < nodes[myInter0Rank].linkedChains.size(); ++chainId) {
        unsigned int auxIndex = nodes[myInter0Rank].linkedChains[chainId].initialPositionIndexInPreviousChain;
        prevChain.getPositionValues(auxIndex,auxInitialPosition); // FIX ME

        //unsigned int auxNumPositions = nodes[myInter0Rank].linkedChains[chainId].numberOfPositions;

        uqMarkovChainSGClass<P_V,P_M> mcSeqGenerator(m_options.m_prefix.c_str(),
                                                     currRv,
                                                     auxInitialPosition,
                                                     unifiedCovMatrix);

        mcSeqGenerator.generateSequence(currChain, // FIX ME: linked chains
                                        &currTargetValues);
      }
    }

    delete unifiedCovMatrix;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << currLevel+1
                              << ", m_options.m_levelOptions[currLevel]->m_rawChainSize = " << m_options.m_levelOptions[currLevel]->m_rawChainSize
                              << std::endl;
    }
  }

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

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
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
