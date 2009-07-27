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
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>

/*! A templated class that represents a Markov Chain generator. 'SG' stands for 'Sequence Generator'.
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
  {
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
  }

  //***********************************************************
  // Take care of remaining levels
  //***********************************************************
  unsigned int currLevel = 0;
  while ((currLevel    < (m_options.m_maxNumberOfLevels-1)) &&
         (currExponent < 1                                )) {
    currLevel++;

    // Step 1 of 6: save [chain and corresponding target pdf values] from previous level
    double prevExponent = currExponent;

    uqSequenceOfVectorsClass<P_V,P_M> prevChain(m_sourceRv.imageSet().vectorSpace(),
                                                0,
                                                m_options.m_prefix+"prev_chain");
    //prevChain = currChain; FIX ME

    uqScalarSequenceClass<double> prevTargetValues(m_env,0);
    //prevTargetValues = currTargetValues; FIX ME

    // Step 2 of 6: create [currExponent and sequence of weights] for current level
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

    // Step 3 of 6: create covariance matrix for current level
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
    *unifiedCovMatrix = *subCovMatrix; // FIX ME
    delete subCovMatrix;

    // Step 4 of 6: create *unified* discrete RV for current level
    std::vector<double> unifiedWeightStdVector(weightSequence.unifiedSequenceSize(),0.);
    // FIX ME: do MPI stuff...
    // FIX ME: uqDiscreteRVClass<P_V,P_M> discreteRv(m_env,unifiedWeightStdVector);

    std::vector<unsigned int> unifiedIndexCounters(weightSequence.unifiedSequenceSize(),0);
    for (unsigned int i = 0; i < unifiedIndexCounters.size(); ++i) {
      // FIX ME: unsigned int index = discreteRv.sample();
      // FIX ME: unifiedIndexCounters[index] += 1;
    }

    // FIX ME: do MPI stuff...

    // FIX ME: load balancing
  
    // Step 5 of 6: create vector RV for current level
    uqPoweredJointPdfClass<P_V,P_M> currPdf(m_options.m_prefix.c_str(),
                                            m_sourceRv.pdf(),
                                            currExponent);

    uqGenericVectorRVClass<P_V,P_M> currRv(m_options.m_prefix.c_str(),
                                           m_sourceRv.pdf().domainSet());

    currRv.setPdf(currPdf);

    // Step 6 of 6: sample vector RV of current level
    uqMarkovChainSGClass<P_V,P_M> mcSeqGenerator(m_options.m_prefix.c_str(),
                                                 currRv,
                                                 m_initialPosition, // FIX ME
                                                 unifiedCovMatrix);

    // FIX ME: linked chains

    mcSeqGenerator.generateSequence(currChain,
                                    &currTargetValues);

    delete unifiedCovMatrix;
  }

  UQ_FATAL_TEST_MACRO((currExponent < 1),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                      "exponent has not achieved value '1' even after maximum number of leves");

  //***********************************************************
  // Prepare to return
  //***********************************************************
  //workingChain                                  = currChain; FIX ME
  //if (workingTargetValues) *workingTargetValues = currTargetValues; FIX ME

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
