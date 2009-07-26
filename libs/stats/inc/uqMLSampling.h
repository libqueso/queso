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
  uqBaseJointPdfClass      <P_V,P_M>* currPdf          = NULL;
  uqGenericVectorRVClass   <P_V,P_M>* currRv           = NULL;
  uqMarkovChainSGClass     <P_V,P_M>* mcSeqGenerator   = NULL;
  P_M*                                unifiedCovMatrix = NULL;

  double                              prevExponent = -1.;
  double                              currExponent = m_options.m_initialExponent;

  uqSequenceOfVectorsClass<P_V,P_M>  prevChain(m_sourceRv.imageSet().vectorSpace(),
                                               0,
                                               m_options.m_prefix+"prev_chain");
  uqSequenceOfVectorsClass<P_V,P_M>  currChain(m_sourceRv.imageSet().vectorSpace(),
                                               0,
                                               m_options.m_prefix+"curr_chain");

  uqScalarSequenceClass    <double>   prevTargetValues(m_env,0);
  uqScalarSequenceClass    <double>   currTargetValues(m_env,0);

  uqScalarSequenceClass    <double>   weightSequence  (m_env,0);

  //***********************************************************
  // Actual loop
  //***********************************************************
  unsigned int currLevel = 0;
  do {
    // Step 1 of 7: save [chain and corresponding target pdf values] from previous level
    if (currLevel > 0) {
      prevExponent = currExponent;

      prevChain.clear();
      //prevChain = currChain; FIX ME
      currChain.clear();

      prevTargetValues.clear();
      //prevTargetValues = currTargetValues; FIX ME
      currTargetValues.clear();

      weightSequence.clear();
      weightSequence.resizeSequence(prevTargetValues.subSequenceSize());

      if (currLevel > 1) delete unifiedCovMatrix;
      delete mcSeqGenerator;
      delete currRv;
      delete currPdf;
    }

    // Step 2 of 7: loop until [currExponent and vector of weights] are set for current level
    if (currLevel > 0) {
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
        double auxRatio = effectiveSampleSize/(double) weightSequence.subSequenceSize();
        testResult = (auxRatio >= m_options.m_levelOptions[currLevel]->m_minEffectiveSizeRatio);
        currAttempt++;
      } while ((currAttempt < m_options.m_levelOptions[currLevel]->m_maxNumberOfAttempts) &&
               (testResult == false));

      UQ_FATAL_TEST_MACRO((testResult == false),
                          m_env.fullRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "test for next exponent failed even after maximum number of attempts");
    }

    // Step 3 of 7: create vector RV for current level
    currPdf = new uqPoweredJointPdfClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                  m_sourceRv.pdf(),
                                                  currExponent);

    currRv = new uqGenericVectorRVClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                 m_sourceRv.pdf().domainSet());
    currRv->setPdf(*currPdf);

    // Step 4 of 7: create covariance matrix for current level
    if (currLevel > 0) {
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

      unifiedCovMatrix = m_sourceRv.imageSet().vectorSpace().newMatrix();
      *unifiedCovMatrix = *subCovMatrix; // FIX ME
      delete subCovMatrix;
    }

    // Step 5 of 7: create discrete RV for current level
    if (currLevel > 0) {
    }

    // Step 6 of 7: sample vector RV of current level
    if (currLevel == 0) {
      mcSeqGenerator = new uqMarkovChainSGClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                         *currRv,
                                                         m_initialPosition,
                                                         m_initialProposalCovMatrix);
    }
    else {
      mcSeqGenerator = new uqMarkovChainSGClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                         *currRv,
                                                         m_initialPosition,
                                                         unifiedCovMatrix);
    }

    mcSeqGenerator->generateSequence(currChain,
                                     &currTargetValues);

    // Step 7 of 7: prepare for next cycle in the loop
    currLevel++;
  } while ((currLevel    < m_options.m_maxNumberOfLevels) &&
           (currExponent < 1                            ));

  //***********************************************************
  // Prepare to return
  //***********************************************************
  //workingChain                                  = currChain; FIX ME
  //if (workingTargetValues) *workingTargetValues = currTargetValues; FIX ME

  prevChain.clear();
  currChain.clear();
  prevTargetValues.clear();
  currTargetValues.clear();
  weightSequence.clear();
  delete unifiedCovMatrix;
  delete mcSeqGenerator;
  delete currRv;
  delete currPdf;

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
