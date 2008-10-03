/* uq/libs/mcmc/inc/uqMarkovChainSG2.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_MAC_SG2_H__
#define __UQ_MAC_SG2_H__

template <class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences(
  const P_M&                          proposalCovMatrix,
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences()..."
              << std::endl;
  }

  P_V valuesOf1stPosition(m_initialPosition);
  int iRC = UQ_OK_RC;
#if 0
  unsigned int chainSumId = 0;
  std::vector<const P_V*> chainSum(0);//,NULL);
  if (m_avgChainCompute.size() > 0) {
    // It is expected that all participating chains will have the same size.
    // The code will check this assumption.
    chainSum.resize(m_chainSizes[0],NULL);
  }
#endif
  for (unsigned int chainId = 0; chainId < 1/*m_chainSizes.size()*/; ++chainId) { // FIXME: make a definitive change
    char tmpChainId[10];
    sprintf(tmpChainId,"%d",chainId);
    //std::string prefixName = m_prefix; // + "c" + tmpChainId + "_"; FIXME: make a definitive change
    workingChain.setName(m_prefix + "chain");

    if (m_chainType == UQ_MAC_SG_WHITE_NOISE_CHAIN_TYPE) {
      //****************************************************
      // Just generate white noise
      //****************************************************
      generateWhiteNoiseChain(m_chainSizes[chainId],
                              workingChain);
    }
    else if (m_chainType == UQ_MAC_SG_UNIFORM_CHAIN_TYPE) {
      //****************************************************
      // Just generate uniform    
      //****************************************************
      generateUniformChain(m_chainSizes[chainId],
                           workingChain);
    }
    else {
      //****************************************************
      // Initialize variables before chain loop
      //****************************************************
      if (chainId > 0) {
        workingChain.getPositionValues(workingChain.sequenceSize()-1,valuesOf1stPosition);
        resetChainAndRelatedInfo();
      }

      //****************************************************
      // Initialize m_lowerCholProposalCovMatrices[0]
      // Initialize m_proposalCovMatrices[0]
      //****************************************************
      iRC = prepareForNextChain(proposalCovMatrix);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences()",
                        "improper prepareForNextChain() return");

      //****************************************************
      // Generate chain
      //****************************************************
      intGenerateSequence(valuesOf1stPosition,
                          workingChain,
                          m_chainSizes[chainId]);
    }

    //****************************************************
    // Open file      
    //****************************************************
    std::ofstream* ofs = NULL;
    if (m_chainOutputFileNames[chainId] == UQ_MAC_SG_FILENAME_FOR_NO_OUTPUT_FILE) {
      if (m_env.rank() == 0) {
        std::cout << "No output file opened for chain loop id = " << chainId
                  << std::endl;
      }
    }
    else {
      if (m_env.rank() == 0) {
        std::cout << "Opening output file '"  << m_chainOutputFileNames[chainId]
                  << "' for chain loop id = " << chainId
                  << " ..."
                  << std::endl;
      }

      // Open file
      ofs = new std::ofstream(m_chainOutputFileNames[chainId].c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
      if ((ofs            == NULL ) ||
          (ofs->is_open() == false)) {
        delete ofs;
        ofs = new std::ofstream(m_chainOutputFileNames[chainId].c_str(), std::ofstream::out | std::ofstream::trunc);
      }

      UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
                          m_env.rank(),
                          "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences()",
                          "failed to open file");
    }
  
    //****************************************************
    // Eventually:
    // --> write chain
    // --> compute statistics on it
    //****************************************************
    if (m_chainWrite && ofs) {
      workingChain.printContents(*ofs);

      // Write likelihoodValues and alphaValues, if they were requested by user
      iRC = writeInfo(workingChain,
                      *ofs);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences()",
                        "improper writeInfo() return");
    }

    if (m_chainComputeStats) {
      workingChain.computeStatistics(*m_chainStatisticalOptions,
                                     ofs);
    }

    //****************************************************
    // Eventually:
    // --> generate an unique chain
    // --> write it
    // --> compute statistics on it
    //****************************************************
    if (m_uniqueChainGenerate) {
      // Select only the unique positions
      workingChain.select(m_idsOfUniquePositions);
      //chainVectorPositionIteratorTypedef positionIterator = m_uniqueChain1.begin();
      //std::advance(positionIterator,uniquePos);
      //m_uniqueChain1.erase(positionIterator,m_uniqueChain1.end());
      //UQ_FATAL_TEST_MACRO((uniquePos != m_uniqueChain1.size()),
      //                    m_env.rank(),
      //                    "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences()",
      //                    "uniquePos != m_uniqueChain1.size()");

      // Write unique chain
      workingChain.setName(m_prefix + "uniqueChain");
      if (m_uniqueChainWrite && ofs) {
        workingChain.printContents(*ofs);
      }

      // Compute statistics
      if (m_uniqueChainComputeStats) {
        workingChain.computeStatistics(*m_uniqueChainStatisticalOptions,
                                       ofs);
      }
    }

    //****************************************************
    // Eventually:
    // --> filter the (maybe unique) chain
    // --> write it
    // --> compute statistics on it
    //****************************************************
    if (m_filteredChainGenerate) {
      // Compute filter parameters
      unsigned int filterInitialPos = (unsigned int) (m_filteredChainDiscardedPortion * (double) workingChain.sequenceSize());
      unsigned int filterSpacing    = m_filteredChainLag;
      if (filterSpacing == 0) {
        workingChain.computeFilterParams(*m_filteredChainStatisticalOptions,
                                         ofs,
                                         filterInitialPos,
                                         filterSpacing);
      }

      // Filter positions from the converged portion of the chain
      workingChain.filter(filterInitialPos,
                          filterSpacing);

      // Write filtered chain
      workingChain.setName(m_prefix + "filteredChain");
      if (m_filteredChainWrite && ofs) {
        workingChain.printContents(*ofs);
      }

      // Compute statistics
      if (m_filteredChainComputeStats) {
        workingChain.computeStatistics(*m_filteredChainStatisticalOptions,
                                       ofs);
      }
    }

    //****************************************************
    // Eventually:
    // --> compute an average chain
    // --> write it
    // --> compute statistics on it
    //****************************************************
#if 0
    if (m_avgChainCompute.size() > chainSumId) {
      // Update workingChainSum

      // Check if it is time to compute an average
      if ((chainId+1) == m_avgChainCompute[chainSumId]) {
        // Compute the average

        // Write the computed average
        
        // Compute statistics on the average

        // Prepare for eventual next chain average
        chainSumId++;
      }
    }
#endif

    //****************************************************
    // Close file      
    //****************************************************
    if (ofs) {
      // Close file
      ofs->close();

      if (m_env.rank() == 0) {
        std::cout << "Closed output file '"   << m_chainOutputFileNames[chainId]
                  << "' for chain loop id = " << chainId
                  << " ..."
                  << std::endl;
      }
    }
    if (m_env.rank() == 0) {
      std::cout << std::endl;
    }
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqMarkovChainSGClass<P_V,P_M>::intGenerateSequences()"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::generateWhiteNoiseChain(
  unsigned int                        chainSize,
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain)
{
  if (m_env.rank() == 0) {
    std::cout << "Starting the generation of white noise chain " << workingChain.name()
              << ", with "                                       << chainSize
              << " positions ..."
              << std::endl;
  }

  int iRC;
  struct timeval timevalTmp;
  double tmpRunTime;

  tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  workingChain.resizeSequence(chainSize); 

  P_V meanVec  (m_vectorSpace.zeroVector());
  P_V stdDevVec(m_vectorSpace.zeroVector());
  meanVec.cwSet(0.);
  stdDevVec.cwSet(1.);
  workingChain.setGaussian(m_env.rng(),meanVec,stdDevVec);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);

  if (m_env.rank() == 0) {
    std::cout << "Finished the generation of white noise chain " << workingChain.name()
              << ", with "                                       << workingChain.sequenceSize()
              << " positions";
  }

  if (m_env.rank() == 0) {
    std::cout << "Chain generation took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::generateUniformChain(
  unsigned int                        chainSize,
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain)
{
  if (m_env.rank() == 0) {
    std::cout << "Starting the generation of uniform chain " << workingChain.name()
              << ", with "                                   << chainSize
              << " positions ..."
              << std::endl;
  }

  int iRC;
  struct timeval timevalTmp;
  double tmpRunTime;

  tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  workingChain.resizeSequence(chainSize); 

  P_V aVec(m_vectorSpace.zeroVector());
  P_V bVec(m_vectorSpace.zeroVector());
  aVec.cwSet(0.);
  bVec.cwSet(1.);
  workingChain.setUniform(m_env.rng(),aVec,bVec);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);

  if (m_env.rank() == 0) {
    std::cout << "Finished the generation of white noise chain " << workingChain.name()
              << ", with "                                       << workingChain.sequenceSize()
              << " positions";
  }

  if (m_env.rank() == 0) {
    std::cout << "Chain generation took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence(
  const P_V&                          valuesOf1stPosition,
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  unsigned int                        chainSize)
{
  if (m_env.rank() == 0) {
    std::cout << "Starting the generation of Markov chain " << workingChain.name()
              << ", with "                                  << chainSize
              << " positions..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalChain;
  struct timeval timevalCandidate;
  struct timeval timevalTargetD;
  struct timeval timevalMhAlpha;
  struct timeval timevalDrAlpha;
  struct timeval timevalDR;
  struct timeval timevalAM;

  double candidateRunTime = 0;
  double targetDRunTime   = 0;
  double mhAlphaRunTime   = 0;
  double drAlphaRunTime   = 0;
  double drRunTime        = 0;
  double amRunTime        = 0;

  iRC = gettimeofday(&timevalChain, NULL);
  bool outOfDomainBounds = m_targetDensity.outOfDomainBounds(valuesOf1stPosition);
  UQ_FATAL_TEST_MACRO(outOfDomainBounds,
                      m_env.rank(),
                      "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()",
                      "initial position should not be out of bound");
  if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalTargetD, NULL);
  double logTarget = -0.5 * m_targetDensity.minus2LnDensity(valuesOf1stPosition);
  if (m_chainMeasureRunTimes) targetDRunTime += uqMiscGetEllapsedSeconds(&timevalTargetD);
  uqMarkovChainPositionClass<P_V> currentPosition(m_env,
                                            valuesOf1stPosition,
                                            outOfDomainBounds,
                                            logTarget);

  P_V gaussianVector(m_vectorSpace.zeroVector());
  P_V tmpVecValues(m_vectorSpace.zeroVector());
  uqMarkovChainPositionClass<P_V> currentCandidate(m_env);

  //****************************************************
  // Begin chain loop from positionId = 1
  //****************************************************
  workingChain.resizeSequence(chainSize); 
  if (m_uniqueChainGenerate) m_idsOfUniquePositions.resize(chainSize,0); 
  if (m_chainGenerateExtra) {
    m_logTargets.resize (chainSize,0.);
    m_alphaQuotients.resize(chainSize,0.);
  }

  unsigned int uniquePos = 0;
  workingChain.setPositionValues(0,currentPosition.vecValues());
  if (m_uniqueChainGenerate) m_idsOfUniquePositions[uniquePos++] = 0;
  if (m_chainGenerateExtra) {
    m_logTargets [0] = currentPosition.logTarget();
    m_alphaQuotients[0] = 1.;
  }

  for (unsigned int positionId = 1; positionId < workingChain.sequenceSize(); ++positionId) {
    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
                << ": beginning chain position of id = " << positionId
                << ", m_maxNumExtraStages =  "           << m_maxNumExtraStages
                << std::endl;
    }
    unsigned int stageId = 0;

    //****************************************************
    // Loop: generate new position
    //****************************************************
    if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
    gaussianVector.cwSetGaussian(m_env.rng(),0.,1.);
    tmpVecValues = currentPosition.vecValues() + *(m_lowerCholProposalCovMatrices[stageId]) * gaussianVector;
    if (m_chainMeasureRunTimes) candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

    outOfDomainBounds = m_targetDensity.outOfDomainBounds(tmpVecValues);
    if (outOfDomainBounds) {
      m_numOutOfBounds++;
      logTarget = -INFINITY;
    }
    else {
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalTargetD, NULL);
      logTarget = -0.5 * m_targetDensity.minus2LnDensity(tmpVecValues);
      if (m_chainMeasureRunTimes) targetDRunTime += uqMiscGetEllapsedSeconds(&timevalTargetD);
    }
    currentCandidate.set(tmpVecValues,
                         outOfDomainBounds,
                         logTarget);

    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      std::cout << "\n-----------------------------------------------------------\n"
                << std::endl;
    }
    bool accept = false;
    if (outOfDomainBounds) {
      if (m_chainGenerateExtra) {
        m_alphaQuotients[positionId] = 0.;
      }
    }
    else {
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalMhAlpha, NULL);
      double alpha = 0.;
      if (m_chainGenerateExtra) {
        alpha = this->alpha(currentPosition,currentCandidate,&m_alphaQuotients[positionId]);
      }
      else {
        alpha = this->alpha(currentPosition,currentCandidate,NULL);
      }
      if (m_chainMeasureRunTimes) mhAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalMhAlpha);
      if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
        std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
                  << ": for chain position of id = " << positionId
                  << ", alpha = " << alpha
                  << std::endl;
      }
      accept = acceptAlpha(alpha);
    }
    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      if (m_env.rank() == 0) std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
                                       << ": for chain position of id = " << positionId
                                       << " contents of currentCandidate.vecValues() are:"
                                       << std::endl;
      std::cout << currentCandidate.vecValues();
      if (m_env.rank() == 0) std::cout << std::endl;

      if (m_env.rank() == 0) std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
                                       << ": for chain position of id = " << positionId
                                       << ", outOfDomainBounds = "        << outOfDomainBounds
                                       << "\n"
                                       << "\n curLogTarget  = "           << currentPosition.logTarget()
                                       << "\n"
                                       << "\n canLogTarget  = "           << currentCandidate.logTarget()
                                       << "\n"
                                       << "\n accept = "                  << accept
                                       << std::endl;
    }
    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      std::cout << "\n-----------------------------------------------------------\n"
                << std::endl;
    }

    //****************************************************
    // Loop: delayed rejection
    //****************************************************
    std::vector<uqMarkovChainPositionClass<P_V>*> drPositions(stageId+2,NULL);
    if ((accept == false) && (outOfDomainBounds == false) && (m_maxNumExtraStages > 0)) {
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalDR, NULL);

      drPositions[0] = new uqMarkovChainPositionClass<P_V>(currentPosition);
      drPositions[1] = new uqMarkovChainPositionClass<P_V>(currentCandidate);

      while ((accept == false) && (stageId < m_maxNumExtraStages)) {
        stageId++;

        if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
        gaussianVector.cwSetGaussian(m_env.rng(),0.,1.);
        tmpVecValues = currentPosition.vecValues() + *(m_lowerCholProposalCovMatrices[stageId]) * gaussianVector;
        if (m_chainMeasureRunTimes) candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

        outOfDomainBounds = m_targetDensity.outOfDomainBounds(tmpVecValues);
        if (outOfDomainBounds) {
          logTarget = -INFINITY;
        }
        else {
          if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalTargetD, NULL);
          logTarget = -0.5 * m_targetDensity.minus2LnDensity(tmpVecValues);
          if (m_chainMeasureRunTimes) targetDRunTime += uqMiscGetEllapsedSeconds(&timevalTargetD);
        }
        currentCandidate.set(tmpVecValues,
                             outOfDomainBounds,
                             logTarget);

        drPositions.push_back(new uqMarkovChainPositionClass<P_V>(currentCandidate));
        if (outOfDomainBounds == false) {
          if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalDrAlpha, NULL);
          double alpha = this->alpha(drPositions);
          if (m_chainMeasureRunTimes) drAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalDrAlpha);
#if 0 // For debug only
          if (m_env.rank() == 0) std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
                                           << ": for chain position of id = " << positionId
                                           << " and stageId = " << stageId
                                           << ", alpha = " << alpha
                                           << std::endl;
#endif
          accept = acceptAlpha(alpha);
        }
      } // while

      if (m_chainMeasureRunTimes) drRunTime += uqMiscGetEllapsedSeconds(&timevalDR);
    } // end of 'delayed rejection' logic

    for (unsigned int i = 0; i < drPositions.size(); ++i) {
      if (drPositions[i]) delete drPositions[i];
    }

    //****************************************************
    // Loop: update chain
    //****************************************************
    if (accept) {
      workingChain.setPositionValues(positionId,currentCandidate.vecValues());
      if (m_uniqueChainGenerate) m_idsOfUniquePositions[uniquePos++] = positionId;
      currentPosition = currentCandidate;
    }
    else {
      workingChain.setPositionValues(positionId,currentPosition.vecValues());
      m_numRejections++;
    }

    if (m_chainGenerateExtra) {
      m_logTargets[positionId] = currentPosition.logTarget();
    }

#ifdef UQ_MAC_SG_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
    if (m_likelihoodObjComputesMisfits) {
      if (m_observableSpace.shouldVariancesBeUpdated()) {
        L_V misfitVec    (currentPosition.misfitVector()           );
        L_V numbersOfObs (m_observableSpace.numbersOfObservations());
        L_V varAccuracies(m_observableSpace.varianceAccuracies()   );
        L_V priorVars    (m_observableSpace.priorVariances()       );
        for (unsigned int i = 0; i < misfitVarianceVector.size(); ++i) {
          double term1 = 0.5*( varAccuracies[i] + numbersOfObs[i]             );
          double term2 =  2./( varAccuracies[i] * priorVars[i] + misfitVec[i] );
          misfitVarianceVector[i] = 1./uqMiscGammar(term1,term2,m_env.rng());
          //if (m_env.rank() == 0) {
          //  std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
          //            << ": for chain position of id = "     << positionId
          //            << ", numbersOfObs = "                 << numbersOfObs
          //            << ", varAccuracies = "                << varAccuracies
          //            << ", priorVars = "                    << priorVars
          //            << ", (*m_misfitChain[positionId]) = " << (*m_misfitChain[positionId])
          //            << ", term1 = "                        << term1
          //            << ", term2 = "                        << term2
          //            << std::endl;
          //}
        }
        //if (m_env.rank() == 0) {
        //  std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
        //            << ": for chain position of id = "        << positionId
        //            << ", misfitVarianceVector changed from " << *(m_misfitVarianceChain[positionId])
        //            << " to "                                 << misfitVarianceVector
        //            << std::endl;
        //}
      }
      if (m_chainGenerateExtra) {
        m_misfitVarianceChain[positionId] = m_observableSpace.newVector(misfitVarianceVector);
      }
    }
#endif

    //****************************************************
    // Loop: adaptive Metropolis (adaptation of covariance matrix)
    //****************************************************
    if ((m_initialNonAdaptInterval > 0) &&
        (m_adaptInterval           > 0)) {
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalAM, NULL);

      // Now might be the moment to adapt
      unsigned int idOfFirstPositionInSubChain = 0;
      uqSequenceOfVectorsClass<P_V,P_M> subChain(m_vectorSpace,0,m_prefix+"subChain");

      // Check if now is indeed the moment to adapt
      if (positionId < m_initialNonAdaptInterval) {
        // Do nothing
      }
      else if (positionId == m_initialNonAdaptInterval) {
        idOfFirstPositionInSubChain = 0;
        subChain.resizeSequence(m_initialNonAdaptInterval+1);
        m_lastMean             = m_vectorSpace.newVector();
        m_lastAdaptedCovMatrix = m_vectorSpace.newMatrix();
      }
      else {
        unsigned int interval = positionId - m_initialNonAdaptInterval;
        if ((interval % m_adaptInterval) == 0) {
          idOfFirstPositionInSubChain = positionId - m_adaptInterval;
          subChain.resizeSequence(m_adaptInterval);
        }
      }

      // If now is indeed the moment to adapt, then do it!
      if (subChain.sequenceSize() > 0) {
        P_V transporterVec(m_vectorSpace.zeroVector());
        for (unsigned int i = 0; i < subChain.sequenceSize(); ++i) {
          workingChain.getPositionValues(idOfFirstPositionInSubChain+i,transporterVec);
          subChain.setPositionValues(i,transporterVec);
        }
        updateCovMatrix(subChain,
			idOfFirstPositionInSubChain,
                        m_lastChainSize,
                        *m_lastMean,
                        *m_lastAdaptedCovMatrix);

        bool tmpCholIsPositiveDefinite = false;
        P_M tmpChol(*m_lastAdaptedCovMatrix);
        //if (m_env.rank() == 0) {
        //  std::cout << "DRAM: chainId = " << chainId
        //            << ", positionId = "  << positionId
        //            << ": 'am' calling first tmpChol.chol()"
        //            << std::endl;
        //}
        iRC = tmpChol.chol();
        //if (m_env.rank() == 0) {
        //  std::cout << "DRAM: chainId = " << chainId
        //            << ", positionId = "  << positionId
        //            << ": 'am' got first tmpChol.chol() with iRC = " << iRC
        //            << std::endl;
        //}
        if (iRC) {
          UQ_FATAL_TEST_MACRO(iRC != UQ_MATRIX_IS_NOT_POS_DEFINITE_RC,
                              m_env.rank(),
                              "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()",
                              "invalid iRC returned from first chol()");
          // Matrix is not positive definite
          P_M* tmpDiag = m_vectorSpace.newDiagMatrix(m_epsilon);
          tmpChol = *m_lastAdaptedCovMatrix + *tmpDiag;
          delete tmpDiag;
          //if (m_env.rank() == 0) {
          //  std::cout << "DRAM: chainId = " << chainId
          //            << ", positionId = "  << positionId
          //            << ": 'am' calling second tmpChol.chol()"
          //            << std::endl;
          //}
          iRC = tmpChol.chol();
          //if (m_env.rank() == 0) {
          //  std::cout << "DRAM: chainId = " << chainId
          //            << ", positionId = "  << positionId
          //            << ": 'am' got second tmpChol.chol() with iRC = " << iRC
          //            << std::endl;
          //}
          if (iRC) {
            UQ_FATAL_TEST_MACRO(iRC != UQ_MATRIX_IS_NOT_POS_DEFINITE_RC,
                                m_env.rank(),
                                "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()",
                                "invalid iRC returned from second chol()");
            // Do nothing
          }
          else {
            tmpCholIsPositiveDefinite = true;
          }
        }
        else {
          tmpCholIsPositiveDefinite = true;
        }
        if (tmpCholIsPositiveDefinite) {
          *(m_lowerCholProposalCovMatrices[0]) = tmpChol;
          *(m_lowerCholProposalCovMatrices[0]) *= sqrt(m_eta);
          m_lowerCholProposalCovMatrices[0]->zeroUpper(false);

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
          UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                            m_env.rank(),
                            "uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()",
                            "need to code the update of m_upperCholProposalPrecMatrices");
#endif

          if (m_maxNumExtraStages > 0) updateCovMatrices();
        }

        //for (unsigned int i = 0; i < subChain.sequenceSize(); ++i) {
        //  if (subChain[i]) delete subChain[i];
        //}
      }

      if (m_chainMeasureRunTimes) amRunTime += uqMiscGetEllapsedSeconds(&timevalAM);
    } // End of 'adaptive Metropolis' logic

    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      std::cout << "In uqMarkovChainSGClass<P_V,P_M>::intGenerateSequence()"
                << ": finishing chain position of id = " << positionId
                << std::endl;
    }

    if ((m_chainDisplayPeriod                     > 0) && 
        (((positionId+1) % m_chainDisplayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << positionId+1
                  << " positions"
                  << std::endl;
      }
    }
  } // end chain loop

  //****************************************************
  // Print basic information about the chain
  //****************************************************
  m_chainRunTime += uqMiscGetEllapsedSeconds(&timevalChain);
  if (m_env.rank() == 0) {
    std::cout << "Finished the generation of Markov chain " << workingChain.name()
              << ", with "                                  << workingChain.sequenceSize()
              << " positions";
    if (m_uniqueChainGenerate) {
      std::cout << " and " << uniquePos
                << " 'unique' positions (i.e., not counting repetitions due to rejections)";
    }
    std::cout << "\nSome information about this chain:"
              << "\n  Chain run time       = " << m_chainRunTime
              << " seconds";
    if (m_chainMeasureRunTimes) {
      std::cout << "\n\n Breaking of the chain run time:\n";
      std::cout << "\n  Candidate run time   = " << candidateRunTime
                << " seconds ("                  << 100.*candidateRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  Target d. run time   = " << targetDRunTime
                << " seconds ("                  << 100.*targetDRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  Mh alpha run time    = " << mhAlphaRunTime
                << " seconds ("                  << 100.*mhAlphaRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  Dr alpha run time    = " << drAlphaRunTime
                << " seconds ("                  << 100.*drAlphaRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n----------------------   --------------";
      double sumRunTime = candidateRunTime + targetDRunTime + mhAlphaRunTime + drAlphaRunTime;
      std::cout << "\n  Sum                  = " << sumRunTime
                << " seconds ("                  << 100.*sumRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n\n Other run times:";
      std::cout << "\n  DR run time          = " << drRunTime
                << " seconds ("                  << 100.*drRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  AM run time          = " << amRunTime
                << " seconds ("                  << 100.*amRunTime/m_chainRunTime
                << "%)";
    }
    std::cout << "\n  Rejection percentage = " << 100. * (double) m_numRejections/(double) workingChain.sequenceSize()
              << " %";
    std::cout << "\n   Outbound percentage = " << 100. * (double) m_numOutOfBounds/(double) workingChain.sequenceSize()
              << " %";
    std::cout << std::endl;
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************

  return;
}

template <class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::updateCovMatrix(
  const uqBaseVectorSequenceClass<P_V,P_M>& subChain,
  unsigned int                          idOfFirstPositionInSubChain,
  double&                               lastChainSize,
  P_V&                                  lastMean,
  P_M&                                  lastAdaptedCovMatrix)
{
  double doubleSubChainSize = (double) subChain.sequenceSize();
  if (lastChainSize == 0) {
    UQ_FATAL_TEST_MACRO(subChain.sequenceSize() < 2,
                        m_env.rank(),
                        "uqMarkovChainSGClass<P_V,P_M>::updateCovMatrix()",
                        "'subChain.sequenceSize()' should be >= 2");

    subChain.mean(0,subChain.sequenceSize(),lastMean);

    P_V tmpVec(m_vectorSpace.zeroVector());
    lastAdaptedCovMatrix = -doubleSubChainSize * matrixProduct(lastMean,lastMean);
    for (unsigned int i = 0; i < subChain.sequenceSize(); ++i) {
      subChain.getPositionValues(i,tmpVec);
      lastAdaptedCovMatrix += matrixProduct(tmpVec,tmpVec);
    }
    lastAdaptedCovMatrix /= (doubleSubChainSize - 1.); // That is why subChain size must be >= 2
  }
  else {
    UQ_FATAL_TEST_MACRO(subChain.sequenceSize() < 1,
                        m_env.rank(),
                        "uqMarkovChainSGClass<P_V,P_M>::updateCovMatrix()",
                        "'subChain.sequenceSize()' should be >= 1");

    UQ_FATAL_TEST_MACRO(idOfFirstPositionInSubChain < 1,
                        m_env.rank(),
                        "uqMarkovChainSGClass<P_V,P_M>::updateCovMatrix()",
                        "'idOfFirstPositionInSubChain' should be >= 1");

    P_V tmpVec (m_vectorSpace.zeroVector());
    P_V diffVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 0; i < subChain.sequenceSize(); ++i) {
      double doubleCurrentId  = (double) (idOfFirstPositionInSubChain+i);
      subChain.getPositionValues(i,tmpVec);
      diffVec = tmpVec - lastMean;

      double ratio1         = (1. - 1./doubleCurrentId); // That is why idOfFirstPositionInSubChain must be >= 1
      double ratio2         = (1./(1.+doubleCurrentId));
      lastAdaptedCovMatrix  = ratio1 * lastAdaptedCovMatrix + ratio2 * matrixProduct(diffVec,diffVec);
      lastMean             += ratio2 * diffVec;
    } 
  }
  lastChainSize += doubleSubChainSize;

  return;
}
#endif // __UQ_MAC_SG2_H__

