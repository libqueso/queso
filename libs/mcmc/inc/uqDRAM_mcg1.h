/* uq/libs/mcmc/inc/uqDRAM_mcg1.h
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

#ifndef __UQ_DRAM_MCG1_H__
#define __UQ_DRAM_MCG1_H__

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1(
  const M* proposalCovMatrix,
  const M* mahalanobisMatrix,
  bool     applyMahalanobisInvert)
{
  //if (m_env.rank() == 0) std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1()..."
  //                                 << std::endl;

  V valuesOf1stPosition(m_paramInitials);
  int iRC = UQ_OK_RC;
  for (unsigned int chainId = 0; chainId < m_sizesOfChains.size(); ++chainId) {
    if (m_generateWhiteNoise) {
      //****************************************************
      // Just generate white noise
      //****************************************************
      iRC = generateWhiteNoise1(chainId);
    }
    else {
      //****************************************************
      // Initialize variables before chain1 loop
      //****************************************************
      if (chainId > 0) {
        valuesOf1stPosition = *(m_chain1[m_chain1.size()-1]);
        resetChainAndRelatedInfo();
      }

      //****************************************************
      // Initialize m_lowerCholProposalCovMatrices[0]
      // Initialize m_proposalCovMatrices[0]
      //****************************************************
      iRC = prepareForNextChain(proposalCovMatrix);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1()",
                        "improper prepareForNextChain() return");

      //****************************************************
      // Generate chain
      //****************************************************
      iRC = generateChain1(chainId,
                           valuesOf1stPosition,
                           proposalCovMatrix,
                           mahalanobisMatrix,
                           applyMahalanobisInvert);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1()",
                        "improper generateChain1() return");
    }

    //****************************************************
    // Open file      
    //****************************************************
    std::ofstream* ofs = NULL;
    if (m_namesOfOutputFiles[chainId] == ".") {
      if (m_env.rank() == 0) {
        std::cout << "No output file opened for chain1 of id " << chainId
                  << std::endl;
      }
    }
    else {
      if (m_env.rank() == 0) {
        std::cout << "Opening output file for chain1 of id " << chainId
                  << ", with name '"                         << m_namesOfOutputFiles[chainId]
                  << "'..."
                  << std::endl;
      }

      // Open file
      ofs = new std::ofstream(m_namesOfOutputFiles[chainId].c_str(), std::ofstream::out | std::ofstream::trunc);
      UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
                          m_env.rank(),
                          "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()",
                          "failed to open file");
    }
  
    //****************************************************
    // Print statistics on the chain
    //****************************************************
    if (m_env.rank() == 0) {
      std::cout << "\n"
                << "\n-----------------------------------------------------"
                << "\n Statistics on chain:"
                << "\n-----------------------------------------------------"
                << "\n"
                << std::endl;
    }
    computeStatistics(m_chain1,m_chain2,ofs);
#if 0
    if (m_env.rank() == 0) {
      std::cout << "\n"
                << "\n-----------------------------------------------------"
                << "\n Statistics on 'unique' chain:"
                << "\n-----------------------------------------------------"
                << "\n"
                << std::endl;
    }
    computeStatistics1(m_uniqueChain1);
#endif

    //****************************************************
    // Write chain1 out
    //****************************************************
    if (ofs) {
      iRC = writeChain(m_chain1,
                       m_chain2,
                       *ofs,
                       mahalanobisMatrix,
                       applyMahalanobisInvert);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1()",
                        "improper writeChain() return");
    }

    //****************************************************
    // Close file      
    //****************************************************
    if (ofs) {
      // Close file
      ofs->close();

      if (m_env.rank() == 0) {
        std::cout << "Closed output file for chain1 of id " << chainId
                  << std::endl;
      }
    }
    if (m_env.rank() == 0) {
      std::cout << std::endl;
    }
  }

  //if (m_env.rank() == 0) std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1()"
  //                                 << std::endl;

  return;
}

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::generateWhiteNoise1(unsigned int chainId)
{
  if (m_env.rank() == 0) {
    std::cout << "Generating white noise for chain1 of id " << chainId
              << ", with "                                  << m_sizesOfChains[chainId]
              << " positions..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  double tmpRunTime;

  V gaussianVector(m_paramSpace.zeroVector());

  tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  m_chain1.resize(m_sizesOfChains[chainId],NULL); 
  std::vector<const V*>(m_chain1).swap(m_chain1);
  for (unsigned int positionId = 0; positionId < m_chain1.size(); ++positionId) {
    gaussianVector.cwSetGaussian(m_env.rng(),0.,1.);
    m_chain1[positionId] = new V(gaussianVector);

    if ((m_chainDisplayPeriod                     > 0) && 
        (((positionId+1) % m_chainDisplayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << positionId+1
                  << " positions"
                  << std::endl;
      }
    }
  }
  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain1 generation took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return iRC;
}

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1(
  unsigned int chainId,
  const V&     valuesOf1stPosition,
  const M*     proposalCovMatrix,
  const M*     mahalanobisMatrix,
  bool         applyMahalanobisInvert)
{
  if (m_env.rank() == 0) {
    std::cout << "Generating chain1 of id " << chainId
              << ", with "                  << m_sizesOfChains[chainId]
              << " positions..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalChain;
  struct timeval timevalCandidate;
  struct timeval timevalPrior;
  struct timeval timevalLH;
  struct timeval timevalMhAlpha;
  struct timeval timevalDrAlpha;
  struct timeval timevalDR;
  struct timeval timevalAM;

  iRC = gettimeofday(&timevalChain, NULL);

  bool   outOfBounds = m_paramSpace.outOfBounds(valuesOf1stPosition);
  UQ_FATAL_TEST_MACRO(outOfBounds,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()",
                      "paramInitials should not be out of bound");
  if (m_measureRunTimes) iRC = gettimeofday(&timevalPrior, NULL);
  double m2lPrior            = m_m2lPriorProbDensity_Obj.minus2LnDensity(valuesOf1stPosition);
  if (m_measureRunTimes) m_priorRunTime += uqMiscGetEllapsedSeconds(&timevalPrior);
  V*     m2lLikelihoodVector = m_observableSpace.newVector();
  V*     misfitVector        = m_observableSpace.newVector();
  V      misfitVarianceVector(m_observableSpace.priorVariances());

  if (m_measureRunTimes) iRC = gettimeofday(&timevalLH, NULL);
  if (m_likelihoodObjComputesMisfits) {
    m_m2lLikelihoodFunction_Obj.computeMisfits(valuesOf1stPosition, *misfitVector);
    *m2lLikelihoodVector = *misfitVector/misfitVarianceVector;
  }
  else {
    m_m2lLikelihoodFunction_Obj.computeMinus2LnLikelihoods(valuesOf1stPosition, *m2lLikelihoodVector);
  }
  if (m_measureRunTimes) m_lhRunTime += uqMiscGetEllapsedSeconds(&timevalLH);
  double m2lLikelihoodScalar  = m2lLikelihoodVector->sumOfComponents();
  double logPosterior = -0.5 * ( m2lPrior + m2lLikelihoodScalar );
  uqChainPositionClass<V> currentPosition(m_env,
                                          valuesOf1stPosition,
                                          outOfBounds,
                                          m2lPrior,
                                          *misfitVector,
                                          misfitVarianceVector,
                                          *m2lLikelihoodVector,
                                          logPosterior);

  V* gaussianVector = m_paramSpace.newVector();
  V* tmpParamValues = m_paramSpace.newVector();
  uqChainPositionClass<V> currentCandidate(m_env);

  //****************************************************
  // Begin chain1 loop from positionId = 1
  //****************************************************
  m_chain1.resize(m_sizesOfChains[chainId],NULL); 
  if (m_generateUniqueChain) m_uniqueChain1.resize(m_sizesOfChains[chainId],NULL); 
  if (m_generateExtraChains) {
    if (m_likelihoodObjComputesMisfits) {
      m_misfitChain.resize        (m_sizesOfChains[chainId],NULL); 
      m_misfitVarianceChain.resize(m_sizesOfChains[chainId],NULL); 
    }
    m_m2lLikelihoodChain.resize(m_sizesOfChains[chainId],NULL); 
    m_alphaQuotients.resize    (m_sizesOfChains[chainId],0.);
  }

  m_chain1[0] = new V(currentPosition.paramValues());
  if (m_generateUniqueChain) m_uniqueChain1[m_uniqueChain1Pos++] = new V(*(m_chain1[0]));
  if (m_generateExtraChains) {
    if (m_likelihoodObjComputesMisfits) {
      m_misfitChain        [0] = m_observableSpace.newVector(*misfitVector);
      m_misfitVarianceChain[0] = m_observableSpace.newVector(misfitVarianceVector);
    }
    m_m2lLikelihoodChain[0] = m_observableSpace.newVector(currentPosition.m2lLikelihoodVector());
    m_alphaQuotients    [0] = 1.;
  }

  for (unsigned int positionId = 1; positionId < m_chain1.size(); ++positionId) {
    //if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
    //                                 << ": beginning chain1 position of id = " << positionId
    //                                 << ", m_maxNumberOfExtraStages =  "       << m_maxNumberOfExtraStages
    //                                 << std::endl;
    unsigned int stageId = 0;

    //****************************************************
    // Loop: generate new parameters
    //****************************************************
    if (m_measureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
    gaussianVector->cwSetGaussian(m_env.rng(),0.,1.);
    *tmpParamValues = currentPosition.paramValues() + *(m_lowerCholProposalCovMatrices[stageId]) * *gaussianVector;
    if (m_measureRunTimes) m_candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

    outOfBounds     = m_paramSpace.outOfBounds(*tmpParamValues);
    if (outOfBounds) {
      m_numOutOfBounds++;
      m2lPrior      = 0.;
      m2lLikelihoodVector->cwSet(INFINITY);
      logPosterior  = -INFINITY;
    }
    else {
      if (m_measureRunTimes) iRC = gettimeofday(&timevalPrior, NULL);
      m2lPrior      = m_m2lPriorProbDensity_Obj.minus2LnDensity(*tmpParamValues);
      if (m_measureRunTimes) m_priorRunTime += uqMiscGetEllapsedSeconds(&timevalPrior);
      if (m_measureRunTimes) iRC = gettimeofday(&timevalLH, NULL);
      if (m_likelihoodObjComputesMisfits) {
        m_m2lLikelihoodFunction_Obj.computeMisfits(*tmpParamValues, *misfitVector);
        *m2lLikelihoodVector = *misfitVector/misfitVarianceVector;
      }
      else {
        m_m2lLikelihoodFunction_Obj.computeMinus2LnLikelihoods(*tmpParamValues, *m2lLikelihoodVector);
      }
      if (m_measureRunTimes) m_lhRunTime += uqMiscGetEllapsedSeconds(&timevalLH);
      m2lLikelihoodScalar = m2lLikelihoodVector->sumOfComponents();
      logPosterior  = -0.5 * ( m2lPrior + m2lLikelihoodScalar );
    }
    currentCandidate.set(*tmpParamValues,
                         outOfBounds,
                         m2lPrior,
                         *misfitVector,
                         misfitVarianceVector,
                         *m2lLikelihoodVector,
                         logPosterior);

    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      std::cout << "\n-----------------------------------------------------------\n"
                << std::endl;
    }
    bool accept = false;
    if (outOfBounds) {
      if (m_generateExtraChains) {
        m_alphaQuotients[positionId] = 0.;
      }
    }
    else {
      if (m_measureRunTimes) iRC = gettimeofday(&timevalMhAlpha, NULL);
      double alpha = 0.;
      if (m_generateExtraChains) {
        alpha = this->alpha(currentPosition,currentCandidate,&m_alphaQuotients[positionId]);
      }
      else {
        alpha = this->alpha(currentPosition,currentCandidate,NULL);
      }
      if (m_measureRunTimes) m_mhAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalMhAlpha);
      if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
        std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
                  << ": for chain1 position of id = " << positionId
                  << ", alpha = " << alpha
                  << std::endl;
      }
      accept = acceptAlpha(alpha);
    }
    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
                                       << ": for chain1 position of id = " << positionId
                                       << " contents of currentCandidate.paramValues() are:"
                                       << std::endl;
      std::cout << currentCandidate.paramValues();
      if (m_env.rank() == 0) std::cout << std::endl;

      if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
                                       << ": for chain1 position of id = " << positionId
                                       << ", outOfBounds = "               << outOfBounds
                                       << "\n"
                                       << "\n curM2lPrior = "             << currentPosition.m2lPrior()
                                       << "\n curMisfitVector = "         << currentPosition.misfitVector()
                                       << "\n curMisfitVarianceVector = " << currentPosition.misfitVarianceVector()
                                       << "\n curM2lLikelihoodVector = "  << currentPosition.m2lLikelihoodVector()
                                       << "\n curLogPosterior = "         << currentPosition.logPosterior()
                                       << "\n"
                                       << "\n canM2lPrior = "             << currentCandidate.m2lPrior()
                                       << "\n canMisfitVector = "         << currentCandidate.misfitVector()
                                       << "\n canMisfitVarianceVector = " << currentCandidate.misfitVarianceVector()
                                       << "\n canM2lLikelihoodVector = "  << currentCandidate.m2lLikelihoodVector()
                                       << "\n canLogPosterior = "         << currentCandidate.logPosterior()
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
    std::vector<uqChainPositionClass<V>*> drPositions(stageId+2,NULL);
    if ((accept == false) && (outOfBounds == false) && (m_maxNumberOfExtraStages > 0)) {
      if (m_measureRunTimes) iRC = gettimeofday(&timevalDR, NULL);

      drPositions[0] = new uqChainPositionClass<V>(currentPosition);
      drPositions[1] = new uqChainPositionClass<V>(currentCandidate);

      while ((accept == false) && (stageId < m_maxNumberOfExtraStages)) {
        stageId++;

        if (m_measureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
        gaussianVector->cwSetGaussian(m_env.rng(),0.,1.);
        *tmpParamValues = currentPosition.paramValues() + *(m_lowerCholProposalCovMatrices[stageId]) * *gaussianVector;
        if (m_measureRunTimes) m_candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

        outOfBounds   = m_paramSpace.outOfBounds(*tmpParamValues);
        if (outOfBounds) {
          m2lPrior      = 0.;
          m2lLikelihoodVector->cwSet(INFINITY);
          logPosterior  = -INFINITY;
        }
        else {
          if (m_measureRunTimes) iRC = gettimeofday(&timevalPrior, NULL);
          m2lPrior      = m_m2lPriorProbDensity_Obj.minus2LnDensity(*tmpParamValues);
          if (m_measureRunTimes) m_priorRunTime += uqMiscGetEllapsedSeconds(&timevalPrior);
          if (m_measureRunTimes) iRC = gettimeofday(&timevalLH, NULL);
          if (m_likelihoodObjComputesMisfits) {
            m_m2lLikelihoodFunction_Obj.computeMisfits(*tmpParamValues, *misfitVector);
            *m2lLikelihoodVector = *misfitVector/misfitVarianceVector;
          }
          else {
            m_m2lLikelihoodFunction_Obj.computeMinus2LnLikelihoods(*tmpParamValues, *m2lLikelihoodVector);
          }
          if (m_measureRunTimes) m_lhRunTime += uqMiscGetEllapsedSeconds(&timevalLH);
          m2lLikelihoodScalar = m2lLikelihoodVector->sumOfComponents();
          logPosterior  = -0.5 * ( m2lPrior + m2lLikelihoodScalar );
        }
        currentCandidate.set(*tmpParamValues,
                             outOfBounds,
                             m2lPrior,
                             *misfitVector,
                             misfitVarianceVector,
                             *m2lLikelihoodVector,
                             logPosterior);

        drPositions.push_back(new uqChainPositionClass<V>(currentCandidate));
        if (outOfBounds == false) {
          if (m_measureRunTimes) iRC = gettimeofday(&timevalDrAlpha, NULL);
          double alpha = this->alpha(drPositions);
          if (m_measureRunTimes) m_drAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalDrAlpha);
#if 0 // For debug only
          if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
                                           << ": for chain1 position of id = " << positionId
                                           << " and stageId = " << stageId
                                           << ", alpha = " << alpha
                                           << std::endl;
#endif
          accept = acceptAlpha(alpha);
        }
      } // while

      if (m_measureRunTimes) m_drRunTime += uqMiscGetEllapsedSeconds(&timevalDR);
    } // end of 'delayed rejection' logic

    for (unsigned int i = 0; i < drPositions.size(); ++i) {
      if (drPositions[i]) delete drPositions[i];
    }

    //****************************************************
    // Loop: update chain
    //****************************************************
    if (accept) {
      m_chain1[positionId] = new V(currentCandidate.paramValues());
      if (m_generateUniqueChain) m_uniqueChain1[m_uniqueChain1Pos++] = new V(*(m_chain1[positionId]));
      if (m_generateExtraChains) {
        if (m_likelihoodObjComputesMisfits) {
          m_misfitChain[positionId] = m_observableSpace.newVector(currentCandidate.misfitVector());
          //m_misfitVarianceChain[positionId] is updated below, after the update of 'misfitVarianceVector'
        }
        m_m2lLikelihoodChain[positionId] = m_observableSpace.newVector(currentCandidate.m2lLikelihoodVector());
      }
      currentPosition = currentCandidate;
    }
    else {
      m_chain1[positionId] = new V(currentPosition.paramValues());
      if (m_generateExtraChains) {
        if (m_likelihoodObjComputesMisfits) {
          m_misfitChain[positionId] = m_observableSpace.newVector(currentPosition.misfitVector());
          //m_misfitVarianceChain[positionId] is updated below, after the update of 'misfitVarianceVector'
        }
        m_m2lLikelihoodChain[positionId] = m_observableSpace.newVector(currentPosition.m2lLikelihoodVector());
      }
      m_numRejections++;
    }

    if (m_likelihoodObjComputesMisfits) {
      if (m_observableSpace.shouldVariancesBeUpdated()) {
        V misfitVec    (currentPosition.misfitVector()           );
        V numbersOfObs (m_observableSpace.numbersOfObservations());
        V varAccuracies(m_observableSpace.varianceAccuracies()   );
        V priorVars    (m_observableSpace.priorVariances()       );
        for (unsigned int i = 0; i < misfitVarianceVector.size(); ++i) {
          double term1 = 0.5*( varAccuracies[i] + numbersOfObs[i]             );
          double term2 =  2./( varAccuracies[i] * priorVars[i] + misfitVec[i] );
          misfitVarianceVector[i] = 1./uqMiscGammar(term1,term2,m_env.rng());
          //if (m_env.rank() == 0) {
          //  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
          //            << ": for chain1 position of id = "    << positionId
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
        //  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()"
        //            << ": for chain1 position of id = "       << positionId
        //            << ", misfitVarianceVector changed from " << *(m_misfitVarianceChain[positionId])
        //            << " to "                                 << misfitVarianceVector
        //            << std::endl;
        //}
      }
      if (m_generateExtraChains) {
        m_misfitVarianceChain[positionId] = m_observableSpace.newVector(misfitVarianceVector);
      }
    }

    //****************************************************
    // Loop: adaptive Metropolis (adaptation of covariance matrix)
    //****************************************************
    if ((m_initialNonAdaptInterval > 0) &&
        (m_adaptInterval           > 0)) {
      if (m_measureRunTimes) iRC = gettimeofday(&timevalAM, NULL);

      // Now might be the moment to adapt
      unsigned int idOfFirstPositionInSubChain = 0;
      std::vector<V*> subChain1(0);//,NULL);

      // Check if now is indeed the moment to adapt
      if (positionId < m_initialNonAdaptInterval) {
        // Do nothing
      }
      else if (positionId == m_initialNonAdaptInterval) {
        idOfFirstPositionInSubChain = 0;
        subChain1.resize(m_initialNonAdaptInterval+1,NULL);
        m_lastMean             = m_paramSpace.newVector();
        m_lastAdaptedCovMatrix = m_paramSpace.newMatrix();
      }
      else {
        unsigned int interval = positionId - m_initialNonAdaptInterval;
        if ((interval % m_adaptInterval) == 0) {
          idOfFirstPositionInSubChain = positionId - m_adaptInterval;
          subChain1.resize(m_adaptInterval,NULL);
        }
      }

      // If now is indeed the moment to adapt, then do it!
      if (subChain1.size() > 0) {
        for (unsigned int i = 0; i < subChain1.size(); ++i) {
          subChain1[i] = new V(*(m_chain1[idOfFirstPositionInSubChain+i]));
        }
        updateCovMatrix1(subChain1,
			 idOfFirstPositionInSubChain,
                         m_lastChainSize,
                         *m_lastMean,
                         *m_lastAdaptedCovMatrix);

        bool tmpCholIsPositiveDefinite = false;
        M tmpChol(*m_lastAdaptedCovMatrix);
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
                              "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()",
                              "invalid iRC returned from first chol()");
          // Matrix is not positive definite
          M* tmpDiag = m_paramSpace.newDiagMatrix(m_epsilon);
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
                                "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()",
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
                            "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()",
                            "need to code the update of m_upperCholProposalPrecMatrices");
#endif

          if (m_maxNumberOfExtraStages > 0) updateCovMatrices();
        }

        for (unsigned int i = 0; i < subChain1.size(); ++i) {
          if (subChain1[i]) delete subChain1[i];
        }
      }

      if (m_measureRunTimes) m_amRunTime += uqMiscGetEllapsedSeconds(&timevalAM);
    } // End of 'adaptive Metropolis' logic

    if ((m_chainDisplayPeriod                     > 0) && 
        (((positionId+1) % m_chainDisplayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << positionId+1
                  << " positions"
                  << std::endl;
      }
    }
  } // end chain1 loop

  if (m_generateUniqueChain) {
    chainPositionIteratorTypedef positionIterator = m_uniqueChain1.begin();
    std::advance(positionIterator,m_uniqueChain1Pos);
    m_uniqueChain1.erase(positionIterator,m_uniqueChain1.end());
    UQ_FATAL_TEST_MACRO((m_uniqueChain1Pos != m_uniqueChain1.size()),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain1()",
                        "m_uniqueChain1Pos != m_uniqueChain1.size()");
  }

  //****************************************************
  // Print basic information about the chain
  //****************************************************
  m_chainRunTime += uqMiscGetEllapsedSeconds(&timevalChain);
  if (m_env.rank() == 0) {
    if (m_generateUniqueChain) {
      std::cout << "Finished generating the chain1 of id " << chainId
                << ", with "                               << m_uniqueChain1Pos
                << " 'unique' positions (i.e., not counting repetitions due to rejections).";
    }
    else {
      std::cout << "Finished generating the chain1 of id " << chainId
                << ", with "                               << m_chain1.size()
                << " positions.";
    }
    std::cout << "\nSome information about this chain:"
              << "\n  Chain1 run time       = " << m_chainRunTime
              << " seconds";
    if (m_measureRunTimes) {
      std::cout << "\n\n Breaking of the chain1 run time:\n";
      std::cout << "\n  Candidate run time   = " << m_candidateRunTime
                << " seconds ("                  << 100.*m_candidateRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  Prior run time       = " << m_priorRunTime
                << " seconds ("                  << 100.*m_priorRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  LH run time          = " << m_lhRunTime
                << " seconds ("                  << 100.*m_lhRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  Mh alpha run time    = " << m_mhAlphaRunTime
                << " seconds ("                  << 100.*m_mhAlphaRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  Dr alpha run time    = " << m_drAlphaRunTime
                << " seconds ("                  << 100.*m_drAlphaRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n----------------------   --------------";
      double sumRunTime = m_candidateRunTime + m_priorRunTime + m_lhRunTime + m_mhAlphaRunTime + m_drAlphaRunTime;
      std::cout << "\n  Sum                  = " << sumRunTime
                << " seconds ("                  << 100.*sumRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n\n Other run times:";
      std::cout << "\n  DR run time          = " << m_drRunTime
                << " seconds ("                  << 100.*m_drRunTime/m_chainRunTime
                << "%)";
      std::cout << "\n  AM run time          = " << m_amRunTime
                << " seconds ("                  << 100.*m_amRunTime/m_chainRunTime
                << "%)";
    }
    std::cout << "\n  Rejection percentage = " << 100. * (double) m_numRejections/(double) m_chain1.size()
              << " %";
    std::cout << "\n   Outbound percentage = " << 100. * (double) m_numOutOfBounds/(double) m_chain1.size()
              << " %";
    std::cout << std::endl;
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************
  if (gaussianVector           ) delete gaussianVector;
  if (misfitVector             ) delete misfitVector;
  if (m2lLikelihoodVector      ) delete m2lLikelihoodVector;
  if (tmpParamValues           ) delete tmpParamValues;

  return iRC;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix1(
  const std::vector<V*>& subChain1,
  unsigned int           idOfFirstPositionInSubChain,
  double&                lastChainSize,
  V&                     lastMean,
  M&                     lastAdaptedCovMatrix)
{
  double doubleSubChainSize = (double) subChain1.size();
  if (lastChainSize == 0) {
    UQ_FATAL_TEST_MACRO(subChain1.size() < 2,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix1()",
                        "'subChain1.size()' should be >= 2");

    lastMean.cwSet(0.);
    double ratio = 1./doubleSubChainSize;
    for (unsigned int i = 0; i < subChain1.size(); ++i) {
      lastMean += ratio * *(subChain1[i]);
    }

    lastAdaptedCovMatrix = -doubleSubChainSize * matrixProduct(lastMean,lastMean);
    for (unsigned int i = 0; i < subChain1.size(); ++i) {
      lastAdaptedCovMatrix += matrixProduct(*(subChain1[i]),*(subChain1[i]));
    }
    lastAdaptedCovMatrix /= (doubleSubChainSize - 1.); // That is why subChain1 size must be >= 2
  }
  else {
    UQ_FATAL_TEST_MACRO(subChain1.size() < 1,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix1()",
                        "'subChain1.size()' should be >= 1");

    UQ_FATAL_TEST_MACRO(idOfFirstPositionInSubChain < 1,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix1()",
                        "'idOfFirstPositionInSubChain' should be >= 1");

    for (unsigned int i = 0; i < subChain1.size(); ++i) {
      double doubleCurrentId  = (double) (idOfFirstPositionInSubChain+i);
      V diffVec(*(subChain1[i]) - lastMean);

      double ratio1         = (1. - 1./doubleCurrentId); // That is why idOfFirstPositionInSubChain must be >= 1
      double ratio2         = (1./(1.+doubleCurrentId));
      lastAdaptedCovMatrix  = ratio1 * lastAdaptedCovMatrix + ratio2 * matrixProduct(diffVec,diffVec);
      lastMean             += ratio2 * diffVec;
    } 
  }
  lastChainSize += doubleSubChainSize;

  return;
}
#endif // __UQ_DRAM_MCG1_H__

