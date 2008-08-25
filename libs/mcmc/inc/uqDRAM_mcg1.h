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
    computeStatistics1(m_chain1,ofs);
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
      iRC = writeChain1(*ofs,
                        mahalanobisMatrix,
                        applyMahalanobisInvert);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains1()",
                        "improper writeChain1() return");
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
uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics1(
  const std::vector<const V*>& chain1,
  std::ofstream* passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  double tmpRunTime;

  // Set initial positions for the computation of chain statistics
  std::vector<unsigned int> initialPosForStatistics(m_finalPercentsForStats.size(),0);
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    initialPosForStatistics[i] = (unsigned int) ((1. - m_finalPercentsForStats[i]*0.01) * (double) chain1.size());
  }
  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics1(): initial positions for statistics =";
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    std::cout << " " << initialPosForStatistics[i];
  }
  std::cout << std::endl;

  //****************************************************
  // Compute mean, sample std, population std
  //****************************************************
  tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  V chainMean(m_paramSpace.zeroVector());
  uqVectorSequenceMean(chain1,
                       0,
                       chain1.size(),
                       chainMean);

  V chainSampleVariance(m_paramSpace.zeroVector());
  uqVectorSequenceSampleVariance(chain1,
                                 0,
                                 chain1.size(),
                                 chainMean,
                                 chainSampleVariance);

  V chainPopulationVariance(m_paramSpace.zeroVector());
  uqVectorSequencePopulationVariance(chain1,
                                     0,
                                     chain1.size(),
                                     chainMean,
                                     chainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain1 statistics took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  if (m_env.rank() == 0) {
    std::cout << "\nMean, sample std, population std"
              << std::endl;
    char line[512];
    sprintf(line,"%s%4s%s%9s%s%9s%s",
	    "Parameter",
            " ",
            "Mean",
            " ",
            "SampleStd",
            " ",
            "Popul.Std");
    std::cout << line;

    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e",
              m_paramSpace.parameter(i).name().c_str(),
              " ",
	      chainMean[i],
              " ",
              sqrt(chainSampleVariance[i]),
              " ",
              sqrt(chainPopulationVariance[i]));
      std::cout << line;
    }
    std::cout << std::endl;
  }

  //****************************************************
  // Compute variance of sample mean through the 'batch means method' (BMM)
  //****************************************************

  if ((m_runBMM                          ) &&
      (initialPosForStatistics.size() > 0)) {
    std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics1(): lengths for batchs in BMM =";
    for (unsigned int i = 0; i < m_lengthsForBMM.size(); ++i) {
      std::cout << " " << m_lengthsForBMM[i];
    }
    std::cout << std::endl;
  }

  if ((m_runBMM                          ) &&
      (initialPosForStatistics.size() > 0) &&
      (m_lengthsForBMM.size()         > 0)) { 
    if (m_env.rank() == 0) {
      std::cout << "\nComputing variance of sample mean through BMM..."
                << std::endl;
    }
    uq2dArrayOfStuff<V> _2dArrayOfBMM(initialPosForStatistics.size(),m_lengthsForBMM.size());
    for (unsigned int i = 0; i < _2dArrayOfBMM.numRows(); ++i) {
      for (unsigned int j = 0; j < _2dArrayOfBMM.numCols(); ++j) {
        _2dArrayOfBMM.setLocation(i,j,m_paramSpace.newVector());
      }
    }
    uqVectorSequenceBMM(chain1,
                        initialPosForStatistics,
                        m_lengthsForBMM,
                        _2dArrayOfBMM);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        std::cout << "\nEstimated variances of sample mean, through batch means method, for subchain beggining at position " << initialPosForStatistics[initialPosId]
                  << " (each column corresponds to a batch length)"
                  << std::endl;

        char line[512];
        sprintf(line,"%s",
                "Parameter");
        std::cout << line;
        for (unsigned int batchLengthId = 0; batchLengthId < m_lengthsForBMM.size(); batchLengthId++) {
          sprintf(line,"%10s%3d",
                  " ",
                  m_lengthsForBMM[batchLengthId]);
          std::cout << line;
        }

        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          sprintf(line,"\n%9.9s",
                  m_paramSpace.parameter(i).name().c_str());
          std::cout << line;
          for (unsigned int batchLengthId = 0; batchLengthId < m_lengthsForBMM.size(); batchLengthId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
                    _2dArrayOfBMM(initialPosId,batchLengthId)[i]);
            std::cout << line;
          }
        }
        std::cout << std::endl;
      }
    }
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain1 at zero frequency
  //****************************************************
  if ((m_computePSDs                     ) &&
      (initialPosForStatistics.size() > 0) &&
      (m_numBlocksForPSD.size()       > 0)) { 
    if (m_env.rank() == 0) {
      std::cout << "\nComputing PSD at frequency zero..."
                << std::endl;
    }
    uq2dArrayOfStuff<V> _2dArrayOfPSDAtZero(initialPosForStatistics.size(),m_numBlocksForPSD.size());
    for (unsigned int i = 0; i < _2dArrayOfPSDAtZero.numRows(); ++i) {
      for (unsigned int j = 0; j < _2dArrayOfPSDAtZero.numCols(); ++j) {
        _2dArrayOfPSDAtZero.setLocation(i,j,m_paramSpace.newVector());
      }
    }
    uqVectorSequencePSD(chain1,
                        initialPosForStatistics,
                        m_numBlocksForPSD,
                        m_hopSizeRatioForPSD,
                        _2dArrayOfPSDAtZero);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        double sizeForPSD = chain1.size() - initialPosForStatistics[initialPosId];
        std::cout << "\nEstimated variances of sample mean, through psd (fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                  << ", so effective data size = " << sizeForPSD
                  << " (each column corresponds to a number of blocks)"
                  << std::endl;

        char line[512];
        sprintf(line,"%s",
                "Parameter");
        std::cout << line;
        for (unsigned int numBlocksId = 0; numBlocksId < m_numBlocksForPSD.size(); numBlocksId++) {
          sprintf(line,"%10s%3d",
                  " ",
                  m_numBlocksForPSD[numBlocksId]);
          std::cout << line;
        }

        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          sprintf(line,"\n%9.9s",
                  m_paramSpace.parameter(i).name().c_str());
          std::cout << line;
          for (unsigned int numBlocksId = 0; numBlocksId < m_numBlocksForPSD.size(); numBlocksId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
                    M_PI*_2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]/sizeForPSD); // CHECK
            std::cout << line;
          }
        }
        std::cout << std::endl;
      }
    }
  }

  //****************************************************
  // Compute Geweke
  //****************************************************
  if ((m_computeGewekeCoefs              ) &&
      (initialPosForStatistics.size() > 0)) {
    if (m_env.rank() == 0) {
      std::cout << "\nComputing Geweke coefficients..."
                << std::endl;
    }

    std::vector<V*> vectorOfGeweke(initialPosForStatistics.size(),NULL);
    uqVectorSequenceGeweke(chain1,
                           initialPosForStatistics,
                           m_ratioNaForGeweke,
                           m_ratioNbForGeweke,
                           vectorOfGeweke);

    if (m_env.rank() == 0) {
      std::cout << "\nComputed Geweke coefficients with 10% and 50% percentages"
                  << " (each column corresponds to a different initial position on the full chain)"
                  << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        sprintf(line,"%10s%3d",
                " ",
                initialPosForStatistics[initialPosId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%9.9s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  (*(vectorOfGeweke[initialPosId]))[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  //****************************************************
  // Compute autocorrelation coefficients
  //****************************************************

  // Set lags for the computation of chain autocorrelations
  std::vector<unsigned int> lagsForCorrs(m_numberOfLagsForCorrs,1);
  if ((m_computeCorrelations             ) &&
      (initialPosForStatistics.size() > 0)) {
    for (unsigned int i = 1; i < lagsForCorrs.size(); ++i) {
      lagsForCorrs[i] = m_secondLagForCorrs + (i-1)*m_lagSpacingForCorrs;
    }

    if ((m_printCorrs) && (m_env.rank() == 0)) {
      std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics1(): lags for autocorrelation =";
      for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
        std::cout << " " << lagsForCorrs[i];
      }
      std::cout << std::endl;
    }
  }

  if ((m_computeCorrelations             ) &&
      (initialPosForStatistics.size() > 0) &&
      (lagsForCorrs.size()            > 0)) { 
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);

    if (m_env.rank() == 0) {
      std::cout << "\nComputing autocorrelation coefficients..."
                << std::endl;
    }
    uq2dArrayOfStuff<V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
    for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
      for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
        _2dArrayOfAutoCorrs.setLocation(i,j,m_paramSpace.newVector());
      }
    }
    uqVectorSequenceAutoCorrelations(chain1,
                                     initialPosForStatistics,
                                     lagsForCorrs,
                                     _2dArrayOfAutoCorrs);

    V estimatedVarianceOfSampleMean(m_paramSpace.zeroVector());
    estimatedVarianceOfSampleMean.cwSet(1.);
    for (unsigned int i = 0; i < 1; ++i) {
      for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
        double ratio = ((double) j)/((double) m_chain1.size());
        estimatedVarianceOfSampleMean += 2.*(1.-ratio)*_2dArrayOfAutoCorrs(i,j);
      }
    }
    estimatedVarianceOfSampleMean *= chainSampleVariance;
    estimatedVarianceOfSampleMean /= (double) m_chain1.size();
    if (m_env.rank() == 0) {
      bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
      estimatedVarianceOfSampleMean.setPrintHorizontally(false);
      std::cout << "\nVariance of sample mean, estimated through autocorrelation:\n"
                << estimatedVarianceOfSampleMean
                << std::endl;
      estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);
    }

    if ((m_printCorrs) && (m_env.rank() == 0)) {
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        std::cout << "\nEstimated autocorrelation coefficients, for subchain beggining at position " << initialPosForStatistics[initialPosId]
                  << " (each column corresponds to a different lag)"
                  << std::endl;

        char line[512];
        sprintf(line,"%s",
	        "Parameter");
        std::cout << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%10s%3d",
                  " ",
                  lagsForCorrs[lagId]);
          std::cout << line;
        }

        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          sprintf(line,"\n%9.9s",
                  m_paramSpace.parameter(i).name().c_str());
          std::cout << line;
          for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
	            _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
            std::cout << line;
          }
        }
        std::cout << std::endl;
      }
    }

    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.rank() == 0) {
      std::cout << "Chain1 autocorrelation took " << tmpRunTime
                << " seconds"
                << std::endl;
    }

    // Write autocorrelations
    if (m_writeCorrs && passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << "queso_" << m_prefix << "corrs_lags = zeros(" << 1
          << ","                                           << lagsForCorrs.size()
          << ");"
          << std::endl;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        ofs << "queso_" << m_prefix << "corrs_lags(" << 1
            << ","                                   << lagId+1
            << ") = "                                << lagsForCorrs[lagId]
            << ";"
            << std::endl;
      }

      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        ofs << "queso_" << m_prefix << "corrs_ip" << initialPosForStatistics[initialPosId] << " = zeros(" << m_paramSpace.dim()
            << ","      << lagsForCorrs.size()
            << ");"
            << std::endl;
        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
            ofs << "queso_" << m_prefix << "corrs_ip" << initialPosForStatistics[initialPosId] << "(" << i+1
                << ","                                                                                << lagId+1
                << ") = "                                                                             << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
                << ";"
                << std::endl;
          }
        }
      }
    } 
  }

  //****************************************************
  // Compute MIN and MAX: for histograms and KDE
  //****************************************************
  unsigned int initialPosForUncorrelation = m_initialPosForUncorrelation;
  if (initialPosForUncorrelation >= chain1.size()) initialPosForUncorrelation = chain1.size() - 1;
  uqVectorSequenceMinMax(chain1,
                         initialPosForUncorrelation,
                         *m_minPositionsForStatistics,
                         *m_maxPositionsForStatistics);

  //****************************************************
  // Compute histograms
  //****************************************************
  if ((m_computeHistograms               ) &&
      (m_numberOfInternalBinsForHists > 0)) {
    if (m_env.rank() == 0) {
      std::cout << "\nComputing histograms..."
                << std::endl;
    }
    for (unsigned int i = 0; i < m_maxPositionsForStatistics->size(); ++i) {
      (*m_maxPositionsForStatistics)[i] *= (1. + 1.e-15);
    }

    m_centersForAllHistogramBins.resize(m_numberOfInternalBinsForHists+2,NULL);
    m_histogramBinsForAllParams.resize (m_numberOfInternalBinsForHists+2,NULL);
    uqVectorSequenceHistogram(chain1,
                              initialPosForUncorrelation,
                              m_spacingForUncorrelation,
                              *m_minPositionsForStatistics,
                              *m_maxPositionsForStatistics,
                              m_centersForAllHistogramBins,
                              m_histogramBinsForAllParams);
  }

  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((m_computeKDEs                     ) &&
      (m_numberOfEvaluationPosForKDEs > 0)) {
    if (m_env.rank() == 0) {
      std::cout << "\nComputing KDE..."
                << std::endl;
    }

    m_evaluationPositionsForKDEs.resize(m_numberOfEvaluationPosForKDEs,NULL);
    uqMiscComputePositionsBetweenMinMax(*m_minPositionsForStatistics,
                                        *m_maxPositionsForStatistics,
                                        m_evaluationPositionsForKDEs);

    V iqrs(*(chain1[0]));
    uqVectorSequenceInterQuantileRange(chain1,
                                       initialPosForUncorrelation,
                                       m_spacingForUncorrelation,
                                       iqrs);
    if (m_env.rank() == 0) {
      std::cout << "\nComputed min values, max values and inter quantile ranges, for subchain beggining at position " << initialPosForUncorrelation
                  << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;

      sprintf(line,"%9s%s%9s%s%9s%s",
              " ",
              "min",
              " ",
              "max",
              " ",
              "iqr");
      std::cout << line;

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%8.8s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;

        sprintf(line,"%2s%11.4e%2s%11.4e%2s%11.4e",
                " ",
                (*m_minPositionsForStatistics)[i],
                " ",
                (*m_maxPositionsForStatistics)[i],
                " ",
                iqrs[i]);
        std::cout << line;
      }
      std::cout << std::endl;
    }

    uqVectorSequenceScalesForKDE(chain1,
                                 initialPosForUncorrelation,
                                 m_spacingForUncorrelation,
                                 iqrs,
                                 *m_scalesForKDEs);

    m_densityValuesFromGaussianKDE.resize(m_numberOfEvaluationPosForKDEs,NULL);
    uqVectorSequenceGaussianKDE(chain1,
                                initialPosForUncorrelation,
                                m_spacingForUncorrelation,
                                m_evaluationPositionsForKDEs,
                                *m_scalesForKDEs,
                                m_densityValuesFromGaussianKDE);
  }

  return;
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

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::writeChain1(
  std::ofstream& ofs,
  const M*       mahalanobisMatrix,
  bool           applyMahalanobisInvert) const
{
  int iRC = UQ_OK_RC;

  // Write m_chain1

  ofs << "queso_" << m_prefix << "chain = [";
  for (unsigned int i = 0; i < m_chain1.size(); ++i) {
    ofs << *(m_chain1[i])
        << std::endl;
  }
  ofs << "];\n";


  if (m_generateExtraChains) {
    if (m_likelihoodObjComputesMisfits) {
      // Write m_misfitChain
      ofs << "queso_" << m_prefix << "misfitChain = [";
      for (unsigned int i = 0; i < m_misfitChain.size(); ++i) {
        ofs << *(m_misfitChain[i])
            << std::endl;
      }
      ofs << "];\n";

      // Write m_misfitVarianceChain
      ofs << "queso_" << m_prefix << "misfitVarianceChain = [";
      for (unsigned int i = 0; i < m_misfitVarianceChain.size(); ++i) {
        ofs << *(m_misfitVarianceChain[i])
            << std::endl;
      }
      ofs << "];\n";
    }

    // Write m_m2lLikelihoodChain
    ofs << "queso_" << m_prefix << "m2lLikelihoodChain = [";
    for (unsigned int i = 0; i < m_m2lLikelihoodChain.size(); ++i) {
      ofs << *(m_m2lLikelihoodChain[i])
          << std::endl;
    }
    ofs << "];\n";

    // Write m_alphaQuotients
    ofs << "queso_" << m_prefix << "alphaQuotients = [";
    for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
      ofs << m_alphaQuotients[i]
          << std::endl;
    }
    ofs << "];\n";
  }

  // Write names of parameters
  ofs << "queso_" << m_prefix << "paramNames = {";
  m_paramSpace.printParameterNames(ofs,false);
  ofs << "};\n";

  // Write mahalanobis distances
  if (mahalanobisMatrix != NULL) {
    V* diffVec = m_paramSpace.newVector();
    ofs << "queso_" << m_prefix << "d = [";
    if (applyMahalanobisInvert) {
      for (unsigned int i = 0; i < m_chain1.size(); ++i) {
        *diffVec = *(m_chain1[i]) - *(m_chain1[0]);
        ofs << scalarProduct(*diffVec, mahalanobisMatrix->invertMultiply(*diffVec))
            << std::endl;
      }
    }
    else {
      for (unsigned int i = 0; i < m_chain1.size(); ++i) {
        *diffVec = *(m_chain1[i]) - *(m_chain1[0]);
        ofs << scalarProduct(*diffVec, *mahalanobisMatrix * *diffVec)
            << std::endl;
      }
    }
    ofs << "];\n";
    delete diffVec;
  }

  // Write prior mean values
  ofs << "queso_" << m_prefix << "priorMeanValues = ["
      << m_paramSpace.priorMuValues()
      << "];\n";

  // Write prior sigma values
  ofs << "queso_" << m_prefix << "priorSigmaValues = ["
      << m_paramSpace.priorSigmaValues()
      << "];\n";

#if 0
  ofs << "queso_" << m_prefix << "results.prior = [queso_priorMeanValues',queso_priorSigmaValues'];\n";
#endif

  // Write param lower bounds
  ofs << "queso_" << m_prefix << "minValues = ["
      << m_paramSpace.minValues()
      << "];\n";

  // Write param upper bounds
  ofs << "queso_" << m_prefix << "maxValues = ["
      << m_paramSpace.maxValues()
      << "];\n";

#if 0
  ofs << "queso_" << m_prefix << "results.limits = [queso_low',queso_upp'];\n";

  // Write out data for mcmcpred.m
  ofs << "queso_" << m_prefix << "results.parind = ["; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << i+1
        << std::endl;
  }
  ofs << "];\n";

  ofs << "queso_" << m_prefix << "results.local = [\n"; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << " 0";
    //<< std::endl;
  }
  ofs << "];\n";

  bool savedVectorPrintState = m_chain1[m_chain1.size()-1]->getPrintHorizontally();
  m_chain1[m_chain1.size()-1]->setPrintHorizontally(false);
  ofs << "queso_" << m_prefix << "results.theta = ["
      << *(m_chain1[m_chain1.size()-1])
      << "];\n";
  m_chain1[m_chain1.size()-1]->setPrintHorizontally(savedVectorPrintState);
  
  ofs << "queso_" << m_prefix << "results.nbatch = 1.;\n"; // FIXME

  if (mahalanobisMatrix != NULL) {
    // Write covMatrix
    ofs << "queso_" << m_prefix << "mahalanobisMatrix = ["
        << *mahalanobisMatrix
        << "];\n";
  }
#endif

  // Write number of rejections
  ofs << "queso_" << m_prefix << "rejected = "  << (double) m_numRejections/(double) (m_chain1.size()-1)
      << ";\n"
      << std::endl;

  // Write chain1 run time
  ofs << "queso_" << m_prefix << "runTime = "  << m_chainRunTime
      << ";\n"
      << std::endl;

  // Write number of outbounds
  ofs << "queso_" << m_prefix << "outbounds = " << (double) m_numOutOfBounds/(double) m_chain1.size()
      << ";\n"
      << std::endl;

  // Write histograms
  if ((m_computeHistograms               ) &&
      (m_numberOfInternalBinsForHists > 0)) {
    // plot(queso_centersOfHistBins(1,:)',queso_histBins(1,:)','r-');
    ofs << "queso_" << m_prefix << "centersOfHistBins = zeros(" << m_paramSpace.dim()
        << ","                                                  << m_centersForAllHistogramBins.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_centersForAllHistogramBins.size(); ++j) {
         ofs << "queso_" << m_prefix << "centersOfHistBins(" << i+1
             << ","                                          << j+1
             << ") = "                                       << (*(m_centersForAllHistogramBins[j]))[i]
             << ";"
             << std::endl;
      }
    }

    ofs << "queso_" << m_prefix << "histBins = zeros(" << m_paramSpace.dim()
        << ","                       << m_histogramBinsForAllParams.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_histogramBinsForAllParams.size(); ++j) {
         ofs << "queso_" << m_prefix << "histBins(" << i+1
             << ","                                 << j+1
             << ") = "                              << (*(m_histogramBinsForAllParams[j]))[i]
             << ";"
             << std::endl;
      }
    }
  }

  // Write estimations of probability densities
  if ((m_computeKDEs                     ) &&
      (m_numberOfEvaluationPosForKDEs > 0)) {
    // hold
    // plot(queso_evalPosForKDE(1,:)',7*queso_densFromGaussianKDE(1,:)','r-');
    ofs << "queso_" << m_prefix << "evalPosForKDE = zeros(" << m_paramSpace.dim()
        << ","                                              << m_evaluationPositionsForKDEs.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_evaluationPositionsForKDEs.size(); ++j) {
        ofs << "queso_" << m_prefix << "evalPosForKDE(" << i+1
            << ","                                      << j+1
            << ") = "                                   << (*(m_evaluationPositionsForKDEs[j]))[i]
            << ";"
            << std::endl;
      }
    }

    ofs << "queso_" << m_prefix << "scalesForKDE = zeros(" << m_paramSpace.dim()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      ofs << "queso_" << m_prefix << "scalesForKDE(" << i+1
          << ") = "                                  << (*m_scalesForKDEs)[i]
          << ";"
          << std::endl;
    }

    ofs << "queso_" << m_prefix << "densFromGaussianKDE = zeros(" << m_paramSpace.dim()
        << ","                                                    << m_densityValuesFromGaussianKDE.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_densityValuesFromGaussianKDE.size(); ++j) {
        ofs << "queso_" << m_prefix << "densFromGaussianKDE(" << i+1
            << ","                                            << j+1
            << ") = "                                         << (*(m_densityValuesFromGaussianKDE[j]))[i]
            << ";"
            << std::endl;
      }
    }
  }

  return iRC;
}
#endif // __UQ_DRAM_MCG1_H__

