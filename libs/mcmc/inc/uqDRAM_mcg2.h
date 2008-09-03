/* uq/libs/mcmc/inc/uqDRAM_mcg2.h
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

#ifndef __UQ_DRAM_MCG2_H__
#define __UQ_DRAM_MCG2_H__

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains2(
  const M* proposalCovMatrix,
  const M* mahalanobisMatrix,
  bool     applyMahalanobisInvert)
{
  //if (m_env.rank() == 0) std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains2()..."
  //                                 << std::endl;

  V valuesOf1stPosition(m_paramInitials);
  int iRC = UQ_OK_RC;

  unsigned int chainSumId = 0;
  uqArrayOfSequencesClass<V> chain2Sum(0,m_paramSpace.zeroVector());
  if (m_avgChainCompute.size() > 0) {
    // It is expected that all participating chains will have the same size.
    // The code will check this assumption.
    chain2Sum.resizeSequence(m_chainSizes[0]);
  }

  for (unsigned int chainId = 0; chainId < m_chainSizes.size(); ++chainId) {
    char tmpChainId[10];
    sprintf(tmpChainId,"%d",chainId);
    std::string chainName = "queso_" + m_prefix + tmpChainId + "_chain";
    std::string prefixName = "queso_" + m_prefix + tmpChainId;

    //****************************************************
    // Open file      
    //****************************************************
    std::ofstream* ofs = NULL;
    if (m_chainOutputFileNames[chainId] == UQ_MCMC_NAME_FOR_NO_OUTPUT_FILE) {
      if (m_env.rank() == 0) {
        std::cout << "No output file opened for chain loop id = " << chainId
                  << std::endl;
      }
    }
    else {
      if (m_env.rank() == 0) {
        std::cout << "Opening output file '"  << m_chainOutputFileNames[chainId]
                  << "' for chain loop id = " << chainId
                  << "'..."
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
                          "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains2()",
                          "failed to open file");
    }
  
    if (m_chainType == UQ_MCMC_WHITE_NOISE_CHAIN_TYPE) {
      //****************************************************
      // Just generate white noise
      //****************************************************
      iRC = generateWhiteNoise2(m_chainSizes[chainId],
                                chainName);
    }
    else {
      //****************************************************
      // Initialize variables before chain2 loop
      //****************************************************
      if (chainId > 0) {
        m_chain2.getPositionValues(m_chain2.sequenceSize()-1,valuesOf1stPosition);
        resetChainAndRelatedInfo();
      }

      //****************************************************
      // Initialize m_lowerCholProposalCovMatrices[0]
      // Initialize m_proposalCovMatrices[0]
      //****************************************************
      iRC = prepareForNextChain(proposalCovMatrix);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains2()",
                        "improper prepareForNextChain() return");

      //****************************************************
      // Generate chain
      //****************************************************
      iRC = generateChain2(m_chainSizes[chainId],
                           valuesOf1stPosition,
                           proposalCovMatrix,
                           mahalanobisMatrix,
                           applyMahalanobisInvert,
                           chainName,
                           prefixName,
                           ofs);
      UQ_FATAL_RC_MACRO(iRC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains2()",
                        "improper generateChain2() return");
    }

    //****************************************************
    // Write chain2 out
    //****************************************************
    if (m_chainWrite && ofs) {
      char tmpChainId[10];
      sprintf(tmpChainId,"%d",chainId);
      const std::string& name = m_prefix + tmpChainId + "_chain";
      m_chain2.write(name,*ofs);
    }

    //****************************************************
    // Compute statistics on the chain
    //****************************************************
    if (m_chainComputeStatistics) {
      computeStatistics(m_chain1,
                        m_chain2,
                        chainName,
                        ofs);
    }

    //****************************************************
    // Eventually:
    // --> generate an unique chain
    // --> write it
    // --> compute its statistics
    //****************************************************
    if (m_uniqueChainGenerate) {
      //m_uniqueChain2.erasePositions(m_uniqueChain2Pos,m_uniqueChain2.sequenceSize());

      //computeStatistics(m_uniqueChain2);
    }

    //****************************************************
    // Eventually:
    // --> compute an average chain
    // --> write it
    // --> compute its statistics
    //****************************************************
    if (m_avgChainCompute.size() > chainSumId) {
      // Update m_chain2Sum

      // Check if it is time to compute an average
      if ((chainId+1) == m_avgChainCompute[chainSumId]) {
        // Compute the average

        // Write the computed average
        
        // Compute statistics on the average

        // Prepare for eventual next chain average
        chainSumId++;
      }
    }

    //****************************************************
    // Close file      
    //****************************************************
    if (ofs) {
      // Close file
      ofs->close();

      if (m_env.rank() == 0) {
        std::cout << "Closed output file for chain loop id = " << chainId
                  << std::endl;
      }
    }
    if (m_env.rank() == 0) {
      std::cout << std::endl;
    }
  }

  //if (m_env.rank() == 0) std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains2()"
  //                                 << std::endl;

  return;
}

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::generateWhiteNoise2(
  unsigned int       chainSize,
  const std::string& chainName)
{
  if (m_env.rank() == 0) {
    std::cout << "Generating white noise for chain " << chainName
              << ", with "                           << chainSize
              << " positions..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  double tmpRunTime;

  tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  m_chain2.resizeSequence(chainSize); 
  //double resizeTime = uqMiscGetEllapsedSeconds(&timevalTmp);
  //if (m_env.rank() == 0) {
  //  std::cout << "Chain2 resize took " << resizeTime
  //            << " seconds"
  //            << std::endl;
  //}

#if 1
  V meanVec  (m_paramSpace.zeroVector());
  V stdDevVec(m_paramSpace.zeroVector());
  meanVec.cwSet(0.);
  stdDevVec.cwSet(1.);
  m_chain2.setGaussian(m_env.rng(),meanVec,stdDevVec);
#else
  V gaussianVector(m_paramSpace.zeroVector());
  for (unsigned int positionId = 0; positionId < m_chain2.sequenceSize(); ++positionId) {
    gaussianVector.cwSetGaussian(m_env.rng(),0.,1.);
    m_chain2.setPositionValues(positionId,gaussianVector);

    if ((m_chainDisplayPeriod                     > 0) && 
        (((positionId+1) % m_chainDisplayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << positionId+1
                  << " positions"
                  << std::endl;
      }
    }
  }
#endif
  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain2 generation took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return iRC;
}

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2(
  unsigned int       chainSize,
  const V&           valuesOf1stPosition,
  const M*           proposalCovMatrix,
  const M*           mahalanobisMatrix,
  bool               applyMahalanobisInvert,
  const std::string& chainName,
  const std::string& prefixName,
  std::ofstream*     passedOfs)
{
  if (m_env.rank() == 0) {
    std::cout << "Generating chain " << chainName
              << ", with "           << chainSize
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
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()",
                      "paramInitials should not be out of bound");
  if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalPrior, NULL);
  double m2lPrior            = m_m2lPriorProbDensity_Obj.minus2LnDensity(valuesOf1stPosition);
  if (m_chainMeasureRunTimes) m_priorRunTime += uqMiscGetEllapsedSeconds(&timevalPrior);
  V     m2lLikelihoodVector (m_observableSpace.zeroVector());
  V     misfitVector        (m_observableSpace.zeroVector());
  V     misfitVarianceVector(m_observableSpace.priorVariances());

  if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalLH, NULL);
  if (m_likelihoodObjComputesMisfits) {
    m_m2lLikelihoodFunction_Obj.computeMisfits(valuesOf1stPosition, misfitVector);
    m2lLikelihoodVector = misfitVector/misfitVarianceVector;
  }
  else {
    m_m2lLikelihoodFunction_Obj.computeMinus2LnLikelihoods(valuesOf1stPosition, m2lLikelihoodVector);
  }
  if (m_chainMeasureRunTimes) m_lhRunTime += uqMiscGetEllapsedSeconds(&timevalLH);
  double m2lLikelihoodScalar  = m2lLikelihoodVector.sumOfComponents();
  double logPosterior = -0.5 * ( m2lPrior + m2lLikelihoodScalar );
  uqChainPositionClass<V> currentPosition(m_env,
                                          valuesOf1stPosition,
                                          outOfBounds,
                                          m2lPrior,
                                          misfitVector,
                                          misfitVarianceVector,
                                          m2lLikelihoodVector,
                                          logPosterior);

  V gaussianVector(m_paramSpace.zeroVector());
  V tmpParamValues(m_paramSpace.zeroVector());
  uqChainPositionClass<V> currentCandidate(m_env);

  //****************************************************
  // Begin chain2 loop from positionId = 1
  //****************************************************
  m_chain2.resizeSequence(chainSize); 
  if (m_uniqueChainGenerate) m_idsOfUniquePositions.resize(chainSize,0);
  if (m_chainGenerateExtra) {
    if (m_likelihoodObjComputesMisfits) {
      m_misfitChain.resize        (chainSize,NULL); 
      m_misfitVarianceChain.resize(chainSize,NULL); 
    }
    m_m2lLikelihoodChain.resize(chainSize,NULL); 
    m_alphaQuotients.resize    (chainSize,0.);
  }

  unsigned int uniquePos = 0;
  m_chain2.setPositionValues(0,currentPosition.paramValues());
  if (m_uniqueChainGenerate) m_idsOfUniquePositions[uniquePos++] = 0;
  if (m_chainGenerateExtra) {
    if (m_likelihoodObjComputesMisfits) {
      m_misfitChain        [0] = m_observableSpace.newVector(misfitVector);
      m_misfitVarianceChain[0] = m_observableSpace.newVector(misfitVarianceVector);
    }
    m_m2lLikelihoodChain[0] = m_observableSpace.newVector(currentPosition.m2lLikelihoodVector());
    m_alphaQuotients    [0] = 1.;
  }

  for (unsigned int positionId = 1; positionId < m_chain2.sequenceSize(); ++positionId) {
    //if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
    //                                 << ": beginning chain2 position of id = " << positionId
    //                                 << ", m_maxNumberOfExtraStages =  "       << m_maxNumberOfExtraStages
    //                                 << std::endl;
    unsigned int stageId = 0;

    //****************************************************
    // Loop: generate new parameters
    //****************************************************
    if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
    gaussianVector.cwSetGaussian(m_env.rng(),0.,1.);
    tmpParamValues = currentPosition.paramValues() + *(m_lowerCholProposalCovMatrices[stageId]) * gaussianVector;
    if (m_chainMeasureRunTimes) m_candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

    outOfBounds     = m_paramSpace.outOfBounds(tmpParamValues);
    if (outOfBounds) {
      m_numOutOfBounds++;
      m2lPrior      = 0.;
      m2lLikelihoodVector.cwSet(INFINITY);
      logPosterior  = -INFINITY;
    }
    else {
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalPrior, NULL);
      m2lPrior      = m_m2lPriorProbDensity_Obj.minus2LnDensity(tmpParamValues);
      if (m_chainMeasureRunTimes) m_priorRunTime += uqMiscGetEllapsedSeconds(&timevalPrior);
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalLH, NULL);
      if (m_likelihoodObjComputesMisfits) {
        m_m2lLikelihoodFunction_Obj.computeMisfits(tmpParamValues, misfitVector);
        m2lLikelihoodVector = misfitVector/misfitVarianceVector;
      }
      else {
        m_m2lLikelihoodFunction_Obj.computeMinus2LnLikelihoods(tmpParamValues, m2lLikelihoodVector);
      }
      if (m_chainMeasureRunTimes) m_lhRunTime += uqMiscGetEllapsedSeconds(&timevalLH);
      m2lLikelihoodScalar = m2lLikelihoodVector.sumOfComponents();
      logPosterior  = -0.5 * ( m2lPrior + m2lLikelihoodScalar );
    }
    currentCandidate.set(tmpParamValues,
                         outOfBounds,
                         m2lPrior,
                         misfitVector,
                         misfitVarianceVector,
                         m2lLikelihoodVector,
                         logPosterior);

    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      std::cout << "\n-----------------------------------------------------------\n"
                << std::endl;
    }
    bool accept = false;
    if (outOfBounds) {
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
      if (m_chainMeasureRunTimes) m_mhAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalMhAlpha);
      if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
        std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
                  << ": for chain2 position of id = " << positionId
                  << ", alpha = " << alpha
                  << std::endl;
      }
      accept = acceptAlpha(alpha);
    }
    if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
      if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
                                       << ": for chain2 position of id = " << positionId
                                       << " contents of currentCandidate.paramValues() are:"
                                       << std::endl;
      std::cout << currentCandidate.paramValues();
      if (m_env.rank() == 0) std::cout << std::endl;

      if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
                                       << ": for chain2 position of id = " << positionId
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
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalDR, NULL);

      drPositions[0] = new uqChainPositionClass<V>(currentPosition);
      drPositions[1] = new uqChainPositionClass<V>(currentCandidate);

      while ((accept == false) && (stageId < m_maxNumberOfExtraStages)) {
        stageId++;

        if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
        gaussianVector.cwSetGaussian(m_env.rng(),0.,1.);
        tmpParamValues = currentPosition.paramValues() + *(m_lowerCholProposalCovMatrices[stageId]) * gaussianVector;
        if (m_chainMeasureRunTimes) m_candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

        outOfBounds   = m_paramSpace.outOfBounds(tmpParamValues);
        if (outOfBounds) {
          m2lPrior      = 0.;
          m2lLikelihoodVector.cwSet(INFINITY);
          logPosterior  = -INFINITY;
        }
        else {
          if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalPrior, NULL);
          m2lPrior      = m_m2lPriorProbDensity_Obj.minus2LnDensity(tmpParamValues);
          if (m_chainMeasureRunTimes) m_priorRunTime += uqMiscGetEllapsedSeconds(&timevalPrior);
          if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalLH, NULL);
          if (m_likelihoodObjComputesMisfits) {
            m_m2lLikelihoodFunction_Obj.computeMisfits(tmpParamValues, misfitVector);
            m2lLikelihoodVector = misfitVector/misfitVarianceVector;
          }
          else {
            m_m2lLikelihoodFunction_Obj.computeMinus2LnLikelihoods(tmpParamValues, m2lLikelihoodVector);
          }
          if (m_chainMeasureRunTimes) m_lhRunTime += uqMiscGetEllapsedSeconds(&timevalLH);
          m2lLikelihoodScalar = m2lLikelihoodVector.sumOfComponents();
          logPosterior  = -0.5 * ( m2lPrior + m2lLikelihoodScalar );
        }
        currentCandidate.set(tmpParamValues,
                             outOfBounds,
                             m2lPrior,
                             misfitVector,
                             misfitVarianceVector,
                             m2lLikelihoodVector,
                             logPosterior);

        drPositions.push_back(new uqChainPositionClass<V>(currentCandidate));
        if (outOfBounds == false) {
          if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalDrAlpha, NULL);
          double alpha = this->alpha(drPositions);
          if (m_chainMeasureRunTimes) m_drAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalDrAlpha);
#if 0 // For debug only
          if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
                                           << ": for chain2 position of id = " << positionId
                                           << " and stageId = " << stageId
                                           << ", alpha = " << alpha
                                           << std::endl;
#endif
          accept = acceptAlpha(alpha);
        }
      } // while

      if (m_chainMeasureRunTimes) m_drRunTime += uqMiscGetEllapsedSeconds(&timevalDR);
    } // end of 'delayed rejection' logic

    for (unsigned int i = 0; i < drPositions.size(); ++i) {
      if (drPositions[i]) delete drPositions[i];
    }

    //****************************************************
    // Loop: update chain
    //****************************************************
    if (accept) {
      m_chain2.setPositionValues(positionId,currentCandidate.paramValues());
      if (m_uniqueChainGenerate) m_idsOfUniquePositions[uniquePos++] = positionId;
      if (m_chainGenerateExtra) {
        if (m_likelihoodObjComputesMisfits) {
          m_misfitChain[positionId] = m_observableSpace.newVector(currentCandidate.misfitVector());
          //m_misfitVarianceChain[positionId] is updated below, after the update of 'misfitVarianceVector'
        }
        m_m2lLikelihoodChain[positionId] = m_observableSpace.newVector(currentCandidate.m2lLikelihoodVector());
      }
      currentPosition = currentCandidate;
    }
    else {
      m_chain2.setPositionValues(positionId,currentPosition.paramValues());
      if (m_chainGenerateExtra) {
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
          double term1 = 0.5*( varAccuracies[i] + numbersOfObs[i]                );
          double term2 =  2./( varAccuracies[i] * priorVars[i] + misfitVector[i] );
          misfitVarianceVector[i] = 1./uqMiscGammar(term1,term2,m_env.rng());
          //if (m_env.rank() == 0) {
          //  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
          //            << ": for chain2 position of id = "    << positionId
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
        //  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()"
        //            << ": for chain2 position of id = "       << positionId
        //            << ", misfitVarianceVector changed from " << *(m_misfitVarianceChain[positionId])
        //            << " to "                                 << misfitVarianceVector
        //            << std::endl;
        //}
      }
      if (m_chainGenerateExtra) {
        m_misfitVarianceChain[positionId] = m_observableSpace.newVector(misfitVarianceVector);
      }
    }

    //****************************************************
    // Loop: adaptive Metropolis (adaptation of covariance matrix)
    //****************************************************
    if ((m_initialNonAdaptInterval > 0) &&
        (m_adaptInterval           > 0)) {
      if (m_chainMeasureRunTimes) iRC = gettimeofday(&timevalAM, NULL);

      // Now might be the moment to adapt
      unsigned int idOfFirstPositionInSubChain = 0;
      uqArrayOfSequencesClass<V> subChain2(0,m_paramSpace.zeroVector());

      // Check if now is indeed the moment to adapt
      if (positionId < m_initialNonAdaptInterval) {
        // Do nothing
      }
      else if (positionId == m_initialNonAdaptInterval) {
        idOfFirstPositionInSubChain = 0;
        subChain2.resizeSequence(m_initialNonAdaptInterval+1);
        m_lastMean             = m_paramSpace.newVector();
        m_lastAdaptedCovMatrix = m_paramSpace.newMatrix();
      }
      else {
        unsigned int interval = positionId - m_initialNonAdaptInterval;
        if ((interval % m_adaptInterval) == 0) {
          idOfFirstPositionInSubChain = positionId - m_adaptInterval;
          subChain2.resizeSequence(m_adaptInterval);
        }
      }

      // If now is indeed the moment to adapt, then do it!
      if (subChain2.sequenceSize() > 0) {
        V transporterVec(m_paramSpace.zeroVector());
        for (unsigned int i = 0; i < subChain2.sequenceSize(); ++i) {
          m_chain2.getPositionValues(idOfFirstPositionInSubChain+i,transporterVec);
          subChain2.setPositionValues(i,transporterVec);
        }
        updateCovMatrix2(subChain2,
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
                              "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()",
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
                                "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()",
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
                            "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()",
                            "need to code the update of m_upperCholProposalPrecMatrices");
#endif

          if (m_maxNumberOfExtraStages > 0) updateCovMatrices();
        }
      }

      if (m_chainMeasureRunTimes) m_amRunTime += uqMiscGetEllapsedSeconds(&timevalAM);
    } // End of 'adaptive Metropolis' logic

    if ((m_chainDisplayPeriod                     > 0) && 
        (((positionId+1) % m_chainDisplayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << positionId+1
                  << " positions"
                  << std::endl;
      }
    }
  } // end chain2 loop

  //****************************************************
  // Print basic information about the chain
  //****************************************************
  m_chainRunTime += uqMiscGetEllapsedSeconds(&timevalChain);
  if (m_env.rank() == 0) {
    std::cout << "Finished generating the chain " << chainName
              << ", with "                        << m_chain2.sequenceSize()
              << " positions.";
    if (m_uniqueChainGenerate) {
      std::cout << " and " << uniquePos
                << " 'unique' positions (i.e., not counting repetitions due to rejections)";
    }
    std::cout << "\nSome information about this chain:"
              << "\n  Chain2 run time       = " << m_chainRunTime
              << " seconds";
    if (m_chainMeasureRunTimes) {
      std::cout << "\n\n Breaking of the chain2 run time:\n";
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
    std::cout << "\n  Rejection percentage = " << 100. * (double) m_numRejections/(double) m_chain2.sequenceSize()
              << " %";
    std::cout << "\n   Outbound percentage = " << 100. * (double) m_numOutOfBounds/(double) m_chain2.sequenceSize()
              << " %";
    std::cout << std::endl;
  }

  if (m_chainWrite && passedOfs) {
    iRC = writeInfo(m_chain1,
                    m_chain2,
                    chainName,
                    prefixName,
                    *passedOfs,
                    mahalanobisMatrix,
                    applyMahalanobisInvert);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain2()",
                      "improper writeInfo() return");
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************

  return iRC;
}
#if 0
template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics2(
  const uqArrayOfSequencesClass<V>& chain2,
  std::ofstream* passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  double tmpRunTime;

  // Set initial positions for the computation of chain statistics
  std::vector<unsigned int> initialPosForStatistics(m_finalPercentsForStats.size(),0);
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    initialPosForStatistics[i] = (unsigned int) ((1. - m_finalPercentsForStats[i]*0.01) * (double) chain2.sequenceSize());
  }
  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics2(): initial positions for statistics =";
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
  chain2.mean(0,
              chain2.sequenceSize(),
              chainMean);

  V chainSampleVariance(m_paramSpace.zeroVector());
  chain2.sampleVariance(0,
                        chain2.sequenceSize(),
                        chainMean,
                        chainSampleVariance);

  V chainPopulationVariance(m_paramSpace.zeroVector());
  chain2.populationVariance(0,
                            chain2.sequenceSize(),
                            chainMean,
                            chainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "\nChain2 statistics took " << tmpRunTime
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
    std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics2(): lengths for batchs in BMM =";
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
    chain2.bmm(initialPosForStatistics,
               m_lengthsForBMM,
               _2dArrayOfBMM);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        std::cout << "\nEstimated covariances of sample mean, through batch means method, for subchain beggining at position " << initialPosForStatistics[initialPosId]
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
  // Compute power spectral density (PSD) of chain2 at zero frequency
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
    chain2.psd(initialPosForStatistics,
               m_numBlocksForPSD,
               m_hopSizeRatioForPSD,
               _2dArrayOfPSDAtZero);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        double sizeForPSD = chain2.sequenceSize() - initialPosForStatistics[initialPosId];
        std::cout << "\nEstimated covariances of sample mean, through psd (fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
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
    chain2.geweke(initialPosForStatistics,
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
      std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::computeStatistics2(): lags for autocorrelation =";
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
    chain2.autoCorrelations(initialPosForStatistics,
                            lagsForCorrs,
                            _2dArrayOfAutoCorrs);

    V estimatedVarianceOfSampleMean(m_paramSpace.zeroVector());
    estimatedVarianceOfSampleMean.cwSet(1.);
    for (unsigned int i = 0; i < 1; ++i) {
      for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
        double ratio = ((double) j)/((double) m_chain2.sequenceSize());
        estimatedVarianceOfSampleMean += 2.*(1.-ratio)*_2dArrayOfAutoCorrs(i,j);
      }
    }
    estimatedVarianceOfSampleMean *= chainSampleVariance;
    estimatedVarianceOfSampleMean /= (double) m_chain2.sequenceSize();
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
      std::cout << "Chain2 autocorrelation took " << tmpRunTime
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
  if (initialPosForUncorrelation >= chain2.sequenceSize()) initialPosForUncorrelation = chain2.sequenceSize() - 1;
  chain2.minMax(initialPosForUncorrelation,
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
    chain2.histogram(initialPosForUncorrelation,
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

    V iqrs(m_paramSpace.zeroVector());
    chain2.interQuantileRange(initialPosForUncorrelation,
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

    chain2.scalesForKDE(initialPosForUncorrelation,
                        m_spacingForUncorrelation,
                        iqrs,
                        *m_scalesForKDEs);

    m_densityValuesFromGaussianKDE.resize(m_numberOfEvaluationPosForKDEs,NULL);
    chain2.gaussianKDE(initialPosForUncorrelation,
                       m_spacingForUncorrelation,
                       m_evaluationPositionsForKDEs,
                       *m_scalesForKDEs,
                       m_densityValuesFromGaussianKDE);
  }

  return;
}
#endif
template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix2(
  const uqArrayOfSequencesClass<V>& subChain2,
  unsigned int                      idOfFirstPositionInSubChain,
  double&                           lastChainSize,
  V&                                lastMean,
  M&                                lastAdaptedCovMatrix)
{
  double doubleSubChainSize = (double) subChain2.sequenceSize();
  if (lastChainSize == 0) {
    UQ_FATAL_TEST_MACRO(subChain2.sequenceSize() < 2,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix2()",
                        "'subChain2.sequenceSize()' should be >= 2");

    subChain2.mean(0,subChain2.sequenceSize(),lastMean);

    V tmpVec(m_paramSpace.zeroVector());
    lastAdaptedCovMatrix = -doubleSubChainSize * matrixProduct(lastMean,lastMean);
    for (unsigned int i = 0; i < subChain2.sequenceSize(); ++i) {
      subChain2.getPositionValues(i,tmpVec);
      lastAdaptedCovMatrix += matrixProduct(tmpVec,tmpVec);
    }
    lastAdaptedCovMatrix /= (doubleSubChainSize - 1.); // That is why subChain2 size must be >= 2
  }
  else {
    UQ_FATAL_TEST_MACRO(subChain2.sequenceSize() < 1,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix2()",
                        "'subChain2.sequenceSize()' should be >= 1");

    UQ_FATAL_TEST_MACRO(idOfFirstPositionInSubChain < 1,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix2()",
                        "'idOfFirstPositionInSubChain' should be >= 1");

    V tmpVec (m_paramSpace.zeroVector());
    V diffVec(m_paramSpace.zeroVector());
    for (unsigned int i = 0; i < subChain2.sequenceSize(); ++i) {
      double doubleCurrentId  = (double) (idOfFirstPositionInSubChain+i);
      subChain2.getPositionValues(i,tmpVec);
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
#if 0
template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::writeChain2(
  std::ofstream& ofs,
  const M*       mahalanobisMatrix,
  bool           applyMahalanobisInvert) const
{
  int iRC = UQ_OK_RC;

  // Write m_chain2
  V tmpVec(m_paramSpace.zeroVector());
  ofs << "queso_" << m_prefix << "chain = [";
  for (unsigned int i = 0; i < m_chain2.sequenceSize(); ++i) {
    m_chain2.getPositionValues(i,tmpVec);
    ofs << tmpVec
        << std::endl;
  }
  ofs << "];\n";

  if (m_chainGenerateExtra) {
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
    V diffVec(m_paramSpace.zeroVector());
    V vec0   (m_paramSpace.zeroVector());
    m_chain2.getPositionValues(0,vec0);
    ofs << "queso_" << m_prefix << "d = [";
    if (applyMahalanobisInvert) {
      for (unsigned int i = 0; i < m_chain2.sequenceSize(); ++i) {
        m_chain2.getPositionValues(i,tmpVec);
        diffVec = tmpVec - vec0;
        ofs << scalarProduct(diffVec, mahalanobisMatrix->invertMultiply(diffVec))
            << std::endl;
      }
    }
    else {
      for (unsigned int i = 0; i < m_chain2.sequenceSize(); ++i) {
        m_chain2.getPositionValues(i,tmpVec);
        diffVec = tmpVec - vec0;
        ofs << scalarProduct(diffVec, *mahalanobisMatrix * diffVec)
            << std::endl;
      }
    }
    ofs << "];\n";
  }

  // Write prior mean values
  ofs << "queso_" << m_prefix << "priorMeanValues = ["
      << m_paramSpace.priorMuValues()
      << "];\n";

  // Write prior sigma values
  ofs << "queso_" << m_prefix << "priorSigmaValues = ["
      << m_paramSpace.priorSigmaValues()
      << "];\n";

  // Write param lower bounds
  ofs << "queso_" << m_prefix << "minValues = ["
      << m_paramSpace.minValues()
      << "];\n";

  // Write param upper bounds
  ofs << "queso_" << m_prefix << "maxValues = ["
      << m_paramSpace.maxValues()
      << "];\n";

  // Write number of rejections
  ofs << "queso_" << m_prefix << "rejected = "  << (double) m_numRejections/(double) (m_chain2.sequenceSize()-1)
      << ";\n"
      << std::endl;

  // Write chain2 run time
  ofs << "queso_" << m_prefix << "runTime = "  << m_chainRunTime
      << ";\n"
      << std::endl;

  // Write number of outbounds
  ofs << "queso_" << m_prefix << "outbounds = " << (double) m_numOutOfBounds/(double) m_chain2.sequenceSize()
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
#endif
#endif // __UQ_DRAM_MCG2_H__
