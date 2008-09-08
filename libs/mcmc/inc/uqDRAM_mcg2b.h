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
      computeStatistics(m_chain2,
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
    iRC = writeInfo(m_chain2,
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

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix2(
  const uqChainBaseClass<V>& subChain2,
  unsigned int               idOfFirstPositionInSubChain,
  double&                    lastChainSize,
  V&                         lastMean,
  M&                         lastAdaptedCovMatrix)
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
#endif // __UQ_DRAM_MCG2_H__
