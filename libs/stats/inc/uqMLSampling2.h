/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyrightg (C) 2008 The PECOS Development Team
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

#ifndef __UQ_MULTI_LEVEL_SAMPLING2_H__
#define __UQ_MULTI_LEVEL_SAMPLING2_H__

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step8(
  const uqSequenceOfVectorsClass<P_V,P_M>& prevChain,                         // input
  unsigned int                             indexOfFirstWeight,                // input
  unsigned int                             indexOfLastWeight,                 // input
  const std::vector<double>&               unifiedWeightStdVectorAtProc0Only, // input
  const uqScalarSequenceClass<double>&     weightSequence,                    // input
  double                                   prevEta,                           // input
  const uqGenericVectorRVClass<P_V,P_M>&   currRv,                            // input
  uqMLSamplingLevelOptionsClass*           currOptions,                       // input (changed temporarily internally)
  P_M&                                     unifiedCovMatrix,                  // input/output
  double&                                  currEta)                           // output
{
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 8 of 9"
                                << std::endl;
      }

      double beforeEta           = prevEta;
      double beforeRejectionRate = 0.;               // To be updated
      bool   beforeRejectionRateIsBelowRange = true; // To be updated

      double nowEta           = prevEta;
      double nowRejectionRate = 0.;               // To be computed
      bool   nowRejectionRateIsBelowRange = true; // To be computed

      std::vector<double> etas(2,0.);
      etas[0] = beforeEta;
      etas[1] = 1.;

      std::vector<double> rejs(2,0.);
      rejs[0] = 0.; // To be computed
      rejs[1] = 0.; // To be computed

      unsigned int nowAttempt = 0;
      bool testResult = false;
      double meanRejectionRate = .5*(currOptions->m_minRejectionRate + currOptions->m_maxRejectionRate);
      bool useMiddlePointLogicForEta = false;
      P_M nowCovMatrix(unifiedCovMatrix);
#if 0 // KAUST, to check
      std::vector<double> unifiedWeightStdVectorAtProc0Only(0);
      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);
#endif
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ": entering loop for assessing rejection rate"
                                  << ", with nowAttempt = "  << nowAttempt
                                  << ", nowRejectionRate = " << nowRejectionRate
                                  << std::endl;
        }
        nowCovMatrix = unifiedCovMatrix;

        if (nowRejectionRate < currOptions->m_minRejectionRate) {
          nowRejectionRateIsBelowRange = true;
        }
        else if (nowRejectionRate > currOptions->m_maxRejectionRate) {
          nowRejectionRateIsBelowRange = false;
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                              "nowRejectionRate should be out of the requested range at this point of the logic");
        }

        if (m_env.inter0Rank() >= 0) { // KAUST
          if (nowAttempt > 0) {
            if (useMiddlePointLogicForEta == false) {
              if (nowAttempt == 1) {
                // Ok, keep useMiddlePointLogicForEta = false
              }
              else if ((beforeRejectionRateIsBelowRange == true) &&
                       (nowRejectionRateIsBelowRange    == true)) {
                // Ok
              }
              else if ((beforeRejectionRateIsBelowRange == false) &&
                       (nowRejectionRateIsBelowRange    == false)) {
                // Ok
              }
              else if ((beforeRejectionRateIsBelowRange == true ) &&
                       (nowRejectionRateIsBelowRange    == false)) {
                useMiddlePointLogicForEta = true;

                // This is the first time the middle point logic will be used below
                etas[0] = std::min(beforeEta,nowEta);
                etas[1] = std::max(beforeEta,nowEta);

                if (etas[0] == beforeEta) {
                  rejs[0] = beforeRejectionRate;
                  rejs[1] = nowRejectionRate;
                }
                else {
                  rejs[0] = nowRejectionRate;
                  rejs[1] = beforeRejectionRate;
                }
              }
              else if ((beforeRejectionRateIsBelowRange == false) &&
                       (nowRejectionRateIsBelowRange    == true )) {
                useMiddlePointLogicForEta = true;

                // This is the first time the middle point logic will be used below
                etas[0] = std::min(beforeEta,nowEta);
                etas[1] = std::max(beforeEta,nowEta);
              }
              else {
                UQ_FATAL_TEST_MACRO(true,
                                    m_env.fullRank(),
                                    "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                                    "before and now range flags are inconsistent");
              }
            } // if (useMiddlePointLogicForEta == false)

            beforeEta                       = nowEta;
            beforeRejectionRate             = nowRejectionRate;
            beforeRejectionRateIsBelowRange = nowRejectionRateIsBelowRange;
            if (useMiddlePointLogicForEta == false) {
              if (beforeRejectionRateIsBelowRange) nowEta *= 4.;
              else                                 nowEta /= 4.;
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                        << ", level " << m_currLevel+LEVEL_REF_ID
                                        << ", step "  << m_currStep
                                        << ": in loop for assessing rejection rate"
                                        << ", with nowAttempt = "  << nowAttempt
                                        << ", useMiddlePointLogicForEta = false"
                                        << ", nowEta just updated to value (to be tested) " << nowEta
                                        << std::endl;
              }
            }
            else {
              if (nowRejectionRate > meanRejectionRate) {
                if (rejs[0] > meanRejectionRate) {
                  etas[0] = nowEta;
                  etas[1] = etas[1];
                }
                else {
                  etas[0] = etas[0];
                  etas[1] = nowEta;
                }
              }
              else {
                if (rejs[0] < meanRejectionRate) {
                  etas[0] = nowEta;
                  etas[1] = etas[1];
                }
                else {
                  etas[0] = etas[0];
                  etas[1] = nowEta;
                }
              }
              nowEta = .5*(etas[0] + etas[1]);
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                        << ", level " << m_currLevel+LEVEL_REF_ID
                                        << ", step "  << m_currStep
                                        << ": in loop for assessing rejection rate"
                                        << ", with nowAttempt = " << nowAttempt
                                        << ", useMiddlePointLogicForEta = true"
                                        << ", nowEta just updated to value (to be tested) " << nowEta
                                        << ", etas[0] = " << etas[0]
                                        << ", etas[1] = " << etas[1]
                                        << std::endl;
              }
            }
          } // if (nowAttempt > 0)
        } // if (m_env.inter0Rank() >= 0) // KAUST

        nowCovMatrix *= nowEta;

        unsigned int originalSubNumSamples = 1 + (unsigned int) ( (1.-meanRejectionRate)/meanRejectionRate/currOptions->m_covRejectionRate/currOptions->m_covRejectionRate );

        if (m_env.inter0Rank() >= 0) { // KAUST
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                    << ", level " << m_currLevel+LEVEL_REF_ID
                                    << ", step "  << m_currStep
                                    << ": in loop for assessing rejection rate"
                                    << ", about to sample "     << originalSubNumSamples << " indexes"
                                    << ", meanRejectionRate = " << meanRejectionRate
                                    << ", covRejectionRate = "  << currOptions->m_covRejectionRate
                                    << std::endl;
          }
        } // KAUST

        std::vector<unsigned int> nowUnifiedIndexCountersAtProc0Only(0); // It will be resized by 'sampleIndexesAtProc0()' below
        if (m_env.inter0Rank() >= 0) { // KAUST
          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          sampleIndexesAtProc0(tmpUnifiedNumSamples,                // input
                               unifiedWeightStdVectorAtProc0Only,   // input
                               nowUnifiedIndexCountersAtProc0Only); // output

          unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
          if (m_env.inter0Rank() == 0) {
            UQ_FATAL_TEST_MACRO(nowUnifiedIndexCountersAtProc0Only.size() != auxUnifiedSize,
                                m_env.fullRank(),
                                "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                                "wrong output from sampleIndexesAtProc0() in step 8");
          }

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                    << ", level " << m_currLevel+LEVEL_REF_ID
                                    << ", step "  << m_currStep
                                    << ": in loop for assessing rejection rate"
                                    << ", about to distribute sampled assessment indexes"
                                    << std::endl;
          }
        } // KAUST

        uqBalancedLinkedChainsPerNodeStruct<P_V> nowBalLinkControl;
        uqUnbalancedLinkedChainsPerNodeStruct    nowUnbLinkControl; // KAUST

        unsigned int Np = 0;
        if (m_env.inter0Rank() == 0) { // Yes, '== 0'
          Np = (unsigned int) m_env.inter0Comm().NumProc();
        }
        std::vector<uqExchangeInfoStruct> exchangeStdVec(0);

        // All processors should call this routine in order to have the same decision value
        bool useBalancedChains = decideOnBalancedChains(currOptions,                        // input
                                                        indexOfFirstWeight,                 // input
                                                        indexOfLastWeight,                  // input
                                                        nowUnifiedIndexCountersAtProc0Only, // input
                                                        exchangeStdVec);                    // output

        if (m_env.inter0Rank() >= 0) { // KAUST
          if (useBalancedChains) {
            prepareBalLinkedChains(currOptions,                        // input
                                   prevChain,                          // input
                                   exchangeStdVec,                     // input/output
                                   nowBalLinkControl);                 // output
          }
          else {
            prepareUnbLinkedChains(indexOfFirstWeight,                 // input
                                   indexOfLastWeight,                  // input
                                   nowUnifiedIndexCountersAtProc0Only, // input
                                   nowUnbLinkControl);                 // output
          }
        } // KAUST

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ": in loop for assessing rejection rate"
                                  << ", about to generate assessment chain"
                                  << std::endl;
        }

        uqSequenceOfVectorsClass<P_V,P_M> nowChain(m_vectorSpace,
                                                   0,
                                                   m_options.m_prefix+"now_chain");
        double       nowRunTime    = 0.;
        unsigned int nowRejections = 0;

        // KAUST: all nodes should call here
        bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
        unsigned int savedRawChainSize          = currOptions->m_rawChainSize; // Ok to use rawChainSize
        bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
        bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;
        unsigned int savedDrMaxNumExtraStages   = currOptions->m_drMaxNumExtraStages;
        unsigned int savedAmAdaptInterval       = currOptions->m_amAdaptInterval;

        currOptions->m_totallyMute           = true;
        currOptions->m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
        currOptions->m_rawChainComputeStats  = false;
        currOptions->m_filteredChainGenerate = false;
        currOptions->m_drMaxNumExtraStages   = 0;
        currOptions->m_amAdaptInterval       = 0;

        // KAUST: all nodes in 'subComm' should call here, important
        if (useBalancedChains) {
          generateBalLinkedChains(*currOptions,       // input, only m_rawChainSize changes
                                  nowCovMatrix,       // input
                                  currRv,             // input
                                  nowBalLinkControl,  // input // Round Rock
                                  nowChain,           // output 
                                  nowRunTime,         // output
                                  nowRejections,      // output
                                  NULL,               // output
                                  NULL);              // output
        }
        else {
          generateUnbLinkedChains(*currOptions,       // input, only m_rawChainSize changes
                                  nowCovMatrix,       // input
                                  currRv,             // input
                                  nowUnbLinkControl,  // input // Round Rock
                                  indexOfFirstWeight, // input // Round Rock
                                  prevChain,          // input // Round Rock
                                  nowChain,           // output 
                                  nowRunTime,         // output
                                  nowRejections,      // output
                                  NULL,               // output
                                  NULL);              // output
        }

        // KAUST: all nodes should call here
        currOptions->m_totallyMute           = savedTotallyMute;
        currOptions->m_rawChainSize          = savedRawChainSize;
        currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
        currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME
        currOptions->m_drMaxNumExtraStages   = savedDrMaxNumExtraStages;
        currOptions->m_amAdaptInterval       = savedAmAdaptInterval;

        for (unsigned int i = 0; i < nowBalLinkControl.balLinkedChains.size(); ++i) {
          UQ_FATAL_TEST_MACRO(nowBalLinkControl.balLinkedChains[i].initialPosition == NULL,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                              "Initial position pointer in step 8 should not be NULL");
          delete nowBalLinkControl.balLinkedChains[i].initialPosition;
          nowBalLinkControl.balLinkedChains[i].initialPosition = NULL;
        }
        nowBalLinkControl.balLinkedChains.clear();

        if (m_env.inter0Rank() >= 0) { // KAUST
          // If only one cov matrix is used, then the rejection should be assessed among all inter0Comm nodes // KAUST3
          unsigned int nowUnifiedRejections = 0;
          int mpiRC = MPI_Allreduce((void *) &nowRejections, (void *) &nowUnifiedRejections, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                              "failed MPI_Allreduce() for now rejections");

#if 0 // Round Rock 2009 12 29
          unsigned int tmpUnifiedNumSamples = 0;
          mpiRC = MPI_Allreduce((void *) &tmpSubNumSamples, (void *) &tmpUnifiedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                              "failed MPI_Allreduce() for num samples in step 8");
#endif

          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          nowRejectionRate = ((double) nowUnifiedRejections) / ((double) tmpUnifiedNumSamples);

          //bool aux1 = (nowRejectionRate == meanRejectionRate);
          bool aux2 = (nowRejectionRate >= currOptions->m_minRejectionRate)
                      &&
                      (nowRejectionRate <= currOptions->m_maxRejectionRate);
          testResult = aux2;

          // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
          uqMiscCheckForSameValueInAllNodes(testResult,
                                            0.,
                                            m_env.inter0Comm(),
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8(), step 8, testResult");
        } // if (m_env.inter0Rank() >= 0) { // KAUST

        // KAUST: all nodes in 'subComm' should have the same 'testResult'
        unsigned int tmpUint = (unsigned int) testResult;
        int mpiRC = MPI_Bcast((void *) &tmpUint, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8()",
                            "failed MPI_Bcast() for testResult");
        testResult = (bool) tmpUint;

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                  << ", level "              << m_currLevel+LEVEL_REF_ID
                                  << ", step "               << m_currStep
                                  << ": in loop for assessing rejection rate"
                                  << ", nowAttempt = "       << nowAttempt
                                  << ", beforeEta = "        << beforeEta
                                  << ", etas[0] = "          << etas[0]
                                  << ", nowEta = "           << nowEta
                                  << ", etas[1] = "          << etas[1]
                                  << ", minRejectionRate = " << currOptions->m_minRejectionRate
                                  << ", nowRejectionRate = " << nowRejectionRate
                                  << ", maxRejectionRate = " << currOptions->m_maxRejectionRate
                                  << std::endl;
        }
        nowAttempt++;

        if (m_env.inter0Rank() >= 0) { // KAUST
          // Make sure all nodes in 'inter0Comm' have the same value of 'nowEta'
          uqMiscCheckForSameValueInAllNodes(nowEta,
                                            1.e-16,
                                            m_env.inter0Comm(),
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence_Step8(), step 8, testResult");
        }
      } while (testResult == false);
      currEta = nowEta;
      if (currEta != 1.) {
        unifiedCovMatrix *= currEta;
      }

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step8()"
                                << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                << ", step "                                   << m_currStep
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", currEta = "                              << currEta
                                << ", assessed rejection rate = "              << nowRejectionRate
                                << std::endl;
      }

  return;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING2_H__


