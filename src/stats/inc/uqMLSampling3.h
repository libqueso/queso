//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id:$
//
//--------------------------------------------------------------------------

#ifndef __UQ_MULTI_LEVEL_SAMPLING3_H__
#define __UQ_MULTI_LEVEL_SAMPLING3_H__

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::sampleIndexes_proc0(
  unsigned int               unifiedRequestedNumSamples,        // input
  const std::vector<double>& unifiedWeightStdVectorAtProc0Only, // input
  std::vector<unsigned int>& unifiedIndexCountersAtProc0Only)   // output
{
  if (m_env.inter0Rank() != 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::sampleIndexes_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": unifiedRequestedNumSamples = "               << unifiedRequestedNumSamples
                            << ", unifiedWeightStdVectorAtProc0Only.size() = " << unifiedWeightStdVectorAtProc0Only.size()
                            << std::endl;
  }

#if 0 // For debug only
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::sampleIndexes_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ":"
                            << std::endl;
    unsigned int numZeros = 0;
    for (unsigned int i = 0; i < unifiedWeightStdVectorAtProc0Only.size(); ++i) {
      *m_env.subDisplayFile() << "  unifiedWeightStdVectorAtProc0Only[" << i
                              << "] = " << unifiedWeightStdVectorAtProc0Only[i]
                              << std::endl;
      if (unifiedWeightStdVectorAtProc0Only[i] == 0.) numZeros++;
    }
    *m_env.subDisplayFile() << "Number of zeros in unifiedWeightStdVectorAtProc0Only = " << numZeros
                            << std::endl;
  }
#endif
 
  if (m_env.inter0Rank() == 0) {
    unsigned int resizeSize = unifiedWeightStdVectorAtProc0Only.size();
    unifiedIndexCountersAtProc0Only.resize(resizeSize,0);

    // Generate 'unifiedRequestedNumSamples' samples from 'tmpFD'
    uqFiniteDistributionClass tmpFd(m_env,
                                    "",
                                    unifiedWeightStdVectorAtProc0Only);
    for (unsigned int i = 0; i < unifiedRequestedNumSamples; ++i) {
      unsigned int index = tmpFd.sample();
      unifiedIndexCountersAtProc0Only[index] += 1;
    }
  }

  return;
}

template <class P_V,class P_M>
bool
uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all(
  const uqMLSamplingLevelOptionsClass* currOptions,                     // input
  unsigned int                         indexOfFirstWeight,              // input
  unsigned int                         indexOfLastWeight,               // input
  const std::vector<unsigned int>&     unifiedIndexCountersAtProc0Only, // input
  std::vector<uqExchangeInfoStruct>&   exchangeStdVec)                  // output
{
  bool result = false;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": indexOfFirstWeight = " << indexOfFirstWeight
                            << ", indexOfLastWeight = "  << indexOfLastWeight
                            << std::endl;
  }

  if ((currOptions->m_loadBalanceAlgorithmId > 0) &&
      (m_env.numSubEnvironments()            > 1)) { // Cannot use 'm_env.inter0Comm().NumProc()' because it might happen that not all nodes at this point of the code belong to 'inter0Comm'
    unsigned int Np = 0;
    if (m_env.inter0Rank() >= 0) { // Yes, '>= 0'
      Np = (unsigned int) m_env.inter0Comm().NumProc();
    }
    std::vector<unsigned int> allFirstIndexes(Np,0); // '0' is already the correct value for recvcnts[0]
    std::vector<unsigned int> allLastIndexes(Np,0);  // '0' is NOT the correct value for recvcnts[0]

    if (m_env.inter0Rank() >= 0) { // Yes, '>= 0'
      //////////////////////////////////////////////////////////////////////////
      // Gather information at proc 0: number of chains and positions per node
      //////////////////////////////////////////////////////////////////////////
      //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
      //                void *recvbuf, int recvcount, MPI_Datatype recvtype, 
      //                int root, MPI_Comm comm )
      unsigned int auxUInt = indexOfFirstWeight;
      int mpiRC = MPI_Gather((void *) &auxUInt, 1, MPI_UNSIGNED, (void *) &allFirstIndexes[0], (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                          "failed MPI_Gather() for first indexes");

      if (m_env.inter0Rank() == 0) {
        UQ_FATAL_TEST_MACRO(allFirstIndexes[0] != indexOfFirstWeight,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                            "failed MPI_Gather() result for first indexes, at proc 0");
      }

      auxUInt = indexOfLastWeight;
      mpiRC = MPI_Gather((void *) &auxUInt, 1, MPI_UNSIGNED, (void *) &allLastIndexes[0], (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                          "failed MPI_Gather() for last indexes");

      if (m_env.inter0Rank() == 0) { // Yes, '== 0'
        //allLastIndexes[0] = indexOfLastWeight; // FIX ME: really necessary????
        UQ_FATAL_TEST_MACRO(allLastIndexes[0] != indexOfLastWeight,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                            "failed MPI_Gather() result for last indexes, at proc 0");
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Proc 0 prepares information to decide if load balancing is needed
    //////////////////////////////////////////////////////////////////////////
    if (m_env.inter0Rank() == 0) {
      *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ": original distribution of unified indexes in 'inter0Comm' is as follows"
                              << std::endl;
      for (unsigned int r = 0; r < Np; ++r) {
        *m_env.subDisplayFile() << "  allFirstIndexes[" << r << "] = " << allFirstIndexes[r]
                                << "  allLastIndexes["  << r << "] = " << allLastIndexes[r]
                                << std::endl;
      }
      for (unsigned int r = 0; r < (Np-1); ++r) { // Yes, '-1'
        UQ_FATAL_TEST_MACRO(allFirstIndexes[r+1] != (allLastIndexes[r]+1),
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                            "wrong indexes");
      }

      for (unsigned int r = 0; r < (Np-1); ++r) { // Yes, '-1'
        UQ_FATAL_TEST_MACRO(allFirstIndexes[r+1] != (allLastIndexes[r]+1),
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                            "wrong indexes");
      }

      std::vector<unsigned int> origNumChainsPerNode   (Np,0);
      std::vector<unsigned int> origNumPositionsPerNode(Np,0);
      int r = 0;
      for (unsigned int i = 0; i < unifiedIndexCountersAtProc0Only.size(); ++i) {
        if ((allFirstIndexes[r] <= i) && // FIX ME: not a robust logic
            (i <= allLastIndexes[r] )) {
          // Ok
        }
        else {
          r++;
          if ((r < (int) Np           ) &&
              (allFirstIndexes[r] <= i) && 
              (i <= allLastIndexes[r] )) {
            // Ok
          }
          else {
	    std::cerr << "In uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()"
                      << ": i = " << i
                      << ", r = " << r
                      << ", allFirstIndexes[r] = " << allFirstIndexes[r]
                      << ", allLastIndexes[r] = "  << allLastIndexes[r]
                      << std::endl;
            UQ_FATAL_TEST_MACRO(true,
                                m_env.worldRank(),
                                "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                                "wrong indexes or 'r' got too large");
          }
        }
        if (unifiedIndexCountersAtProc0Only[i] != 0) {
          origNumChainsPerNode   [r] += 1;
          origNumPositionsPerNode[r] += unifiedIndexCountersAtProc0Only[i];

          uqExchangeInfoStruct auxInfo;
          auxInfo.originalNodeOfInitialPosition  = r;
          auxInfo.originalIndexOfInitialPosition = i - allFirstIndexes[r];
          auxInfo.finalNodeOfInitialPosition     = -1; // Yes, '-1' for now, important
          auxInfo.numberOfPositions              = unifiedIndexCountersAtProc0Only[i];
          exchangeStdVec.push_back(auxInfo);
        }
        // FIX ME: swap trick to save memory
      }

      // Check if number of procs is too large
      unsigned int totalNumberOfChains = 0;
      for (unsigned int r = 0; r < Np; ++r) {
        totalNumberOfChains += origNumChainsPerNode[r];
      }
      *m_env.subDisplayFile() << "  KEY"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", Np = "  << Np
                              << ", totalNumberOfChains = " << totalNumberOfChains
                              << std::endl;

      // Check if ratio max/min justifies optimization
      unsigned int origMinPosPerNode  = *std::min_element(origNumPositionsPerNode.begin(), origNumPositionsPerNode.end());
      unsigned int origMaxPosPerNode  = *std::max_element(origNumPositionsPerNode.begin(), origNumPositionsPerNode.end());
      for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
        *m_env.subDisplayFile() << "  KEY"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", origNumChainsPerNode["     << nodeId << "] = " << origNumChainsPerNode[nodeId]
                                << ", origNumPositionsPerNode["  << nodeId << "] = " << origNumPositionsPerNode[nodeId]
                                << std::endl;
      }
      double origRatioOfPosPerNode = ((double) origMaxPosPerNode ) / ((double) origMinPosPerNode);
      *m_env.subDisplayFile() << "  KEY"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", origRatioOfPosPerNode = "      << origRatioOfPosPerNode
                              << ", option loadBalanceTreshold = " << currOptions->m_loadBalanceTreshold
                              << std::endl;

      // At this point, only proc 0 is running...
      // Set boolean 'result' for good
      if ((                   Np < totalNumberOfChains               ) &&
          (origRatioOfPosPerNode > currOptions->m_loadBalanceTreshold)) {
        result = true;
      }
    }
  }

  m_env.fullComm().Barrier();
  unsigned int tmpValue = result;
  int mpiRC = MPI_Bcast((void *) &tmpValue, (int) 1, MPI_UNSIGNED, 0, m_env.fullComm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()",
                      "failed MPI_Bcast() for 'result'");
  if (m_env.inter0Rank() != 0) result = tmpValue;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::decideOnBalancedChains_all()"
                            << ", level "    << m_currLevel+LEVEL_REF_ID
                            << ", step "     << m_currStep
                            << ": result = " << result
                            << std::endl;
  }

  return result;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0( // EXTRA FOR LOAD BALANCE
  const uqMLSamplingLevelOptionsClass*      currOptions,         // input
  const uqSequenceOfVectorsClass<P_V,P_M>&  prevChain,           // input
  std::vector<uqExchangeInfoStruct>&        exchangeStdVec,      // input/output
  uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl) // output
{
  if (m_env.inter0Rank() < 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << std::endl;
  }

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  if (m_env.inter0Rank() == 0) {
    switch (currOptions->m_loadBalanceAlgorithmId) {
      case 2:
        justBalance_proc0(currOptions,     // input
                          exchangeStdVec); // input/output
      break;

      case 1:
      default:
#ifdef QUESO_HAS_GLPK
        // Get final node responsible for a linked chain by solving BIP at node zero only
        solveBIP_proc0(exchangeStdVec); // input/output
#else
#endif
      break;
    }
  } // if (m_env.inter0Rank() == 0)

  m_env.inter0Comm().Barrier();

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 now broadcasts the information on 'exchangeStdVec'
  //////////////////////////////////////////////////////////////////////////
  unsigned int exchangeStdVecSize = exchangeStdVec.size();
  int mpiRC = MPI_Bcast((void *) &exchangeStdVecSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed MPI_Bcast() for exchangeStdVec size");
  if (m_env.inter0Rank() > 0) exchangeStdVec.resize(exchangeStdVecSize);

  mpiRC = MPI_Bcast((void *) &exchangeStdVec[0], (int) (exchangeStdVecSize*sizeof(uqExchangeInfoStruct)), MPI_CHAR, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed MPI_Bcast() for exchangeStdVec data");

  //////////////////////////////////////////////////////////////////////////
  // All "management" nodes update 'finalNumChainsPerNode' and 'finalNumPostionsPerNode'
  //////////////////////////////////////////////////////////////////////////
  std::vector<unsigned int> finalNumChainsPerNode   (Np,0);
  std::vector<unsigned int> finalNumPositionsPerNode(Np,0);
  unsigned int Nc = exchangeStdVec.size();
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    unsigned int nodeId = exchangeStdVec[chainId].finalNodeOfInitialPosition;
    finalNumChainsPerNode   [nodeId] += 1;
    finalNumPositionsPerNode[nodeId] += exchangeStdVec[chainId].numberOfPositions;
  }

  //////////////////////////////////////////////////////////////////////////
  // Sanity check
  //////////////////////////////////////////////////////////////////////////
  unsigned int finalMinPosPerNode = *std::min_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  unsigned int finalMaxPosPerNode = *std::max_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  double finalRatioOfPosPerNode = ((double) finalMaxPosPerNode) / ((double)finalMinPosPerNode);
  //std::cout << m_env.worldRank() << ", finalRatioOfPosPerNode = " << finalRatioOfPosPerNode << std::endl;

  std::vector<double> auxBuf(1,0.);
  double minRatio = 0.;
  auxBuf[0] = finalRatioOfPosPerNode;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minRatio, (int) auxBuf.size(), MPI_DOUBLE, MPI_MIN, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed MPI_Allreduce() for min");
  //std::cout << m_env.worldRank() << ", minRatio = " << minRatio << std::endl;
  UQ_FATAL_TEST_MACRO(minRatio != finalRatioOfPosPerNode,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed minRatio sanity check");

  double maxRatio = 0.;
  auxBuf[0] = finalRatioOfPosPerNode;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxRatio, (int) auxBuf.size(), MPI_DOUBLE, MPI_MAX, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed MPI_Allreduce() for max");
  //std::cout << m_env.worldRank() << ", maxRatio = " << maxRatio << std::endl;
  UQ_FATAL_TEST_MACRO(maxRatio != finalRatioOfPosPerNode,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed maxRatio sanity check");

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 now broadcasts the information on 'finalNumChainsPerNode'
  //////////////////////////////////////////////////////////////////////////
  unsigned int finalNumChainsPerNodeSize = finalNumChainsPerNode.size();
  mpiRC = MPI_Bcast((void *) &finalNumChainsPerNodeSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed MPI_Bcast() for finalNumChainsPerNode size");
  if (m_env.inter0Rank() > 0) finalNumChainsPerNode.resize(finalNumChainsPerNodeSize);

  mpiRC = MPI_Bcast((void *) &finalNumChainsPerNode[0], (int) finalNumChainsPerNodeSize, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm()); // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareBalLinkedChains_inter0()",
                      "failed MPI_Bcast() for finalNumChainsPerNode data");

  //////////////////////////////////////////////////////////////////////////
  // Mpi exchange information between nodes and properly populate
  // balancedLinkControl.linkedChains at each node
  //////////////////////////////////////////////////////////////////////////
  mpiExchangePositions_inter0(prevChain,
                              exchangeStdVec,
                              finalNumChainsPerNode,
                              finalNumPositionsPerNode, // It is already valid at all "management" nodes (not only at node 0) because of the sanity check above
                              balancedLinkControl);

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0(
  unsigned int                           indexOfFirstWeight,              // input
  unsigned int                           indexOfLastWeight,               // input
  const std::vector<unsigned int>&       unifiedIndexCountersAtProc0Only, // input
  uqUnbalancedLinkedChainsPerNodeStruct& unbalancedLinkControl)           // output
{
  if (m_env.inter0Rank() < 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": indexOfFirstWeight = " << indexOfFirstWeight
                            << ", indexOfLastWeight = "  << indexOfLastWeight
                            << std::endl;
  }

  unsigned int              subNumSamples = 0;
  std::vector<unsigned int> unifiedIndexCountersAtAllProcs(0);

  // All nodes in 'inter0Comm' should resize to the same size // KAUST3
  unsigned int resizeSize = unifiedIndexCountersAtProc0Only.size();
  int mpiRC = MPI_Bcast((void *) &resizeSize, (int) 1, MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "failed MPI_Bcast() for resizeSize");
  unifiedIndexCountersAtAllProcs.resize(resizeSize,0);

  if (m_env.inter0Rank() == 0) unifiedIndexCountersAtAllProcs = unifiedIndexCountersAtProc0Only;

  // Broadcast index counters to all nodes
  mpiRC = MPI_Bcast((void *) &unifiedIndexCountersAtAllProcs[0], (int) unifiedIndexCountersAtAllProcs.size(), MPI_UNSIGNED, 0, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "failed MPI_Bcast() for unified index counters");
#if 0 // Use allgatherv ??? for subNumSamples instead
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ":"
                            << std::endl;
    for (int r = 0; r < m_env.inter0Comm().NumProc(); ++r) {
      *m_env.subDisplayFile() << "  unifiedIndexCountersAtAllProcs[" << r << "] = " << unifiedIndexCountersAtAllProcs[r]
                              << std::endl;
    }
  }
#endif
  //for (unsigned int i = 0; i < unifiedIndexCountersAtAllProcs.size(); ++i) {
  //  *m_env.subDisplayFile() << "unifiedIndexCountersAtAllProcs[" << i
  //                          << "] = " << unifiedIndexCountersAtAllProcs[i]
  //                          << std::endl;
  //}

  // Use 'indexOfFirstWeight' and 'indexOfLastWeight' in order to update 'subNumSamples'
  UQ_FATAL_TEST_MACRO(indexOfFirstWeight >= unifiedIndexCountersAtAllProcs.size(),
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "invalid indexOfFirstWeight");
  UQ_FATAL_TEST_MACRO(indexOfLastWeight >= unifiedIndexCountersAtAllProcs.size(),
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "invalid indexOfLastWeight");
  subNumSamples = 0;
  for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
    subNumSamples += unifiedIndexCountersAtAllProcs[i];
  }

  std::vector<unsigned int> auxBuf(1,0);

  unsigned int minModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "failed MPI_Allreduce() for min");

  unsigned int maxModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "failed MPI_Allreduce() for max");

  unsigned int sumModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumModifiedSubNumSamples, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "failed MPI_Allreduce() for sum");

  //UQ_FATAL_TEST_MACRO(unifiedRequestedNumSamples != sumModifiedSubNumSamples,
  //                    m_env.worldRank(),
  //                    "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
  //                    "invalid state");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                                   << m_currLevel+LEVEL_REF_ID
                            << ", step "                                    << m_currStep
                            << ": subNumSamples = "                         << subNumSamples
                            << ", unifiedIndexCountersAtAllProcs.size() = " << unifiedIndexCountersAtAllProcs.size()
                            << std::endl;
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                      << m_currLevel+LEVEL_REF_ID
                            << ", step "                       << m_currStep
                            << ": minModifiedSubNumSamples = " << minModifiedSubNumSamples
                            << ", avgModifiedSubNumSamples = " << ((double) sumModifiedSubNumSamples)/((double) m_env.inter0Comm().NumProc())
                            << ", maxModifiedSubNumSamples = " << maxModifiedSubNumSamples
                            << std::endl;
  }

  unsigned int numberOfPositionsToGuaranteeForNode = subNumSamples;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                                 << m_currLevel+LEVEL_REF_ID
                            << ", step "                                  << m_currStep
                            << ": numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                            << std::endl;
  }
  for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
//for (unsigned int i = 0; i < unifiedIndexCountersAtAllProcs.size(); ++i) { // KAUST4: important
    while (unifiedIndexCountersAtAllProcs[i] != 0) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 30)) {
        *m_env.subDisplayFile() << ", numberOfPositionsToGuaranteeForNode = " << numberOfPositionsToGuaranteeForNode
                                << ", unifiedIndexCountersAtAllProcs["        << i
                                << "] = "                                     << unifiedIndexCountersAtAllProcs[i]
                                << std::endl;
      }
      if (unifiedIndexCountersAtAllProcs[i] < numberOfPositionsToGuaranteeForNode) {
        uqUnbalancedLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = unifiedIndexCountersAtAllProcs[i];
        unbalancedLinkControl.unbLinkedChains.push_back(auxControl);

        numberOfPositionsToGuaranteeForNode -= unifiedIndexCountersAtAllProcs[i];
        unifiedIndexCountersAtAllProcs[i] = 0;
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] > 0                                   )) {
      //else { // KAUST4
        uqUnbalancedLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = numberOfPositionsToGuaranteeForNode;
        unbalancedLinkControl.unbLinkedChains.push_back(auxControl);

        unifiedIndexCountersAtAllProcs[i] -= numberOfPositionsToGuaranteeForNode;
        numberOfPositionsToGuaranteeForNode = 0;
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] == 0                                  )) {
        // Ok
      }
      else {
        UQ_FATAL_TEST_MACRO(true, // KAUST4
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                            "should never get here");
      }
    }
  }
  UQ_FATAL_TEST_MACRO(numberOfPositionsToGuaranteeForNode != 0, // subNumSamples, // KAUST4
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                      "numberOfPositionsToGuaranteeForNode exited loop with wrong value");
  // FIX ME: swap trick to save memory

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving uqMLSamplingClass<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                                          << m_currLevel+LEVEL_REF_ID
                            << ", step "                                           << m_currStep
                            << ": unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all( // EXTRA FOR LOAD BALANCE
  uqMLSamplingLevelOptionsClass&                  inputOptions,            // input, only m_rawChainSize changes
  const P_M&                                      unifiedCovMatrix,        // input
  const uqGenericVectorRVClass  <P_V,P_M>&        rv,                      // input
  const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,     // input // Round Rock
  uqSequenceOfVectorsClass      <P_V,P_M>&        workingChain,            // output
  double&                                         cumulativeRunTime,       // output
  unsigned int&                                   cumulativeRejections,    // output
  uqScalarSequenceClass         <double>*         currLogLikelihoodValues, // output
  uqScalarSequenceClass         <double>*         currLogTargetValues)     // output
{
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()"
                            << ": balancedLinkControl.balLinkedChains.size() = " << balancedLinkControl.balLinkedChains.size()
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = balancedLinkControl.balLinkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  int mpiRC = MPI_Bcast((void *) &chainIdMax, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important // LOAD BALANCE
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
                      "failed MPI_Bcast() for chainIdMax");

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfPositions = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfPositions += balancedLinkControl.balLinkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm()); // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
                        "failed MPI_Allreduce() for min");

    unsigned int maxNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm()); // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
                        "failed MPI_Allreduce() for max");

    unsigned int sumNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm()); // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
                        "failed MPI_Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()"
                              << ", level "               << m_currLevel+LEVEL_REF_ID
                              << ", step "                << m_currStep
                              << ": chainIdMax = "        << chainIdMax
                              << ", numberOfPositions = " << numberOfPositions
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()"
                              << ", level "                  << m_currLevel+LEVEL_REF_ID
                              << ", step "                   << m_currStep
                              << ": minNumberOfPositions = " << minNumberOfPositions
                              << ", avgNumberOfPositions = " << ((double) sumNumberOfPositions)/((double) m_env.inter0Comm().NumProc())
                              << ", maxNumberOfPositions = " << maxNumberOfPositions
                              << std::endl;
    }
  }
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      // aqui 4
      auxInitialPosition = *(balancedLinkControl.balLinkedChains[chainId].initialPosition); // Round Rock
      tmpChainSize = balancedLinkControl.balLinkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99/*2*/)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()"
                                << ", level "          << m_currLevel+LEVEL_REF_ID
                                << ", step "           << m_currStep
                                << ": chainId = "      << chainId
                                << ", tmpChainSize = " << tmpChainSize
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm().Comm()); // Yes, 'subComm', important // KAUST
#if 0 // For debug only
    for (int r = 0; r < m_env.subComm().NumProc(); ++r) {
      if (r == m_env.subComm().MyPID()) {
	std::cout << "Vector 'auxInitialPosition at rank " << r
                  << " has contents "                      << auxInitialPosition
                  << std::endl;
      }
      m_env.subComm().Barrier();
    }
    sleep(1);
#endif

    // KAUST: all nodes in 'subComm' should have the same 'tmpChainSize'
    mpiRC = MPI_Bcast((void *) &tmpChainSize, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important // LOAD BALANCE
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
                        "failed MPI_Bcast() for tmpChainSize");

    inputOptions.m_rawChainSize = tmpChainSize;
    uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    uqScalarSequenceClass<double> tmpLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> tmpLogTargetValues    (m_env,0,"");

    // KAUST: all nodes should call here
    uqMetropolisHastingsSGClass<P_V,P_M> mcSeqGenerator(inputOptions,
                                                        rv,
                                                        auxInitialPosition,
                                                        &unifiedCovMatrix);

    // KAUST: all nodes should call here
    mcSeqGenerator.generateSequence(tmpChain,
                                    &tmpLogLikelihoodValues, // likelihood is IMPORTANT
                                    &tmpLogTargetValues);
    uqMHRawChainInfoStruct mcRawInfo;
    mcSeqGenerator.getRawChainInfo(mcRawInfo);
    cumulativeRunTime    += mcRawInfo.runTime;
    cumulativeRejections += mcRawInfo.numRejections;

    if (m_env.inter0Rank() >= 0) {
      if ((m_env.subDisplayFile()       ) &&
          (m_env.displayVerbosity() >= 0)) { // detailed output debug
        for (unsigned int i = 0; i < tmpLogLikelihoodValues.subSequenceSize(); ++i) {
          *m_env.subDisplayFile() << "tmpLogLikelihoodValues[" << i << "] = " << tmpLogLikelihoodValues[i]
                                  << ", tmpLogTargetValues["   << i << "] = " << tmpLogTargetValues[i]
                                  << std::endl;
        }
      }
        
      if ((m_env.subDisplayFile()             ) &&
          (m_env.displayVerbosity()   >= 0    ) &&
          (inputOptions.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                                << ", level "               << m_currLevel+LEVEL_REF_ID
                                << ", step "                << m_currStep
                                << ", chainId = "           << chainId
                                << ": finished generating " << tmpChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      // KAUST5: what if workingChain ends up with different size in different nodes? Important
      workingChain.append              (tmpChain,              1,tmpChain.subSequenceSize()-1              ); // IMPORTANT: '1' in order to discard initial position
      if (currLogLikelihoodValues) {
        currLogLikelihoodValues->append(tmpLogLikelihoodValues,1,tmpLogLikelihoodValues.subSequenceSize()-1); // IMPORTANT: '1' in order to discard initial position
        if ((m_env.subDisplayFile()        ) &&
            (m_env.displayVerbosity() >= 99) &&
            (chainId == 0                  )) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                                  << ", level "     << m_currLevel+LEVEL_REF_ID
                                  << ", step "      << m_currStep
                                  << ", chainId = " << chainId
                                  << ", tmpLogLikelihoodValues.subSequenceSize() = " << tmpLogLikelihoodValues.subSequenceSize()
                                  << ", tmpLogLikelihoodValues[0] = "                << tmpLogLikelihoodValues[0]
                                  << ", tmpLogLikelihoodValues[1] = "                << tmpLogLikelihoodValues[1]
                                  << ", currLogLikelihoodValues[0] = "               << (*currLogLikelihoodValues)[0]
                                  << std::endl;
        }
      }
      if (currLogTargetValues) {
        currLogTargetValues->append    (tmpLogTargetValues,    1,tmpLogTargetValues.subSequenceSize()-1    ); // IMPORTANT: '1' in order to discard initial position
      }
    }
  } // for

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()"
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all(
  uqMLSamplingLevelOptionsClass&               inputOptions,            // input, only m_rawChainSize changes
  const P_M&                                   unifiedCovMatrix,        // input
  const uqGenericVectorRVClass  <P_V,P_M>&     rv,                      // input
  const uqUnbalancedLinkedChainsPerNodeStruct& unbalancedLinkControl,   // input // Round Rock
  unsigned int                                 indexOfFirstWeight,      // input // Round Rock
  const uqSequenceOfVectorsClass<P_V,P_M>&     prevChain,               // input // Round Rock
  uqSequenceOfVectorsClass      <P_V,P_M>&     workingChain,            // output
  double&                                      cumulativeRunTime,       // output
  unsigned int&                                cumulativeRejections,    // output
  uqScalarSequenceClass         <double>*      currLogLikelihoodValues, // output
  uqScalarSequenceClass         <double>*      currLogTargetValues)     // output
{
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()"
                            << ": unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                            << ", indexOfFirstWeight = "                           << indexOfFirstWeight
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = unbalancedLinkControl.unbLinkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  int mpiRC = MPI_Bcast((void *) &chainIdMax, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()",
                      "failed MPI_Bcast() for chainIdMax");

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfPositions = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfPositions += unbalancedLinkControl.unbLinkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &minNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MIN, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()",
                        "failed MPI_Allreduce() for min");

    unsigned int maxNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &maxNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_MAX, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()",
                        "failed MPI_Allreduce() for max");

    unsigned int sumNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    mpiRC = MPI_Allreduce((void *) &auxBuf[0], (void *) &sumNumberOfPositions, (int) auxBuf.size(), MPI_UNSIGNED, MPI_SUM, m_env.inter0Comm().Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()",
                        "failed MPI_Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()"
                              << ", level "               << m_currLevel+LEVEL_REF_ID
                              << ", step "                << m_currStep
                              << ": chainIdMax = "        << chainIdMax
                              << ", numberOfPositions = " << numberOfPositions
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()"
                              << ", level "                  << m_currLevel+LEVEL_REF_ID
                              << ", step "                   << m_currStep
                              << ": minNumberOfPositions = " << minNumberOfPositions
                              << ", avgNumberOfPositions = " << ((double) sumNumberOfPositions)/((double) m_env.inter0Comm().NumProc())
                              << ", maxNumberOfPositions = " << maxNumberOfPositions
                              << std::endl;
    }
  }
  if ((m_debugExponent == 1.) && 
      (m_currStep      == 10)) {
    //m_env.setExceptionalCircunstance(true);
  }
  unsigned int cumulativeNumPositions = 0;
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      unsigned int auxIndex = unbalancedLinkControl.unbLinkedChains[chainId].initialPositionIndexInPreviousChain - indexOfFirstWeight; // KAUST4 // Round Rock
      prevChain.getPositionValues(auxIndex,auxInitialPosition); // Round Rock
      tmpChainSize = unbalancedLinkControl.unbLinkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99/*2*/)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()"
                                << ", level "          << m_currLevel+LEVEL_REF_ID
                                << ", step "           << m_currStep
                                << ": chainId = "      << chainId
                                << ", tmpChainSize = " << tmpChainSize
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm().Comm()); // Yes, 'subComm', important // KAUST
#if 0 // For debug only
    for (int r = 0; r < m_env.subComm().NumProc(); ++r) {
      if (r == m_env.subComm().MyPID()) {
	std::cout << "Vector 'auxInitialPosition at rank " << r
                  << " has contents "                      << auxInitialPosition
                  << std::endl;
      }
      m_env.subComm().Barrier();
    }
    sleep(1);
#endif

    // KAUST: all nodes in 'subComm' should have the same 'tmpChainSize'
    mpiRC = MPI_Bcast((void *) &tmpChainSize, (int) 1, MPI_UNSIGNED, 0, m_env.subComm().Comm()); // Yes, 'subComm', important
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()",
                        "failed MPI_Bcast() for tmpChainSize");

    inputOptions.m_rawChainSize = tmpChainSize;
    uqSequenceOfVectorsClass<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    uqScalarSequenceClass<double> tmpLogLikelihoodValues(m_env,0,"");
    uqScalarSequenceClass<double> tmpLogTargetValues    (m_env,0,"");

    // KAUST: all nodes should call here
    uqMetropolisHastingsSGClass<P_V,P_M> mcSeqGenerator(inputOptions,
                                                        rv,
                                                        auxInitialPosition,
                                                        &unifiedCovMatrix);

    // KAUST: all nodes should call here
    mcSeqGenerator.generateSequence(tmpChain,
                                    &tmpLogLikelihoodValues, // likelihood is IMPORTANT
                                    &tmpLogTargetValues);
    uqMHRawChainInfoStruct mcRawInfo;
    mcSeqGenerator.getRawChainInfo(mcRawInfo);
    cumulativeRunTime    += mcRawInfo.runTime;
    cumulativeRejections += mcRawInfo.numRejections;

    if (m_env.inter0Rank() >= 0) {
      if (m_env.exceptionalCircunstance()) {
        if ((m_env.subDisplayFile()       ) &&
            (m_env.displayVerbosity() >= 0)) { // detailed output debug
          P_V tmpVec(m_vectorSpace.zeroVector());
          for (unsigned int i = 0; i < tmpLogLikelihoodValues.subSequenceSize(); ++i) {
            tmpChain.getPositionValues(i,tmpVec);
            *m_env.subDisplayFile() << "DEBUG finalChain[" << cumulativeNumPositions+i << "] "
                                    << "= tmpChain["               << i << "] = " << tmpVec 
                                    << ", tmpLogLikelihoodValues[" << i << "] = " << tmpLogLikelihoodValues[i]
                                    << ", tmpLogTargetValues["     << i << "] = " << tmpLogTargetValues[i]
                                    << std::endl;
          }
        }
      } // exceptional
    
      cumulativeNumPositions += tmpChainSize;
      if (cumulativeNumPositions > 100) m_env.setExceptionalCircunstance(false);
        
      if ((m_env.subDisplayFile()             ) &&
          (m_env.displayVerbosity()   >= 0    ) &&
          (inputOptions.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                                << ", level "               << m_currLevel+LEVEL_REF_ID
                                << ", step "                << m_currStep
                                << ", chainId = "           << chainId
                                << ": finished generating " << tmpChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;
      }

      // KAUST5: what if workingChain ends up with different size in different nodes? Important
      workingChain.append              (tmpChain,              1,tmpChain.subSequenceSize()-1              ); // IMPORTANT: '1' in order to discard initial position
      if (currLogLikelihoodValues) {
        currLogLikelihoodValues->append(tmpLogLikelihoodValues,1,tmpLogLikelihoodValues.subSequenceSize()-1); // IMPORTANT: '1' in order to discard initial position
        if ((m_env.subDisplayFile()        ) &&
            (m_env.displayVerbosity() >= 99) &&
            (chainId == 0                  )) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                                  << ", level "     << m_currLevel+LEVEL_REF_ID
                                  << ", step "      << m_currStep
                                  << ", chainId = " << chainId
                                  << ", tmpLogLikelihoodValues.subSequenceSize() = " << tmpLogLikelihoodValues.subSequenceSize()
                                  << ", tmpLogLikelihoodValues[0] = "                << tmpLogLikelihoodValues[0]
                                  << ", tmpLogLikelihoodValues[1] = "                << tmpLogLikelihoodValues[1]
                                  << ", currLogLikelihoodValues[0] = "               << (*currLogLikelihoodValues)[0]
                                  << std::endl;
        }
      }
      if (currLogTargetValues) {
        currLogTargetValues->append    (tmpLogTargetValues,    1,tmpLogTargetValues.subSequenceSize()-1    ); // IMPORTANT: '1' in order to discard initial position
      }
    }
  } // for

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()"
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  return;
}

#ifdef QUESO_HAS_GLPK
template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::solveBIP_proc0( // EXTRA FOR LOAD BALANCE 
  std::vector<uqExchangeInfoStruct>& exchangeStdVec) // input/output
{
  if (m_env.inter0Rank() != 0) return;

  int iRC = UQ_OK_RC;
  struct timeval timevalBIP;
  iRC = gettimeofday(&timevalBIP, NULL);

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  unsigned int Nc = exchangeStdVec.size();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": Np = "  << Np
                            << ", Nc = "  << Nc
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Instantiate BIP
  //////////////////////////////////////////////////////////////////////////
  glp_prob *lp; 
  lp = glp_create_prob(); 
  glp_set_prob_name(lp, "sample"); 

  //////////////////////////////////////////////////////////////////////////
  // Set rows and colums of BIP constraint matrix
  //////////////////////////////////////////////////////////////////////////
  unsigned int m = Nc+Np-1;
  unsigned int n = Nc*Np;
  unsigned int ne = Nc*Np + 2*Nc*(Np -1);

  glp_add_rows(lp, m); // Not 'm+1'
  for (int i = 1; i <= (int) Nc; ++i) {
    glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0); 
    glp_set_row_name(lp, i, ""); 
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0); 
    glp_set_row_name(lp, i, ""); 
  }
 
  glp_add_cols(lp, n); // Not 'n+1'
  for (int j = 1; j <= (int) n; ++j) {
    //glp_set_col_kind(lp, j, GLP_BV);
    glp_set_col_kind(lp, j, GLP_IV);
    glp_set_col_bnds(lp, j, GLP_DB, 0.0, 1.0); 
    glp_set_col_name(lp, j, ""); 
  }

  glp_set_obj_dir(lp, GLP_MIN); 
  for (int chainId = 0; chainId <= (int) (Nc-1); ++chainId) {
    glp_set_obj_coef(lp, (chainId*Np)+1, exchangeStdVec[chainId].numberOfPositions); 
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished setting BIP rows and cols"
                            << ", m = "  << m
                            << ", n = "  << n
                            << ", ne = " << ne
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Load constraint matrix
  //////////////////////////////////////////////////////////////////////////
  std::vector<int   > iVec(ne+1,0);
  std::vector<int   > jVec(ne+1,0);
  std::vector<double> aVec(ne+1,0.);
  int coefId = 1; // Yes, '1'
  for (int i = 1; i <= (int) Nc; ++i) {
    for (int j = 1; j <= (int) Np; ++j) {
      iVec[coefId] = i;
      jVec[coefId] = (i-1)*Np + j;
      aVec[coefId] = 1.;
      coefId++;
    }
  }
  for (int i = 1; i <= (int) (Np-1); ++i) {
    for (int j = 1; j <= (int) Nc; ++j) {
      iVec[coefId] = Nc+i;
      jVec[coefId] = (j-1)*Np + i;
      aVec[coefId] = -((double) exchangeStdVec[j-1].numberOfPositions);
      coefId++;

      iVec[coefId] = Nc+i;
      jVec[coefId] = (j-1)*Np + i + 1;
      aVec[coefId] = exchangeStdVec[j-1].numberOfPositions;
      coefId++;
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished setting BIP constraint matrix"
                            << ", ne = "     << ne
                            << ", coefId = " << coefId
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(coefId != (int) (ne+1),
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "invalid final coefId");

  glp_load_matrix(lp, ne, &iVec[0], &jVec[0], &aVec[0]);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished loading BIP constraint matrix"
                            << ", glp_get_num_rows(lp) = " << glp_get_num_rows(lp)
                            << ", glp_get_num_cols(lp) = " << glp_get_num_cols(lp)
                            << ", glp_get_num_nz(lp) = "   << glp_get_num_nz(lp)
                            << ", glp_get_num_int(lp) = "  << glp_get_num_int(lp)
                            << ", glp_get_num_bin(lp) = "  << glp_get_num_bin(lp)
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Check BIP before solving it
  //////////////////////////////////////////////////////////////////////////
  UQ_FATAL_TEST_MACRO(glp_get_num_rows(lp) != (int) m, // Not 'm+1'
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "invalid number of rows");

  UQ_FATAL_TEST_MACRO(glp_get_num_cols(lp) != (int) n, // Not 'n+1'
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "invalid number of columnss");

  UQ_FATAL_TEST_MACRO(glp_get_num_nz(lp) != (int) ne,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "invalid number of nonzero constraint coefficients");

  UQ_FATAL_TEST_MACRO(glp_get_num_int(lp) != (int) n, // ????
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "invalid number of integer structural variables");

  UQ_FATAL_TEST_MACRO(glp_get_num_bin(lp) != (int) n,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "invalid number of binary structural variables");

  //////////////////////////////////////////////////////////////////////////
  // Set initial state
  //////////////////////////////////////////////////////////////////////////
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      int j = chainId*Np + nodeId + 1;
      if (nodeId == 0) {
        glp_set_col_stat(lp, j, GLP_BS);
      }
      else {
        glp_set_col_stat(lp, j, GLP_BS);
      }
    }
  }
#if 0
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      int j = chainId*Np + nodeId + 1;
      int initialState = glp_mip_col_val(lp, j);
      if (nodeId == 0) {
        UQ_FATAL_TEST_MACRO(initialState != 1,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                            "for nodeId = 0, initial state should be '1'");
      }
      else {
        UQ_FATAL_TEST_MACRO(initialState != 0,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                            "for nodeId > 0, initial state should be '0'");
      }
    }
  }
#endif
  for (int i = 1; i <= (int) Nc; ++i) {
    glp_set_row_stat(lp, i, GLP_NS); 
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    glp_set_row_stat(lp, i, GLP_BS); 
  }

  //glp_write_mps(lp, GLP_MPS_DECK, NULL, "nada.fixed_mps");
  //glp_write_mps(lp, GLP_MPS_FILE, NULL, "nada.free_mps" );
  //glp_write_lp (lp, NULL, "nada.cplex");

  //////////////////////////////////////////////////////////////////////////
  // Solve BIP
  //////////////////////////////////////////////////////////////////////////
  BIP_routine_struct BIP_routine_info;
  BIP_routine_info.env = &m_env;
  BIP_routine_info.currLevel = m_currLevel;

  glp_iocp BIP_params;
  glp_init_iocp(&BIP_params);
  BIP_params.presolve = GLP_ON;
  // aqui 2
  //BIP_params.binarize = GLP_ON;
  //BIP_params.cb_func = BIP_routine; 
  //BIP_params.cb_info = (void *) (&BIP_routine_info);
  int BIP_rc = glp_intopt(lp, &BIP_params);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished solving BIP"
                            << ", BIP_rc = " << BIP_rc
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(BIP_rc != 0,
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "BIP returned rc != 0");

  //////////////////////////////////////////////////////////////////////////
  // Check BIP status after solution
  //////////////////////////////////////////////////////////////////////////
  int BIP_Status = glp_mip_status(lp);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": BIP_Status = " << BIP_Status
                            << std::endl;
  }

  switch (BIP_Status) {
    case GLP_OPT:
      // Ok 
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": BIP solution is optimal"
                                << std::endl;
      }
    break;

    case GLP_FEAS:
      // Ok 
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": BIP solution is guaranteed to be 'only' feasible"
                                << std::endl;
      }
    break;

    default:
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                          "BIP has an undefined solution or has no solution");
    break;
  }

  for (int i = 1; i <= (int) Nc; ++i) {
    UQ_FATAL_TEST_MACRO(glp_mip_row_val(lp, i) != 1,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                        "row should have value 1 at solution");
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    UQ_FATAL_TEST_MACRO(glp_mip_row_val(lp, i) > 0,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                        "row should have value 0 or should be negative at solution");
  }

  //////////////////////////////////////////////////////////////////////////
  // Prepare output information, needed to for MPI distribution afterwards
  //////////////////////////////////////////////////////////////////////////
  std::vector<unsigned int> finalNumChainsPerNode   (Np,0);
  std::vector<unsigned int> finalNumPositionsPerNode(Np,0);
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      int j = chainId*Np + nodeId + 1;
      if (glp_mip_col_val(lp, j) == 0) {
        // Do nothing
      }
      else if (glp_mip_col_val(lp, j) == 1) {
        UQ_FATAL_TEST_MACRO(exchangeStdVec[chainId].finalNodeOfInitialPosition != -1, // important
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                            "chain has already been taken care of");
        exchangeStdVec[chainId].finalNodeOfInitialPosition = nodeId;
        finalNumChainsPerNode   [nodeId] += 1;
        finalNumPositionsPerNode[nodeId] += exchangeStdVec[chainId].numberOfPositions;
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                            "control variable should be either '0' or '1'");
      }
    }
  }

  unsigned int finalMinPosPerNode = *std::min_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  unsigned int finalMaxPosPerNode = *std::max_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished preparing output information"
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Printout solution information
  //////////////////////////////////////////////////////////////////////////
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": solution gives the following redistribution"
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY In uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", finalNumChainsPerNode["    << nodeId << "] = " << finalNumChainsPerNode[nodeId]
                              << ", finalNumPositionsPerNode[" << nodeId << "] = " << finalNumPositionsPerNode[nodeId]
                              << std::endl;
    }
    *m_env.subDisplayFile() << "  KEY In uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", finalRatioOfPosPerNode = " << ((double) finalMaxPosPerNode) / ((double)finalMinPosPerNode)
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Make sanity checks
  //////////////////////////////////////////////////////////////////////////
  UQ_FATAL_TEST_MACRO(glp_mip_obj_val(lp) != (double) finalNumPositionsPerNode[0],
                      m_env.worldRank(),
                      "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                      "Invalid objective value");

  for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, '1'
    UQ_FATAL_TEST_MACRO(finalNumPositionsPerNode[nodeId-1] < finalNumPositionsPerNode[nodeId],
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                        "Next node should have a number of positions equal or less than the current node");
  }

  for (int i = (int) (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    unsigned int nodeId = i - Nc;
    int diff = ((int) finalNumPositionsPerNode[nodeId]) - ((int) finalNumPositionsPerNode[nodeId-1]);
    UQ_FATAL_TEST_MACRO(glp_mip_row_val(lp, i) != diff,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::solveBIP_proc0()",
                        "wrong state");
  }
    
  //////////////////////////////////////////////////////////////////////////
  // Free memory and return
  //////////////////////////////////////////////////////////////////////////
  glp_delete_prob(lp); 

  double bipRunTime = uqMiscGetEllapsedSeconds(&timevalBIP);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << bipRunTime << " seconds"
                            << std::endl;
  }

  return;
}
#endif // QUESO_HAS_GLPK

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::justBalance_proc0(
  const uqMLSamplingLevelOptionsClass* currOptions,    // input
  std::vector<uqExchangeInfoStruct>&   exchangeStdVec) // input/output
{
  if (m_env.inter0Rank() != 0) return;

  int iRC = UQ_OK_RC;
  struct timeval timevalBal;
  iRC = gettimeofday(&timevalBal, NULL);

  unsigned int Np = m_env.numSubEnvironments();
  unsigned int Nc = exchangeStdVec.size();

  std::vector<uqExchangeInfoStruct> currExchangeStdVec(Nc);
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    currExchangeStdVec[chainId] = exchangeStdVec[chainId];
    currExchangeStdVec[chainId].finalNodeOfInitialPosition = currExchangeStdVec[chainId].originalNodeOfInitialPosition; // final = original
  }

  //////////////////////////////////////////////////////////////////////////
  // Compute original ratio of positions per node
  //////////////////////////////////////////////////////////////////////////
  unsigned int iterIdMax = 0;
  std::vector<unsigned int> currNumChainsPerNode   (Np,0);
  std::vector<unsigned int> currNumPositionsPerNode(Np,0);
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    unsigned int nodeId = currExchangeStdVec[chainId].finalNodeOfInitialPosition; // Yes, 'final'
    currNumChainsPerNode   [nodeId] += 1;
    currNumPositionsPerNode[nodeId] += currExchangeStdVec[chainId].numberOfPositions;
    iterIdMax                       += currExchangeStdVec[chainId].numberOfPositions;
  }
  unsigned int currMinPosPerNode = *std::min_element(currNumPositionsPerNode.begin(), currNumPositionsPerNode.end());
  unsigned int currMaxPosPerNode = *std::max_element(currNumPositionsPerNode.begin(), currNumPositionsPerNode.end());
  double currRatioOfPosPerNode = ((double) currMaxPosPerNode ) / ((double) currMinPosPerNode);

  //////////////////////////////////////////////////////////////////////////
  // Loop
  //////////////////////////////////////////////////////////////////////////
  //iterIdMax /= 2;
  int iterId = -1;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::justBalance_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", iter "  << iterId
                            << ", currRatioOfPosPerNode = " << currRatioOfPosPerNode
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY In uqMLSamplingClass<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", iter "  << iterId
                              << ", currNumChainsPerNode["    << nodeId << "] = " << currNumChainsPerNode[nodeId]
                              << ", currNumPositionsPerNode[" << nodeId << "] = " << currNumPositionsPerNode[nodeId]
                              << std::endl;
    }
  }

  std::vector<std::vector<double> > vectorOfChainSizesPerNode(Np);
  while ((iterId                < (int) iterIdMax                   ) &&
         (currRatioOfPosPerNode > currOptions->m_loadBalanceTreshold)) {
    iterId++;

    //////////////////////////////////////////////////////////////////////////
    // Initialize information
    //////////////////////////////////////////////////////////////////////////
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      vectorOfChainSizesPerNode[nodeId].clear(); // make sure vectors have size 0
    }
    for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
      unsigned int nodeId = currExchangeStdVec[chainId].finalNodeOfInitialPosition; // Yes, 'final'
      vectorOfChainSizesPerNode[nodeId].push_back(currExchangeStdVec[chainId].numberOfPositions);
    }
    // FIX ME: swap to save memory
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      std::sort(vectorOfChainSizesPerNode[nodeId].begin(), vectorOfChainSizesPerNode[nodeId].end());
      UQ_FATAL_TEST_MACRO(vectorOfChainSizesPerNode[nodeId].size() != currNumChainsPerNode[nodeId],
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::justBalance_proc0()",
                          "inconsistent number of chains in node");
    }
   
    //////////////////////////////////////////////////////////////////////////
    // Find [node with most postions], [node with least positions] and [number of positions to move]
    //////////////////////////////////////////////////////////////////////////
    unsigned int currBiggestAmountOfPositionsPerNode  = currNumPositionsPerNode[0];
    unsigned int currSmallestAmountOfPositionsPerNode = currNumPositionsPerNode[0];
    unsigned int currNodeWithMostPositions = 0;
    unsigned int currNodeWithLeastPositions = 0;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      if (currNumPositionsPerNode[nodeId] > currBiggestAmountOfPositionsPerNode) {
        currBiggestAmountOfPositionsPerNode = currNumPositionsPerNode[nodeId];
        currNodeWithMostPositions = nodeId;
      }
      if (currNumPositionsPerNode[nodeId] < currSmallestAmountOfPositionsPerNode) {
        currSmallestAmountOfPositionsPerNode = currNumPositionsPerNode[nodeId];
        currNodeWithLeastPositions = nodeId;
      }
    }

    UQ_FATAL_TEST_MACRO(currMinPosPerNode != currNumPositionsPerNode[currNodeWithLeastPositions],
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::justBalance_proc0()",
                        "inconsistent currMinPosPerNode");

    UQ_FATAL_TEST_MACRO(currMaxPosPerNode != currNumPositionsPerNode[currNodeWithMostPositions],
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::justBalance_proc0()",
                        "inconsistent currMaxPosPerNode");

    unsigned int numberOfPositionsToMove = vectorOfChainSizesPerNode[currNodeWithMostPositions][0];

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", iter "  << iterId
                              << ", before update"
                              << ", node w/ most pos is "
                              << currNodeWithMostPositions  << "(cs=" << currNumChainsPerNode[currNodeWithMostPositions ] << ", ps=" << currNumPositionsPerNode[currNodeWithMostPositions ] << ")"
                              << ", node w/ least pos is "
                              << currNodeWithLeastPositions << "(cs=" << currNumChainsPerNode[currNodeWithLeastPositions] << ", ps=" << currNumPositionsPerNode[currNodeWithLeastPositions] << ")"
                              << ", number of pos to move = " << numberOfPositionsToMove
                              << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////
    // Update 'final' fields in the two nodes
    //////////////////////////////////////////////////////////////////////////
    std::vector<uqExchangeInfoStruct> newExchangeStdVec(Nc);
    for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
      newExchangeStdVec[chainId] = currExchangeStdVec[chainId];
    }

    for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
      if ((newExchangeStdVec[chainId].finalNodeOfInitialPosition == (int) currNodeWithMostPositions) &&
          (newExchangeStdVec[chainId].numberOfPositions          == numberOfPositionsToMove        )) {
        newExchangeStdVec[chainId].finalNodeOfInitialPosition = currNodeWithLeastPositions;
        break; // exit 'for'
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Compute new ratio of positions per node
    //////////////////////////////////////////////////////////////////////////
    std::vector<unsigned int> newNumChainsPerNode   (Np,0);
    std::vector<unsigned int> newNumPositionsPerNode(Np,0);
    for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
      unsigned int nodeId = newExchangeStdVec[chainId].finalNodeOfInitialPosition; // Yes, 'final'
      newNumChainsPerNode   [nodeId] += 1;
      newNumPositionsPerNode[nodeId] += newExchangeStdVec[chainId].numberOfPositions;
    }

    unsigned int newBiggestAmountOfPositionsPerNode  = newNumPositionsPerNode[0];
    unsigned int newSmallestAmountOfPositionsPerNode = newNumPositionsPerNode[0];
    unsigned int newNodeWithMostPositions = 0;
    unsigned int newNodeWithLeastPositions = 0;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      if (newNumPositionsPerNode[nodeId] > newBiggestAmountOfPositionsPerNode) {
        newBiggestAmountOfPositionsPerNode = newNumPositionsPerNode[nodeId];
        newNodeWithMostPositions = nodeId;
      }
      if (newNumPositionsPerNode[nodeId] < newSmallestAmountOfPositionsPerNode) {
        newSmallestAmountOfPositionsPerNode = newNumPositionsPerNode[nodeId];
        newNodeWithLeastPositions = nodeId;
      }
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", iter "  << iterId
                              << ", after update"
                              << ", node w/ most pos is "
                              << newNodeWithMostPositions  << "(cs=" << newNumChainsPerNode[newNodeWithMostPositions ] << ", ps=" << newNumPositionsPerNode[newNodeWithMostPositions ] << ")"
                              << ", node w/ least pos is "
                              << newNodeWithLeastPositions << "(cs=" << newNumChainsPerNode[newNodeWithLeastPositions] << ", ps=" << newNumPositionsPerNode[newNodeWithLeastPositions] << ")"
                              << std::endl;
    }

    unsigned int newMinPosPerNode = *std::min_element(newNumPositionsPerNode.begin(), newNumPositionsPerNode.end());
    unsigned int newMaxPosPerNode = *std::max_element(newNumPositionsPerNode.begin(), newNumPositionsPerNode.end());
    double newRatioOfPosPerNode = ((double) newMaxPosPerNode ) / ((double) newMinPosPerNode);

    UQ_FATAL_TEST_MACRO(newMinPosPerNode != newNumPositionsPerNode[newNodeWithLeastPositions],
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::justBalance_proc0()",
                        "inconsistent newMinPosPerNode");

    UQ_FATAL_TEST_MACRO(newMaxPosPerNode != newNumPositionsPerNode[newNodeWithMostPositions],
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::justBalance_proc0()",
                        "inconsistent newMaxPosPerNode");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", iter "  << iterId
                              << ", newMaxPosPerNode = "     << newMaxPosPerNode
                              << ", newMinPosPerNode = "     << newMinPosPerNode
                              << ", newRatioOfPosPerNode = " << newRatioOfPosPerNode
                              << std::endl;
      for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
        *m_env.subDisplayFile() << "  KEY In uqMLSamplingClass<P_V,P_M>::justBalance_proc0()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", iter "  << iterId
                                << ", newNumChainsPerNode["    << nodeId << "] = " << newNumChainsPerNode   [nodeId]
                                << ", newNumPositionsPerNode[" << nodeId << "] = " << newNumPositionsPerNode[nodeId]
                                << std::endl;
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // See if we need to exit 'while'
    //////////////////////////////////////////////////////////////////////////
    if (newRatioOfPosPerNode > currRatioOfPosPerNode) {
      break; // exit 'while'
    }

    //////////////////////////////////////////////////////////////////////////
    // Prepare for next iteration
    //////////////////////////////////////////////////////////////////////////
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      currNumChainsPerNode   [nodeId] = 0;
      currNumPositionsPerNode[nodeId] = 0;
    }
    currRatioOfPosPerNode = newRatioOfPosPerNode;
    for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
      currExchangeStdVec[chainId] = newExchangeStdVec[chainId];
      unsigned int nodeId = currExchangeStdVec[chainId].finalNodeOfInitialPosition; // Yes, 'final'
      currNumChainsPerNode   [nodeId] += 1;
      currNumPositionsPerNode[nodeId] += currExchangeStdVec[chainId].numberOfPositions;
    }
    currMinPosPerNode = *std::min_element(currNumPositionsPerNode.begin(), currNumPositionsPerNode.end());
    currMaxPosPerNode = *std::max_element(currNumPositionsPerNode.begin(), currNumPositionsPerNode.end());
    currRatioOfPosPerNode = ((double) currMaxPosPerNode ) / ((double) currMinPosPerNode);
  }

  //////////////////////////////////////////////////////////////////////////
  // Prepare output information
  //////////////////////////////////////////////////////////////////////////
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    exchangeStdVec[chainId].finalNodeOfInitialPosition = currExchangeStdVec[chainId].finalNodeOfInitialPosition; // Yes, 'final' = 'final'
  }

  //////////////////////////////////////////////////////////////////////////
  // Printout solution information
  //////////////////////////////////////////////////////////////////////////
  std::vector<unsigned int> finalNumChainsPerNode   (Np,0);
  std::vector<unsigned int> finalNumPositionsPerNode(Np,0);
  for (unsigned int chainId = 0; chainId < Nc; ++chainId) {
    unsigned int nodeId = exchangeStdVec[chainId].finalNodeOfInitialPosition; // Yes, 'final'
    finalNumChainsPerNode   [nodeId] += 1;
    finalNumPositionsPerNode[nodeId] += exchangeStdVec[chainId].numberOfPositions;
  }
  unsigned int finalMinPosPerNode = *std::min_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  unsigned int finalMaxPosPerNode = *std::max_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  double finalRatioOfPosPerNode = ((double) finalMaxPosPerNode ) / ((double) finalMinPosPerNode);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::justBalance_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": solution gives the following redistribution"
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY In uqMLSamplingClass<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", finalNumChainsPerNode["    << nodeId << "] = " << finalNumChainsPerNode[nodeId]
                              << ", finalNumPositionsPerNode[" << nodeId << "] = " << finalNumPositionsPerNode[nodeId]
                              << std::endl;
    }
    *m_env.subDisplayFile() << "  KEY In uqMLSamplingClass<P_V,P_M>::justBalance_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", finalRatioOfPosPerNode = " << finalRatioOfPosPerNode
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Measure time
  //////////////////////////////////////////////////////////////////////////
  double balRunTime = uqMiscGetEllapsedSeconds(&timevalBal);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::justBalance_proc0()"
                            << ", level "                   << m_currLevel+LEVEL_REF_ID
                            << ", step "                    << m_currStep
                            << ", iterId = "                << iterId
                            << ", currRatioOfPosPerNode = " << currRatioOfPosPerNode
                            << ", after " << balRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0( // EXTRA FOR LOAD BALANCE
  const uqSequenceOfVectorsClass<P_V,P_M>&  prevChain,                // input
  const std::vector<uqExchangeInfoStruct>&  exchangeStdVec,           // input
  const std::vector<unsigned int>&          finalNumChainsPerNode,    // input
  const std::vector<unsigned int>&          finalNumPositionsPerNode, // input
  uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl)      // output
{
  if (m_env.inter0Rank() < 0) return;

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  unsigned int Nc = exchangeStdVec.size();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSampling<P_V,P_M>::mpiExchangePositions_inter0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": Np = " << Np
                            << ", Nc = " << Nc
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Each node performs:
  // --> a 'gatherv' for collecting all necessary initial positions from other nodes
  // --> a 'gatherv' for collecting all necessary chain lenghts from other nodes
  //////////////////////////////////////////////////////////////////////////
  for (unsigned int r = 0; r < Np; ++r) {
    //////////////////////////////////////////////////////////////////////////
    // Prepare some counters
    //////////////////////////////////////////////////////////////////////////
    unsigned int              numberOfInitialPositionsNodeRAlreadyHas = 0;
    std::vector<unsigned int> numberOfInitialPositionsNodeRHasToReceiveFromNode(Np,0);
    std::vector<unsigned int> indexesOfInitialPositionsNodeRHasToReceiveFromMe(0);

    unsigned int              sumOfChainLenghtsNodeRAlreadyHas = 0;
    std::vector<unsigned int> chainLenghtsNodeRHasToInherit(0);

    for (unsigned int i = 0; i < Nc; ++i) {
      if (exchangeStdVec[i].finalNodeOfInitialPosition == (int) r) {
        if (exchangeStdVec[i].originalNodeOfInitialPosition == (int) r) {
          numberOfInitialPositionsNodeRAlreadyHas++;
          sumOfChainLenghtsNodeRAlreadyHas += exchangeStdVec[i].numberOfPositions;
        }
        else {
          numberOfInitialPositionsNodeRHasToReceiveFromNode[exchangeStdVec[i].originalNodeOfInitialPosition]++;
          chainLenghtsNodeRHasToInherit.push_back(exchangeStdVec[i].numberOfPositions);
          if (m_env.inter0Rank() == exchangeStdVec[i].originalNodeOfInitialPosition) {
            indexesOfInitialPositionsNodeRHasToReceiveFromMe.push_back(exchangeStdVec[i].originalIndexOfInitialPosition);
          }
        }
      }
    }

    unsigned int totalNumberOfInitialPositionsNodeRHasToReceive = 0;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      totalNumberOfInitialPositionsNodeRHasToReceive += numberOfInitialPositionsNodeRHasToReceiveFromNode[nodeId];
    }

    unsigned int totalNumberOfChainLenghtsNodeRHasToInherit = chainLenghtsNodeRHasToInherit.size();
    unsigned int totalSumOfChainLenghtsNodeRHasToInherit = 0;
    for (unsigned int i = 0; i < totalNumberOfChainLenghtsNodeRHasToInherit; ++i) {
      totalSumOfChainLenghtsNodeRHasToInherit += chainLenghtsNodeRHasToInherit[i];
    }

    //////////////////////////////////////////////////////////////////////////
    // Printout important information
    //////////////////////////////////////////////////////////////////////////
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ": r = "                                              << r
                              << ", finalNumChainsPerNode[r] = "                       << finalNumChainsPerNode[r]
                              << ", totalNumberOfInitialPositionsNodeRHasToReceive = " << totalNumberOfInitialPositionsNodeRHasToReceive
                              << ", numberOfInitialPositionsNodeRAlreadyHas = "        << numberOfInitialPositionsNodeRAlreadyHas
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ": r = "                                       << r
                              << ", finalNumPositionsPerNode[r] = "             << finalNumPositionsPerNode[r]
                              << ", totalSumOfChainLenghtsNodeRHasToInherit = " << totalSumOfChainLenghtsNodeRHasToInherit
                              << ", sumOfChainLenghtsNodeRAlreadyHas = "        << sumOfChainLenghtsNodeRAlreadyHas
                              << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////
    // Make sanity checks
    //////////////////////////////////////////////////////////////////////////
    UQ_FATAL_TEST_MACRO(indexesOfInitialPositionsNodeRHasToReceiveFromMe.size() != numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()],
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()",
                        "inconsistent number of initial positions to send to node 'r'");

    UQ_FATAL_TEST_MACRO(finalNumChainsPerNode[r] != (totalNumberOfInitialPositionsNodeRHasToReceive + numberOfInitialPositionsNodeRAlreadyHas),
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()",
                        "inconsistent number of chains in node 'r'");

    UQ_FATAL_TEST_MACRO(finalNumPositionsPerNode[r] != (totalSumOfChainLenghtsNodeRHasToInherit + sumOfChainLenghtsNodeRAlreadyHas),
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()",
                        "inconsistent sum of chain lenghts in node 'r'");

    UQ_FATAL_TEST_MACRO(totalNumberOfInitialPositionsNodeRHasToReceive != totalNumberOfChainLenghtsNodeRHasToInherit,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()",
                        "inconsistent on total number of initial positions to receive in node 'r'");

    // Optimize use of memory (FIX ME: don't need to use swap here ????)
    indexesOfInitialPositionsNodeRHasToReceiveFromMe.resize(numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()]);
    chainLenghtsNodeRHasToInherit.resize                   (totalSumOfChainLenghtsNodeRHasToInherit);

    //////////////////////////////////////////////////////////////////////////
    // Prepare counters and buffers for gatherv of initial positions
    //////////////////////////////////////////////////////////////////////////
    unsigned int dimSize = m_vectorSpace.dimLocal();
    P_V auxInitialPosition(m_vectorSpace.zeroVector());
    std::vector<double> sendbuf(0);
    unsigned int sendcnt = 0;
    if (m_env.inter0Rank() != (int) r) {
      sendcnt = numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()] * dimSize;
      sendbuf.resize(sendcnt);
      for (unsigned int i = 0; i < numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()]; ++i) {
        unsigned int auxIndex = indexesOfInitialPositionsNodeRHasToReceiveFromMe[i];
        prevChain.getPositionValues(auxIndex,auxInitialPosition);
        for (unsigned int j = 0; j < dimSize; ++j) {
          sendbuf[i*dimSize + j] = auxInitialPosition[j];
        }
      }
    }

    std::vector<double> recvbuf(0);
    std::vector<int> recvcnts(Np,0); // '0' is already the correct value for recvcnts[r]
    if (m_env.inter0Rank() == (int) r) {
      recvbuf.resize(totalNumberOfInitialPositionsNodeRHasToReceive * dimSize);
      for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) { // Yes, from '0' on (for 'r', numberOf...ToReceiveFromNode[r] = 0 anyway)
        recvcnts[nodeId] = numberOfInitialPositionsNodeRHasToReceiveFromNode[nodeId/*m_env.inter0Rank()*/] * dimSize;
      }
    }

    std::vector<int> displs(Np,0);
    for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, from '1' on
      displs[nodeId] = displs[nodeId-1] + recvcnts[nodeId-1];
    }

    //int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
    //                void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, 
    //                int root, MPI_Comm comm )
    int mpiRC = 0;
#if 0
    if (m_env.inter0Rank() == r) {
      mpiRC = MPI_Gatherv(MPI_IN_PLACE, (int) sendcnt, MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], MPI_DOUBLE, r, m_env.inter0Comm().Comm()); // LOAD BALANCE
    }
    else {
      mpiRC = MPI_Gatherv((void *) &sendbuf[0], (int) sendcnt, MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], MPI_DOUBLE, r, m_env.inter0Comm().Comm()); // LOAD BALANCE
    }
#else
    mpiRC = MPI_Gatherv((void *) &sendbuf[0], (int) sendcnt, MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], MPI_DOUBLE, r, m_env.inter0Comm().Comm()); // LOAD BALANCE
#endif
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::mpiExchangePositions_inter0()",
                        "failed MPI_Gatherv()");

    //////////////////////////////////////////////////////////////////////////
    // Make sanity checks
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // Transfer data from 'recvbuf' to 'balancedLinkControl'
    // Remember that finalNumChainsPerNode[r] = (totalNumberOfInitialPositionsNodeRHasToReceive + numberOfInitialPositionsNodeRAlreadyHas)
    // Remember that totalNumberOfInitialPositionsNodeRHasToReceive = totalNumberOfChainLenghtsNodeRHasToInherit
    //////////////////////////////////////////////////////////////////////////
    if (m_env.inter0Rank() == (int) r) {
      balancedLinkControl.balLinkedChains.resize(finalNumChainsPerNode[r]);
      unsigned int auxIndex = 0;

      for (unsigned int i = 0; i < Nc; ++i) {
        if ((exchangeStdVec[i].finalNodeOfInitialPosition    == (int) r) &&
            (exchangeStdVec[i].originalNodeOfInitialPosition == (int) r)) {
          prevChain.getPositionValues(exchangeStdVec[i].originalIndexOfInitialPosition,auxInitialPosition);
          balancedLinkControl.balLinkedChains[auxIndex].initialPosition = new P_V(auxInitialPosition);
          balancedLinkControl.balLinkedChains[auxIndex].numberOfPositions = exchangeStdVec[i].numberOfPositions;
          auxIndex++;
	}
      }

      for (unsigned int i = 0; i < totalNumberOfInitialPositionsNodeRHasToReceive; ++i) {
        for (unsigned int j = 0; j < dimSize; ++j) {
          auxInitialPosition[j] = recvbuf[i*dimSize + j];
        }
        balancedLinkControl.balLinkedChains[auxIndex].initialPosition = new P_V(auxInitialPosition);
        balancedLinkControl.balLinkedChains[auxIndex].numberOfPositions = chainLenghtsNodeRHasToInherit[i]; // aqui 3
        auxIndex++;
      }
    }

    m_env.inter0Comm().Barrier();
  } // for 'r'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::mpiExchangePositions_inter0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << std::endl;
  }

  return;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING3_H__
