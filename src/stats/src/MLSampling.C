//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <algorithm>

#include <unistd.h> // sleep

#include <queso/MLSampling.h>
#include <queso/InstantiateIntersection.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BayesianJointPdf.h>

namespace QUESO {

#ifdef QUESO_HAS_GLPK

void BIP_routine(glp_tree *tree, void *info)
{
  const BaseEnvironment& env = *(((BIP_routine_struct *) info)->env);
  unsigned int currLevel            =   ((BIP_routine_struct *) info)->currLevel;

  int reason = glp_ios_reason(tree);

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 1)) {
    *env.subDisplayFile() << "In BIP_routine()"
                          << ", level " << currLevel+LEVEL_REF_ID
                          << ": glp_ios_reason() = " << reason
                          << std::endl;
  }
  std::cout << "In BIP_routine: reason = " << reason << std::endl;

  switch (reason) {
    case GLP_IROWGEN: // 0x01  /* request for row generation       */
      sleep(1);
    break;

    case GLP_IBINGO:  // 0x02  /* better integer solution found    */
      sleep(1);
    break;

    case GLP_IHEUR:   // 0x03  /* request for heuristic solution   */
      // Do nothing
    break;

    case GLP_ICUTGEN: // 0x04  /* request for cut generation       */
      // Do nothing
    break;

    case GLP_IBRANCH: // 0x05  /* request for branching            */
      // Do nothing
    break;

    case GLP_ISELECT: // 0x06  /* request for subproblem selection */
      // Do nothing
    break;

    case GLP_IPREPRO: // 0x07  /* request for preprocessing        */
      // Do nothing
    break;

    default:
      queso_error_msg("invalid glp_ios_readon");
    break;
  }

  return;
}

#endif // QUESO_HAS_GLPK

template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::sampleIndexes_proc0(
  unsigned int               unifiedRequestedNumSamples,        // input
  const std::vector<double>& unifiedWeightStdVectorAtProc0Only, // input
  std::vector<unsigned int>& unifiedIndexCountersAtProc0Only)   // output
{
  if (m_env.inter0Rank() != 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::sampleIndexes_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": unifiedRequestedNumSamples = "               << unifiedRequestedNumSamples
                            << ", unifiedWeightStdVectorAtProc0Only.size() = " << unifiedWeightStdVectorAtProc0Only.size()
                            << std::endl;
  }

#if 0 // For debug only
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::sampleIndexes_proc0()"
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
    FiniteDistribution tmpFd(m_env,
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
MLSampling<P_V,P_M>::decideOnBalancedChains_all(
  const MLSamplingLevelOptions* currOptions,                     // input
  unsigned int                         indexOfFirstWeight,              // input
  unsigned int                         indexOfLastWeight,               // input
  const std::vector<unsigned int>&     unifiedIndexCountersAtProc0Only, // input
  std::vector<ExchangeInfoStruct>&   exchangeStdVec)                  // output
{
  bool result = false;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::decideOnBalancedChains_all()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": indexOfFirstWeight = " << indexOfFirstWeight
                            << ", indexOfLastWeight = "  << indexOfLastWeight
                            << std::endl;
  }

  if (true) {
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
      unsigned int auxUInt = indexOfFirstWeight;
      m_env.inter0Comm().template Gather<unsigned int>(&auxUInt, 1, &allFirstIndexes[0], (int) 1, 0, // LOAD BALANCE
                                "MLSampling<P_V,P_M>::decideOnBalancedChains_all()",
                                "failed MPI.Gather() for first indexes");

      if (m_env.inter0Rank() == 0) {
        queso_require_equal_to_msg(allFirstIndexes[0], indexOfFirstWeight, "failed MPI.Gather() result for first indexes, at proc 0");
      }

      auxUInt = indexOfLastWeight;
      m_env.inter0Comm().template Gather<unsigned int>(&auxUInt, 1, &allLastIndexes[0], (int) 1, 0, // LOAD BALANCE
                                "MLSampling<P_V,P_M>::decideOnBalancedChains_all()",
                                "failed MPI.Gather() for last indexes");

      if (m_env.inter0Rank() == 0) { // Yes, '== 0'
        //allLastIndexes[0] = indexOfLastWeight; // FIX ME: really necessary????
        queso_require_equal_to_msg(allLastIndexes[0], indexOfLastWeight, "failed MPI.Gather() result for last indexes, at proc 0");
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Proc 0 prepares information to decide if load balancing is needed
    //////////////////////////////////////////////////////////////////////////
    if (m_env.inter0Rank() == 0) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::decideOnBalancedChains_all()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": original distribution of unified indexes in 'inter0Comm' is as follows"
                                << std::endl;
        for (unsigned int r = 0; r < Np; ++r) {
          *m_env.subDisplayFile() << "  allFirstIndexes[" << r << "] = " << allFirstIndexes[r]
                                  << "  allLastIndexes["  << r << "] = " << allLastIndexes[r]
                                  << std::endl;
        }
      }
      for (unsigned int r = 0; r < (Np-1); ++r) { // Yes, '-1'
        queso_require_equal_to_msg(allFirstIndexes[r+1], (allLastIndexes[r]+1), "wrong indexes");
      }

      for (unsigned int r = 0; r < (Np-1); ++r) { // Yes, '-1'
        queso_require_equal_to_msg(allFirstIndexes[r+1], (allLastIndexes[r]+1), "wrong indexes");
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
      std::cerr << "In MLSampling<P_V,P_M>::decideOnBalancedChains_all()"
                      << ": i = " << i
                      << ", r = " << r
                      << ", allFirstIndexes[r] = " << allFirstIndexes[r]
                      << ", allLastIndexes[r] = "  << allLastIndexes[r]
                      << std::endl;
            queso_error_msg("wrong indexes or 'r' got too large");
          }
        }
        if (unifiedIndexCountersAtProc0Only[i] != 0) {
          origNumChainsPerNode   [r] += 1;
          origNumPositionsPerNode[r] += unifiedIndexCountersAtProc0Only[i];

          ExchangeInfoStruct auxInfo;
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
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "  KEY"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", Np = "  << Np
                                << ", totalNumberOfChains = " << totalNumberOfChains
                                << std::endl;
      }

      // Check if ratio max/min justifies optimization
      unsigned int origMinPosPerNode  = *std::min_element(origNumPositionsPerNode.begin(), origNumPositionsPerNode.end());
      unsigned int origMaxPosPerNode  = *std::max_element(origNumPositionsPerNode.begin(), origNumPositionsPerNode.end());
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
          *m_env.subDisplayFile() << "  KEY"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ", origNumChainsPerNode["     << nodeId << "] = " << origNumChainsPerNode[nodeId]
                                  << ", origNumPositionsPerNode["  << nodeId << "] = " << origNumPositionsPerNode[nodeId]
                                  << std::endl;
        }
      }
      double origRatioOfPosPerNode = ((double) origMaxPosPerNode ) / ((double) origMinPosPerNode);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "  KEY"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", origRatioOfPosPerNode = "      << origRatioOfPosPerNode
                                << ", option loadBalanceTreshold = " << currOptions->m_loadBalanceTreshold
                                << std::endl;
      }

      // At this point, only proc 0 is running...
      // Set boolean 'result' for good
      if ((currOptions->m_loadBalanceAlgorithmId > 0                                 ) &&
          (m_env.numSubEnvironments()            > 1                                 ) && // Cannot use 'm_env.inter0Comm().NumProc()': not all nodes at this point of the code belong to 'inter0Comm'
          (Np                                    < totalNumberOfChains               ) &&
          (origRatioOfPosPerNode                 > currOptions->m_loadBalanceTreshold)) {
        result = true;
      }
    } // if (m_env.inter0Rank() == 0)
  } // if (true)

  m_env.fullComm().Barrier();
  unsigned int tmpValue = result;
  m_env.fullComm().Bcast((void *) &tmpValue, (int) 1, RawValue_MPI_UNSIGNED, 0, // LOAD BALANCE
                         "MLSampling<P_V,P_M>::decideOnBalancedChains_all()",
                         "failed MPI.Bcast() for 'result'");
  if (m_env.inter0Rank() != 0) result = tmpValue;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::decideOnBalancedChains_all()"
                            << ", level "    << m_currLevel+LEVEL_REF_ID
                            << ", step "     << m_currStep
                            << ": result = " << result
                            << std::endl;
  }

  return result;
}

template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0( // EXTRA FOR LOAD BALANCE
  const MLSamplingLevelOptions*      currOptions,         // input
  const SequenceOfVectors<P_V,P_M>&  prevChain,           // input
  double                                  prevExponent,                       // input
  double                                  currExponent,                       // input
  const ScalarSequence<double>&           prevLogLikelihoodValues,            // input
  const ScalarSequence<double>&           prevLogTargetValues,                // input
  std::vector<ExchangeInfoStruct>&        exchangeStdVec,      // input/output
  BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl) // output
{
  if (m_env.inter0Rank() < 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()"
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
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "WARNING in MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()"
                                  << ": algorithm id '" << currOptions->m_loadBalanceAlgorithmId
                                  << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                                  << ". Code will therefore process the algorithm id '" << 2
                                  << "' instead..."
                                  << std::endl;
        }
        if (m_env.subRank() == 0) {
          std::cerr << "WARNING in MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()"
                    << ": algorithm id '" << currOptions->m_loadBalanceAlgorithmId
                    << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                    << ". Code will therefore process the algorithm id '" << 2
                    << "' instead..."
                    << std::endl;
        }
        justBalance_proc0(currOptions,     // input
                          exchangeStdVec); // input/output
#endif
      break;
    }
  } // if (m_env.inter0Rank() == 0)

  m_env.inter0Comm().Barrier();

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 now broadcasts the information on 'exchangeStdVec'
  //////////////////////////////////////////////////////////////////////////
  unsigned int exchangeStdVecSize = exchangeStdVec.size();
  m_env.inter0Comm().Bcast((void *) &exchangeStdVecSize, (int) 1, RawValue_MPI_UNSIGNED, 0, // LOAD BALANCE
                           "MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()",
                           "failed MPI.Bcast() for exchangeStdVec size");
  if (m_env.inter0Rank() > 0) exchangeStdVec.resize(exchangeStdVecSize);

  m_env.inter0Comm().Bcast((void *) &exchangeStdVec[0], (int) (exchangeStdVecSize*sizeof(ExchangeInfoStruct)), RawValue_MPI_CHAR, 0, // LOAD BALANCE
                           "MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()",
                           "failed MPI.Bcast() for exchangeStdVec data");

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
  m_env.inter0Comm().template Allreduce<double>(&auxBuf[0], &minRatio, (int) auxBuf.size(), RawValue_MPI_MIN, // LOAD BALANCE
                               "MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()",
                               "failed MPI.Allreduce() for min");
  //std::cout << m_env.worldRank() << ", minRatio = " << minRatio << std::endl;
  queso_require_equal_to_msg(minRatio, finalRatioOfPosPerNode, "failed minRatio sanity check");

  double maxRatio = 0.;
  auxBuf[0] = finalRatioOfPosPerNode;
  m_env.inter0Comm().template Allreduce<double>(&auxBuf[0], &maxRatio, (int) auxBuf.size(), RawValue_MPI_MAX, // LOAD BALANCE
                               "MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()",
                               "failed MPI.Allreduce() for max");
  //std::cout << m_env.worldRank() << ", maxRatio = " << maxRatio << std::endl;
  queso_require_equal_to_msg(maxRatio, finalRatioOfPosPerNode, "failed maxRatio sanity check");

  //////////////////////////////////////////////////////////////////////////
  // Proc 0 now broadcasts the information on 'finalNumChainsPerNode'
  //////////////////////////////////////////////////////////////////////////
  unsigned int finalNumChainsPerNodeSize = finalNumChainsPerNode.size();
  m_env.inter0Comm().Bcast((void *) &finalNumChainsPerNodeSize, (int) 1, RawValue_MPI_UNSIGNED, 0, // LOAD BALANCE
                           "MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()",
                           "failed MPI.Bcast() for finalNumChainsPerNode size");
  if (m_env.inter0Rank() > 0) finalNumChainsPerNode.resize(finalNumChainsPerNodeSize);

  m_env.inter0Comm().Bcast((void *) &finalNumChainsPerNode[0], (int) finalNumChainsPerNodeSize, RawValue_MPI_UNSIGNED, 0, // LOAD BALANCE
                           "MLSampling<P_V,P_M>::prepareBalLinkedChains_inter0()",
                           "failed MPI.Bcast() for finalNumChainsPerNode data");

  //////////////////////////////////////////////////////////////////////////
  // Mpi exchange information between nodes and properly populate
  // balancedLinkControl.linkedChains at each node
  //////////////////////////////////////////////////////////////////////////
  mpiExchangePositions_inter0(prevChain,
                              prevExponent,
                              currExponent,
                              prevLogLikelihoodValues,
                              prevLogTargetValues,
                              exchangeStdVec,
                              finalNumChainsPerNode,
                              finalNumPositionsPerNode, // It is already valid at all "management" nodes (not only at node 0) because of the sanity check above
                              balancedLinkControl);

  return;
}

template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0(
  unsigned int                           indexOfFirstWeight,              // input
  unsigned int                           indexOfLastWeight,               // input
  const std::vector<unsigned int>&       unifiedIndexCountersAtProc0Only, // input
  UnbalancedLinkedChainsPerNodeStruct& unbalancedLinkControl)           // output
{
  if (m_env.inter0Rank() < 0) return;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
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
  m_env.inter0Comm().Bcast((void *) &resizeSize, (int) 1, RawValue_MPI_UNSIGNED, 0,
                           "MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                           "failed MPI.Bcast() for resizeSize");
  unifiedIndexCountersAtAllProcs.resize(resizeSize,0);

  if (m_env.inter0Rank() == 0) unifiedIndexCountersAtAllProcs = unifiedIndexCountersAtProc0Only;

  // Broadcast index counters to all nodes
  m_env.inter0Comm().Bcast((void *) &unifiedIndexCountersAtAllProcs[0], (int) unifiedIndexCountersAtAllProcs.size(), RawValue_MPI_UNSIGNED, 0,
                           "MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                           "failed MPI.Bcast() for unified index counters");
#if 0 // Use allgatherv ??? for subNumSamples instead
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
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
  queso_require_less_msg(indexOfFirstWeight, unifiedIndexCountersAtAllProcs.size(), "invalid indexOfFirstWeight");
  queso_require_less_msg(indexOfLastWeight, unifiedIndexCountersAtAllProcs.size(), "invalid indexOfLastWeight");
  subNumSamples = 0;
  for (unsigned int i = indexOfFirstWeight; i <= indexOfLastWeight; ++i) {
    subNumSamples += unifiedIndexCountersAtAllProcs[i];
  }

  std::vector<unsigned int> auxBuf(1,0);

  unsigned int minModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &minModifiedSubNumSamples, (int) auxBuf.size(), RawValue_MPI_MIN,
                               "MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                               "failed MPI.Allreduce() for min");

  unsigned int maxModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &maxModifiedSubNumSamples, (int) auxBuf.size(), RawValue_MPI_MAX,
                               "MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                               "failed MPI.Allreduce() for max");

  unsigned int sumModifiedSubNumSamples = 0;
  auxBuf[0] = subNumSamples;
  m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &sumModifiedSubNumSamples, (int) auxBuf.size(), RawValue_MPI_SUM,
                               "MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()",
                               "failed MPI.Allreduce() for sum");


  //                    m_env.worldRank(),
  //                    "MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()",
  //                    "invalid state");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                                   << m_currLevel+LEVEL_REF_ID
                            << ", step "                                    << m_currStep
                            << ": subNumSamples = "                         << subNumSamples
                            << ", unifiedIndexCountersAtAllProcs.size() = " << unifiedIndexCountersAtAllProcs.size()
                            << std::endl;
    *m_env.subDisplayFile() << "KEY Leaving MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                      << m_currLevel+LEVEL_REF_ID
                            << ", step "                       << m_currStep
                            << ": minModifiedSubNumSamples = " << minModifiedSubNumSamples
                            << ", avgModifiedSubNumSamples = " << ((double) sumModifiedSubNumSamples)/((double) m_env.inter0Comm().NumProc())
                            << ", maxModifiedSubNumSamples = " << maxModifiedSubNumSamples
                            << std::endl;
  }

  unsigned int numberOfPositionsToGuaranteeForNode = subNumSamples;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
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
        UnbalancedLinkedChainControlStruct auxControl;
        auxControl.initialPositionIndexInPreviousChain = i;
        auxControl.numberOfPositions = unifiedIndexCountersAtAllProcs[i];
        unbalancedLinkControl.unbLinkedChains.push_back(auxControl);

        numberOfPositionsToGuaranteeForNode -= unifiedIndexCountersAtAllProcs[i];
        unifiedIndexCountersAtAllProcs[i] = 0;
      }
      else if ((unifiedIndexCountersAtAllProcs[i] == numberOfPositionsToGuaranteeForNode) &&
               (unifiedIndexCountersAtAllProcs[i] > 0                                   )) {
      //else { // KAUST4
        UnbalancedLinkedChainControlStruct auxControl;
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
        queso_error_msg("should never get here");
      }
    }
  }
  queso_require_equal_to_msg(numberOfPositionsToGuaranteeForNode, 0, "numberOfPositionsToGuaranteeForNode exited loop with wrong value");
  // FIX ME: swap trick to save memory

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY Leaving MLSampling<P_V,P_M>::prepareUnbLinkedChains_inter0()"
                            << ", level "                                          << m_currLevel+LEVEL_REF_ID
                            << ", step "                                           << m_currStep
                            << ": unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateBalLinkedChains_all( // EXTRA FOR LOAD BALANCE
  MLSamplingLevelOptions&                  inputOptions,            // input, only m_rawChainSize changes
  const P_M&                                      unifiedCovMatrix,        // input
  const GenericVectorRV  <P_V,P_M>&        rv,                      // input
  const BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,     // input // Round Rock
  SequenceOfVectors      <P_V,P_M>&        workingChain,            // output
  double&                                         cumulativeRunTime,       // output
  unsigned int&                                   cumulativeRejections,    // output
  ScalarSequence         <double>*         currLogLikelihoodValues, // output
  ScalarSequence         <double>*         currLogTargetValues)     // output
{
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                            << ": balancedLinkControl.balLinkedChains.size() = " << balancedLinkControl.balLinkedChains.size()
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());
  double auxInitialLogPrior;
  double auxInitialLogLikelihood;

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = balancedLinkControl.balLinkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  m_env.subComm().Bcast((void *) &chainIdMax, (int) 1, RawValue_MPI_UNSIGNED, 0, // Yes, 'subComm', important // LOAD BALANCE
                        "MLSampling<P_V,P_M>::generateBalLinkedChains_all()",
                        "failed MPI.Bcast() for chainIdMax");

  struct timeval timevalEntering;
  int iRC = 0;
  iRC = gettimeofday(&timevalEntering, NULL);
  if (iRC) {}; // just to remove compiler warning

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfPositions = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfPositions += balancedLinkControl.balLinkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &minNumberOfPositions, (int) auxBuf.size(), RawValue_MPI_MIN, // LOAD BALANCE
                                 "MLSampling<P_V,P_M>::generateBalLinkedChains_all()",
                                 "failed MPI.Allreduce() for min");

    unsigned int maxNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &maxNumberOfPositions, (int) auxBuf.size(), RawValue_MPI_MAX, // LOAD BALANCE
                                 "MLSampling<P_V,P_M>::generateBalLinkedChains_all()",
                                 "failed MPI.Allreduce() for max");

    unsigned int sumNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &sumNumberOfPositions, (int) auxBuf.size(), RawValue_MPI_SUM, // LOAD BALANCE
                                 "MLSampling<P_V,P_M>::generateBalLinkedChains_all()",
                                 "failed MPI.Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                              << ", level "               << m_currLevel+LEVEL_REF_ID
                              << ", step "                << m_currStep
                              << ": chainIdMax = "        << chainIdMax
                              << ", numberOfPositions = " << numberOfPositions
                              << ", at "                  << ctime(&timevalEntering.tv_sec)
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                              << ", level "                  << m_currLevel+LEVEL_REF_ID
                              << ", step "                   << m_currStep
                              << ": minNumberOfPositions = " << minNumberOfPositions
                              << ", avgNumberOfPositions = " << ((double) sumNumberOfPositions)/((double) m_env.inter0Comm().NumProc())
                              << ", maxNumberOfPositions = " << maxNumberOfPositions
                              << std::endl;
    }

    // 2013-02-23: print sizes, and expected final size
  }
  if ((m_debugExponent == 1.) &&
      (m_currStep      == 10)) {
    //m_env.setExceptionalCircumstance(true);
  }
  unsigned int cumulativeNumPositions = 0;
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      // aqui 4
      auxInitialPosition = *(balancedLinkControl.balLinkedChains[chainId].initialPosition); // Round Rock
      auxInitialLogPrior = balancedLinkControl.balLinkedChains[chainId].initialLogPrior;
      auxInitialLogLikelihood = balancedLinkControl.balLinkedChains[chainId].initialLogLikelihood;
      tmpChainSize = balancedLinkControl.balLinkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()       ) &&
          (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                                << ", level "            << m_currLevel+LEVEL_REF_ID
                                << ", step "             << m_currStep
                                << ", chainId = "        << chainId
                                << " < "                 << chainIdMax
                                << ": begin generating " << tmpChainSize
                                << " chain positions"
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm()); // Yes, 'subComm', important // KAUST

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
    m_env.subComm().Bcast((void *) &tmpChainSize, (int) 1, RawValue_MPI_UNSIGNED, 0, // Yes, 'subComm', important // LOAD BALANCE
                          "MLSampling<P_V,P_M>::generateBalLinkedChains_all()",
                          "failed MPI.Bcast() for tmpChainSize");

    inputOptions.m_rawChainSize = tmpChainSize;
    SequenceOfVectors<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    ScalarSequence<double> tmpLogLikelihoodValues(m_env,0,"");
    ScalarSequence<double> tmpLogTargetValues    (m_env,0,"");

    MHRawChainInfoStruct mcRawInfo;
    if (inputOptions.m_initialPositionUsePreviousLevelLikelihood) {  // ml_likelihood_caching
      m_env.subComm().Bcast((void *) &auxInitialLogPrior, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm', important
          "MLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
          "failed MPI.Bcast() for auxInitialLogPrior");
      m_env.subComm().Bcast((void *) &auxInitialLogLikelihood, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm', important
          "MLSamplingClass<P_V,P_M>::generateBalLinkedChains_all()",
          "failed MPI.Bcast() for auxInitialLogLikelihood");

      MetropolisHastingsSG<P_V,P_M> mcSeqGenerator(inputOptions,
                                                   rv,
                                                   auxInitialPosition, // KEY new: pass logPrior and logLikelihood
                                                   auxInitialLogPrior,
                                                   auxInitialLogLikelihood,
                                                   &unifiedCovMatrix);
      mcSeqGenerator.generateSequence(tmpChain,
                                      &tmpLogLikelihoodValues, // likelihood is IMPORTANT
                                      &tmpLogTargetValues);
      mcSeqGenerator.getRawChainInfo(mcRawInfo);
    }
    else {
      MetropolisHastingsSG<P_V,P_M> mcSeqGenerator(inputOptions,
          rv,
          auxInitialPosition,
          &unifiedCovMatrix);
      mcSeqGenerator.generateSequence(tmpChain,
          &tmpLogLikelihoodValues, // likelihood is IMPORTANT
          &tmpLogTargetValues);
      mcSeqGenerator.getRawChainInfo(mcRawInfo);
    }

    cumulativeRunTime    += mcRawInfo.runTime;
    cumulativeRejections += mcRawInfo.numRejections;

    if (m_env.inter0Rank() >= 0) {
      if (m_env.exceptionalCircumstance()) {
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
      if (cumulativeNumPositions > 100) m_env.setExceptionalCircumstance(false);

      if ((m_env.subDisplayFile()       ) &&
          (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                                << ", level "               << m_currLevel+LEVEL_REF_ID
                                << ", step "                << m_currStep
                                << ", chainId = "           << chainId
                                << " < "                    << chainIdMax
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
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
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
      // 2013-02-23: print size just appended
    }
  } // for 'chainId'

  // 2013-02-23: print final size

  struct timeval timevalBarrier;
  iRC = gettimeofday(&timevalBarrier, NULL);
  if (iRC) {}; // just to remove compiler warning
  double loopTime = MiscGetEllapsedSeconds(&timevalEntering);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": ended chain loop after " << loopTime << " seconds"
                            << ", calling fullComm().Barrier() at " << ctime(&timevalBarrier.tv_sec)
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  struct timeval timevalLeaving;
  iRC = gettimeofday(&timevalLeaving, NULL);
  if (iRC) {}; // just to remove compiler warning
  double barrierTime = MiscGetEllapsedSeconds(&timevalBarrier);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateBalLinkedChains_all()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": after " << barrierTime << " seconds in fullComm().Barrier()"
                            << ", at " << ctime(&timevalLeaving.tv_sec)
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateUnbLinkedChains_all(
  MLSamplingLevelOptions&               inputOptions,            // input, only m_rawChainSize changes
  const P_M&                                   unifiedCovMatrix,        // input
  const GenericVectorRV  <P_V,P_M>&     rv,                      // input
  const UnbalancedLinkedChainsPerNodeStruct& unbalancedLinkControl,   // input // Round Rock
  unsigned int                                 indexOfFirstWeight,      // input // Round Rock
  const SequenceOfVectors<P_V,P_M>&     prevChain,               // input // Round Rock
  double                                prevExponent,            // input
  double                                currExponent,            // input
  const ScalarSequence<double>&         prevLogLikelihoodValues, // input
  const ScalarSequence<double>&         prevLogTargetValues,     // input
  SequenceOfVectors      <P_V,P_M>&     workingChain,            // output
  double&                                      cumulativeRunTime,       // output
  unsigned int&                                cumulativeRejections,    // output
  ScalarSequence         <double>*      currLogLikelihoodValues, // output
  ScalarSequence         <double>*      currLogTargetValues)     // output
{
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                            << ": unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                            << ", indexOfFirstWeight = "                           << indexOfFirstWeight
                            << std::endl;
  }

  P_V auxInitialPosition(m_vectorSpace.zeroVector());
  double auxInitialLogPrior;
  double auxInitialLogLikelihood;

  unsigned int chainIdMax = 0;
  if (m_env.inter0Rank() >= 0) {
    chainIdMax = unbalancedLinkControl.unbLinkedChains.size();
  }
  // KAUST: all nodes in 'subComm' should have the same 'chainIdMax'
  m_env.subComm().Bcast((void *) &chainIdMax, (int) 1, RawValue_MPI_UNSIGNED, 0, // Yes, 'subComm', important
                        "MLSampling<P_V,P_M>::generateUnbLinkedChains_all()",
                        "failed MPI.Bcast() for chainIdMax");

  struct timeval timevalEntering;
  int iRC = 0;
  iRC = gettimeofday(&timevalEntering, NULL);
  if (iRC) {}; // just to remove compiler warning

  if (m_env.inter0Rank() >= 0) {
    unsigned int numberOfPositions = 0;
    for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
      numberOfPositions += unbalancedLinkControl.unbLinkedChains[chainId].numberOfPositions;
    }

    std::vector<unsigned int> auxBuf(1,0);

    unsigned int minNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &minNumberOfPositions, (int) auxBuf.size(), RawValue_MPI_MIN,
                                 "MLSampling<P_V,P_M>::generateUnbLinkedChains_all()",
                                 "failed MPI.Allreduce() for min");

    unsigned int maxNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &maxNumberOfPositions, (int) auxBuf.size(), RawValue_MPI_MAX,
                                 "MLSampling<P_V,P_M>::generateUnbLinkedChains_all()",
                                 "failed MPI.Allreduce() for max");

    unsigned int sumNumberOfPositions = 0;
    auxBuf[0] = numberOfPositions;
    m_env.inter0Comm().template Allreduce<unsigned int>(&auxBuf[0], &sumNumberOfPositions, (int) auxBuf.size(), RawValue_MPI_SUM,
                                 "MLSampling<P_V,P_M>::generateUnbLinkedChains_all()",
                                 "failed MPI.Allreduce() for sum");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                              << ", level "               << m_currLevel+LEVEL_REF_ID
                              << ", step "                << m_currStep
                              << ": chainIdMax = "        << chainIdMax
                              << ", numberOfPositions = " << numberOfPositions
                              << ", at "                  << ctime(&timevalEntering.tv_sec)
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
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
    //m_env.setExceptionalCircumstance(true);
  }
  double expRatio = currExponent;
  if (prevExponent > 0.0) {
    expRatio /= prevExponent;
  }
  unsigned int cumulativeNumPositions = 0;
  for (unsigned int chainId = 0; chainId < chainIdMax; ++chainId) {
    unsigned int tmpChainSize = 0;
    if (m_env.inter0Rank() >= 0) {
      unsigned int auxIndex = unbalancedLinkControl.unbLinkedChains[chainId].initialPositionIndexInPreviousChain - indexOfFirstWeight; // KAUST4 // Round Rock
      prevChain.getPositionValues(auxIndex,auxInitialPosition); // Round Rock
      auxInitialLogPrior = prevLogTargetValues[auxIndex] - prevLogLikelihoodValues[auxIndex];
      auxInitialLogLikelihood = expRatio * prevLogLikelihoodValues[auxIndex];
      tmpChainSize = unbalancedLinkControl.unbLinkedChains[chainId].numberOfPositions+1; // IMPORTANT: '+1' in order to discard initial position afterwards
      if ((m_env.subDisplayFile()       ) &&
          (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                                << ", level "            << m_currLevel+LEVEL_REF_ID
                                << ", step "             << m_currStep
                                << ", chainId = "        << chainId
                                << " < "                 << chainIdMax
                                << ": begin generating " << tmpChainSize
                                << " chain positions"
                                << std::endl;
      }
    }
    auxInitialPosition.mpiBcast(0, m_env.subComm()); // Yes, 'subComm', important // KAUST

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
    m_env.subComm().Bcast((void *) &tmpChainSize, (int) 1, RawValue_MPI_UNSIGNED, 0, // Yes, 'subComm', important
                          "MLSampling<P_V,P_M>::generateUnbLinkedChains_all()",
                          "failed MPI.Bcast() for tmpChainSize");

    inputOptions.m_rawChainSize = tmpChainSize;
    SequenceOfVectors<P_V,P_M> tmpChain(m_vectorSpace,
                                               0,
                                               m_options.m_prefix+"tmp_chain");
    ScalarSequence<double> tmpLogLikelihoodValues(m_env,0,"");
    ScalarSequence<double> tmpLogTargetValues    (m_env,0,"");

    // KAUST: all nodes should call here
    MHRawChainInfoStruct mcRawInfo;
    if (inputOptions.m_initialPositionUsePreviousLevelLikelihood) {  // ml_likelihood_caching
      m_env.subComm().Bcast((void *) &auxInitialLogPrior, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm', important
          "MLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all()",
          "failed MPI.Bcast() for auxInitialLogPrior");
      m_env.subComm().Bcast((void *) &auxInitialLogLikelihood, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm', important
          "MLSamplingClass<P_V,P_M>::generateUnbLinkedChains_all",
          "failed MPI.Bcast() for auxInitialLogLikelihood");
      MetropolisHastingsSG<P_V,P_M> mcSeqGenerator(inputOptions,
          rv,
          auxInitialPosition, // KEY new: pass logPrior and logLikelihood
          auxInitialLogPrior,
          auxInitialLogLikelihood,
          &unifiedCovMatrix);
      mcSeqGenerator.generateSequence(tmpChain,
          &tmpLogLikelihoodValues, // likelihood is IMPORTANT
          &tmpLogTargetValues);
      mcSeqGenerator.getRawChainInfo(mcRawInfo);
    }
    else {
      MetropolisHastingsSG<P_V,P_M> mcSeqGenerator(inputOptions,
          rv,
          auxInitialPosition,
          &unifiedCovMatrix);
      mcSeqGenerator.generateSequence(tmpChain,
          &tmpLogLikelihoodValues, // likelihood is IMPORTANT
          &tmpLogTargetValues);
      mcSeqGenerator.getRawChainInfo(mcRawInfo);
    }

    cumulativeRunTime    += mcRawInfo.runTime;
    cumulativeRejections += mcRawInfo.numRejections;

    if (m_env.inter0Rank() >= 0) {
      if (m_env.exceptionalCircumstance()) {
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
      if (cumulativeNumPositions > 100) m_env.setExceptionalCircumstance(false);

      if ((m_env.subDisplayFile()       ) &&
          (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                                << ", level "               << m_currLevel+LEVEL_REF_ID
                                << ", step "                << m_currStep
                                << ", chainId = "           << chainId
                                << " < "                    << chainIdMax
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
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
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
  } // for 'chainId'

  struct timeval timevalBarrier;
  iRC = gettimeofday(&timevalBarrier, NULL);
  if (iRC) {}; // just to remove compiler warning
  double loopTime = MiscGetEllapsedSeconds(&timevalEntering);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": ended chain loop after " << loopTime << " seconds"
                            << ", calling fullComm().Barrier() at " << ctime(&timevalBarrier.tv_sec)
                            << std::endl;
  }

  m_env.fullComm().Barrier(); // KAUST4

  struct timeval timevalLeaving;
  iRC = gettimeofday(&timevalLeaving, NULL);
  if (iRC) {}; // just to remove compiler warning
  double barrierTime = MiscGetEllapsedSeconds(&timevalBarrier);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateUnbLinkedChains_all()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": after " << barrierTime << " seconds in fullComm().Barrier()"
                            << ", at " << ctime(&timevalLeaving.tv_sec)
                            << std::endl;
  }

  return;
}

#ifdef QUESO_HAS_GLPK
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::solveBIP_proc0( // EXTRA FOR LOAD BALANCE
  std::vector<ExchangeInfoStruct>& exchangeStdVec) // input/output
{
  if (m_env.inter0Rank() != 0) return;

  int iRC = UQ_OK_RC;
  struct timeval timevalBIP;
  iRC = gettimeofday(&timevalBIP, NULL);
  if (iRC) {}; // just to remove compiler warning

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  unsigned int Nc = exchangeStdVec.size();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::solveBIP_proc0()"
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
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
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
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished setting BIP constraint matrix"
                            << ", ne = "     << ne
                            << ", coefId = " << coefId
                            << std::endl;
  }

  queso_require_equal_to_msg(coefId, (int) (ne+1), "invalid final coefId");

  glp_load_matrix(lp, ne, &iVec[0], &jVec[0], &aVec[0]);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
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
  queso_require_equal_to_msg(glp_get_num_rows(lp), (int) m, "invalid number of rows");

  queso_require_equal_to_msg(glp_get_num_cols(lp), (int) n, "invalid number of columnss");

  queso_require_equal_to_msg(glp_get_num_nz(lp), (int) ne, "invalid number of nonzero constraint coefficients");

  queso_require_equal_to_msg(glp_get_num_int(lp), (int) n, "invalid number of integer structural variables");

  queso_require_equal_to_msg(glp_get_num_bin(lp), (int) n, "invalid number of binary structural variables");

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
        queso_require_equal_to_msg(initialState, 1, "for nodeId = 0, initial state should be '1'");
      }
      else {
        queso_require_equal_to_msg(initialState, 0, "for nodeId > 0, initial state should be '0'");
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
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished solving BIP"
                            << ", BIP_rc = " << BIP_rc
                            << std::endl;
  }

  queso_require_equal_to_msg(BIP_rc, 0, "BIP returned rc != 0");

  //////////////////////////////////////////////////////////////////////////
  // Check BIP status after solution
  //////////////////////////////////////////////////////////////////////////
  int BIP_Status = glp_mip_status(lp);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": BIP_Status = " << BIP_Status
                            << std::endl;
  }

  switch (BIP_Status) {
    case GLP_OPT:
      // Ok
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": BIP solution is optimal"
                                << std::endl;
      }
    break;

    case GLP_FEAS:
      // Ok
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": BIP solution is guaranteed to be 'only' feasible"
                                << std::endl;
      }
    break;

    default:
      queso_error_msg("BIP has an undefined solution or has no solution");
    break;
  }

  for (int i = 1; i <= (int) Nc; ++i) {
    queso_require_equal_to_msg(glp_mip_row_val(lp, i), 1, "row should have value 1 at solution");
  }
  for (int i = (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    queso_require_less_equal_msg(glp_mip_row_val(lp, i), 0, "row should have value 0 or should be negative at solution");
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
        queso_require_equal_to_msg(exchangeStdVec[chainId].finalNodeOfInitialPosition, -1, "chain has already been taken care of");
        exchangeStdVec[chainId].finalNodeOfInitialPosition = nodeId;
        finalNumChainsPerNode   [nodeId] += 1;
        finalNumPositionsPerNode[nodeId] += exchangeStdVec[chainId].numberOfPositions;
      }
      else {
        queso_error_msg("control variable should be either '0' or '1'");
      }
    }
  }

  unsigned int finalMinPosPerNode = *std::min_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());
  unsigned int finalMaxPosPerNode = *std::max_element(finalNumPositionsPerNode.begin(), finalNumPositionsPerNode.end());

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": finished preparing output information"
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Printout solution information
  //////////////////////////////////////////////////////////////////////////
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": solution gives the following redistribution"
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY In MLSampling<P_V,P_M>::solveBIP_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", finalNumChainsPerNode["    << nodeId << "] = " << finalNumChainsPerNode[nodeId]
                              << ", finalNumPositionsPerNode[" << nodeId << "] = " << finalNumPositionsPerNode[nodeId]
                              << std::endl;
    }
    *m_env.subDisplayFile() << "  KEY In MLSampling<P_V,P_M>::solveBIP_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", finalRatioOfPosPerNode = " << ((double) finalMaxPosPerNode) / ((double)finalMinPosPerNode)
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Make sanity checks
  //////////////////////////////////////////////////////////////////////////
  queso_require_equal_to_msg(glp_mip_obj_val(lp), (double) finalNumPositionsPerNode[0], "Invalid objective value");

  for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, '1'
    queso_require_greater_equal_msg(finalNumPositionsPerNode[nodeId-1], finalNumPositionsPerNode[nodeId], "Next node should have a number of positions equal or less than the current node");
  }

  for (int i = (int) (Nc+1); i <= (int) (Nc+Np-1); ++i) {
    unsigned int nodeId = i - Nc;
    int diff = ((int) finalNumPositionsPerNode[nodeId]) - ((int) finalNumPositionsPerNode[nodeId-1]);
    queso_require_equal_to_msg(glp_mip_row_val(lp, i), diff, "wrong state");
  }

  //////////////////////////////////////////////////////////////////////////
  // Free memory and return
  //////////////////////////////////////////////////////////////////////////
  glp_delete_prob(lp);

  double bipRunTime = MiscGetEllapsedSeconds(&timevalBIP);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::solveBIP_proc0()"
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
MLSampling<P_V,P_M>::justBalance_proc0(
  const MLSamplingLevelOptions* currOptions,    // input
  std::vector<ExchangeInfoStruct>&   exchangeStdVec) // input/output
{
  if (m_env.inter0Rank() != 0) return;

  int iRC = UQ_OK_RC;
  struct timeval timevalBal;
  iRC = gettimeofday(&timevalBal, NULL);
  if (iRC) {}; // just to remove compiler warning

  unsigned int Np = m_env.numSubEnvironments();
  unsigned int Nc = exchangeStdVec.size();

  std::vector<ExchangeInfoStruct> currExchangeStdVec(Nc);
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
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::justBalance_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", iter "  << iterId
                            << ", currRatioOfPosPerNode = " << currRatioOfPosPerNode
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY In MLSampling<P_V,P_M>::justBalance_proc0()"
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
      queso_require_equal_to_msg(vectorOfChainSizesPerNode[nodeId].size(), currNumChainsPerNode[nodeId], "inconsistent number of chains in node");
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

    queso_require_equal_to_msg(currMinPosPerNode, currNumPositionsPerNode[currNodeWithLeastPositions], "inconsistent currMinPosPerNode");

    queso_require_equal_to_msg(currMaxPosPerNode, currNumPositionsPerNode[currNodeWithMostPositions], "inconsistent currMaxPosPerNode");

    unsigned int numberOfPositionsToMove = vectorOfChainSizesPerNode[currNodeWithMostPositions][0];

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::justBalance_proc0()"
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
    std::vector<ExchangeInfoStruct> newExchangeStdVec(Nc);
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
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::justBalance_proc0()"
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

    queso_require_equal_to_msg(newMinPosPerNode, newNumPositionsPerNode[newNodeWithLeastPositions], "inconsistent newMinPosPerNode");

    queso_require_equal_to_msg(newMaxPosPerNode, newNumPositionsPerNode[newNodeWithMostPositions], "inconsistent newMaxPosPerNode");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", iter "  << iterId
                              << ", newMaxPosPerNode = "     << newMaxPosPerNode
                              << ", newMinPosPerNode = "     << newMinPosPerNode
                              << ", newRatioOfPosPerNode = " << newRatioOfPosPerNode
                              << std::endl;
      for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
        *m_env.subDisplayFile() << "  KEY In MLSampling<P_V,P_M>::justBalance_proc0()"
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
    *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::justBalance_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": solution gives the following redistribution"
                            << std::endl;
    for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) {
      *m_env.subDisplayFile() << "  KEY In MLSampling<P_V,P_M>::justBalance_proc0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", finalNumChainsPerNode["    << nodeId << "] = " << finalNumChainsPerNode[nodeId]
                              << ", finalNumPositionsPerNode[" << nodeId << "] = " << finalNumPositionsPerNode[nodeId]
                              << std::endl;
    }
    *m_env.subDisplayFile() << "  KEY In MLSampling<P_V,P_M>::justBalance_proc0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", finalRatioOfPosPerNode = " << finalRatioOfPosPerNode
                            << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////
  // Measure time
  //////////////////////////////////////////////////////////////////////////
  double balRunTime = MiscGetEllapsedSeconds(&timevalBal);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::justBalance_proc0()"
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
MLSampling<P_V,P_M>::mpiExchangePositions_inter0( // EXTRA FOR LOAD BALANCE
  const SequenceOfVectors<P_V,P_M>&  prevChain,                // input
  double                             prevExponent,             // input
  double                             currExponent,             // input
  const ScalarSequence<double>&      prevLogLikelihoodValues,  // input
  const ScalarSequence<double>&      prevLogTargetValues,      // input
  const std::vector<ExchangeInfoStruct>&  exchangeStdVec,           // input
  const std::vector<unsigned int>&          finalNumChainsPerNode,    // input
  const std::vector<unsigned int>&          finalNumPositionsPerNode, // input
  BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl)      // output
{
  if (m_env.inter0Rank() < 0) {
    return;
  }

  unsigned int Np = (unsigned int) m_env.inter0Comm().NumProc();
  unsigned int Nc = exchangeStdVec.size();

  double expRatio = currExponent;
  if (prevExponent > 0.0) {
    expRatio /= prevExponent;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::mpiExchangePositions_inter0()"
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
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::mpiExchangePositions_inter0()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ": r = "                                              << r
                              << ", finalNumChainsPerNode[r] = "                       << finalNumChainsPerNode[r]
                              << ", totalNumberOfInitialPositionsNodeRHasToReceive = " << totalNumberOfInitialPositionsNodeRHasToReceive
                              << ", numberOfInitialPositionsNodeRAlreadyHas = "        << numberOfInitialPositionsNodeRAlreadyHas
                              << std::endl;
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::mpiExchangePositions_inter0()"
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
    queso_require_equal_to_msg(indexesOfInitialPositionsNodeRHasToReceiveFromMe.size(), numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()], "inconsistent number of initial positions to send to node 'r'");

    queso_require_equal_to_msg(finalNumChainsPerNode[r], (totalNumberOfInitialPositionsNodeRHasToReceive + numberOfInitialPositionsNodeRAlreadyHas), "inconsistent number of chains in node 'r'");

    queso_require_equal_to_msg(finalNumPositionsPerNode[r], (totalSumOfChainLenghtsNodeRHasToInherit + sumOfChainLenghtsNodeRAlreadyHas), "inconsistent sum of chain lenghts in node 'r'");

    queso_require_equal_to_msg(totalNumberOfInitialPositionsNodeRHasToReceive, totalNumberOfChainLenghtsNodeRHasToInherit, "inconsistent on total number of initial positions to receive in node 'r'");

    // Optimize use of memory (FIX ME: don't need to use swap here ????)
    indexesOfInitialPositionsNodeRHasToReceiveFromMe.resize(numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()]);
    chainLenghtsNodeRHasToInherit.resize                   (totalSumOfChainLenghtsNodeRHasToInherit);

    //////////////////////////////////////////////////////////////////////////
    // Prepare counters and buffers for gatherv of initial positions
    //////////////////////////////////////////////////////////////////////////
    unsigned int dimSize = m_vectorSpace.dimLocal();
    unsigned int nValuesPerInitialPosition = dimSize + 2;
    P_V auxInitialPosition(m_vectorSpace.zeroVector());
    std::vector<double> sendbuf(0);
    unsigned int sendcnt = 0;
    if (m_env.inter0Rank() != (int) r) {
      sendcnt = numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()] * nValuesPerInitialPosition;
      sendbuf.resize(sendcnt);
      for (unsigned int i = 0; i < numberOfInitialPositionsNodeRHasToReceiveFromNode[m_env.inter0Rank()]; ++i) {
        unsigned int auxIndex = indexesOfInitialPositionsNodeRHasToReceiveFromMe[i];
        prevChain.getPositionValues(auxIndex,auxInitialPosition);
        for (unsigned int j = 0; j < dimSize; ++j) {
          sendbuf[i*nValuesPerInitialPosition + j] = auxInitialPosition[j];
        }
      sendbuf[i*nValuesPerInitialPosition + dimSize]     = prevLogLikelihoodValues[auxIndex];
      sendbuf[i*nValuesPerInitialPosition + dimSize + 1] = prevLogTargetValues[auxIndex];
      }
    }

    std::vector<double> recvbuf(0);
    std::vector<int> recvcnts(Np,0); // '0' is already the correct value for recvcnts[r]
    if (m_env.inter0Rank() == (int) r) {
      recvbuf.resize(totalNumberOfInitialPositionsNodeRHasToReceive * nValuesPerInitialPosition);
      for (unsigned int nodeId = 0; nodeId < Np; ++nodeId) { // Yes, from '0' on (for 'r', numberOf...ToReceiveFromNode[r] = 0 anyway)
        recvcnts[nodeId] = numberOfInitialPositionsNodeRHasToReceiveFromNode[nodeId/*m_env.inter0Rank()*/] * nValuesPerInitialPosition;
      }
    }

    std::vector<int> displs(Np,0);
    for (unsigned int nodeId = 1; nodeId < Np; ++nodeId) { // Yes, from '1' on
      displs[nodeId] = displs[nodeId-1] + recvcnts[nodeId-1];
    }

#if 0
    if (m_env.inter0Rank() == r) {
      m_env.inter0Comm().Gatherv(RawValue_MPI_IN_PLACE, (int) sendcnt, RawValue_MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], RawValue_MPI_DOUBLE, r, // LOAD BALANCE
                                 "MLSampling<P_V,P_M>::mpiExchangePositions_inter0(1)",
                                 "failed MPI.Gatherv()");
    }
    else {
      m_env.inter0Comm().Gatherv((void *) &sendbuf[0], (int) sendcnt, RawValue_MPI_DOUBLE, (void *) &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0], RawValue_MPI_DOUBLE, r, // LOAD BALANCE
                                 "MLSampling<P_V,P_M>::mpiExchangePositions_inter0(2)",
                                 "failed MPI.Gatherv()");
    }
#else
    m_env.inter0Comm().template Gatherv<double>(&sendbuf[0], (int) sendcnt,
        &recvbuf[0], (int *) &recvcnts[0], (int *) &displs[0],
        r, // LOAD BALANCE
        "MLSampling<P_V,P_M>::mpiExchangePositions_inter0()",
        "failed MPI.Gatherv()");
#endif

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
          unsigned int originalIndex = exchangeStdVec[i].originalIndexOfInitialPosition;
          prevChain.getPositionValues(originalIndex, auxInitialPosition);
          balancedLinkControl.balLinkedChains[auxIndex].initialPosition = new P_V(auxInitialPosition);
          balancedLinkControl.balLinkedChains[auxIndex].initialLogPrior = prevLogTargetValues[originalIndex] - prevLogLikelihoodValues[originalIndex];
          balancedLinkControl.balLinkedChains[auxIndex].initialLogLikelihood = expRatio*prevLogLikelihoodValues[originalIndex];
          balancedLinkControl.balLinkedChains[auxIndex].numberOfPositions = exchangeStdVec[i].numberOfPositions;
          auxIndex++;
        }
      }

      for (unsigned int i = 0; i < totalNumberOfInitialPositionsNodeRHasToReceive; ++i) {
        for (unsigned int j = 0; j < dimSize; ++j) {
          auxInitialPosition[j] = recvbuf[i*nValuesPerInitialPosition + j];
        }
        balancedLinkControl.balLinkedChains[auxIndex].initialPosition = new P_V(auxInitialPosition);
        double prevLogLikelihood = recvbuf[i*nValuesPerInitialPosition + dimSize];
        double prevLogTarget     = recvbuf[i*nValuesPerInitialPosition + dimSize + 1];
        balancedLinkControl.balLinkedChains[auxIndex].initialLogPrior = prevLogTarget - prevLogLikelihood;
        balancedLinkControl.balLinkedChains[auxIndex].initialLogLikelihood = expRatio*prevLogLikelihood;
        balancedLinkControl.balLinkedChains[auxIndex].numberOfPositions = chainLenghtsNodeRHasToInherit[i]; // aqui 3
        auxIndex++;
      }
    }

    m_env.inter0Comm().Barrier();
  } // for 'r'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::mpiExchangePositions_inter0()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << std::endl;
  }

  return;
}

// Statistical/private methods-----------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::checkpointML(
  double                                   currExponent,            // input
  double                                   currEta,                 // input
  const SequenceOfVectors<P_V,P_M>& currChain,               // input
  const ScalarSequence<double>&     currLogLikelihoodValues, // input
  const ScalarSequence<double>&     currLogTargetValues)     // input
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING initiating at level " << m_currLevel
                            << "\n" << std::endl;
  }

  //******************************************************************************
  // Write 'control' file without 'level' spefication in name
  //******************************************************************************
  unsigned int quantity1 = currChain.unifiedSequenceSize();
  unsigned int quantity2 = currLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
  unsigned int quantity3 = currLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
  if (m_env.inter0Rank() >= 0) {
    queso_require_equal_to_msg(m_logEvidenceFactors.size(), m_currLevel, "number of evidence factors is not consistent");
    queso_require_equal_to_msg(quantity1, quantity2, "quantity2 is not consistent");
    queso_require_equal_to_msg(quantity1, quantity3, "quantity3 is not consistent");
  }

  if (m_env.fullRank() == 0) {
    std::ofstream* ofsVar = new std::ofstream((m_options.m_restartOutput_baseNameForFiles + "Control.txt").c_str(),
                                              std::ofstream::out | std::ofstream::trunc);
    *ofsVar << m_currLevel               << std::endl  // 1
            << m_vectorSpace.dimGlobal() << std::endl  // 2
            << currExponent              << std::endl  // 3
            << currEta                   << std::endl  // 4
            << quantity1                 << std::endl; // 5
    unsigned int savedPrecision = ofsVar->precision();
    ofsVar->precision(16);
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      *ofsVar << m_logEvidenceFactors[i] << std::endl;
    }
    ofsVar->precision(savedPrecision);
    *ofsVar << "COMPLETE"                << std::endl; // 6 = ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA

    delete ofsVar;
  }
  m_env.fullComm().Barrier();

  //******************************************************************************
  // Write three 'data' files
  //******************************************************************************
  char levelSufix[256];
  sprintf(levelSufix,"%d",m_currLevel+LEVEL_REF_ID); // Yes, '+0'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING chain at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currChain.unifiedWriteContents(m_options.m_restartOutput_baseNameForFiles + "Chain_l" + levelSufix,
                                 m_options.m_restartOutput_fileType);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING like at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogLikelihoodValues.unifiedWriteContents(m_options.m_restartOutput_baseNameForFiles + "LogLike_l" + levelSufix,
                                               m_options.m_restartOutput_fileType);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING target at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogTargetValues.unifiedWriteContents(m_options.m_restartOutput_baseNameForFiles + "LogTarget_l" + levelSufix,
                                           m_options.m_restartOutput_fileType);
  m_env.fullComm().Barrier();

  //******************************************************************************
  // Write 'control' file *with* 'level' spefication in name
  //******************************************************************************
  if (m_env.fullRank() == 0) {
    std::ofstream* ofsVar = new std::ofstream((m_options.m_restartOutput_baseNameForFiles + "Control_l" + levelSufix + ".txt").c_str(),
                                              std::ofstream::out | std::ofstream::trunc);
    *ofsVar << m_currLevel               << std::endl  // 1
            << m_vectorSpace.dimGlobal() << std::endl  // 2
            << currExponent              << std::endl  // 3
            << currEta                   << std::endl  // 4
            << quantity1                 << std::endl; // 5
    unsigned int savedPrecision = ofsVar->precision();
    ofsVar->precision(16);
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      *ofsVar << m_logEvidenceFactors[i] << std::endl;
    }
    ofsVar->precision(savedPrecision);
    *ofsVar << "COMPLETE"                << std::endl; // 6 = ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA

    delete ofsVar;
  }
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING done at level " << m_currLevel
                            << "\n" << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::restartML(
  double&                            currExponent,            // output
  double&                            currEta,                 // output
  SequenceOfVectors<P_V,P_M>& currChain,               // output
  ScalarSequence<double>&     currLogLikelihoodValues, // output
  ScalarSequence<double>&     currLogTargetValues)     // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING initiating at level " << m_currLevel
                            << "\n" << std::endl;
  }

  //******************************************************************************
  // Read 'control' file
  //******************************************************************************
  unsigned int vectorSpaceDim  = 0;
  unsigned int quantity1       = 0;
  std::string  checkingString("");
  if (m_env.fullRank() == 0) {
    std::ifstream* ifsVar = new std::ifstream((m_options.m_restartInput_baseNameForFiles + "Control.txt").c_str(),
                                              std::ifstream::in);

    //******************************************************************************
    // Determine number of lines
    //******************************************************************************
    unsigned int numLines = std::count(std::istreambuf_iterator<char>(*ifsVar),
                                       std::istreambuf_iterator<char>(),
                                       '\n');
    ifsVar->seekg(0,std::ios_base::beg);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "Restart input file has " << numLines
                              << " lines"
                              << std::endl;
    }

    //******************************************************************************
    // Read all values
    //******************************************************************************
    *ifsVar >> m_currLevel; // 1
    queso_require_equal_to_msg(numLines, (ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA + m_currLevel), "number of lines read is different than pre-established number of lines in control file");

    m_logEvidenceFactors.clear();
    m_logEvidenceFactors.resize(m_currLevel,0.);
    *ifsVar >> vectorSpaceDim  // 2
            >> currExponent    // 3
            >> currEta         // 4
            >> quantity1;      // 5
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      *ifsVar >> m_logEvidenceFactors[i];
    }
    *ifsVar >> checkingString; // 6 = ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA
    queso_require_equal_to_msg(checkingString, std::string("COMPLETE"), std::string("control txt input file is not complete"));

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "Restart input file has the following information:"
                              << "\n m_currLevel = "      << m_currLevel
                              << "\n vectorSpaceDim = "   << vectorSpaceDim
                              << "\n currExponent = "     << currExponent
                              << "\n currEta = "          << currEta
                              << "\n quantity1 = "        << quantity1;
      for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
        *m_env.subDisplayFile() << "\n [" << i << "] = " << m_logEvidenceFactors[i];
      }
      *m_env.subDisplayFile() << std::endl;
    }

#if 0 // For debug only
    std::string tmpString;
    for (unsigned int i = 0; i < 2; ++i) {
      *ifsVar >> tmpString;
      std::cout << "Just read '" << tmpString << "'" << std::endl;
    }
    while ((lineId < numLines) && (ifsVar->eof() == false)) {
    }
    ifsVar->ignore(maxCharsPerLine,'\n');
#endif

    delete ifsVar;
  } // if (m_env.fullRank() == 0)
  m_env.fullComm().Barrier();

  //******************************************************************************
  // MPI_Bcast 'm_currLevel'
  //******************************************************************************
  unsigned int tmpUint = (unsigned int) m_currLevel;
  m_env.fullComm().Bcast((void *) &tmpUint, (int) 1, RawValue_MPI_UNSIGNED, 0, // Yes, 'fullComm'
                         "MLSampling<P_V,P_M>::restartML()",
                         "failed MPI.Bcast() for m_currLevel");
  if (m_env.fullRank() != 0) {
    m_currLevel = tmpUint;
  }

  //******************************************************************************
  // MPI_Bcast the rest of the information just read
  //******************************************************************************
  std::vector<double> tmpData(ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA-1+m_currLevel,0.);
  if (m_env.fullRank() == 0) {
    tmpData[0] = vectorSpaceDim;
    tmpData[1] = currExponent;
    tmpData[2] = currEta;
    tmpData[3] = quantity1;
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      tmpData[4+i] = m_logEvidenceFactors[i];
    }
  }
  else {
    m_logEvidenceFactors.clear();
    m_logEvidenceFactors.resize(m_currLevel,0.);
  }
  m_env.fullComm().Bcast((void *) &tmpData[0], (int) tmpData.size(), RawValue_MPI_DOUBLE, 0, // Yes, 'fullComm'
                         "MLSampling<P_V,P_M>::restartML()",
                         "failed MPI.Bcast() for rest of information read from input file");
  if (m_env.fullRank() != 0) {
    vectorSpaceDim = tmpData[0];
    currExponent   = tmpData[1];
    currEta        = tmpData[2];
    quantity1      = tmpData[3];
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      m_logEvidenceFactors[i] = tmpData[4+i];
    }
  }

  //******************************************************************************
  // Process read data in all MPI nodes now
  //******************************************************************************
  queso_require_equal_to_msg(vectorSpaceDim, m_vectorSpace.dimGlobal(), "read vector space dimension is not consistent");
  queso_require_msg(!((currExponent < 0.) || (currExponent > 1.)), "read currExponent is not consistent");
  queso_require_equal_to_msg((quantity1 % m_env.numSubEnvironments()), 0, "read size of chain should be a multiple of the number of subenvironments");
  unsigned int subSequenceSize = 0;
  subSequenceSize = ((double) quantity1) / ((double) m_env.numSubEnvironments());

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Restart input file has the following information"
                            << ": subSequenceSize = " << subSequenceSize
                            << std::endl;
  }

  //******************************************************************************
  // Read three 'data' files
  //******************************************************************************
  char levelSufix[256];
  sprintf(levelSufix,"%d",m_currLevel+LEVEL_REF_ID); // Yes, '+0'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING chain at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currChain.unifiedReadContents(m_options.m_restartInput_baseNameForFiles + "Chain_l" + levelSufix,
                                m_options.m_restartInput_fileType,
                                subSequenceSize);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING like at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogLikelihoodValues.unifiedReadContents(m_options.m_restartInput_baseNameForFiles + "LogLike_l" + levelSufix,
                                              m_options.m_restartInput_fileType,
                                              subSequenceSize);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING target at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogTargetValues.unifiedReadContents(m_options.m_restartInput_baseNameForFiles + "LogTarget_l" + levelSufix,
                                          m_options.m_restartInput_fileType,
                                          subSequenceSize);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING done at level " << m_currLevel
                            << "\n" << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Level0_all(
  const MLSamplingLevelOptions& currOptions,                // input
  unsigned int&                        unifiedRequestedNumSamples, // output
  SequenceOfVectors<P_V,P_M>&   currChain,                  // output
  ScalarSequence<double>&       currLogLikelihoodValues,    // output
  ScalarSequence<double>&       currLogTargetValues)        // output
{
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level "              << m_currLevel+LEVEL_REF_ID
                              << ", currOptions.m_rawChainSize = " << currOptions.m_rawChainSize // Ok to use rawChainSize
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);
    if (iRC) {}; // just to remove compiler warning

    if (m_env.inter0Rank() >= 0) {
      unsigned int tmpSize = currOptions.m_rawChainSize;
      m_env.inter0Comm().template Allreduce<unsigned int>(&tmpSize, &unifiedRequestedNumSamples, (int) 1, RawValue_MPI_SUM,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for requested num samples in level 0");
    }
    else {
      unifiedRequestedNumSamples = currOptions.m_rawChainSize;
    }

    currChain.setName              (currOptions.m_prefix + "rawChain"        );
    currLogLikelihoodValues.setName(currOptions.m_prefix + "rawLogLikelihood");
    currLogTargetValues.setName    (currOptions.m_prefix + "rawLogTarget"    );

    currChain.resizeSequence              (currOptions.m_rawChainSize); // Ok to use rawChainSize
    currLogLikelihoodValues.resizeSequence(currOptions.m_rawChainSize); // Ok to use rawChainSize
    currLogTargetValues.resizeSequence    (currOptions.m_rawChainSize); // Ok to use rawChainSize

    P_V auxVec(m_vectorSpace.zeroVector());
    ScalarFunctionSynchronizer<P_V,P_M> likelihoodSynchronizer(m_likelihoodFunction,auxVec); // prudencio 2010-08-01
    for (unsigned int i = 0; i < currChain.subSequenceSize(); ++i) {
      //std::cout << "In QUESO: before prior realizer with i = " << i << std::endl;
      bool outOfSupport = true;
      do {

        m_priorRv.realizer().realization(auxVec);  // gpmsa2
        if (m_numDisabledParameters > 0) { // gpmsa2
          unsigned int disabledCounter = 0;
          for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
            if (m_parameterEnabledStatus[paramId] == false) {
              auxVec[paramId] = currOptions.m_initialValuesOfDisabledParameters[disabledCounter];
              disabledCounter++;
            }
          }
        }
        auxVec.mpiBcast(0, m_env.subComm()); // prudencio 2010-08-01

        outOfSupport = !(m_targetDomain->contains(auxVec));
      } while (outOfSupport); // prudenci 2011-Oct-04

      currChain.setPositionValues(i,auxVec);
      // KAUST: all nodes should call likelihood
#if 1 // prudencio 2010-08-01
      currLogLikelihoodValues[i] = likelihoodSynchronizer.callFunction(&auxVec,NULL,NULL,NULL,NULL,NULL,NULL); // likelihood is important
#else
      currLogLikelihoodValues[i] = m_likelihoodFunction.lnValue(auxVec,NULL,NULL,NULL,NULL); // likelihood is important
#endif
      currLogTargetValues[i]     = m_priorRv.pdf().lnValue(auxVec,NULL,NULL,NULL,NULL) + currLogLikelihoodValues[i];
      //std::cout << "In QUESO: currLogTargetValues[" << i << "] = " << currLogTargetValues[i] << std::endl;
    }

    if (m_env.inter0Rank() >= 0) { // KAUST
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      if (currOptions.m_rawChainComputeStats) {
        FilePtrSetStruct filePtrSet;
        m_env.openOutputFile(currOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                             currOptions.m_dataOutputAllowedSet,
                             false,
                             filePtrSet);

        //m_env.syncPrintDebugMsg("At level 0, calling computeStatistics for chain",1,10,m_env.inter0Comm()); // output debug
        currChain.computeStatistics(*currOptions.m_rawChainStatisticalOptionsObj,
                                    filePtrSet.ofsVar);

        m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
      }
      // Compute MLE and MAP
      // rr0
#endif
      if (currOptions.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
        currChain.unifiedWriteContents              (currOptions.m_rawChainDataOutputFileName,
                                                     currOptions.m_rawChainDataOutputFileType);
        currLogLikelihoodValues.unifiedWriteContents(currOptions.m_rawChainDataOutputFileName,
                                                     currOptions.m_rawChainDataOutputFileType);
        currLogTargetValues.unifiedWriteContents    (currOptions.m_rawChainDataOutputFileName,
                                                     currOptions.m_rawChainDataOutputFileType);
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;

        //unsigned int numZeros = 0;
        //for (unsigned int i = 0; i < currTargetValues.subSequenceSize(); ++i) {
        //  *m_env.subDisplayFile() << "currTargetValues[" << i
        //                          << "] = " << currTargetValues[i]
        //                          << std::endl;
        //  if (currTargetValues[i] == 0.) numZeros++;
        //}
        //*m_env.subDisplayFile() << "Number of zeros in currTargetValues = " << numZeros
        //                        << std::endl;
      }

      if (currOptions.m_filteredChainGenerate) {
        // todo
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        if (currOptions.m_filteredChainComputeStats) {
          // todo

          //currChain.computeStatistics(*currOptions.m_filteredChainStatisticalOptionsObj,
          //                            filePtrSet.ofsVar);
        }
        // Compute MLE and MAP
        // rr0
#endif
      }

    } // KAUST

    queso_require_equal_to_msg(currChain.subSequenceSize(), currOptions.m_rawChainSize, "currChain (first one) has been generated with invalid size");

    double levelRunTime = MiscGetEllapsedSeconds(&timevalLevel);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << m_currLevel+LEVEL_REF_ID
                              << ", total level time = " << levelRunTime << " seconds"
                              << std::endl;
    }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step01_inter0(
  const MLSamplingLevelOptions* currOptions,                // input
  unsigned int&                        unifiedRequestedNumSamples) // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 1 of 11"
                                << std::endl;
      }

      unsigned int tmpSize = currOptions->m_rawChainSize;
      // This computed 'unifiedRequestedNumSamples' needs to be recomputed only at the last
      // level, when 'currOptions' is replaced by 'lastLevelOptions' (see step 3 of 11)
      m_env.inter0Comm().template Allreduce<unsigned int>(&tmpSize, &unifiedRequestedNumSamples, (int) 1, RawValue_MPI_SUM,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for requested num samples in step 1");

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "KEY In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", currOptions->m_rawChainSize = " << currOptions->m_rawChainSize // Ok to use rawChainSize
                                << std::endl;
      }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step02_inter0(
  const MLSamplingLevelOptions* currOptions,             // input
  SequenceOfVectors<P_V,P_M>&   currChain,               // input/output
  ScalarSequence<double>&       currLogLikelihoodValues, // input/output
  ScalarSequence<double>&       currLogTargetValues,     // input/output
  SequenceOfVectors<P_V,P_M>&   prevChain,               // output
  ScalarSequence<double>&       prevLogLikelihoodValues, // output
  ScalarSequence<double>&       prevLogTargetValues,     // output
  unsigned int&                        indexOfFirstWeight,      // output
  unsigned int&                        indexOfLastWeight)       // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 2 of 11"
                                << std::endl;
      }

      prevChain = currChain;
      currChain.clear();
      currChain.setName(currOptions->m_prefix + "rawChain");

      prevLogLikelihoodValues = currLogLikelihoodValues; // likelihood is important
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", prevLogLikelihoodValues[0] = " << prevLogLikelihoodValues[0]
                                << std::endl;
      }
      prevLogTargetValues     = currLogTargetValues;

      currLogLikelihoodValues.clear();
      currLogLikelihoodValues.setName(currOptions->m_prefix + "rawLogLikelihood");

      currLogTargetValues.clear();
      currLogTargetValues.setName(currOptions->m_prefix + "rawLogTarget");

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        P_V prevPosition(m_vectorSpace.zeroVector());
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ":"
                                << std::endl;
        for (unsigned int i = 0; i < prevChain.subSequenceSize(); ++i) {
          prevChain.getPositionValues(i,prevPosition);
          *m_env.subDisplayFile() << "  prevChain[" << i
                                  << "] = " << prevPosition
                                  << ", prevLogLikelihoodValues[" << i
                                  << "] = " << prevLogLikelihoodValues[i]
                                  << ", prevLogTargetValues[" << i
                                  << "] = " << prevLogTargetValues[i]
                                  << std::endl;
        }
      }
#endif

      unsigned int quantity1 = prevChain.unifiedSequenceSize();
      unsigned int quantity2 = currChain.unifiedSequenceSize();
      unsigned int quantity3 = prevLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity4 = currLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity5 = prevLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity6 = currLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": prevChain.unifiedSequenceSize() = " << quantity1
                                << ", currChain.unifiedSequenceSize() = " << quantity2
                                << ", prevLogLikelihoodValues.unifiedSequenceSize() = " << quantity3
                                << ", currLogLikelihoodValues.unifiedSequenceSize() = " << quantity4
                                << ", prevLogTargetValues.unifiedSequenceSize() = " << quantity5
                                << ", currLogTargetValues.unifiedSequenceSize() = " << quantity6
                                << std::endl;
      }

      queso_require_equal_to_msg(prevChain.subSequenceSize(), prevLogLikelihoodValues.subSequenceSize(), "different sizes between previous chain and previous sequence of likelihood values");

      queso_require_equal_to_msg(prevChain.subSequenceSize(), prevLogTargetValues.subSequenceSize(), "different sizes between previous chain and previous sequence of target values");

      // Set 'indexOfFirstWeight' and 'indexOfLastWeight' // KAUST
      indexOfFirstWeight = 0;
      indexOfLastWeight  = indexOfFirstWeight + prevChain.subSequenceSize()-1;
      {
        //std::cout << "m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc() << std::endl;
        int r = m_env.inter0Rank();
        //std::cout << "r = " << r << std::endl;
        m_env.inter0Comm().Barrier();
        unsigned int auxUint = 0;
        if (r > 0) {
          RawType_MPI_Status status;
    //std::cout << "Rank " << r << " is entering MPI_Recv()" << std::endl;
          m_env.inter0Comm().Recv((void*) &auxUint, 1, RawValue_MPI_UNSIGNED, r-1, r-1, &status,
                                  "MLSampling<P_V,P_M>::generateSequence()",
                                  "failed MPI.Recv()");
    //std::cout << "Rank " << r << " received auxUint = " << auxUint << std::endl;
          indexOfFirstWeight = auxUint;
          indexOfLastWeight = indexOfFirstWeight + prevChain.subSequenceSize()-1;
        }
        if (r < (m_env.inter0Comm().NumProc()-1)) {
          auxUint = indexOfLastWeight + 1;
    //std::cout << "Rank " << r << " is sending auxUint = " << auxUint << std::endl;
          m_env.inter0Comm().Send((void*) &auxUint, 1, RawValue_MPI_UNSIGNED, r+1, r,
                                  "MLSampling<P_V,P_M>::generateSequence()",
                                  "failed MPI.Send()");
    //std::cout << "Rank " << r << " sent auxUint = " << auxUint << std::endl;
        }
        m_env.inter0Comm().Barrier();
      }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step03_inter0(
  const MLSamplingLevelOptions* currOptions,             // input
  const ScalarSequence<double>& prevLogLikelihoodValues, // input
  double                               prevExponent,            // input
  double                               failedExponent,          // input // gpmsa1
  double&                              currExponent,            // output
  ScalarSequence<double>&       weightSequence)          // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", failedExponent = " << failedExponent // gpmsa1
                                << ": beginning step 3 of 11"
                                << std::endl;
      }

      std::vector<double> exponents(2,0.);
      exponents[0] = prevExponent;
      exponents[1] = 1.;

      double nowExponent = 1.; // Try '1.' right away
      double nowEffectiveSizeRatio = 0.; // To be computed

      unsigned int nowAttempt = 0;
      bool testResult = false;
      double meanEffectiveSizeRatio = .5*(currOptions->m_minEffectiveSizeRatio + currOptions->m_maxEffectiveSizeRatio);
      ScalarSequence<double> omegaLnDiffSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");

      double nowUnifiedEvidenceLnFactor = 0.;
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ", failedExponent = " << failedExponent // gpmsa1
                                  << ": entering loop for computing next exponent"
                                  << ", with nowAttempt = " << nowAttempt
                                  << std::endl;
        }

        if (failedExponent > 0.) { // gpmsa1
          nowExponent = .5*(prevExponent+failedExponent);
        }
        else {
          if (nowAttempt > 0) {
            if (nowEffectiveSizeRatio > meanEffectiveSizeRatio) {
              exponents[0] = nowExponent;
            }
            else {
              exponents[1] = nowExponent;
            }
            nowExponent = .5*(exponents[0] + exponents[1]);
          }
        }
        double auxExponent = nowExponent;
        if (prevExponent != 0.) {
          auxExponent /= prevExponent;
          auxExponent -= 1.;
        }
        double subWeightRatioSum     = 0.;
        double unifiedWeightRatioSum = 0.;

        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] = prevLogLikelihoodValues[i]*auxExponent; // likelihood is important
        }

#if 1 // prudenci-2012-07-06
      //double unifiedOmegaLnMin = omegaLnDiffSequence.unifiedMinPlain(m_vectorSpace.numOfProcsForStorage() == 1);
        double unifiedOmegaLnMax = omegaLnDiffSequence.unifiedMaxPlain(m_vectorSpace.numOfProcsForStorage() == 1);
#else
        double unifiedOmegaLnMin = 0.;
        double unifiedOmegaLnMax = 0.;
        omegaLnDiffSequence.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1, // KAUST3
                                               0,
                                               omegaLnDiffSequence.subSequenceSize(),
                                               unifiedOmegaLnMin,
                                               unifiedOmegaLnMax);
#endif
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] -= unifiedOmegaLnMax;
          weightSequence[i] = exp(omegaLnDiffSequence[i]);
          subWeightRatioSum += weightSequence[i];
#if 0 // For debug only
          if ((m_currLevel == 1) && (nowAttempt == 6))  {
            if (m_env.subDisplayFile() && (m_env.displayVerbosity() >= 99)) {
              *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                      << ", level "                        << m_currLevel+LEVEL_REF_ID
                                      << ", step "                         << m_currStep
                                      << ", i = "                          << i
                                      << ", prevLogLikelihoodValues[i] = " << prevLogLikelihoodValues[i]
                                      << ", omegaLnDiffSequence[i] = "     << omegaLnDiffSequence[i]
                                      << ", weightSequence[i] = "          << weightSequence[i]
    //<< ", subWeightRatioSum = "          << subWeightRatioSum
                                      << std::endl;
            }
          }
#endif
        }
        m_env.inter0Comm().template Allreduce<double>(&subWeightRatioSum, &unifiedWeightRatioSum, (int) 1, RawValue_MPI_SUM,
                                     "MLSampling<P_V,P_M>::generateSequence()",
                                     "failed MPI.Allreduce() for weight ratio sum");

        unsigned int auxQuantity = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
        nowUnifiedEvidenceLnFactor = log(unifiedWeightRatioSum) + unifiedOmegaLnMax - log(auxQuantity);

        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] /= unifiedWeightRatioSum;
          effectiveSampleSize += weightSequence[i]*weightSequence[i];
          //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          //  *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
          //                          << ", level "                 << m_currLevel+LEVEL_REF_ID
          //                          << ", step "                  << m_currStep
          //                          << ": i = "                   << i
          //                          << ", effectiveSampleSize = " << effectiveSampleSize
          //                          << std::endl;
          //}
        }

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                  << ", step "                                   << m_currStep
                                  << ": nowAttempt = "                           << nowAttempt
                                  << ", prevExponent = "                         << prevExponent
                                  << ", exponents[0] = "                         << exponents[0]
                                  << ", nowExponent = "                          << nowExponent
                                  << ", exponents[1] = "                         << exponents[1]
                                  << ", subWeightRatioSum = "                    << subWeightRatioSum
                                  << ", unifiedWeightRatioSum = "                << unifiedWeightRatioSum
                                  << ", unifiedOmegaLnMax = "                    << unifiedOmegaLnMax
                                  << ", weightSequence.unifiedSequenceSize() = " << auxQuantity
                                  << ", nowUnifiedEvidenceLnFactor = "           << nowUnifiedEvidenceLnFactor
                                  << ", effectiveSampleSize = "                  << effectiveSampleSize
                                  << std::endl;
        }

#if 0 // For debug only
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ":"
                                  << std::endl;
        }
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "  weightSequence[" << i
                                    << "] = "              << weightSequence[i]
                                    << std::endl;
          }
        }
#endif

        double subQuantity = effectiveSampleSize;
        effectiveSampleSize = 0.;
        m_env.inter0Comm().template Allreduce<double>(&subQuantity, &effectiveSampleSize, (int) 1, RawValue_MPI_SUM,
                                     "MLSampling<P_V,P_M>::generateSequence()",
                                     "failed MPI.Allreduce() for effective sample size");

        effectiveSampleSize = 1./effectiveSampleSize;
        nowEffectiveSizeRatio = effectiveSampleSize/((double) weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1));
        queso_require_less_equal_msg(nowEffectiveSizeRatio, (1.+1.e-8), "effective sample size ratio cannot be > 1");

        //                    m_env.worldRank(),
        //                    "MLSampling<P_V,P_M>::generateSequence()",
        //                    "effective sample size ratio cannot be < 1");

        if (failedExponent > 0.) { // gpmsa1
          testResult = true;
        }
        else {
          //bool aux1 = (nowEffectiveSizeRatio == meanEffectiveSizeRatio);
          bool aux2 = (nowExponent == 1.                             )
                      &&
                      (nowEffectiveSizeRatio > meanEffectiveSizeRatio);
          bool aux3 = (nowEffectiveSizeRatio >= currOptions->m_minEffectiveSizeRatio)
                      &&
                      (nowEffectiveSizeRatio <= currOptions->m_maxEffectiveSizeRatio);
          testResult = aux2 || aux3;
        }

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "                   << m_currLevel+LEVEL_REF_ID
                                  << ", step "                    << m_currStep
                                  << ": nowAttempt = "            << nowAttempt
                                  << ", prevExponent = "          << prevExponent
                                  << ", failedExponent = "        << failedExponent // gpmsa1
                                  << ", exponents[0] = "          << exponents[0]
                                  << ", nowExponent = "           << nowExponent
                                  << ", exponents[1] = "          << exponents[1]
                                  << ", effectiveSampleSize = "   << effectiveSampleSize
                                  << ", weightSequenceSize = "    << weightSequence.subSequenceSize()
                                  << ", minEffectiveSizeRatio = " << currOptions->m_minEffectiveSizeRatio
                                  << ", nowEffectiveSizeRatio = " << nowEffectiveSizeRatio
                                  << ", maxEffectiveSizeRatio = " << currOptions->m_maxEffectiveSizeRatio
      //<< ", aux2 = "                  << aux2
      //<< ", aux3 = "                  << aux3
                                  << ", testResult = "            << testResult
                                  << std::endl;
        }
        nowAttempt++;

        // Make sure all nodes in 'inter0Comm' have the same value of 'nowExponent'
        if (MiscCheckForSameValueInAllNodes(nowExponent,
                                              0., // kept 'zero' on 2010/03/05
                                              m_env.inter0Comm(),
                                              "MLSampling<P_V,P_M>::generateSequence(), step 3, nowExponent") == false) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "WARNING, In MLSampling<P_V,P_M>::generateSequence()"
                                    << ", level "        << m_currLevel+LEVEL_REF_ID
                                    << ", step "         << m_currStep
                                    << ": nowAttempt = " << nowAttempt
                                    << ", MiscCheck for 'nowExponent' detected a problem"
                                    << std::endl;
          }
        }

        // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
        if (MiscCheckForSameValueInAllNodes(testResult,
                                              0., // kept 'zero' on 2010/03/05
                                              m_env.inter0Comm(),
                                              "MLSampling<P_V,P_M>::generateSequence(), step 3, testResult") == false) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "WARNING, In MLSampling<P_V,P_M>::generateSequence()"
                                    << ", level "        << m_currLevel+LEVEL_REF_ID
                                    << ", step "         << m_currStep
                                    << ": nowAttempt = " << nowAttempt
                                    << ", MiscCheck for 'testResult' detected a problem"
                                    << std::endl;
          }
        }
      } while (testResult == false);
      currExponent = nowExponent;
      if (failedExponent > 0.) { // gpmsa1
        m_logEvidenceFactors[m_logEvidenceFactors.size()-1] = nowUnifiedEvidenceLnFactor;
      }
      else {
        m_logEvidenceFactors.push_back(nowUnifiedEvidenceLnFactor); // restart
      }

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                << ", step "                                   << m_currStep
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", failedExponent = "                       << failedExponent // gpmsa1
                                << ", currExponent = "                         << currExponent
                                << ", effective ratio = "                      << nowEffectiveSizeRatio
                                << ", log(evidence factor) = "                 << m_logEvidenceFactors[m_logEvidenceFactors.size()-1]
                                << ", evidence factor = "                      << exp(m_logEvidenceFactors[m_logEvidenceFactors.size()-1])
                                << std::endl;

        //unsigned int numZeros = 0;
        //for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        //  *m_env.subDisplayFile() << "weightSequence[" << i
        //                          << "] = " << weightSequence[i]
        //                         << std::endl;
        //  if (weightSequence[i] == 0.) numZeros++;
        //}
        //*m_env.subDisplayFile() << "Number of zeros in weightSequence = " << numZeros
        //                        << std::endl;
      }

      // Make sure all nodes in 'inter0Comm' have the same value of 'logEvidenceFactor'
      if (MiscCheckForSameValueInAllNodes(m_logEvidenceFactors[m_logEvidenceFactors.size()-1],
                                            3.0e-16, // changed from 'zero' on 2010/03/03
                                            m_env.inter0Comm(),
                                            "MLSampling<P_V,P_M>::generateSequence(), step 3, logEvidenceFactor") == false) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "WARNING, In MLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "        << m_currLevel+LEVEL_REF_ID
                                  << ", step "         << m_currStep
                                  << ", failedExponent = " << failedExponent // gpmsa1
                                  << ": nowAttempt = " << nowAttempt
                                  << ", MiscCheck for 'logEvidenceFactor' detected a problem"
                                  << std::endl;
        }
      }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", failedExponent = " << failedExponent // gpmsa
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step04_inter0(
  const SequenceOfVectors<P_V,P_M>& prevChain,        // input
  const ScalarSequence<double>&     weightSequence,   // input
  P_M&                                     unifiedCovMatrix) // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 4 of 11"
                                << std::endl;
      }

      P_V auxVec(m_vectorSpace.zeroVector());
      P_V subWeightedMeanVec(m_vectorSpace.zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        subWeightedMeanVec += weightSequence[i]*auxVec;
      }

      // Todd Oliver 2010-09-07: compute weighted mean over all processors
      P_V unifiedWeightedMeanVec(m_vectorSpace.zeroVector());
      if (m_env.inter0Rank() >= 0) {
        subWeightedMeanVec.mpiAllReduce(RawValue_MPI_SUM,m_env.inter0Comm(),unifiedWeightedMeanVec);
      }
      else {
        unifiedWeightedMeanVec = subWeightedMeanVec;
      }

      P_V diffVec(m_vectorSpace.zeroVector());
      P_M subCovMatrix(m_vectorSpace.zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        diffVec = auxVec - unifiedWeightedMeanVec;
        subCovMatrix += weightSequence[i]*matrixProduct(diffVec,diffVec);
      }

      for (unsigned int i = 0; i < unifiedCovMatrix.numRowsLocal(); ++i) { // KAUST5
        for (unsigned int j = 0; j < unifiedCovMatrix.numCols(); ++j) {
          double localValue = subCovMatrix(i,j);
          double sumValue = 0.;
          if (m_env.inter0Rank() >= 0) {
            m_env.inter0Comm().template Allreduce<double>(&localValue, &sumValue, (int) 1, RawValue_MPI_SUM,
                                         "MLSampling<P_V,P_M>::generateSequence()",
                                         "failed MPI.Allreduce() for cov matrix");
          }
          else {
            sumValue = localValue;
          }
          unifiedCovMatrix(i,j) = sumValue;
        }
      }

      if (m_numDisabledParameters > 0) { // gpmsa2
        for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
          if (m_parameterEnabledStatus[paramId] == false) {
            for (unsigned int i = 0; i < m_vectorSpace.dimLocal(); ++i) {
              unifiedCovMatrix(i,paramId) = 0.;
            }
            for (unsigned int j = 0; j < m_vectorSpace.dimLocal(); ++j) {
              unifiedCovMatrix(paramId,j) = 0.;
            }
            unifiedCovMatrix(paramId,paramId) = 1.;
          }
        }
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level "              << m_currLevel+LEVEL_REF_ID
                                << ", step "               << m_currStep
                                << ": unifiedCovMatrix = " << unifiedCovMatrix
                                << std::endl;
      }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step05_inter0(
  unsigned int                         unifiedRequestedNumSamples,        // input
  const ScalarSequence<double>& weightSequence,                    // input
  std::vector<unsigned int>&           unifiedIndexCountersAtProc0Only,   // output
  std::vector<double>&                 unifiedWeightStdVectorAtProc0Only) // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 5 of 11"
                                << std::endl;
      }

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", before weightSequence.getUnifiedContentsAtProc0Only()"
                                << ":"
                                << std::endl;
      }
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << ", weightSequence[" << i
                                  << "] = "              << weightSequence[i]
                                  << std::endl;
        }
      }
#endif

      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", after weightSequence.getUnifiedContentsAtProc0Only()"
                                << ":"
                                << std::endl;
      }
      for (unsigned int i = 0; i < unifiedWeightStdVectorAtProc0Only.size(); ++i) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "  unifiedWeightStdVectorAtProc0Only[" << i
                                  << "] = "                                 << unifiedWeightStdVectorAtProc0Only[i]
                                  << std::endl;
        }
      }
#endif
      sampleIndexes_proc0(unifiedRequestedNumSamples,        // input
                          unifiedWeightStdVectorAtProc0Only, // input
                          unifiedIndexCountersAtProc0Only);  // output

      unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if (m_env.inter0Rank() == 0) {
        queso_require_equal_to_msg(unifiedIndexCountersAtProc0Only.size(), auxUnifiedSize, "wrong output from sampleIndexesAtProc0() in step 5");
      }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step06_all(
  const MLSamplingLevelOptions* currOptions,                     // input
  unsigned int                         indexOfFirstWeight,              // input
  unsigned int                         indexOfLastWeight,               // input
  const std::vector<unsigned int>&     unifiedIndexCountersAtProc0Only, // input
  bool&                                useBalancedChains,               // output
  std::vector<ExchangeInfoStruct>&   exchangeStdVec)                  // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

  useBalancedChains = decideOnBalancedChains_all(currOptions,                     // input
                                                 indexOfFirstWeight,              // input
                                                 indexOfLastWeight,               // input
                                                 unifiedIndexCountersAtProc0Only, // input
                                                 exchangeStdVec);                 // output

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step07_inter0(
  bool                                      useBalancedChains,               // input
  unsigned int                              indexOfFirstWeight,              // input
  unsigned int                              indexOfLastWeight,               // input
  const std::vector<unsigned int>&          unifiedIndexCountersAtProc0Only, // input
  UnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,           // (possible) output
  const MLSamplingLevelOptions*      currOptions,                     // input
  const SequenceOfVectors<P_V,P_M>&  prevChain,                       // input
  double                                  prevExponent,               // input
  double                                  currExponent,               // input
  const ScalarSequence<double>&           prevLogLikelihoodValues,    // input
  const ScalarSequence<double>&           prevLogTargetValues,        // input
  std::vector<ExchangeInfoStruct>&        exchangeStdVec,                  // (possible) input/output
  BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl)             // (possible) output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 7 of 11"
                                << std::endl;
      }

      if (useBalancedChains) {
        prepareBalLinkedChains_inter0(currOptions,                     // input
                                      prevChain,                       // input
                                      prevExponent,                    // input
                                      currExponent,                    // input
                                      prevLogLikelihoodValues,         // input
                                      prevLogTargetValues,             // input
                                      exchangeStdVec,                  // input/output
                                      balancedLinkControl);            // output
      }
      else {
        prepareUnbLinkedChains_inter0(indexOfFirstWeight,              // input
                                      indexOfLastWeight,               // input
                                      unifiedIndexCountersAtProc0Only, // input
                                      unbalancedLinkControl);          // output
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": balancedLinkControl.balLinkedChains.size() = "   << balancedLinkControl.balLinkedChains.size()
                                << ", unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                                << std::endl;
      }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step08_all(
  BayesianJointPdf<P_V,P_M>& currPdf, // input/output
  GenericVectorRV<P_V,P_M>&  currRv)  // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 8 of 11"
                                << std::endl;
      }

      currRv.setPdf(currPdf);

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step09_all(
  const SequenceOfVectors<P_V,P_M>& prevChain,                         // input
  double                                   prevExponent,                      // input
  double                                   currExponent,                      // input
  const ScalarSequence<double>&            prevLogLikelihoodValues,           // input
  const ScalarSequence<double>&            prevLogTargetValues,               // input
  unsigned int                             indexOfFirstWeight,                // input
  unsigned int                             indexOfLastWeight,                 // input
  const std::vector<double>&               unifiedWeightStdVectorAtProc0Only, // input
  const ScalarSequence<double>&     weightSequence,                    // input
  double                                   prevEta,                           // input
  const GenericVectorRV<P_V,P_M>&   currRv,                            // input
  MLSamplingLevelOptions*           currOptions,                       // input (changed temporarily internally)
  P_M&                                     unifiedCovMatrix,                  // input/output
  double&                                  currEta)                           // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

    if (currOptions->m_scaleCovMatrix == false) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": skipping step 9 of 11"
                                << std::endl;
      }
    }
    else {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 9 of 11"
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
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
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
          queso_error_msg("nowRejectionRate should be out of the requested range at this point of the logic");
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
                queso_error_msg("before and now range flags are inconsistent");
              }
            } // if (useMiddlePointLogicForEta == false)

            beforeEta                       = nowEta;
            beforeRejectionRate             = nowRejectionRate;
            beforeRejectionRateIsBelowRange = nowRejectionRateIsBelowRange;
            if (useMiddlePointLogicForEta == false) {
              if (beforeRejectionRateIsBelowRange) nowEta *= 4.;
              else                                 nowEta /= 4.;
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
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
                *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
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

        // prudencio 2010-12-09: logic 'originalSubNumSamples += 1' added because of the difference of results between GNU and INTEL compiled codes
        double       doubSubNumSamples     = (1.-meanRejectionRate)/meanRejectionRate/currOptions->m_covRejectionRate/currOptions->m_covRejectionRate; // e.g. 19.99...; or 20.0; or 20.1; or 20.9
        unsigned int originalSubNumSamples = 1 + (unsigned int) (doubSubNumSamples); // e.g. 20; or 21; or 21; or 21
        double       auxDouble             = (double) originalSubNumSamples; // e.g. 20.0; or 21.0; or 21.0; or 21.0
        if ((auxDouble - doubSubNumSamples) < 1.e-8) { // e.g. 0.00...01; or 1.0; or 0.9; or 0.1
          originalSubNumSamples += 1;
        }

        if (m_env.inter0Rank() >= 0) { // KAUST
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                    << ", level " << m_currLevel+LEVEL_REF_ID
                                    << ", step "  << m_currStep
                                    << ": in loop for assessing rejection rate"
                                    << ", about to sample "     << originalSubNumSamples << " indexes"
                                    << ", meanRejectionRate = " << meanRejectionRate
                                    << ", covRejectionRate = "  << currOptions->m_covRejectionRate
                                    << std::endl;
          }
        } // KAUST

        std::vector<unsigned int> nowUnifiedIndexCountersAtProc0Only(0); // It will be resized by 'sampleIndexes_proc0()' below
        if (m_env.inter0Rank() >= 0) { // KAUST
          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          sampleIndexes_proc0(tmpUnifiedNumSamples,                // input
                              unifiedWeightStdVectorAtProc0Only,   // input
                              nowUnifiedIndexCountersAtProc0Only); // output

          unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
          if (m_env.inter0Rank() == 0) {
            queso_require_equal_to_msg(nowUnifiedIndexCountersAtProc0Only.size(), auxUnifiedSize, "wrong output from sampleIndexesAtProc0() in step 9");
          }

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                    << ", level " << m_currLevel+LEVEL_REF_ID
                                    << ", step "  << m_currStep
                                    << ": in loop for assessing rejection rate"
                                    << ", about to distribute sampled assessment indexes"
                                    << std::endl;
          }
        } // KAUST

        std::vector<ExchangeInfoStruct>        exchangeStdVec(0);
        BalancedLinkedChainsPerNodeStruct<P_V> nowBalLinkControl;
        UnbalancedLinkedChainsPerNodeStruct    nowUnbLinkControl; // KAUST

        // All processors should call this routine in order to have the same decision value
        bool useBalancedChains = decideOnBalancedChains_all(currOptions,                        // input
                                                            indexOfFirstWeight,                 // input
                                                            indexOfLastWeight,                  // input
                                                            nowUnifiedIndexCountersAtProc0Only, // input
                                                            exchangeStdVec);                    // output

        if (m_env.inter0Rank() >= 0) { // KAUST
          if (useBalancedChains) {
            prepareBalLinkedChains_inter0(currOptions,                        // input
                                          prevChain,                          // input
                                          prevExponent,                       // input
                                          currExponent,                       // input
                                          prevLogLikelihoodValues,            // input
                                          prevLogTargetValues,                // input
                                          exchangeStdVec,                     // input/output
                                          nowBalLinkControl);                 // output
          }
          else {
            prepareUnbLinkedChains_inter0(indexOfFirstWeight,                 // input
                                          indexOfLastWeight,                  // input
                                          nowUnifiedIndexCountersAtProc0Only, // input
                                          nowUnbLinkControl);                 // output
          }
        } // KAUST

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ": in loop for assessing rejection rate"
                                  << ", about to generate assessment chain"
                                  << std::endl;
        }

        SequenceOfVectors<P_V,P_M> nowChain(m_vectorSpace,
                                                   0,
                                                   m_options.m_prefix+"now_chain");
        double       nowRunTime    = 0.;
        unsigned int nowRejections = 0;

        // KAUST: all nodes should call here
        bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
        unsigned int savedRawChainSize          = currOptions->m_rawChainSize; // Ok to use rawChainSize
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
#endif
        bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;
        unsigned int savedDrMaxNumExtraStages   = currOptions->m_drMaxNumExtraStages;
        unsigned int savedAmAdaptInterval       = currOptions->m_amAdaptInterval;

        currOptions->m_totallyMute = true;
        if (m_env.displayVerbosity() >= 999999) {
          currOptions->m_totallyMute = false;
        }
        currOptions->m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        currOptions->m_rawChainComputeStats  = false;
#endif
        currOptions->m_filteredChainGenerate = false;
        currOptions->m_drMaxNumExtraStages   = 0;
        currOptions->m_amAdaptInterval       = 0;

        // KAUST: all nodes in 'subComm' should call here, important
        if (useBalancedChains) {
          generateBalLinkedChains_all(*currOptions,       // input, only m_rawChainSize changes
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
          generateUnbLinkedChains_all(*currOptions,       // input, only m_rawChainSize changes
                                      nowCovMatrix,       // input
                                      currRv,             // input
                                      nowUnbLinkControl,  // input // Round Rock
                                      indexOfFirstWeight, // input // Round Rock
                                      prevChain,          // input // Round Rock
                                      prevExponent,             // input
                                      currExponent,             // input
                                      prevLogLikelihoodValues,  // input
                                      prevLogTargetValues,      // input
                                      nowChain,           // output
                                      nowRunTime,         // output
                                      nowRejections,      // output
                                      NULL,               // output
                                      NULL);              // output
        }

        // KAUST: all nodes should call here
        currOptions->m_totallyMute           = savedTotallyMute;
        currOptions->m_rawChainSize          = savedRawChainSize;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
#endif
        currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME
        currOptions->m_drMaxNumExtraStages   = savedDrMaxNumExtraStages;
        currOptions->m_amAdaptInterval       = savedAmAdaptInterval;

        for (unsigned int i = 0; i < nowBalLinkControl.balLinkedChains.size(); ++i) {
          queso_require_msg(nowBalLinkControl.balLinkedChains[i].initialPosition, "Initial position pointer in step 9 should not be NULL");
          delete nowBalLinkControl.balLinkedChains[i].initialPosition;
          nowBalLinkControl.balLinkedChains[i].initialPosition = NULL;
        }
        nowBalLinkControl.balLinkedChains.clear();

        if (m_env.inter0Rank() >= 0) { // KAUST
          // If only one cov matrix is used, then the rejection should be assessed among all inter0Comm nodes // KAUST3
          unsigned int nowUnifiedRejections = 0;
          m_env.inter0Comm().template Allreduce<unsigned int>(&nowRejections, &nowUnifiedRejections, (int) 1, RawValue_MPI_SUM,
                                       "MLSampling<P_V,P_M>::generateSequence_Step09_all()",
                                       "failed MPI.Allreduce() for now rejections");

#if 0 // Round Rock 2009 12 29
          unsigned int tmpUnifiedNumSamples = 0;
          m_env.inter0Comm().Allreduce((void *) &tmpSubNumSamples, (void *) &tmpUnifiedNumSamples, (int) 1, RawValue_MPI_UNSIGNED, RawValue_MPI_SUM,
                                       "MLSampling<P_V,P_M>::generateSequence_Step09_all()",
                                       "failed MPI.Allreduce() for num samples in step 9");
#endif

          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          nowRejectionRate = ((double) nowUnifiedRejections) / ((double) tmpUnifiedNumSamples);

          //bool aux1 = (nowRejectionRate == meanRejectionRate);
          bool aux2 = (nowRejectionRate >= currOptions->m_minRejectionRate)
                      &&
                      (nowRejectionRate <= currOptions->m_maxRejectionRate);
          testResult = aux2;

          // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
          if (MiscCheckForSameValueInAllNodes(testResult,
                                                0., // kept 'zero' on 2010/03/03
                                                m_env.inter0Comm(),
                                                "MLSampling<P_V,P_M>::generateSequence_Step09_all(), step 9, testResult") == false) {
            if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
              *m_env.subDisplayFile() << "WARNING, In MLSampling<P_V,P_M>::generateSequence()"
                                      << ", level "        << m_currLevel+LEVEL_REF_ID
                                      << ", step "         << m_currStep
                                      << ": nowAttempt = " << nowAttempt
                                      << ", MiscCheck for 'testResult' detected a problem"
                                      << std::endl;
            }
    }
        } // if (m_env.inter0Rank() >= 0) { // KAUST

        // KAUST: all nodes in 'subComm' should have the same 'testResult'
        unsigned int tmpUint = (unsigned int) testResult;
        m_env.subComm().Bcast((void *) &tmpUint, (int) 1, RawValue_MPI_UNSIGNED, 0, // Yes, 'subComm', important
                              "MLSampling<P_V,P_M>::generateSequence_Step09_all()",
                              "failed MPI.Bcast() for testResult");
        testResult = (bool) tmpUint;

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
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
          if (MiscCheckForSameValueInAllNodes(nowEta,
                                                1.0e-16, // changed from 'zero' on 2009/11/dd
                                                m_env.inter0Comm(),
                                                "MLSampling<P_V,P_M>::generateSequence_Step09_all(), step 9, nowEta") == false) {
            if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
              *m_env.subDisplayFile() << "WARNING, In MLSampling<P_V,P_M>::generateSequence()"
                                      << ", level "        << m_currLevel+LEVEL_REF_ID
                                      << ", step "         << m_currStep
                                      << ": nowAttempt = " << nowAttempt
                                      << ", MiscCheck for 'nowEta' detected a problem"
                                      << std::endl;
            }
          }
        }
      } while (testResult == false);
      currEta = nowEta;
      if (currEta != 1.) {
        unifiedCovMatrix *= currEta;
        if (m_numDisabledParameters > 0) { // gpmsa2
          for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
            if (m_parameterEnabledStatus[paramId] == false) {
              for (unsigned int i = 0; i < m_vectorSpace.dimLocal(); ++i) {
                unifiedCovMatrix(i,paramId) = 0.;
              }
              for (unsigned int j = 0; j < m_vectorSpace.dimLocal(); ++j) {
                unifiedCovMatrix(paramId,j) = 0.;
              }
              unifiedCovMatrix(paramId,paramId) = 1.;
            }
          }
        }
      }

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                << ", step "                                   << m_currStep
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", currEta = "                              << currEta
                                << ", assessed rejection rate = "              << nowRejectionRate
                                << std::endl;
      }
    }

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step10_all(
  MLSamplingLevelOptions&                  currOptions,                  // input (changed temporarily internally)
  const P_M&                                      unifiedCovMatrix,             // input
  const GenericVectorRV  <P_V,P_M>&        currRv,                       // input
  bool                                            useBalancedChains,            // input
  const UnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,        // input // Round Rock
  unsigned int                                    indexOfFirstWeight,           // input // Round Rock
  const SequenceOfVectors<P_V,P_M>&        prevChain,                    // input // Round Rock
  double                                   prevExponent,                       // input
  double                                   currExponent,                       // input
  const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
  const ScalarSequence<double>&            prevLogTargetValues,                // input
  const BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,          // input // Round Rock
  SequenceOfVectors      <P_V,P_M>&        currChain,                    // output
  double&                                         cumulativeRawChainRunTime,    // output
  unsigned int&                                   cumulativeRawChainRejections, // output
  ScalarSequence         <double>*         currLogLikelihoodValues,      // output
  ScalarSequence         <double>*         currLogTargetValues)          // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 10 of 11"
                                << ", currLogLikelihoodValues = " << currLogLikelihoodValues
                                << std::endl;
      }

      // All nodes should call here
      bool         savedTotallyMute           = currOptions.m_totallyMute; // HERE - ENHANCEMENT
      unsigned int savedRawChainSize          = currOptions.m_rawChainSize; // Ok to use rawChainSize
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      bool         savedRawChainComputeStats  = currOptions.m_rawChainComputeStats;
#endif
      bool         savedFilteredChainGenerate = currOptions.m_filteredChainGenerate;

      currOptions.m_totallyMute = true;
      if (m_env.displayVerbosity() >= 999999) {
        currOptions.m_totallyMute = false;
      }
      currOptions.m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      currOptions.m_rawChainComputeStats  = false;
#endif
      currOptions.m_filteredChainGenerate = false;

      // All nodes should call here
      if (useBalancedChains) {
        generateBalLinkedChains_all(currOptions,                  // input, only m_rawChainSize changes
                                    unifiedCovMatrix,             // input
                                    currRv,                       // input
                                    balancedLinkControl,          // input // Round Rock
                                    currChain,                    // output
                                    cumulativeRawChainRunTime,    // output
                                    cumulativeRawChainRejections, // output
                                    currLogLikelihoodValues,      // output // likelihood is important
                                    currLogTargetValues);         // output
      }
      else {
        generateUnbLinkedChains_all(currOptions,                  // input, only m_rawChainSize changes
                                    unifiedCovMatrix,             // input
                                    currRv,                       // input
                                    unbalancedLinkControl,        // input // Round Rock
                                    indexOfFirstWeight,           // input // Round Rock
                                    prevChain,                    // input // Round Rock
                                    prevExponent,                 // input
                                    currExponent,                 // input
                                    prevLogLikelihoodValues,      // input
                                    prevLogTargetValues,          // input
                                    currChain,                    // output
                                    cumulativeRawChainRunTime,    // output
                                    cumulativeRawChainRejections, // output
                                    currLogLikelihoodValues,      // output // likelihood is important
                                    currLogTargetValues);         // output
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        double tmpValue = INFINITY;
        if (currLogLikelihoodValues) tmpValue = (*currLogLikelihoodValues)[0];
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", after chain generatrion"
                                << ", currLogLikelihoodValues[0] = " << tmpValue
                                << std::endl;
      }

      // All nodes should call here
      currOptions.m_totallyMute           = savedTotallyMute;
      currOptions.m_rawChainSize          = savedRawChainSize;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      currOptions.m_rawChainComputeStats  = savedRawChainComputeStats;
#endif
      currOptions.m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence_Step11_inter0(
  const MLSamplingLevelOptions* currOptions,                  // input
  unsigned int                         unifiedRequestedNumSamples,   // input
  unsigned int                         cumulativeRawChainRejections, // input
  SequenceOfVectors<P_V,P_M>&   currChain,                    // input/output
  ScalarSequence<double>&       currLogLikelihoodValues,      // input/output
  ScalarSequence<double>&       currLogTargetValues,          // input/output
  unsigned int&                        unifiedNumberOfRejections)    // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": beginning step 11 of 11"
                            << std::endl;
  }

  //if (m_env.subComm().MyPID() == 0) std::cout << "Aqui 000" << std::endl;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (currOptions->m_rawChainComputeStats) {
    FilePtrSetStruct filePtrSet;
    m_env.openOutputFile(currOptions->m_dataOutputFileName,
                         UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                         currOptions->m_dataOutputAllowedSet,
                         false,
                         filePtrSet);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) { // output debug
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", calling computeStatistics for raw chain"
                              << ". Ofstream pointer value = " << filePtrSet.ofsVar
                              << ", statistical options are"
                              << "\n" << *currOptions->m_rawChainStatisticalOptionsObj
                              << std::endl;
    }
    //m_env.syncPrintDebugMsg("At step 11, calling computeStatistics for raw chain",1,10,m_env.inter0Comm()); // output debug
    currChain.computeStatistics(*currOptions->m_rawChainStatisticalOptionsObj,
                                filePtrSet.ofsVar);

    m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
  }
  // Compute MLE and MAP
  // rr0
#endif
  if (currOptions->m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
    currChain.unifiedWriteContents(currOptions->m_rawChainDataOutputFileName,
                                   currOptions->m_rawChainDataOutputFileType); // KAUST5
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", before calling currLogLikelihoodValues.unifiedWriteContents()"
                              << ", currLogLikelihoodValues[0] = " << currLogLikelihoodValues[0]
                              << std::endl;
    }
    currLogLikelihoodValues.unifiedWriteContents(currOptions->m_rawChainDataOutputFileName,
                                                 currOptions->m_rawChainDataOutputFileType);
    currLogTargetValues.unifiedWriteContents    (currOptions->m_rawChainDataOutputFileName,
                                                 currOptions->m_rawChainDataOutputFileType);
  }

  if (currOptions->m_filteredChainGenerate) {
    FilePtrSetStruct filePtrSet;
    m_env.openOutputFile(currOptions->m_dataOutputFileName,
                         UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                         currOptions->m_dataOutputAllowedSet,
                         false,
                         filePtrSet);

    unsigned int filterInitialPos = (unsigned int) (currOptions->m_filteredChainDiscardedPortion * (double) currChain.subSequenceSize());
    unsigned int filterSpacing    = currOptions->m_filteredChainLag;
    if (filterSpacing == 0) {
      currChain.computeFilterParams(filePtrSet.ofsVar,
                                    filterInitialPos,
                                    filterSpacing);
    }

    // Filter positions from the converged portion of the chain
    currChain.filter(filterInitialPos,
                     filterSpacing);
    currChain.setName(currOptions->m_prefix + "filtChain");

    currLogLikelihoodValues.filter(filterInitialPos,
                                   filterSpacing);
    currLogLikelihoodValues.setName(currOptions->m_prefix + "filtLogLikelihood");

    currLogTargetValues.filter(filterInitialPos,
                               filterSpacing);
    currLogTargetValues.setName(currOptions->m_prefix + "filtLogTarget");

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    if (currOptions->m_filteredChainComputeStats) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) { // output debug
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence_Step()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", calling computeStatistics for filtered chain"
                                << ". Ofstream pointer value = " << filePtrSet.ofsVar
                                << ", statistical options are"
                                << "\n" << *currOptions->m_rawChainStatisticalOptionsObj
                                << std::endl;
      }

      //m_env.syncPrintDebugMsg("At step 11, calling computeStatistics for filtered chain",1,10,m_env.inter0Comm()); // output debug
      currChain.computeStatistics(*currOptions->m_filteredChainStatisticalOptionsObj,
                                  filePtrSet.ofsVar);

    }
#endif
    // Compute MLE and MAP
    // rr0
    m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);

    if (currOptions->m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
      currChain.unifiedWriteContents              (currOptions->m_filteredChainDataOutputFileName,
                                                   currOptions->m_filteredChainDataOutputFileType);
      currLogLikelihoodValues.unifiedWriteContents(currOptions->m_filteredChainDataOutputFileName,
                                                   currOptions->m_filteredChainDataOutputFileType);
      currLogTargetValues.unifiedWriteContents    (currOptions->m_filteredChainDataOutputFileName,
                                                   currOptions->m_filteredChainDataOutputFileType);
    }
  } // if (currOptions->m_filteredChainGenerate)

  if (currOptions->m_filteredChainGenerate) {
    // Do not check
  }
  else {
    // Check if unified size of generated chain matches the unified requested size // KAUST
    unsigned int tmpSize = currChain.subSequenceSize();
    unsigned int unifiedGeneratedNumSamples = 0;
    m_env.inter0Comm().template Allreduce<unsigned int>(&tmpSize, &unifiedGeneratedNumSamples, (int) 1, RawValue_MPI_SUM,
                                 "MLSampling<P_V,P_M>::generateSequence()",
                                 "failed MPI.Allreduce() for generated num samples in step 11");
    //std::cout << "unifiedGeneratedNumSamples = "   << unifiedGeneratedNumSamples
    //          << ", unifiedRequestedNumSamples = " << unifiedRequestedNumSamples
    //          << std::endl;
    queso_require_equal_to_msg(unifiedGeneratedNumSamples, unifiedRequestedNumSamples, "currChain (linked one) has been generated with invalid size");
  }

  // Compute unified number of rejections
  m_env.inter0Comm().template Allreduce<unsigned int>(&cumulativeRawChainRejections, &unifiedNumberOfRejections, (int) 1, RawValue_MPI_SUM,
                               "MLSampling<P_V,P_M>::generateSequence()",
                               "failed MPI.Allreduce() for number of rejections");

  double stepRunTime = MiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

// Default constructor -----------------------------
template<class P_V,class P_M>
MLSampling<P_V,P_M>::MLSampling(
  const char*                        prefix,
  const BaseVectorRV      <P_V,P_M>& priorRv,
  const BaseScalarFunction<P_V,P_M>& likelihoodFunction)
  :
  m_env               (priorRv.env()),
  m_priorRv           (priorRv),
  m_likelihoodFunction(likelihoodFunction),
  m_vectorSpace       (m_priorRv.imageSet().vectorSpace()),
  m_targetDomain      (InstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet())),
  m_numDisabledParameters (0), // gpmsa2
  m_parameterEnabledStatus(m_vectorSpace.dimLocal(),true), // gpmsa2
  m_options           (m_env,prefix),
  m_currLevel         (0),
  m_currStep          (0),
  m_debugExponent     (0.),
  m_logEvidenceFactors(0),
  m_logEvidence       (0.),
  m_meanLogLikelihood (0.),
  m_eig               (0.)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::constructor()"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::constructor()"
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V,class P_M>
MLSampling<P_V,P_M>::~MLSampling()
{
  m_numDisabledParameters = 0; // gpmsa2
  m_parameterEnabledStatus.clear(); // gpmsa2
  if (m_targetDomain) delete m_targetDomain;
}
// Statistical methods-------------------------------
/* This operation currently implements the PAMSSA algorithm (S. H. Cheung and E. E. Prudencio. Parallel adaptive multilevel
 * sampling algorithms for the Bayesian analysis of mathematical models. International Journal
 * for Uncertainty Quantification, 2(3):215237, 2012.)*/
template <class P_V,class P_M>
void
MLSampling<P_V,P_M>::generateSequence(
  BaseVectorSequence<P_V,P_M>& workingChain,
  ScalarSequence<double>*      workingLogLikelihoodValues,
  ScalarSequence<double>*      workingLogTargetValues)
{
  struct timeval timevalRoutineBegin;
  int iRC = 0;
  iRC = gettimeofday(&timevalRoutineBegin, NULL);
  if (iRC) {}; // just to remove compiler warning

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering MLSampling<P_V,P_M>::generateSequence()"
                            << ", at  "   << ctime(&timevalRoutineBegin.tv_sec)
                            << ", after " << timevalRoutineBegin.tv_sec - m_env.timevalBegin().tv_sec
                            << " seconds from queso environment instatiation..."
                            << std::endl;
  }

  //***********************************************************
  // Declaration of Variables
  //***********************************************************
  double                            currExponent                   = 0.;   // restate
  double                            currEta                        = 1.;   // restate
  unsigned int                      currUnifiedRequestedNumSamples = 0;
  SequenceOfVectors<P_V,P_M> currChain              (m_vectorSpace, // restate
                                                            0,
                                                            m_options.m_prefix+"curr_chain");
  ScalarSequence<double>     currLogLikelihoodValues(m_env,         // restate
                                                            0,
                                                            "");
  ScalarSequence<double>     currLogTargetValues    (m_env,         // restate
                                                            0,
                                                            "");

  bool stopAtEndOfLevel = false;
  char levelPrefix[256];

  //***********************************************************
  // Take care of first level (level '0')
  //***********************************************************
  MLSamplingLevelOptions defaultLevelOptions(m_env,(m_options.m_prefix + "default_").c_str());
  defaultLevelOptions.scanOptionsValues(NULL);

  MLSamplingLevelOptions lastLevelOptions(m_env,(m_options.m_prefix + "last_").c_str());
  lastLevelOptions.scanOptionsValues(&defaultLevelOptions);

  if (m_options.m_restartInput_baseNameForFiles != ".") {
    restartML(currExponent,            // output
              currEta,                 // output
              currChain,               // output
              currLogLikelihoodValues, // output
              currLogTargetValues);    // output

    if (currExponent == 1.) {
      if (lastLevelOptions.m_parameterDisabledSet.size() > 0) { // gpmsa2
        for (std::set<unsigned int>::iterator setIt = lastLevelOptions.m_parameterDisabledSet.begin(); setIt != lastLevelOptions.m_parameterDisabledSet.end(); ++setIt) {
          unsigned int paramId = *setIt;
          if (paramId < m_vectorSpace.dimLocal()) {
            m_numDisabledParameters++;
            m_parameterEnabledStatus[paramId] = false;
          }
        }
      }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      if (lastLevelOptions.m_rawChainComputeStats) {
        FilePtrSetStruct filePtrSet;
        m_env.openOutputFile(lastLevelOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                             lastLevelOptions.m_dataOutputAllowedSet,
                             false,
                             filePtrSet);

        currChain.computeStatistics(*lastLevelOptions.m_rawChainStatisticalOptionsObj,
                                    filePtrSet.ofsVar);

        m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
      }
#endif
      // Compute MLE and MAP
      // rr0

      if (lastLevelOptions.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
        currChain.unifiedWriteContents              (lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType);
        currLogLikelihoodValues.unifiedWriteContents(lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType);
        currLogTargetValues.unifiedWriteContents    (lastLevelOptions.m_rawChainDataOutputFileName,lastLevelOptions.m_rawChainDataOutputFileType);
      }

      if (lastLevelOptions.m_filteredChainGenerate) {
        FilePtrSetStruct filePtrSet;
        m_env.openOutputFile(lastLevelOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                             lastLevelOptions.m_dataOutputAllowedSet,
                             false,
                             filePtrSet);

        unsigned int filterInitialPos = (unsigned int) (lastLevelOptions.m_filteredChainDiscardedPortion * (double) currChain.subSequenceSize());
        unsigned int filterSpacing    = lastLevelOptions.m_filteredChainLag;
        if (filterSpacing == 0) {
          currChain.computeFilterParams(filePtrSet.ofsVar,
                                        filterInitialPos,
                                        filterSpacing);
        }

        // Filter positions from the converged portion of the chain
        currChain.filter(filterInitialPos,
                         filterSpacing);
        currChain.setName(lastLevelOptions.m_prefix + "filtChain");

        currLogLikelihoodValues.filter(filterInitialPos,
                                       filterSpacing);
        currLogLikelihoodValues.setName(lastLevelOptions.m_prefix + "filtLogLikelihood");

        currLogTargetValues.filter(filterInitialPos,
                                   filterSpacing);
        currLogTargetValues.setName(lastLevelOptions.m_prefix + "filtLogTarget");

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        if (lastLevelOptions.m_filteredChainComputeStats) {
          currChain.computeStatistics(*lastLevelOptions.m_filteredChainStatisticalOptionsObj,
                                      filePtrSet.ofsVar);
        }
#endif
        // Compute MLE and MAP
        // rr0
        m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);

        if (lastLevelOptions.m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
          currChain.unifiedWriteContents              (lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
          currLogLikelihoodValues.unifiedWriteContents(lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
          currLogTargetValues.unifiedWriteContents    (lastLevelOptions.m_filteredChainDataOutputFileName,lastLevelOptions.m_filteredChainDataOutputFileType);
        }
      } // if (lastLevelOptions.m_filteredChainGenerate)
    }
  }
  else {
    sprintf(levelPrefix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    MLSamplingLevelOptions currOptions(m_env,(m_options.m_prefix + levelPrefix).c_str());
    currOptions.scanOptionsValues(&defaultLevelOptions);

    if (currOptions.m_parameterDisabledSet.size() > 0) { // gpmsa2
      for (std::set<unsigned int>::iterator setIt = currOptions.m_parameterDisabledSet.begin(); setIt != currOptions.m_parameterDisabledSet.end(); ++setIt) {
        unsigned int paramId = *setIt;
        if (paramId < m_vectorSpace.dimLocal()) {
          m_numDisabledParameters++;
          m_parameterEnabledStatus[paramId] = false;
        }
      }
    }

    generateSequence_Level0_all(currOptions,                    // input
                                currUnifiedRequestedNumSamples, // output
                                currChain,                      // output
                                currLogLikelihoodValues,        // output
                                currLogTargetValues);           // output

    stopAtEndOfLevel = currOptions.m_stopAtEnd;
    bool performCheckpoint = stopAtEndOfLevel;
    if (m_options.m_restartOutput_levelPeriod > 0) {
      performCheckpoint = performCheckpoint || ( ((m_currLevel + 1) % m_options.m_restartOutput_levelPeriod) == 0 );
    }
    if (performCheckpoint) {
      checkpointML(currExponent,            // input
                   currEta,                 // input
                   currChain,               // input
                   currLogLikelihoodValues, // input
                   currLogTargetValues);    // input
    }
  }
  //std::cout << "In QUESO: end of level 0. Exiting on purpose" << std::endl;
  //exit(1);

  double minLogLike = -INFINITY;
  double maxLogLike =  INFINITY;
  currLogLikelihoodValues.subMinMaxExtra(0,
                                         currLogLikelihoodValues.subSequenceSize(),
                                         minLogLike,
                                         maxLogLike);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                            << ": at end of level "  << m_currLevel+LEVEL_REF_ID
                            << ", sub minLogLike = " << minLogLike
                            << ", sub maxLogLike = " << maxLogLike
                            << std::endl;
  }

  m_env.fullComm().Barrier();

  minLogLike = -INFINITY;
  maxLogLike =  INFINITY;
  currLogLikelihoodValues.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                             0,
                                             currLogLikelihoodValues.subSequenceSize(),
                                             minLogLike,
                                             maxLogLike);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                            << ": at end of level "      << m_currLevel+LEVEL_REF_ID
                            << ", unified minLogLike = " << minLogLike
                            << ", unified maxLogLike = " << maxLogLike
                            << std::endl;
  }

  //***********************************************************
  // Take care of next levels
  //***********************************************************
  while ((currExponent     <  1.   ) && // begin level while
         (stopAtEndOfLevel == false)) {
    m_currLevel++; // restate

    struct timeval timevalLevelBegin;
    iRC = 0;
    iRC = gettimeofday(&timevalLevelBegin, NULL);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level " << m_currLevel+LEVEL_REF_ID
                              << ", at  "   << ctime(&timevalLevelBegin.tv_sec)
                              << ", after " << timevalLevelBegin.tv_sec - timevalRoutineBegin.tv_sec
                              << " seconds from entering the routine"
                              << ", after " << timevalLevelBegin.tv_sec - m_env.timevalBegin().tv_sec
                              << " seconds from queso environment instatiation"
                              << std::endl;
    }

    iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);
    double       cumulativeRawChainRunTime    = 0.;
    unsigned int cumulativeRawChainRejections = 0;

    bool   tryExponentEta = true; // gpmsa1
    double failedExponent = 0.;   // gpmsa1
    double failedEta      = 0.;   // gpmsa1

    MLSamplingLevelOptions*            currOptions           = NULL;  // step 1
    SequenceOfVectors<P_V,P_M>*        prevChain             = NULL;  // step 2
    ScalarSequence<double> prevLogLikelihoodValues(m_env, 0, "");
    ScalarSequence<double> prevLogTargetValues    (m_env, 0, "");
    double                             prevExponent          = 0.;    // step 3
    unsigned int                              indexOfFirstWeight    = 0;     // step 2
    unsigned int                              indexOfLastWeight     = 0;     // step 2
    P_M*                                      unifiedCovMatrix      = NULL;  // step 4
    bool                                      useBalancedChains     = false; // step 6
    BalancedLinkedChainsPerNodeStruct<P_V>* balancedLinkControl   = NULL;  // step 7
    UnbalancedLinkedChainsPerNodeStruct*    unbalancedLinkControl = NULL;  // step 7
    BayesianJointPdf<P_V,P_M>*         currPdf               = NULL;  // step 8
    GenericVectorRV<P_V,P_M>*          currRv                = NULL;  // step 8

    unsigned int exponentEtaTriedAmount = 0;
    while (tryExponentEta) { // gpmsa1
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In IMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", beginning 'do-while(tryExponentEta)"
                                << ": failedExponent = " << failedExponent
                                << ", failedEta = "      << failedEta
                                << std::endl;
      }

    //***********************************************************
    // Step 1 of 11: read options
    //***********************************************************
    m_currStep = 1;
    sprintf(levelPrefix,"%d_",m_currLevel+LEVEL_REF_ID); // Yes, '+0'
    currOptions = new MLSamplingLevelOptions(m_env,(m_options.m_prefix + levelPrefix).c_str());
    currOptions->scanOptionsValues(&defaultLevelOptions);

    if (currOptions->m_parameterDisabledSet.size() > 0) { // gpmsa2
      for (std::set<unsigned int>::iterator setIt = currOptions->m_parameterDisabledSet.begin(); setIt != currOptions->m_parameterDisabledSet.end(); ++setIt) {
        unsigned int paramId = *setIt;
        if (paramId < m_vectorSpace.dimLocal()) {
          m_numDisabledParameters++;
          m_parameterEnabledStatus[paramId] = false;
        }
      }
    }

    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step01_inter0(currOptions,                     // input
                                     currUnifiedRequestedNumSamples); // output
    }

    //***********************************************************
    // Step 2 of 11: save [chain and corresponding target pdf values] from previous level
    //***********************************************************
    m_currStep = 2;
    prevExponent                   = currExponent;
    double       prevEta                        = currEta;
    unsigned int prevUnifiedRequestedNumSamples = currUnifiedRequestedNumSamples;
    prevChain = new SequenceOfVectors<P_V,P_M>(m_vectorSpace,
                                                      0,
                                                      m_options.m_prefix+"prev_chain");

    indexOfFirstWeight = 0;
    indexOfLastWeight  = 0;

    //std::cout << "m_env.inter0Rank() = " << m_env.inter0Rank() << std::endl;
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step02_inter0(currOptions,             // input
                                     currChain,               // input/output // restate
                                     currLogLikelihoodValues, // input/output // restate

                                     currLogTargetValues,     // input/output // restate
                                     *prevChain,              // output
                                     prevLogLikelihoodValues, // output
                                     prevLogTargetValues,     // output
                                     indexOfFirstWeight,      // output
                                     indexOfLastWeight);      // output
    }

    //***********************************************************
    // Step 3 of 11: compute [currExponent and sequence of weights] for current level
    //               update 'm_logEvidenceFactors'
    //***********************************************************
    m_currStep = 3;
    ScalarSequence<double> weightSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step03_inter0(currOptions,             // input
                                     prevLogLikelihoodValues, // input
                                     prevExponent,            // input
                                     failedExponent,          // input // gpmsa1
                                     currExponent,            // output
                                     weightSequence);         // output
    }

    // All nodes in 'subComm' should have the same 'currExponent'
    m_env.subComm().Bcast((void *) &currExponent, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm', important
                          "MLSampling<P_V,P_M>::generateSequence()",
                          "failed MPI.Bcast() for currExponent");
    m_debugExponent = currExponent;

    if (currExponent == 1.) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": copying 'last' level options to current options" // In all nodes of 'subComm', important
                                << std::endl;
      }
      delete currOptions;
      currOptions = &lastLevelOptions;

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": after copying 'last' level options to current options, the current options are"
                                << "\n" << *currOptions
                                << std::endl;
      }

      if (m_env.inter0Rank() >= 0) {
        // It is necessary to recompute 'currUnifiedRequestedNumSamples' because
        // 'currOptions' has just been replaced by 'lastLevelOptions'
        unsigned int tmpSize = currOptions->m_rawChainSize;
        m_env.inter0Comm().template Allreduce<unsigned int>(&tmpSize, &currUnifiedRequestedNumSamples, (int) 1, RawValue_MPI_SUM,
                                     "MLSampling<P_V,P_M>::generateSequence()",
                                     "failed MPI.Allreduce() for requested num samples in step 3");
      }
    }

    //***********************************************************
    // Step 4 of 11: create covariance matrix for current level
    //***********************************************************
    m_currStep = 4;
    P_V oneVec(m_vectorSpace.zeroVector());
    oneVec.cwSet(1.);

    unifiedCovMatrix = NULL;
    if (m_env.inter0Rank() >= 0) {
      unifiedCovMatrix = m_vectorSpace.newMatrix();
    }
    else {
      unifiedCovMatrix = new P_M(oneVec);
    }

    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step04_inter0(*prevChain,         // input
                                     weightSequence,     // input
                                     *unifiedCovMatrix); // output
    }

    //***********************************************************
    // Step 5 of 11: create *unified* finite distribution for current level
    //***********************************************************
    m_currStep = 5;
    std::vector<unsigned int> unifiedIndexCountersAtProc0Only(0);
    std::vector<double>       unifiedWeightStdVectorAtProc0Only(0); // KAUST, to check
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step05_inter0(currUnifiedRequestedNumSamples,     // input
                                     weightSequence,                     // input
                                     unifiedIndexCountersAtProc0Only,    // output
                                     unifiedWeightStdVectorAtProc0Only); // output
    }

    //***********************************************************
    // Step 6 of 11: decide on using balanced chains or not
    //***********************************************************
    m_currStep = 6;
    useBalancedChains = false;
    std::vector<ExchangeInfoStruct> exchangeStdVec(0);
    // All processors should call this routine in order to have the same decision value
    generateSequence_Step06_all(currOptions,                     // input
                                indexOfFirstWeight,              // input
                                indexOfLastWeight,               // input
                                unifiedIndexCountersAtProc0Only, // input
                                useBalancedChains,               // output
                                exchangeStdVec);                 // output

    //***********************************************************
    // Step 7 of 11: plan for number of linked chains for each node so that all
    //               nodes generate the closest possible to the same number of positions
    //***********************************************************
    m_currStep = 7;
    balancedLinkControl   = new BalancedLinkedChainsPerNodeStruct<P_V>();
    unbalancedLinkControl = new UnbalancedLinkedChainsPerNodeStruct   ();
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step07_inter0(useBalancedChains,               // input
                                     indexOfFirstWeight,              // input
                                     indexOfLastWeight,               // input
                                     unifiedIndexCountersAtProc0Only, // input
                                     *unbalancedLinkControl,          // (possible) output
                                     currOptions,                     // input
                                     *prevChain,                      // input
                                     prevExponent,                    // input
                                     currExponent,                    // input
                                     prevLogLikelihoodValues,         // input
                                     prevLogTargetValues,             // input
                                     exchangeStdVec,                  // (possible) input/output
                                     *balancedLinkControl);           // (possible) output
    }

    // aqui: prevChain might not be needed anymore. Delete it to save memory

    //***********************************************************
    // Step 8 of 11: create vector RV for current level
    //***********************************************************
    m_currStep = 8;
    currPdf = new BayesianJointPdf<P_V,P_M> (m_options.m_prefix.c_str(),
                                                    m_priorRv.pdf(),
                                                    m_likelihoodFunction,
                                                    currExponent,
                                                    *m_targetDomain);

    currRv = new GenericVectorRV<P_V,P_M> (m_options.m_prefix.c_str(),
                                                  *m_targetDomain);

    // All nodes should set 'currRv'
    generateSequence_Step08_all(*currPdf,
                                *currRv);

    //***********************************************************
    // Step 9 of 11: scale unified covariance matrix until min <= rejection rate <= max
    //***********************************************************
    m_currStep = 9;
    generateSequence_Step09_all(*prevChain,                        // input
				prevExponent,                    // input
				currExponent,                    // input
				prevLogLikelihoodValues,         // input
				prevLogTargetValues,             // input
                                indexOfFirstWeight,                // input
                                indexOfLastWeight,                 // input
                                unifiedWeightStdVectorAtProc0Only, // input
                                weightSequence,                    // input
                                prevEta,                           // input
                                *currRv,                           // input
                                currOptions,                       // input (changed temporarily internally)
                                *unifiedCovMatrix,                 // input/output
                                currEta);                          // output

    tryExponentEta = false; // gpmsa1
    if ((currOptions->m_minAcceptableEta > 0.     ) && // gpmsa1
        (currEta < currOptions->m_minAcceptableEta)) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", preparing to retry ExponentEta"
                                << ": currExponent = " << currExponent
                                << ", currEta = "      << currEta
                                << std::endl;
      }
      exponentEtaTriedAmount++;
      tryExponentEta = true;
      failedExponent = currExponent;
      failedEta      = currEta;

      // "Return" to previous level
      delete currRv;                 // Step 8
      currRv = NULL;                 // Step 8
      delete currPdf;                // Step 8
      currPdf = NULL;                // Step 8

      delete balancedLinkControl;    // Step 7
      balancedLinkControl   = NULL;  // Step 7
      delete unbalancedLinkControl;  // Step 7
      unbalancedLinkControl = NULL;  // Step 7

      delete unifiedCovMatrix;       // Step 4
      unifiedCovMatrix = NULL;       // Step 4

      currExponent                   = prevExponent;
      currEta                        = 1.; // prevEta;
      currUnifiedRequestedNumSamples = prevUnifiedRequestedNumSamples;

      currChain.clear();             // Step 2
      currChain = (*prevChain);      // Step 2
      delete prevChain;              // Step 2
      prevChain = NULL;              // Step 2

      currLogLikelihoodValues        = prevLogLikelihoodValues;
      currLogTargetValues            = prevLogTargetValues;
    }
    } // while (tryExponentEta) // gpmsa1

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) { // gpmsa1
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", exited 'do-while(tryExponentEta)"
                              << ", failedExponent = " << failedExponent
                              << ", failedEta = "      << failedEta
                              << std::endl;
    }

    //***********************************************************
    // Step 10 of 11: sample vector RV of current level
    //***********************************************************
    m_currStep = 10;
    // All nodes should call here
    generateSequence_Step10_all(*currOptions,                 // input (changed temporarily internally)
                                *unifiedCovMatrix,            // input
                                *currRv,                      // input
                                useBalancedChains,            // input
                                *unbalancedLinkControl,       // input // Round Rock
                                indexOfFirstWeight,           // input // Round Rock
                                *prevChain,                   // input // Round Rock
                                prevExponent,                 // input
                                currExponent,                 // input
                                prevLogLikelihoodValues,      // input
                                prevLogTargetValues,          // input
                                *balancedLinkControl,         // input // Round Rock
                                currChain,                    // output
                                cumulativeRawChainRunTime,    // output
                                cumulativeRawChainRejections, // output
                                &currLogLikelihoodValues,     // output // likelihood is important
                                &currLogTargetValues);        // output);

    //***********************************************************
    // Perform checkpoint if necessary
    //***********************************************************
    stopAtEndOfLevel = currOptions->m_stopAtEnd;
    bool performCheckpoint = stopAtEndOfLevel;
    if (m_options.m_restartOutput_levelPeriod > 0) {
      performCheckpoint = performCheckpoint || ( ((m_currLevel + 1) % m_options.m_restartOutput_levelPeriod) == 0 );
      if (currExponent == 1.) {
        performCheckpoint = true;
      }
    }
    if (performCheckpoint) {
      checkpointML(currExponent,            // input
                   currEta,                 // input
                   currChain,               // input
                   currLogLikelihoodValues, // input
                   currLogTargetValues);    // input
    }

    //***********************************************************
    // Just free some memory
    //***********************************************************
    {
      delete unifiedCovMatrix;

      for (unsigned int i = 0; i < balancedLinkControl->balLinkedChains.size(); ++i) {
        queso_require_msg(balancedLinkControl->balLinkedChains[i].initialPosition, "Initial position pointer in step 9 should not be NULL");
        delete balancedLinkControl->balLinkedChains[i].initialPosition;
        balancedLinkControl->balLinkedChains[i].initialPosition = NULL;
      }
      balancedLinkControl->balLinkedChains.clear();
    }

    //***********************************************************
    // Step 11 of 11: filter chain if requested
    //***********************************************************
    m_currStep = 11;
    unsigned int unifiedNumberOfRejections = 0;
    if (m_env.inter0Rank() >= 0) {
      generateSequence_Step11_inter0(currOptions,                      // input
                                     currUnifiedRequestedNumSamples,   // input
                                     cumulativeRawChainRejections,     // input
                                     currChain,                        // input/output
                                     currLogLikelihoodValues,          // input/output
                                     currLogTargetValues,              // input/output
                                     unifiedNumberOfRejections);       // output
    }

    minLogLike = -INFINITY;
    maxLogLike =  INFINITY;
    currLogLikelihoodValues.subMinMaxExtra(0,
                                           currLogLikelihoodValues.subSequenceSize(),
                                           minLogLike,
                                           maxLogLike);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": at end of level "  << m_currLevel+LEVEL_REF_ID
                              << ", sub minLogLike = " << minLogLike
                              << ", sub maxLogLike = " << maxLogLike
                              << std::endl;
    }

    m_env.fullComm().Barrier();

    minLogLike = -INFINITY;
    maxLogLike =  INFINITY;
    currLogLikelihoodValues.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                               0,
                                               currLogLikelihoodValues.subSequenceSize(),
                                               minLogLike,
                                               maxLogLike);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": at end of level "      << m_currLevel+LEVEL_REF_ID
                              << ", unified minLogLike = " << minLogLike
                              << ", unified maxLogLike = " << maxLogLike
                              << std::endl;
    }

    //***********************************************************
    // Prepare to end current level
    //***********************************************************
    double levelRunTime = MiscGetEllapsedSeconds(&timevalLevel);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level "                   << m_currLevel+LEVEL_REF_ID
                              << ", having generated "               << currChain.subSequenceSize()
                              << " chain positions"
                              << ", cumulativeRawChainRunTime = "    << cumulativeRawChainRunTime << " seconds"
                              << ", total level time = "             << levelRunTime              << " seconds"
                              << ", cumulativeRawChainRejections = " << cumulativeRawChainRejections
                              << " (" << 100.*((double) cumulativeRawChainRejections)/((double) currOptions->m_rawChainSize)
                              << "% at this processor)"
                              << " (" << 100.*((double) unifiedNumberOfRejections)/((double) currUnifiedRequestedNumSamples)
                              << "% over all processors)"
                              << ", stopAtEndOfLevel = " << stopAtEndOfLevel
                              << std::endl;
    }

    if (m_env.inter0Rank() >= 0) {
      double minCumulativeRawChainRunTime = 0.;
      m_env.inter0Comm().template Allreduce<double>(&cumulativeRawChainRunTime, &minCumulativeRawChainRunTime, (int) 1, RawValue_MPI_MIN,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for min cumulative raw chain run time");

      double maxCumulativeRawChainRunTime = 0.;
      m_env.inter0Comm().template Allreduce<double>(&cumulativeRawChainRunTime, &maxCumulativeRawChainRunTime, (int) 1, RawValue_MPI_MAX,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for max cumulative raw chain run time");

      double avgCumulativeRawChainRunTime = 0.;
      m_env.inter0Comm().template Allreduce<double>(&cumulativeRawChainRunTime, &avgCumulativeRawChainRunTime, (int) 1, RawValue_MPI_SUM,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for sum cumulative raw chain run time");
      avgCumulativeRawChainRunTime /= ((double) m_env.inter0Comm().NumProc());

      double minLevelRunTime = 0.;
      m_env.inter0Comm().template Allreduce<double>(&levelRunTime, &minLevelRunTime, (int) 1, RawValue_MPI_MIN,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for min level run time");

      double maxLevelRunTime = 0.;
      m_env.inter0Comm().template Allreduce<double>(&levelRunTime, &maxLevelRunTime, (int) 1, RawValue_MPI_MAX,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for max level run time");

      double avgLevelRunTime = 0.;
      m_env.inter0Comm().template Allreduce<double>(&levelRunTime, &avgLevelRunTime, (int) 1, RawValue_MPI_SUM,
                                   "MLSampling<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for sum level run time");
      avgLevelRunTime /= ((double) m_env.inter0Comm().NumProc());

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                                << ", level "               << m_currLevel+LEVEL_REF_ID
                                << ": min cumul seconds = " << minCumulativeRawChainRunTime
                                << ", avg cumul seconds = " << avgCumulativeRawChainRunTime
                                << ", max cumul seconds = " << maxCumulativeRawChainRunTime
                                << ", min level seconds = " << minLevelRunTime
                                << ", avg level seconds = " << avgLevelRunTime
                                << ", max level seconds = " << maxLevelRunTime
                                << std::endl;
      }
    }

    if (currExponent != 1.) delete currOptions;

    struct timeval timevalLevelEnd;
    iRC = 0;
    iRC = gettimeofday(&timevalLevelEnd, NULL);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "Getting at the end of level " << m_currLevel+LEVEL_REF_ID
                              << ", as part of a 'while' on levels"
                              << ", at  "   << ctime(&timevalLevelEnd.tv_sec)
                              << ", after " << timevalLevelEnd.tv_sec - timevalRoutineBegin.tv_sec
                              << " seconds from entering the routine"
                              << ", after " << timevalLevelEnd.tv_sec - m_env.timevalBegin().tv_sec
                              << " seconds from queso environment instatiation"
                              << std::endl;
    }
  } // end of level while


  //                    m_env.worldRank(),
  //                    "MLSampling<P_V,P_M>::generateSequence()",
  //                    "exponent has not achieved value '1' even after exiting level loop");

  //***********************************************************
  // Compute information gain
  // ln( \pi(D|M) ) = E[ln( \pi(D|\theta,M) )] - E[ln( \pi(\theta|D,M) / \pi(\theta|M) )]
  //***********************************************************
  if (m_env.inter0Rank() >= 0) { // KAUST
    queso_require_equal_to_msg(m_currLevel, m_logEvidenceFactors.size(), "invalid m_currLevel at the exit of the level loop");
    m_logEvidence = 0.;
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      m_logEvidence += m_logEvidenceFactors[i];
    }

#if 1 // prudenci-2012-07-06
    m_meanLogLikelihood = currLogLikelihoodValues.unifiedMeanPlain(m_vectorSpace.numOfProcsForStorage() == 1);
#else
    m_meanLogLikelihood = currLogLikelihoodValues.unifiedMeanExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                                                   0,
                                                                   currLogLikelihoodValues.subSequenceSize());
#endif

    m_eig = m_meanLogLikelihood - m_logEvidence;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MLSampling<P_V,P_M>::generateSequence()"
                              << ", log(evidence) = "     << m_logEvidence
                              << ", evidence = "          << exp(m_logEvidence)
                              << ", meanLogLikelihood = " << m_meanLogLikelihood
                              << ", eig = "               << m_eig
                              << std::endl;
    }
  }

  m_env.subComm().Bcast((void *) &m_logEvidence, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm'
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "failed MPI.Bcast() for m_logEvidence");

  m_env.subComm().Bcast((void *) &m_meanLogLikelihood, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm'
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "failed MPI.Bcast() for m_meanLogLikelihood");

  m_env.subComm().Bcast((void *) &m_eig, (int) 1, RawValue_MPI_DOUBLE, 0, // Yes, 'subComm'
                        "MLSampling<P_V,P_M>::generateSequence()",
                        "failed MPI.Bcast() for m_eig");

  //***********************************************************
  // Prepare to return
  //***********************************************************
  workingChain.clear();
  workingChain.resizeSequence(currChain.subSequenceSize());
  P_V auxVec(m_vectorSpace.zeroVector());
  for (unsigned int i = 0; i < workingChain.subSequenceSize(); ++i) {
    if (m_env.inter0Rank() >= 0) {
      currChain.getPositionValues(i,auxVec);
    }
    workingChain.setPositionValues(i,auxVec);
  }

  if (workingLogLikelihoodValues) *workingLogLikelihoodValues = currLogLikelihoodValues;
  if (workingLogTargetValues    ) *workingLogTargetValues     = currLogTargetValues;

  struct timeval timevalRoutineEnd;
  iRC = 0;
  iRC = gettimeofday(&timevalRoutineEnd, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving MLSampling<P_V,P_M>::generateSequence()"
                            << ", at  "   << ctime(&timevalRoutineEnd.tv_sec)
                            << ", after " << timevalRoutineEnd.tv_sec - timevalRoutineBegin.tv_sec
                            << " seconds from entering the routine"
                            << ", after " << timevalRoutineEnd.tv_sec - m_env.timevalBegin().tv_sec
                            << " seconds from queso environment instatiation"
                            << std::endl;
  }

  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const MLSampling<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}

template <class P_V,class P_M>
double MLSampling<P_V,P_M>::logEvidence() const
{
  return m_logEvidence;
}

template <class P_V,class P_M>
double MLSampling<P_V,P_M>::meanLogLikelihood() const
{
  return m_meanLogLikelihood;
}

template <class P_V,class P_M>
double MLSampling<P_V,P_M>::eig() const
{
  return m_eig;
}

}  // End namespace QUESO

template class QUESO::MLSampling<QUESO::GslVector, QUESO::GslMatrix>;
