//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#ifndef UQ_GCM_TOTAL_INFO_H
#define UQ_GCM_TOTAL_INFO_H

#include <queso/GcmSimulationInfo.h>
#include <queso/GcmExperimentInfo.h>
#include <queso/GcmJointInfo.h>
#include <queso/MetropolisHastingsSG.h>
#include <queso/MLSampling.h>
#include <queso/ConcatenationSubset.h>
#include <queso/ConcatenatedVectorRV.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class D_V = GslVector,
         class D_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix,
         class Q_V = GslVector, class Q_M = GslMatrix>
class GcmTotalInfo
{
public:
  GcmTotalInfo(const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s);
  GcmTotalInfo(const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
                      const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>& e);
 ~GcmTotalInfo();

  unsigned int initializeTotalDim(const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s);

  unsigned int initializeTotalDim(const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
                                  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>& e);

  const BaseEnvironment&                           m_env;
        unsigned int                                      m_numConstituents;
        std::vector<const VectorSet   <P_V,P_M>* > m_constitutiveDomains;
        std::vector<const BaseVectorRV<P_V,P_M>* > m_constitutivePriorRvs;
        unsigned int                                      m_totalDim;
        VectorSpace         <P_V,P_M>              m_totalSpace;
        double                                            m_totalDomainVolume;
        ConcatenationSubset <P_V,P_M>              m_totalDomain;
        ConcatenatedVectorRV<P_V,P_M>              m_totalPriorRv;
        GenericVectorRV     <P_V,P_M>              m_totalPostRv;

        P_V                                               m_like_previousTotal;
        P_V                                               m_totalPostMean;
        P_V                                               m_totalPostMedian;
        P_V                                               m_totalPostMode;
        double                                            m_totalPostMaxLnValue;
        P_V                                               m_totalMLE;
        double                                            m_totalLikeMaxLnValue;

        VectorSet           <P_V,P_M>*             m_solutionDomain;
        BaseJointPdf        <P_V,P_M>*             m_solutionPdf;
        BaseVectorRealizer  <P_V,P_M>*             m_solutionRealizer;

        MetropolisHastingsSG<P_V,P_M>*             m_mhSeqGenerator;
        MLSampling          <P_V,P_M>*             m_mlSampler;
        BaseVectorSequence  <P_V,P_M>*             m_chain;

private:
  void commonConstructor();
};

}  // End namespace QUESO

#endif // UQ_GCM_TOTAL_INFO_H
