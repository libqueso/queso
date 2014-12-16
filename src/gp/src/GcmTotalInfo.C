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

#include <queso/GcmTotalInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GcmTotalInfo(
  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s)
  :
  m_env                 (s.m_env),
  m_numConstituents     (4),
  m_constitutiveDomains (m_numConstituents,(const VectorSet   <P_V,P_M>*) NULL),
  m_constitutivePriorRvs(m_numConstituents,(const BaseVectorRV<P_V,P_M>*) NULL),
  m_totalDim            (initializeTotalDim(s)),
  m_totalSpace          (m_env,"total_",m_totalDim,NULL),
  m_totalDomainVolume   (1.),
  m_totalDomain         ("total_",m_totalSpace,m_totalDomainVolume,m_constitutiveDomains),
  m_totalPriorRv        ("total_prior_",m_constitutivePriorRvs,m_totalDomain),
  m_totalPostRv         ("total_post_",m_totalSpace),
  m_like_previousTotal  (m_totalSpace.zeroVector()),
  m_totalPostMean       (m_totalSpace.zeroVector()),
  m_totalPostMedian     (m_totalSpace.zeroVector()),
  m_totalPostMode       (m_totalSpace.zeroVector()),
  m_totalPostMaxLnValue (-INFINITY),
  m_totalMLE            (m_totalSpace.zeroVector()),
  m_totalLikeMaxLnValue (-INFINITY),
  m_solutionDomain      (NULL),
  m_solutionPdf         (NULL),
  m_solutionRealizer    (NULL),
  m_mhSeqGenerator      (NULL),
  m_mlSampler           (NULL),
  m_chain               (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GcmTotalInfo(
  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>& e)
  :
  m_env                 (s.m_env),
  m_numConstituents     (8),
  m_constitutiveDomains (m_numConstituents,(const VectorSet<P_V,P_M>*) NULL),
  m_constitutivePriorRvs(m_numConstituents,(const BaseVectorRV<P_V,P_M>*) NULL),
  m_totalDim            (initializeTotalDim(s,e)),
  m_totalSpace          (m_env, "total_", m_totalDim, NULL),
  m_totalDomainVolume   (1.),
  m_totalDomain         ("total_",m_totalSpace,m_totalDomainVolume,m_constitutiveDomains),
  m_totalPriorRv        ("total_prior_",m_constitutivePriorRvs,m_totalDomain),
  m_totalPostRv         ("total_post_",m_totalSpace),
  m_like_previousTotal  (m_totalSpace.zeroVector()),
  m_totalPostMean       (m_totalSpace.zeroVector()),
  m_totalPostMedian     (m_totalSpace.zeroVector()),
  m_totalPostMode       (m_totalSpace.zeroVector()),
  m_totalPostMaxLnValue (-INFINITY),
  m_totalMLE            (m_totalSpace.zeroVector()),
  m_totalLikeMaxLnValue (-INFINITY),
  m_solutionDomain      (NULL),
  m_solutionPdf         (NULL),
  m_solutionRealizer    (NULL),
  m_mhSeqGenerator      (NULL),
  m_mlSampler           (NULL),
  m_chain               (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~GcmTotalInfo()
{
  if (m_chain) {
    m_chain->clear();
    delete m_chain;
  }
  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
unsigned int
GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::initializeTotalDim(const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s)
{
  m_constitutiveDomains[0] = &(s.m_1lambdaEtaDomain);
  m_constitutiveDomains[1] = &(s.m_2lambdaWDomain);
  m_constitutiveDomains[2] = &(s.m_3rhoWDomain);
  m_constitutiveDomains[3] = &(s.m_4lambdaSDomain);

  m_constitutivePriorRvs[0] = &(s.m_1lambdaEtaPriorRv);
  m_constitutivePriorRvs[1] = &(s.m_2lambdaWPriorRv);
  m_constitutivePriorRvs[2] = &(s.m_3rhoWPriorRv);
  m_constitutivePriorRvs[3] = &(s.m_4lambdaSPriorRv);

  return (s.m_1lambdaEtaDim + s.m_2lambdaWDim + s.m_3rhoWDim + s.m_4lambdaSDim);
}


template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
unsigned int
GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::initializeTotalDim(
  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>& e)
{
  m_constitutiveDomains[0] = &(s.m_1lambdaEtaDomain);
  m_constitutiveDomains[1] = &(s.m_2lambdaWDomain);
  m_constitutiveDomains[2] = &(s.m_3rhoWDomain);
  m_constitutiveDomains[3] = &(s.m_4lambdaSDomain);
  m_constitutiveDomains[4] = &(e.m_5lambdaYDomain);
  m_constitutiveDomains[5] = &(e.m_6lambdaVDomain);
  m_constitutiveDomains[6] = &(e.m_7rhoVDomain);
  m_constitutiveDomains[7] = &(e.m_8thetaPriorRv.imageSet());

  m_constitutivePriorRvs[0] = &(s.m_1lambdaEtaPriorRv);
  m_constitutivePriorRvs[1] = &(s.m_2lambdaWPriorRv);
  m_constitutivePriorRvs[2] = &(s.m_3rhoWPriorRv);
  m_constitutivePriorRvs[3] = &(s.m_4lambdaSPriorRv);
  m_constitutivePriorRvs[4] = &(e.m_5lambdaYPriorRv);
  m_constitutivePriorRvs[5] = &(e.m_6lambdaVPriorRv);
  m_constitutivePriorRvs[6] = &(e.m_7rhoVPriorRv);
  m_constitutivePriorRvs[7] = &(e.m_8thetaPriorRv);

  return (s.m_1lambdaEtaDim + s.m_2lambdaWDim + s.m_3rhoWDim + s.m_4lambdaSDim + e.m_5lambdaYDim + e.m_6lambdaVDim + e.m_7rhoVDim + e.m_8thetaDim);
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor()
{
  m_totalPriorRv.pdf().setNormalizationStyle(0); // CSRI - 2013-Aug-06 - with Laura
  for (unsigned int i = 0; i < m_numConstituents; ++i) {
    m_totalDomainVolume *= m_constitutiveDomains[i]->volume();
  }

  //********************************************************************************
  // Display information
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor()"
                            << "\n  m_totalDim          = "  << m_totalDim
                            << "\n  m_numConstituents   = "  << m_numConstituents
                            << "\n  m_totalDomainVolume = "  << m_totalDomainVolume
                            << std::endl;
  }

  return;
}

}  // End namespace QUESO

template class QUESO::GcmTotalInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
