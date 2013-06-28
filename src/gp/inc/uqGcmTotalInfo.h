//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_GCM_TOTAL_INFO_H__
#define __UQ_GCM_TOTAL_INFO_H__

#include <uqGcmSimulationInfo.h>
#include <uqGcmExperimentInfo.h>
#include <uqGcmJointInfo.h>
#include <uqMetropolisHastingsSG1.h>

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
class uqGcmTotalInfoClass
{
public:
  uqGcmTotalInfoClass(const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s);
  uqGcmTotalInfoClass(const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
                      const uqGcmExperimentInfoClass<S_V,S_M,D_V,D_M,P_V,P_M>& e);
 ~uqGcmTotalInfoClass();

  unsigned int initializeTotalDim(const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s);

  unsigned int initializeTotalDim(const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
                                  const uqGcmExperimentInfoClass<S_V,S_M,D_V,D_M,P_V,P_M>& e);

  const uqBaseEnvironmentClass&                           m_env;
        unsigned int                                      m_numConstituents;
        std::vector<const uqVectorSetClass   <P_V,P_M>* > m_constitutiveDomains;
        std::vector<const uqBaseVectorRVClass<P_V,P_M>* > m_constitutivePriorRvs;
        unsigned int                                      m_totalDim;
        uqVectorSpaceClass         <P_V,P_M>              m_totalSpace;
        double                                            m_totalDomainVolume;
        uqConcatenationSubsetClass <P_V,P_M>              m_totalDomain;
        uqConcatenatedVectorRVClass<P_V,P_M>              m_totalPriorRv;
        uqGenericVectorRVClass     <P_V,P_M>              m_totalPostRv;

        P_V                                               m_like_previousTotal;
        P_V                                               m_totalPostMean;
        P_V                                               m_totalPostMedian;
        P_V                                               m_totalPostMode;
        double                                            m_totalPostMaxLnValue;
        P_V                                               m_totalMLE;
        double                                            m_totalLikeMaxLnValue;

        uqVectorSetClass           <P_V,P_M>*             m_solutionDomain;
        uqBaseJointPdfClass        <P_V,P_M>*             m_solutionPdf;
        uqBaseVectorRealizerClass  <P_V,P_M>*             m_solutionRealizer;

        uqMetropolisHastingsSGClass<P_V,P_M>*             m_mhSeqGenerator;
        uqBaseVectorSequenceClass  <P_V,P_M>*             m_chain;

private:
  void commonConstructor();
};

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::uqGcmTotalInfoClass(
  const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s)
  :
  m_env                 (s.m_env),
  m_numConstituents     (4),
  m_constitutiveDomains (m_numConstituents,(const uqVectorSetClass   <P_V,P_M>*) NULL),
  m_constitutivePriorRvs(m_numConstituents,(const uqBaseVectorRVClass<P_V,P_M>*) NULL),
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
  m_chain               (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::uqGcmTotalInfoClass(
  const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
  const uqGcmExperimentInfoClass<S_V,S_M,D_V,D_M,P_V,P_M>& e)
  :
  m_env                 (s.m_env),
  m_numConstituents     (8),
  m_constitutiveDomains (m_numConstituents,(const uqVectorSetClass<P_V,P_M>*) NULL),
  m_constitutivePriorRvs(m_numConstituents,(const uqBaseVectorRVClass<P_V,P_M>*) NULL),
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
  m_chain               (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~uqGcmTotalInfoClass()
{
  if (m_chain) {
    m_chain->clear();
    delete m_chain;
  }
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
unsigned int
uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::initializeTotalDim(const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s)
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
uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::initializeTotalDim(
  const uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
  const uqGcmExperimentInfoClass<S_V,S_M,D_V,D_M,P_V,P_M>& e)
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
uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor()
{
  m_totalPriorRv.pdf().setNormalizationStyle(1);
  for (unsigned int i = 0; i < m_numConstituents; ++i) {
    m_totalDomainVolume *= m_constitutiveDomains[i]->volume();
  }

  //********************************************************************************
  // Display information
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor()"
                            << "\n  m_totalDim          = "  << m_totalDim
                            << "\n  m_numConstituents   = "  << m_numConstituents
                            << "\n  m_totalDomainVolume = "  << m_totalDomainVolume
                            << std::endl;
  }

  return;
}

#endif // __UQ_GCM_TOTAL_INFO_H__
