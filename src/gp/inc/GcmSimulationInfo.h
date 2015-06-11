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

#ifndef UQ_GCM_SIMULATION_INFO_H
#define UQ_GCM_SIMULATION_INFO_H

#include <queso/SimulationStorage.h>
#include <queso/SimulationModel.h>
#include <queso/VectorRV.h>
#include <queso/GammaVectorRV.h>
#include <queso/BetaVectorRV.h>
#include <queso/GpmsaComputerModelOptions.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class GcmSimulationInfo
{
public:
  GcmSimulationInfo(const GpmsaComputerModelOptions&                  gcmOptionsObj,
                           bool                                                     allOutputsAreScalar,
                           const SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationStorage,
                           const SimulationModel  <S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationModel);
 ~GcmSimulationInfo();

  const BaseEnvironment&                            m_env;
  const SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>& m_simulationStorage;
  const SimulationModel  <S_V,S_M,P_V,P_M,Q_V,Q_M>& m_simulationModel;

        unsigned int                                       m_paper_p_x;
  const std::vector<const S_V* >&                          m_paper_xs_asterisks_standard;
  const std::vector<const P_V* >&                          m_paper_ts_asterisks_standard;
        unsigned int                                       m_paper_p_t;
        unsigned int                                       m_paper_m;
        unsigned int                                       m_paper_n_eta;
        unsigned int                                       m_paper_p_eta;
        VectorSpace         <P_V,P_M>               m_paper_m_space;

        unsigned int                                       m_1lambdaEtaDim; // '1' in paper
        VectorSpace         <P_V,P_M>               m_1lambdaEtaSpace;
        P_V                                                m_1lambdaEtaMins;
        P_V                                                m_1lambdaEtaMaxs;
        BoxSubset           <P_V,P_M>               m_1lambdaEtaDomain;
        P_V                                                m_1lambdaEtaGammaAVec;
        P_V                                                m_1lambdaEtaGammaBVec;
        GammaVectorRV       <P_V,P_M>               m_1lambdaEtaPriorRv;
        P_V                                                m_like_previous1;
        P_V                                                m_tmp_1lambdaEtaVec;

        unsigned int                                       m_2lambdaWDim; // 'p_eta' in paper
        VectorSpace         <P_V,P_M>               m_2lambdaWSpace;
        P_V                                                m_2lambdaWMins;
        P_V                                                m_2lambdaWMaxs;
        BoxSubset           <P_V,P_M>               m_2lambdaWDomain;
        P_V                                                m_2lambdaWGammaAVec;
        P_V                                                m_2lambdaWGammaBVec;
        GammaVectorRV       <P_V,P_M>               m_2lambdaWPriorRv;
        P_V                                                m_like_previous2;
        P_V                                                m_tmp_2lambdaWVec;

        unsigned int                                       m_3rhoWDim; // 'p_eta * (p_x + p_t)' in paper
        VectorSpace         <P_V,P_M>               m_3rhoWSpace;
        P_V                                                m_3rhoWMins;
        P_V                                                m_3rhoWMaxs;
        BoxSubset           <P_V,P_M>               m_3rhoWDomain;
        P_V                                                m_3rhoWBetaAVec;
        P_V                                                m_3rhoWBetaBVec;
        BetaVectorRV        <P_V,P_M>               m_3rhoWPriorRv;
        P_V                                                m_like_previous3;
        P_V                                                m_tmp_3rhoWVec;

        unsigned int                                       m_4lambdaSDim; // 'p_eta' in matlab code
        VectorSpace         <P_V,P_M>               m_4lambdaSSpace;
        P_V                                                m_4lambdaSMins;
        P_V                                                m_4lambdaSMaxs;
        BoxSubset           <P_V,P_M>               m_4lambdaSDomain;
        P_V                                                m_4lambdaSGammaAVec;
        P_V                                                m_4lambdaSGammaBVec;
        GammaVectorRV       <P_V,P_M>               m_4lambdaSPriorRv;
        P_V                                                m_like_previous4;
        P_V                                                m_tmp_4lambdaSVec;

        unsigned int                                       m_eta_size;
        VectorSpace<Q_V,Q_M>                        m_eta_space;

        unsigned int                                       m_w_size;
        VectorSpace<Q_V,Q_M>                        m_w_space;

        VectorSpace<Q_V,Q_M>                        m_unique_w_space;

        Q_V                                                m_Zvec_hat_w;
        VectorSpace<P_V,P_M>                        m_rho_w_space;
        P_V                                                m_tmp_rho_w_vec;
	std::vector<Q_M* >                                 m_Rmat_w_is;       // to be deleted on destructor
	std::vector<Q_M* >                                 m_Smat_w_is;       // to be deleted on destructor
        Q_M                                                m_Smat_w;  // Computed with 'm_simulationModel'
        Q_M                                                m_Smat_w_hat;
	std::vector<Q_M* >                                 m_Rmat_w_hat_w_asterisk_is; // to be deleted on destructor
	std::vector<Q_M* >                                 m_Smat_w_hat_w_asterisk_is; // to be deleted on destructor
        Q_M                                                m_Smat_w_hat_w_asterisk; // Computed with 'm_experimentlModel' and 'm_simulationModel'
        Q_M                                                m_Smat_w_hat_w_asterisk_t;
        Q_M                                                m_Smat_w_asterisk_w_asterisk;

        const Q_M&                                         m_Kmat;      // Equal to 'm_simulationModel.Kmat()'
        const Q_M&                                         m_Kmat_eta;  // Equal to 'm_simulationModel.Kmat_eta()'
        unsigned int                                       m_Kmat_rank;
        Q_M*                                               m_Kt_K;      // to be deleted on destructor
        Q_M*                                               m_Kt_K_inv;  // to be deleted on destructor

        double                                             m_a_eta_modifier;
        double                                             m_b_eta_modifier;

        unsigned int                                       m_predW_counter;
#if 0
        P_V                                                m_predW_samplingRVs_unique_w_meanVec; // todo_rr
        P_M                                                m_predW_samplingRVs_unique_w_covMatrix;
        P_M                                                m_predW_samplingRVs_unique_w_corrMatrix;
#endif
        P_V                                                m_predW_summingRVs_unique_w_meanVec;
        P_M                                                m_predW_summingRVs_mean_of_unique_w_covMatrices;
        P_M                                                m_predW_summingRVs_covMatrix_of_unique_w_means;
        P_M                                                m_predW_summingRVs_corrMatrix_of_unique_w_means;
#if 0
        P_V                                                m_predW_atMean_unique_w_meanVec; // todo_rr
        P_M                                                m_predW_atMean_unique_w_covMatrix;
        P_M                                                m_predW_atMean_unique_w_corrMatrix;

        P_V                                                m_predW_atMedian_unique_w_meanVec; // todo_rr
        P_M                                                m_predW_atMedian_unique_w_covMatrix;
        P_M                                                m_predW_atMedian_unique_w_corrMatrix;

        P_V                                                m_predW_atMode_unique_w_meanVec; // todo_rr
        P_M                                                m_predW_atMode_unique_w_covMatrix;
        P_M                                                m_predW_atMode_unique_w_corrMatrix;

        P_V                                                m_predW_atMLE_unique_w_meanVec; // todo_rr
        P_M                                                m_predW_atMLE_unique_w_covMatrix;
        P_M                                                m_predW_atMLE_unique_w_corrMatrix;
#endif
};

}  // End namespace QUESO

#endif // UQ_GCM_SIMULATION_INFO_H
