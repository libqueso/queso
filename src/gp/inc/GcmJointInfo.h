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

#ifndef UQ_GCM_JOINT_INFO_H
#define UQ_GCM_JOINT_INFO_H

#include <queso/GcmSimulationInfo.h>
#include <queso/GcmExperimentInfo.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class D_V = GslVector,
  class D_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix,
  class Q_V = GslVector, class Q_M = GslMatrix>
class GcmJointInfo
{
public:
  GcmJointInfo(const GpmsaComputerModelOptions&                  gcmOptionsObj,
                      bool                                                     allOutputsAreScalar,
                      const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
                      const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>& e);
 ~GcmJointInfo();

  const BaseEnvironment&                    m_env;

        VectorSpace<Q_V,Q_M>                m_unique_u_space;
        Q_M                                        m_Smat_u_asterisk_u_asterisk;
        unsigned int                               m_u_size;
        VectorSpace<D_V,D_M>                m_u_space;
	std::vector<D_M* >                         m_Rmat_u_is;       // to be deleted on destructor
	std::vector<D_M* >                         m_Smat_u_is;       // to be deleted on destructor
	std::vector<D_M* >                         m_Rmat_uw_is;      // to be deleted on destructor
	std::vector<D_M* >                         m_Smat_uw_is;      // to be deleted on destructor
        D_M                                        m_Smat_uw; // Computed with 'experimentlModel' and 'simulationModel'
        D_M                                        m_Smat_uw_t;

	std::vector<D_M* >                         m_Rmat_u_hat_u_asterisk_is; // to be deleted on destructor
	std::vector<D_M* >                         m_Smat_u_hat_u_asterisk_is; // to be deleted on destructor
        D_M                                        m_Smat_u_hat_u_asterisk;
        D_M                                        m_Smat_u_hat_u_asterisk_t;

	std::vector<D_M* >                         m_Rmat_w_hat_u_asterisk_is; // to be deleted on destructor
	std::vector<D_M* >                         m_Smat_w_hat_u_asterisk_is; // to be deleted on destructor
        D_M                                        m_Smat_w_hat_u_asterisk;
        D_M                                        m_Smat_w_hat_u_asterisk_t;

        unsigned int                               m_vu_size;
        VectorSpace<D_V,D_M>                m_vu_space;

        VectorSpace<D_V,D_M>                m_unique_vu_space;
        unsigned int                               m_predVU_counter;
#if 0
        P_V                                        m_predVU_samplingRVs_unique_vu_meanVec; // todo_rr
        P_M                                        m_predVU_samplingRVs_unique_vu_covMatrix;
        P_M                                        m_predVU_samplingRVs_unique_vu_corrMatrix;
#endif
        P_V                                        m_predVU_summingRVs_unique_vu_meanVec; // todo_rr0
        P_M                                        m_predVU_summingRVs_mean_of_unique_vu_covMatrices;
        P_M                                        m_predVU_summingRVs_covMatrix_of_unique_vu_means;
        P_M                                        m_predVU_summingRVs_corrMatrix_of_unique_vu_means;
#if 0
        P_V                                        m_predVU_atMean_unique_vu_meanVec; // todo_rr
        P_M                                        m_predVU_atMean_unique_vu_covMatrix;
        P_M                                        m_predVU_atMean_unique_vu_corrMatrix;

        P_V                                        m_predVU_atMedian_unique_vu_meanVec; // todo_rr
        P_M                                        m_predVU_atMedian_unique_vu_covMatrix;
        P_M                                        m_predVU_atMedian_unique_vu_corrMatrix;

        P_V                                        m_predVU_atMode_unique_vu_meanVec; // todo_rr
        P_M                                        m_predVU_atMode_unique_vu_covMatrix;
        P_M                                        m_predVU_atMode_unique_vu_corrMatrix;

        P_V                                        m_predVU_atMLE_unique_vu_meanVec; // todo_rr
        P_M                                        m_predVU_atMLE_unique_vu_covMatrix;
        P_M                                        m_predVU_atMLE_unique_vu_corrMatrix;
#endif

        unsigned int                               m_omega_size;
        VectorSpace<D_V,D_M>                m_omega_space;

        D_V                                        m_Zvec_hat_vu;
        D_M                                        m_Smat_u;  // Computed with 'experimentModel'

        D_M*                                       m_Bmat_with_permut;    // to be deleted on destructor
        D_M*                                       m_Bmat_without_permut; // to be deleted on destructor
        unsigned int                               m_Bmat_rank;
        D_M*                                       m_Bwp_t__Wy__Bwp;      // to be deleted on destructor
        D_M*                                       m_Bop_t__Wy__Bop;      // to be deleted on destructor
        D_M*                                       m_Bwp_t__Wy__Bwp__inv; // to be deleted on destructor
        D_M*                                       m_Bop_t__Wy__Bop__inv; // to be deleted on destructor

        double                                     m_a_y_modifier;
        double                                     m_b_y_modifier;
};

}  // End namespace QUESO

#endif // UQ_GCM_JOINT_INFO_H
