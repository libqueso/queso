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

#ifndef UQ_GCM_EXPERIMENT_INFO_H
#define UQ_GCM_EXPERIMENT_INFO_H

#include <queso/ExperimentStorage.h>
#include <queso/ExperimentModel.h>
#include <queso/GpmsaComputerModelOptions.h>
#include <queso/VectorRV.h>
#include <queso/GammaVectorRV.h>
#include <queso/BetaVectorRV.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class D_V = GslVector, class D_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix>
class GcmExperimentInfo
{
public:
  GcmExperimentInfo(const GpmsaComputerModelOptions&          gcmOptionsObj,
                                 bool                                       allOutputsAreScalar,
                           const ExperimentStorage<S_V,S_M,D_V,D_M>& experimentStorage,
                           const ExperimentModel  <S_V,S_M,D_V,D_M>& experimentModel,
                           const BaseVectorRV     <P_V,P_M>&         thetaPriorRv);
 ~GcmExperimentInfo();

  const BaseEnvironment&                           m_env;
  const ExperimentStorage<S_V,S_M,D_V,D_M>&        m_experimentStorage;
  const ExperimentModel  <S_V,S_M,D_V,D_M>&        m_experimentModel;

        unsigned int                                      m_paper_p_x;
        unsigned int                                      m_paper_n;
  const std::vector<const S_V*>&                          m_paper_xs_standard;
	std::vector<unsigned int>                         m_paper_n_ys_transformed;
        unsigned int                                      m_paper_n_y;
        unsigned int                                      m_paper_p_delta;
        unsigned int                                      m_paper_F;
        std::vector<unsigned int>                         m_paper_Gs;
        VectorSpace         <P_V,P_M>              m_paper_n_space;

        unsigned int                                      m_5lambdaYDim; // '1' in paper
        VectorSpace         <P_V,P_M>              m_5lambdaYSpace;
        P_V                                               m_5lambdaYMins;
        P_V                                               m_5lambdaYMaxs;
        BoxSubset           <P_V,P_M>              m_5lambdaYDomain;
        P_V                                               m_5lambdaYGammaAVec;
        P_V                                               m_5lambdaYGammaBVec;
        GammaVectorRV       <P_V,P_M>              m_5lambdaYPriorRv;
        P_V                                               m_like_previous5;
        P_V                                               m_tmp_5lambdaYVec;

        unsigned int                                      m_6lambdaVDim; // 'F' in paper
        VectorSpace         <P_V,P_M>              m_6lambdaVSpace;
        P_V                                               m_6lambdaVMins;
        P_V                                               m_6lambdaVMaxs;
        BoxSubset           <P_V,P_M>              m_6lambdaVDomain;
        P_V                                               m_6lambdaVGammaAVec;
        P_V                                               m_6lambdaVGammaBVec;
        GammaVectorRV       <P_V,P_M>              m_6lambdaVPriorRv;
        P_V                                               m_like_previous6;
        P_V                                               m_tmp_6lambdaVVec;

        unsigned int                                      m_7rhoVDim; // 'F * p_x' in paper
        VectorSpace         <P_V,P_M>              m_7rhoVSpace;
        P_V                                               m_7rhoVMins;
        P_V                                               m_7rhoVMaxs;
        BoxSubset           <P_V,P_M>              m_7rhoVDomain;
        P_V                                               m_7rhoVBetaAVec;
        P_V                                               m_7rhoVBetaBVec;
        BetaVectorRV        <P_V,P_M>              m_7rhoVPriorRv;
        P_V                                               m_like_previous7;
        P_V                                               m_tmp_7rhoVVec;

        unsigned int                                      m_8thetaDim;
        VectorSpace         <P_V,P_M>              m_8thetaSpace;
  const BaseVectorRV        <P_V,P_M>&             m_8thetaPriorRv;
        P_V                                               m_like_previous8;
        P_V                                               m_tmp_8thetaVec;

        unsigned int                                      m_v_size;
        VectorSpace<D_V,D_M>                       m_v_space;
        VectorSpace<D_V,D_M>                       m_unique_v_space;
        VectorSpace<P_V,P_M>                       m_rho_v_space;
        VectorSpace<D_V,D_M>                       m_y_space;

        P_V                                               m_tmp_rho_v_vec;
	std::vector<VectorSpace<D_V,D_M>* >        m_Imat_v_i_spaces; // to be deleted on destructor
	std::vector<D_M* >                                m_Imat_v_is;       // to be deleted on destructor

	std::vector<VectorSpace<D_V,D_M>* >        m_Rmat_v_i_spaces; // to be deleted on destructor
	std::vector<D_M* >                                m_Rmat_v_is;       // to be deleted on destructor

	std::vector<VectorSpace<D_V,D_M>* >        m_Smat_v_i_spaces; // to be deleted on destructor
	std::vector<D_M* >                                m_Smat_v_is;       // to be deleted on destructor
        D_M                                               m_Smat_v; // Computed with 'experimentModel'

	std::vector<D_M* >                                m_Rmat_v_hat_v_asterisk_is; // to be deleted on destructor
	std::vector<D_M* >                                m_Smat_v_hat_v_asterisk_is; // to be deleted on destructor
        D_M                                               m_Smat_v_hat_v_asterisk;
        D_M                                               m_Smat_v_hat_v_asterisk_t;

        D_M*                                              m_PD;                    // to be deleted on destructor
        const D_M*                                        m_Dmat_BlockDiag;        // Equal to '&experimentModel.Dmat_BlockDiag()'
        D_M*                                              m_Dmat_BlockDiag_permut; // to be deleted on destructor
        const D_M*                                        m_Wy;                    // Equal to '&experimentStorage.Wy()'

        D_M                                               m_Smat_v_asterisk_v_asterisk;
};

}  // End namespace QUESO

#endif // UQ_GCM_EXPERIMENT_INFO_H
