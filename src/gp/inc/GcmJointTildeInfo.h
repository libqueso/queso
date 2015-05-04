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

#ifndef UQ_GCM_JOINT_TILDE_INFO_H
#define UQ_GCM_JOINT_TILDE_INFO_H

#include <queso/GpmsaComputerModelOptions.h>
#include <queso/GcmExperimentInfo.h>
#include <queso/GcmJointInfo.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class D_V = GslVector,
  class D_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix,
  class Q_V = GslVector, class Q_M = GslMatrix>
class GcmJointTildeInfo
{
public:
  GcmJointTildeInfo(const GpmsaComputerModelOptions&                     gcmOptionsObj,
                           const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>&    e,
                           const GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& jj);
 ~GcmJointTildeInfo();

  const BaseEnvironment&     m_env;
        D_M                         m_Bmat_tilde;
        unsigned int                m_Bmat_tilde_rank;
        VectorSpace<D_V,D_M> m_vu_tilde_space;
        D_M                         m_Lbmat;
        D_M                         m_Btildet_Wy_Btilde;
        D_M                         m_Btildet_Wy_Btilde_inv;
        D_V                         m_Zvec_tilde_hat_vu;
        unsigned int                m_a_y_modifier_tilde;
        unsigned int                m_b_y_modifier_tilde;
};

}  // End namespace QUESO

#endif // UQ_GCM_JOINT_TILDE_INFO_H
