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

#ifndef UQ_GCM_SIMULATION_TILDE_INFO_H
#define UQ_GCM_SIMULATION_TILDE_INFO_H

#include <queso/SimulationStorage.h>
#include <queso/SimulationModel.h>
#include <queso/VectorRV.h>
#include <queso/GpmsaComputerModelOptions.h>
#include <queso/GcmSimulationInfo.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class GcmSimulationTildeInfo
{
public:
  GcmSimulationTildeInfo(const GpmsaComputerModelOptions&                  gcmOptionsObj,
                                const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s);
 ~GcmSimulationTildeInfo();

  const BaseEnvironment&     m_env;
        Q_M                         m_Kmat_tilde;
        VectorSpace<Q_V,Q_M> m_w_tilde_space;
        Q_M                         m_Lkmat;
        Q_M                         m_Ktildet_Ktilde;
        Q_M                         m_Ktildet_Ktilde_inv;
        Q_V                         m_Zvec_tilde_hat_w;
        unsigned int                m_a_eta_modifier_tilde;
        unsigned int                m_b_eta_modifier_tilde;
};

}  // End namespace QUESO

#endif // UQ_GCM_SIMULATION_TILDE_INFO_H
