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

#ifndef UQ_GCM_Z_INFO_H
#define UQ_GCM_Z_INFO_H

#include <queso/GcmSimulationInfo.h>
#include <queso/GcmExperimentInfo.h>
#include <queso/GcmJointInfo.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
class GcmZInfo
{
public:
  // Case with no experiments
  GcmZInfo(bool                                                             formCMatrix,
                  bool                                                             allOutputsAreScalar,
                  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>&         s);
  // Case with experiments, and with vector outputs
  GcmZInfo(bool                                                             formCMatrix,
                  bool                                                             allOutputsAreScalar,
                  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>&         s,
                  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>&         e,
                  const GcmJointInfo     <S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& jj);
  // Case with experiments, and with scalar outputs
  GcmZInfo(bool                                                             allOutputsAreScalar,
                  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>&         s,
                  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>&         e);
 ~GcmZInfo();

  const BaseEnvironment&     m_env;
        unsigned int                m_z_size;
        VectorSpace<D_V,D_M> m_z_space;
        D_V                         m_Zvec_hat;

        D_M*                        m_Cmat; // to be deleted on destructor
        unsigned int                m_Cmat_rank;

        D_M                         m_tmp_Smat_z;
        D_M                         m_tmp_Smat_extra;
        D_M                         m_tmp_Smat_z_hat;
        D_M                         m_tmp_Smat_z_hat_inv;

private:
  void commonConstructor();
};

}  // End namespace QUESO

#endif // UQ_GCM_Z_INFO_H
