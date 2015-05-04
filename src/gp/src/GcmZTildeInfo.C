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

#include <queso/GcmZTildeInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GcmZTildeInfo(
  const GpmsaComputerModelOptions&                     gcmOptionsObj,
  const GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& jj,
  const GcmZInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>&     z)
  :
  m_env                     (jj.m_env),
  m_Cmat_tilde              (m_env,jj.m_omega_space.map(),z.m_Cmat_rank),
  m_z_tilde_space           (m_env, "z_tilde_", z.m_Cmat_rank, NULL),
  m_Lmat                    (m_env,m_z_tilde_space.map(),z.m_Cmat->numCols()),
  m_Lmat_t                  (m_env,z.m_z_space.map(),z.m_Cmat_rank),
  m_Zvec_tilde_hat          (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_z_tilde        (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_extra_tilde    (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_z_tilde_hat    (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_z_tilde_hat_inv(m_z_tilde_space.zeroVector())
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

    // Naive formation of 'm_Cmat_tilde'
    D_M matU(z.m_Cmat->svdMatU());
    unsigned int uMatRank   = matU.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int uMatRank14 = matU.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                              << ": matU.numRowsLocal() = "  << matU.numRowsLocal()
                              << ", matU.numCols() = "       << matU.numCols()
                              << ", matU.rank(0.,1.e-8) = "  << uMatRank
                              << ", matU.rank(0.,1.e-14) = " << uMatRank14
                              << std::endl;
    }

    if (m_env.checkingLevel() >= 1) {
      D_M matUcheck(z.m_z_space.zeroVector());
      D_V vecI(jj.m_omega_space.zeroVector());
      D_V vecJ(jj.m_omega_space.zeroVector());
      for (unsigned int i = 0; i < matU.numCols(); ++i) {
        matU.getColumn(i,vecI);
        for (unsigned int j = i; j < matU.numCols(); ++j) {
          matU.getColumn(j,vecJ);
          matUcheck(i,j) = scalarProduct(vecI,vecJ);
        }
      }
      matUcheck.setPrintHorizontally(false);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                                << ": m_Cmat_rank = "              << z.m_Cmat_rank
                                << ", matUcheck.numRowsLocal() = " << matUcheck.numRowsLocal()
                                << ", matUcheck.numCols() = "      << matUcheck.numCols()
                                << ", matUcheck =\n"               << matUcheck
                                << std::endl;
      }
    }

    D_V vecJ(jj.m_omega_space.zeroVector());
    for (unsigned int j = 0; j < z.m_Cmat_rank; ++j) {
      matU.getColumn(j,vecJ);
      m_Cmat_tilde.setColumn(j,vecJ);
    }

    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Cmat_tilde.subWriteContents("Ctilde",
                                    "Ctilde2",
                                    "m",
                                    tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(1)"
                              << ": m_Cmat_tilde formed (2)"
                              << std::endl;
    }

    // Naive formation of 'm_Lmat'
    m_Cmat_tilde.svdSolve(*(z.m_Cmat),m_Lmat);

  commonConstructor(z);
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GcmZTildeInfo(
  const GpmsaComputerModelOptions&                               gcmOptionsObj,
  const GcmJointInfo          <S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& jj,
  const GcmZInfo              <S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& z,
  const GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>&         st,
  const GcmJointTildeInfo     <S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& jt)
  :
  m_env                     (jj.m_env),
  m_Cmat_tilde              (m_env,jj.m_omega_space.map(),z.m_Cmat_rank),
  m_z_tilde_space           (m_env, "z_tilde_", z.m_Cmat_rank, NULL),
  m_Lmat                    (m_env,m_z_tilde_space.map(),z.m_Cmat->numCols()),
  m_Lmat_t                  (m_env,z.m_z_space.map(),z.m_Cmat_rank),
  m_Zvec_tilde_hat          (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_z_tilde        (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_extra_tilde    (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_z_tilde_hat    (m_z_tilde_space.zeroVector()),
  m_tmp_Smat_z_tilde_hat_inv(m_z_tilde_space.zeroVector())
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

    m_Cmat_tilde.cwSet(0.);
    m_Cmat_tilde.cwSet(                             0,                        0,jt.m_Bmat_tilde);
    m_Cmat_tilde.cwSet(jt.m_Bmat_tilde.numRowsLocal(),jt.m_Bmat_tilde.numCols(),st.m_Kmat_tilde);
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Cmat_tilde.subWriteContents("Ctilde",
                                    "Ctilde1",
                                    "m",
                                    tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                              << ": m_Cmat_tilde formed (1)"
                              << std::endl;
    }
    m_Lmat.cwSet(0.);
    m_Lmat.cwSet(                        0,                   0,jt.m_Lbmat);
    m_Lmat.cwSet(jt.m_Lbmat.numRowsLocal(),jt.m_Lbmat.numCols(),st.m_Lkmat);
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Lmat.subWriteContents("Lmat",
                              "Lmat",
                              "m",
                              tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor(2)"
                              << ": m_Lmat_tilde formed"
                              << std::endl;
    }

  m_Zvec_tilde_hat.cwSetConcatenated(jt.m_Zvec_tilde_hat_vu,st.m_Zvec_tilde_hat_w); // todo_rrr: not in constructor(1) ???

  commonConstructor(z);
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~GcmZTildeInfo()
{
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor(const GcmZInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& z)
{
    unsigned int cMatTildeRank   = m_Cmat_tilde.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int cMatTildeRank14 = m_Cmat_tilde.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Cmat_tilde.numRowsLocal() = "  << m_Cmat_tilde.numRowsLocal()
                              << ", m_Cmat_tilde.numCols() = "       << m_Cmat_tilde.numCols()
                              << ", m_Cmat_tilde.rank(0.,1.e-8) = "  << cMatTildeRank
                              << ", m_Cmat_tilde.rank(0.,1.e-14) = " << cMatTildeRank14
                              << std::endl;
    }
    queso_require_equal_to_msg(cMatTildeRank, z.m_Cmat_rank, "'m_Cmat_tilde' should have full column rank");

    queso_require_less_msg(m_Lmat.numRowsGlobal(), m_Lmat.numCols(), "'m_Lmat' should be a 'horizontal' rectangular matrix");

    //******************************************************************************
    // Tilde situation: form 'm_Lmat_t'
    //******************************************************************************
    m_Lmat_t.fillWithTranspose(0,0,m_Lmat,true,true);
    unsigned int lMatRank   = m_Lmat_t.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int lMatRank14 = m_Lmat_t.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor()"
                              << ": m_Lmat.numRowsLocal() = "  << m_Lmat.numRowsLocal()
                              << ", m_Lmat.numCols() = "       << m_Lmat.numCols()
                              << ", m_Lmat.rank(0.,1.e-8) = "  << lMatRank
                              << ", m_Lmat.rank(0.,1.e-14) = " << lMatRank14
                              << std::endl;
    }

    queso_require_equal_to_msg(lMatRank, z.m_Cmat_rank, "'m_Lmat' should have full row rank");

    if (m_env.checkingLevel() >= 1) {
      // Check if C == C_tilde * L
      D_M tmpCmat(m_Cmat_tilde * m_Lmat);
      tmpCmat -= *z.m_Cmat;
      double cDiffNorm = tmpCmat.normFrob();
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::commonConstructor()"
                                << ": ||tmpC - C||_2 = " << cDiffNorm
                                << std::endl;
      }
    }

  return;
}

}  // End namespace QUESO

template class QUESO::GcmZTildeInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
