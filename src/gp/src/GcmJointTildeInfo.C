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

#include <queso/GcmJointTildeInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
  GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GcmJointTildeInfo(
  const GpmsaComputerModelOptions&                     gcmOptionsObj,
  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>&    e,
  const GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>& jj)
  :
  m_env                  (jj.m_env),
  m_Bmat_tilde           (m_env,e.m_y_space.map(),jj.m_Bmat_rank),
  m_Bmat_tilde_rank      (0),
  m_vu_tilde_space       (m_env, "vu_tilde_", jj.m_Bmat_rank, NULL),          // rr0 check
  m_Lbmat                (m_env,m_vu_tilde_space.map(),jj.m_Bmat_with_permut->numCols()), // rr0 check
  m_Btildet_Wy_Btilde    (m_vu_tilde_space.zeroVector()),
  m_Btildet_Wy_Btilde_inv(m_vu_tilde_space.zeroVector()),
  m_Zvec_tilde_hat_vu    (m_vu_tilde_space.zeroVector()),
  m_a_y_modifier_tilde   (0),
  m_b_y_modifier_tilde   (0)
{
  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());
    //******************************************************************************
    // Tilde situation: form 'm_Bmat_tilde'
    // Tilde situation: form 'm_vu_tilde_space'
    // Tilde situation: form 'm_Lbmat'
    //******************************************************************************
    if (jj.m_Bmat_with_permut->numRowsGlobal() >= jj.m_Bmat_with_permut->numCols()) {
      D_M matbU(m_env,e.m_y_space.map(),jj.m_vu_size); // same as m_Bmat
      matbU = jj.m_Bmat_with_permut->svdMatU();
      unsigned int buMatRank   = matbU.rank(0.,1.e-8 ); // todo: should be an option
      unsigned int buMatRank14 = matbU.rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": matbU.numRowsLocal() = "  << matbU.numRowsLocal()
                                << ", matbU.numCols() = "       << matbU.numCols()
                                << ", matbU.rank(0.,1.e-8) = "  << buMatRank
                                << ", matbU.rank(0.,1.e-14) = " << buMatRank14
                                << std::endl;
      }

      if (m_env.checkingLevel() >= 1) {
        D_M matbUcheck(jj.m_vu_space.zeroVector());
        D_V vecI(e.m_y_space.zeroVector());
        D_V vecJ(e.m_y_space.zeroVector());
        for (unsigned int i = 0; i < matbU.numCols(); ++i) {
          matbU.getColumn(i,vecI);
          for (unsigned int j = i; j < matbU.numCols(); ++j) {
            matbU.getColumn(j,vecJ);
            matbUcheck(i,j) = scalarProduct(vecI,vecJ);
          }
        }
        matbUcheck.setPrintHorizontally(false);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                  << ": m_Bmat_with_permut->numRowsLocal() = " << jj.m_Bmat_with_permut->numRowsLocal()
                                  << ", m_Bmat_with_permut->numCols() = "      << jj.m_Bmat_with_permut->numCols()
                                  << ", m_Bmat_rank = "                       << jj.m_Bmat_rank
                                  << ", matbU.numRowsLocal() = "              << matbU.numRowsLocal()
                                  << ", matbU.numCols() = "                   << matbU.numCols()
                                  << ", buMatrank(0.,1.e-8) = "               << buMatRank
                                  << ", buMatrank(0.,1.e-14) = "              << buMatRank14
                                  << ", matbUcheck.numRowsLocal() = "         << matbUcheck.numRowsLocal()
                                  << ", matbUcheck.numCols() = "              << matbUcheck.numCols()
                                  << ", matbUcheck =\n"                       << matbUcheck
                                  << std::endl;
        }
      }

      D_V vecbJ(e.m_y_space.zeroVector());
      for (unsigned int j = 0; j < jj.m_Bmat_rank; ++j) {
        matbU.getColumn(j,vecbJ);
        m_Bmat_tilde.setColumn(j,vecbJ);
      }
      if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
        m_Bmat_tilde.subWriteContents("Btilde",
                                       "Btilde1",
                                       "m",
                                       tmpSet);
      }
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Bmat_tilde computed (1)"
                                << std::endl;
      }
    }
    else {
      D_M Bmat_t(m_env,jj.m_vu_space.map(),e.m_paper_n_y); // same as m_Bmat^T
      Bmat_t.fillWithTranspose(0,0,*jj.m_Bmat_with_permut,true,true);

      D_M matbV(e.m_y_space.zeroVector());
      matbV = Bmat_t.svdMatV();
      unsigned int bvMatRank   = matbV.rank(0.,1.e-8); // todo: should be an option
      unsigned int bvMatRank14 = matbV.rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": matbV.numRowsLocal() = "  << matbV.numRowsLocal()
                                << ", matbV.numCols() = "       << matbV.numCols()
                                << ", matbV.rank(0.,1.e-8) = "  << bvMatRank
                                << ", matbV.rank(0.,1.e-14) = " << bvMatRank14
                                << std::endl;
      }

      if (m_env.checkingLevel() >= 1) {
        D_M matbVcheck(e.m_y_space.zeroVector());
        D_V vecI(e.m_y_space.zeroVector());
        D_V vecJ(e.m_y_space.zeroVector());
        for (unsigned int i = 0; i < matbV.numCols(); ++i) {
          matbV.getColumn(i,vecI);
          for (unsigned int j = i; j < matbV.numCols(); ++j) {
            matbV.getColumn(j,vecJ);
            matbVcheck(i,j) = scalarProduct(vecI,vecJ);
          }
        }
        matbVcheck.setPrintHorizontally(false);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                  << ": m_Bmat_with_permut->numRowsLocal() = " << jj.m_Bmat_with_permut->numRowsLocal()
                                  << ", m_Bmat_with_permut->numCols() = "      << jj.m_Bmat_with_permut->numCols()
                                  << ", m_Bmat_rank = "                       << jj.m_Bmat_rank
                                  << ", matbV.numRowsLocal() = "              << matbV.numRowsLocal()
                                  << ", matbV.numCols() = "                   << matbV.numCols()
                                  << ", bvMatrank(0.,1.e-8) = "               << bvMatRank
                                  << ", bvMatrank(0.,1.e-14) = "              << bvMatRank14
                                  << ", matbVcheck.numRowsLocal() = "         << matbVcheck.numRowsLocal()
                                  << ", matbVcheck.numCols() = "              << matbVcheck.numCols()
                                  << ", matbVcheck =\n"                       << matbVcheck
                                  << std::endl;
        }
      }

      D_V vecbJ(e.m_y_space.zeroVector());
      for (unsigned int j = 0; j < jj.m_Bmat_rank; ++j) {
        matbV.getColumn(j,vecbJ);
        m_Bmat_tilde.setColumn(j,vecbJ);
      }
      if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
        m_Bmat_tilde.subWriteContents("Btilde",
                                       "Btilde2",
                                       "m",
                                       tmpSet);
      }
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Bmat_tilde computed (2)"
                                << std::endl;
      }
    }

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished forming 'm_Bmat_tilde'"
                              << ", m_Bmat_tilde.numRowsLocal() = " << m_Bmat_tilde.numRowsLocal()
                              << ", m_Bmat_tilde.numCols() = "      << m_Bmat_tilde.numCols()
                              << std::endl;
    }

    m_Bmat_tilde.svdSolve(*jj.m_Bmat_with_permut,m_Lbmat);
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Lbmat.subWriteContents("Lbmat",
                                "Lbmat",
                                "m",
                                tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Lbmat_tilde computed"
                              << std::endl;
    }

    //********************************************************************************
    // Form 'Btilde^T' matrix
    //********************************************************************************
    D_M Btildet(m_env,m_vu_tilde_space.map(),e.m_paper_n_y);
    Btildet.fillWithTranspose(0,0,m_Bmat_tilde,true,true);

    if (m_env.checkingLevel() >= 1) {
      // Check transpose operation
      D_M Btildett(m_env,e.m_y_space.map(),m_vu_tilde_space.dimGlobal());
      Btildett.fillWithTranspose(0,0,Btildet,true,true);
      Btildett -= m_Bmat_tilde;
      double btDiffNorm = Btildett.normFrob();
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": ||Btildett - Btilde||_2 = " << btDiffNorm
                                << std::endl;
      }
    }

    //********************************************************************************
    // Compute 'Btilde' rank
    //********************************************************************************
    double bTildeRank14 = 0.;
    if (m_Bmat_tilde.numRowsGlobal() >= m_Bmat_tilde.numCols()) {
      m_Bmat_tilde_rank = m_Bmat_tilde.rank(0.,1.e-8 ); // todo: should be an option
      bTildeRank14      = m_Bmat_tilde.rank(0.,1.e-14);
    }
    else {
      m_Bmat_tilde_rank = Btildet.rank(0.,1.e-8 ); // todo: should be an option
      bTildeRank14      = Btildet.rank(0.,1.e-14);
    }

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Bmat_tilde.numRowsLocal() = "  << m_Bmat_tilde.numRowsLocal()
                              << ", m_Bmat_tilde.numCols() = "       << m_Bmat_tilde.numCols()
                              << ", m_Bmat_tilde.rank(0.,1.e-8) = "  << m_Bmat_tilde_rank
                              << ", m_Bmat_tilde.rank(0.,1.e-14) = " << bTildeRank14
                              << std::endl;
    }
    queso_require_equal_to_msg(m_Bmat_tilde_rank, std::min(m_Bmat_tilde.numRowsGlobal(),m_Bmat_tilde.numCols()), "'m_Bmat_tilde' does not have a proper rank");

    //******************************************************************************
    // Tilde situation: compute 'Btilde^T W_y Btilde' matrix, and its inverse
    //******************************************************************************
    m_Btildet_Wy_Btilde = Btildet * (*e.m_Wy * m_Bmat_tilde); // todo: add 1.e-4 to diagonal
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Btildet_Wy_Btilde'"
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Btildet_Wy_Btilde.subWriteContents("Btildet_Wy_Btilde",
                                            "Btildet_Wy_Btilde",
                                            "m",
                                            tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Btildet_Wy_Btilde computed"
                              << std::endl;
    }

    double       btildetWyBtildeLnDeterminant = m_Btildet_Wy_Btilde.lnDeterminant();
    unsigned int btildetWyBtildeRank          = m_Btildet_Wy_Btilde.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int btildetWyBtildeRank14        = m_Btildet_Wy_Btilde.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Btildet_Wy_Btilde.numRowsLocal() = "  << m_Btildet_Wy_Btilde.numRowsLocal()
                              << ", m_Btildet_Wy_Btilde.numCols() = "       << m_Btildet_Wy_Btilde.numCols()
                              << ", m_Btildet_Wy_Btilde.lnDeterminant() = " << btildetWyBtildeLnDeterminant
                              << ": m_Btildet_Wy_Btilde.rank(0.,1.e-8) = "  << btildetWyBtildeRank
                              << ": m_Btildet_Wy_Btilde.rank(0.,1.e-14) = " << btildetWyBtildeRank14
                              << std::endl;
    }

    m_Btildet_Wy_Btilde_inv = m_Btildet_Wy_Btilde.inverse(); // todo: add 1.e-6 to diagonal
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Btildet_Wy_Btilde_inv'"
                              << ", m_Btildet_Wy_Btilde_inv.lnDeterminant() = " << m_Btildet_Wy_Btilde_inv.lnDeterminant()
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Btildet_Wy_Btilde_inv.subWriteContents("Btildet_Wy_Btilde_inv",
                                                "Btildet_Wy_Btilde_inv",
                                                "m",
                                                tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Btildet_Wy_Btilde_inv computed"
                              << std::endl;
    }

    double       btildetWyBtildeInvLnDeterminant = m_Btildet_Wy_Btilde_inv.lnDeterminant();
    unsigned int btildetWyBtildeInvRank          = m_Btildet_Wy_Btilde_inv.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int btildetWyBtildeInvRank14        = m_Btildet_Wy_Btilde_inv.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Btildet_Wy_Btilde_inv.numRowsLocal() = "  << m_Btildet_Wy_Btilde_inv.numRowsLocal()
                              << ", m_Btildet_Wy_Btilde_inv.numCols() = "       << m_Btildet_Wy_Btilde_inv.numCols()
                              << ": m_Btildet_Wy_Btilde_inv.lnDeterminant() = " << btildetWyBtildeInvLnDeterminant
                              << ": m_Btildet_Wy_Btilde_inv.rank(0.,1.e-8) = "  << btildetWyBtildeInvRank
                              << ": m_Btildet_Wy_Btilde_inv.rank(0.,1.e-14) = " << btildetWyBtildeInvRank14
                              << std::endl;
    }

    //********************************************************************************
    // Compute 'tilde' exponent modifiers
    //********************************************************************************
    m_a_y_modifier_tilde   = ((double) (e.m_paper_n_y - m_Bmat_tilde_rank)) / 2.;
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_a_y_modifier_tilde = " << m_a_y_modifier_tilde
                              << std::endl;
    }

    D_V yVec_transformed(e.m_experimentStorage.yVec_transformed());
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Zvec_tilde_hat_vu.sizeLocal() = "        << m_Zvec_tilde_hat_vu.sizeLocal()
                              << ", m_Btildet_Wy_Btilde_inv.numRowsLocal() = " << m_Btildet_Wy_Btilde_inv.numRowsLocal()
                              << ", m_Btildet_Wy_Btilde_inv.numCols() = "      << m_Btildet_Wy_Btilde_inv.numCols()
                              << ", Btildet.numRowsLocal() = "                 << Btildet.numRowsLocal()
                              << ", Btildet.numCols() = "                      << Btildet.numCols()
                              << ", m_Wy->numRowsLocal() = "                   << e.m_Wy->numRowsLocal()
                              << ", m_Wy->numCols() = "                        << e.m_Wy->numCols()
                              << ", yVec_transformed.sizeLocal() = "           << yVec_transformed.sizeLocal()
                              << std::endl;
    }
    m_Zvec_tilde_hat_vu = m_Btildet_Wy_Btilde_inv * (Btildet * (*e.m_Wy * yVec_transformed));
    D_V tmpVec2(yVec_transformed - (m_Bmat_tilde * m_Zvec_tilde_hat_vu));
    tmpVec2 = *e.m_Wy * tmpVec2;
    m_b_y_modifier_tilde = scalarProduct(yVec_transformed,tmpVec2) / 2.;

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_b_y_modifier_tilde = " << m_b_y_modifier_tilde
                              << std::endl;
    }

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'tilde' exponent modifiers"
                              << std::endl;
    }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~GcmJointTildeInfo()
{
}

}  // End namespace QUESO

template class QUESO::GcmJointTildeInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
