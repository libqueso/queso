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

#include <queso/GcmSimulationTildeInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::GcmSimulationTildeInfo(
  const GpmsaComputerModelOptions&                  gcmOptionsObj,
  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s)
  :
  m_env                 (s.m_env),
  m_Kmat_tilde          (m_env,s.m_eta_space.map(),s.m_Kmat_rank),
  m_w_tilde_space       (m_env, "w_tilde_", s.m_Kmat_rank, NULL),         // rr0 check
  m_Lkmat               (m_env,m_w_tilde_space.map(),s.m_Kmat.numCols()), // rr0 check
  m_Ktildet_Ktilde      (m_w_tilde_space.zeroVector()),
  m_Ktildet_Ktilde_inv  (m_w_tilde_space.zeroVector()),
  m_Zvec_tilde_hat_w    (m_w_tilde_space.zeroVector()),
  m_a_eta_modifier_tilde(0),
  m_b_eta_modifier_tilde(0)
{
  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

    //******************************************************************************
    // Tilde situation: form 'm_Kmat_tilde'
    // Tilde situation: form 'm_w_tilde_space'
    // Tilde situation: form 'm_Lkmat'
    //******************************************************************************
    if (s.m_Kmat.numRowsGlobal() >= s.m_Kmat.numCols()) {
      Q_M matkU(s.m_Kmat.svdMatU());
      unsigned int kuMatRank   = matkU.rank(0.,1.e-8 ); // todo: should be an option
      unsigned int kuMatRank14 = matkU.rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": matkU.numRowsLocal() = "  << matkU.numRowsLocal()
                                << ", matkU.numCols() = "       << matkU.numCols()
                                << ", matkU.rank(0.,1.e-8) = "  << kuMatRank
                                << ", matkU.rank(0.,1.e-14) = " << kuMatRank14
                                << std::endl;
      }

      if (m_env.checkingLevel() >= 1) {
        Q_M matkUcheck(s.m_w_space.zeroVector());
        Q_V vecI(s.m_eta_space.zeroVector());
        Q_V vecJ(s.m_eta_space.zeroVector());
        for (unsigned int i = 0; i < matkU.numCols(); ++i) {
          matkU.getColumn(i,vecI);
          for (unsigned int j = i; j < matkU.numCols(); ++j) {
            matkU.getColumn(j,vecJ);
            matkUcheck(i,j) = scalarProduct(vecI,vecJ);
          }
        }
        matkUcheck.setPrintHorizontally(false);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                  << ": matkUcheck.numRowsLocal() = " << matkUcheck.numRowsLocal()
                                  << ", matkUcheck.numCols() = "      << matkUcheck.numCols()
                                  << ", m_Kmat_rank = "               << s.m_Kmat_rank
                                  << ", matkUcheck =\n"               << matkUcheck
                                  << std::endl;
        }
      }

      Q_V veckJ(s.m_eta_space.zeroVector());
      for (unsigned int j = 0; j < s.m_Kmat_rank; ++j) {
        matkU.getColumn(j,veckJ);
        m_Kmat_tilde.setColumn(j,veckJ);
      }
      if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
        m_Kmat_tilde.subWriteContents("Ktilde",
                                      "Ktilde1",
                                      "m",
                                      tmpSet);
      }
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Kmat_tilde computed (1)"
                                << std::endl;
      }
    }
    else {
      queso_error_msg("code for horizontal 'm_Kmat' is not ready yet");
      if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
        m_Kmat_tilde.subWriteContents("Ktilde",
                                      "Ktilde2",
                                      "m",
                                      tmpSet);
      }
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Kmat_tilde computed (2)"
                                << std::endl;
      }
    }

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished forming 'm_Kmat_tilde'"
                              << ", m_Kmat_tilde.numRowsLocal() = " << m_Kmat_tilde.numRowsLocal()
                              << ", m_Kmat_tilde.numCols() = "      << m_Kmat_tilde.numCols()
                              << std::endl;
    }

    m_Kmat_tilde.svdSolve(s.m_Kmat,m_Lkmat);
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Lkmat.subWriteContents("Lkmat",
                               "Lkmat",
                               "m",
                               tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Lkmat_tilde computed"
                              << std::endl;
    }

    //********************************************************************************
    // Form 'Ktilde^T' matrix
    //********************************************************************************
    Q_M Ktildet(m_env,m_w_tilde_space.map(),s.m_eta_size);
    Ktildet.fillWithTranspose(0,0,s.m_Kmat,true,true);

    //********************************************************************************
    // Compute 'Ktilde' rank
    //********************************************************************************
    unsigned int kTildeRank = m_Kmat_tilde.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int kRank14    = m_Kmat_tilde.rank(0.,1.e-14);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Kmat_tilde.numRowsLocal() = "  << m_Kmat_tilde.numRowsLocal()
                              << ", m_Kmat_tilde.numCols() = "       << m_Kmat_tilde.numCols()
                              << ", m_Kmat_tilde.rank(0.,1.e-8) = "  << kTildeRank
                              << ", m_Kmat_tilde.rank(0.,1.e-14) = " << kRank14
                              << std::endl;
    }

    queso_require_equal_to_msg(kTildeRank, std::min(m_Kmat_tilde.numRowsGlobal(),m_Kmat_tilde.numCols()), "'m_Kmat_tilde' does not have a proper rank");

    //******************************************************************************
    // Tilde situation: compute 'Ktilde^T Ktilde' matrix, and its inverse
    //******************************************************************************
    m_Ktildet_Ktilde = Ktildet * m_Kmat_tilde;
    if (m_env.subDisplayFile()) {
      m_Ktildet_Ktilde.setPrintHorizontally(false);
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Ktildet_Ktilde'"
        //<< "\n m_Ktildet_Ktilde =\n" << m_Ktildet_Ktilde
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Ktildet_Ktilde.subWriteContents("Ktildet_Ktilde",
                                        "Ktildet_Ktilde",
                                        "m",
                                        tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Ktildet_Ktilde computed"
                              << std::endl;
    }

    double       ktildetKtildeLnDeterminant = m_Ktildet_Ktilde.lnDeterminant();
    unsigned int ktildetKtildeRank          = m_Ktildet_Ktilde.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int ktildetKtildeRank14        = m_Ktildet_Ktilde.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Ktildet_Ktilde.numRowsLocal() = "  << m_Ktildet_Ktilde.numRowsLocal()
                              << ", m_Ktildet_Ktilde.numCols() = "       << m_Ktildet_Ktilde.numCols()
                              << ", m_Ktildet_Ktilde.lnDeterminant() = " << ktildetKtildeLnDeterminant
                              << ", m_Ktildet_Ktilde.rank(0.,1.e-8) = "  << ktildetKtildeRank
                              << ", m_Ktildet_Ktilde.rank(0.,1.e-14) = " << ktildetKtildeRank14
                              << std::endl;
    }

    m_Ktildet_Ktilde_inv = m_Ktildet_Ktilde.inverse();
    if (m_env.subDisplayFile()) {
      m_Ktildet_Ktilde_inv.setPrintHorizontally(false);
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Ktildet_Ktilde_inv'"
        //<< "\n m_Ktildet_Ktilde_inv =\n" << m_Ktildet_Ktilde_inv
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Ktildet_Ktilde_inv.subWriteContents("Ktildet_Ktilde_inv",
                                            "Ktildet_Ktilde_inv",
                                            "m",
                                            tmpSet);
    }
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Ktildet_Ktilde_inv computed"
                              << std::endl;
    }

    double       ktildetKtildeInvLnDeterminant = m_Ktildet_Ktilde_inv.lnDeterminant();
    unsigned int ktildetKtildeInvRank          = m_Ktildet_Ktilde_inv.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int ktildetKtildeInvRank14        = m_Ktildet_Ktilde_inv.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Ktildet_Ktilde_inv.numRowsLocal() = "  << m_Ktildet_Ktilde_inv.numRowsLocal()
                              << ", m_Ktildet_Ktilde_inv.numCols() = "       << m_Ktildet_Ktilde_inv.numCols()
                              << ", m_Ktildet_Ktilde_inv.lnDeterminant() = " << ktildetKtildeInvLnDeterminant
                              << ", m_Ktildet_Ktilde_inv.rank(0.,1.e-8) = "  << ktildetKtildeInvRank
                              << ", m_Ktildet_Ktilde_inv.rank(0.,1.e-14) = " << ktildetKtildeInvRank14
                              << std::endl;
    }

    //********************************************************************************
    // Compute 'tilde' exponent modifiers
    //********************************************************************************
    m_a_eta_modifier_tilde = ((double) s.m_paper_m) * ((double) (s.m_paper_n_eta - s.m_paper_p_eta)) / 2.;
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_a_eta_modifier_tilde = " << m_a_eta_modifier_tilde
                              << std::endl;
    }

    Q_V etaVec_transformed(s.m_simulationModel.etaVec_transformed("Gp.h.003"));
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Zvec_tilde_hat_w.sizeLocal() = "      << m_Zvec_tilde_hat_w.sizeLocal()
                              << ", m_Ktildet_Ktilde_inv.numRowsLocal() = " << m_Ktildet_Ktilde_inv.numRowsLocal()
                              << ", m_Ktildet_Ktilde_inv.numCols() = "      << m_Ktildet_Ktilde_inv.numCols()
                              << ", Ktildet.numRowsLocal() = "              << Ktildet.numRowsLocal()
                              << ", Ktildet.numCols() = "                   << Ktildet.numCols()
                              << ", etaVec_transformed.sizeLocal() = "      << etaVec_transformed.sizeLocal()
                              << std::endl;
    }
    m_Zvec_tilde_hat_w = m_Ktildet_Ktilde_inv * (Ktildet * etaVec_transformed);
    Q_V tmpVec1(etaVec_transformed - (m_Kmat_tilde * m_Zvec_tilde_hat_w));
    m_b_eta_modifier_tilde = scalarProduct(etaVec_transformed,tmpVec1) / 2.;

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_b_eta_modifier_tilde = " << m_b_eta_modifier_tilde
                              << std::endl;
    }
}

template <class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::~GcmSimulationTildeInfo()
{
}

}  // End namespace QUESO

template class QUESO::GcmSimulationTildeInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
