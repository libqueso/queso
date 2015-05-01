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

#include <queso/GcmJointInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GcmJointInfo(
  const GpmsaComputerModelOptions&                  gcmOptionsObj,
  bool                                                     allOutputsAreScalar,
  const GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>& s,
  const GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>& e)
  :
  m_env                                            (s.m_env),
  m_unique_u_space                                 (m_env, "unique_u_", s.m_paper_p_eta, NULL),
  m_Smat_u_asterisk_u_asterisk                     (m_unique_u_space.zeroVector()),
  m_u_size                                         (e.m_paper_n * s.m_paper_p_eta),
  m_u_space                                        (m_env, "u_", m_u_size, NULL),
  m_Rmat_u_is                                      (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_u_is                                      (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Rmat_uw_is                                     (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_uw_is                                     (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_uw                                        (m_env,  m_u_space.map(),s.m_w_size),
  m_Smat_uw_t                                      (m_env,s.m_w_space.map(),  m_u_size),
  m_Rmat_u_hat_u_asterisk_is                       (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_u_hat_u_asterisk_is                       (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_u_hat_u_asterisk                          (m_env, m_u_space.map(),        s.m_paper_p_eta),
  m_Smat_u_hat_u_asterisk_t                        (m_env, m_unique_u_space.map(), m_u_size),
  m_Rmat_w_hat_u_asterisk_is                       (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_w_hat_u_asterisk_is                       (s.m_paper_p_eta, (D_M*) NULL), // to be deleted on destructor
  m_Smat_w_hat_u_asterisk                          (m_env, s.m_w_space.map(),      s.m_paper_p_eta),
  m_Smat_w_hat_u_asterisk_t                        (m_env, m_unique_u_space.map(), s.m_w_size),
  m_vu_size                                        (e.m_v_size + m_u_size),
  m_vu_space                                       (m_env, "vu_", m_vu_size, NULL),
  m_unique_vu_space                                (m_env, "unique_vu_", e.m_paper_p_delta + s.m_paper_p_eta, NULL),
  m_predVU_counter                                 (MiscUintDebugMessage(0,NULL)),
  m_predVU_summingRVs_unique_vu_meanVec            (m_unique_vu_space.zeroVector()),
  m_predVU_summingRVs_mean_of_unique_vu_covMatrices(m_unique_vu_space.zeroVector()),
  m_predVU_summingRVs_covMatrix_of_unique_vu_means (m_unique_vu_space.zeroVector()),
  m_predVU_summingRVs_corrMatrix_of_unique_vu_means(m_unique_vu_space.zeroVector()),
  m_omega_size                                     (e.m_paper_n_y + s.m_eta_size),
  m_omega_space                                    (m_env, "omega_", m_omega_size, NULL),
  m_Zvec_hat_vu                                    (m_vu_space.zeroVector()),
  m_Smat_u                                         (m_u_space.zeroVector()),
  m_Bmat_with_permut                               (NULL), // to be deleted on destructor
  m_Bmat_without_permut                            (NULL), // to be deleted on destructor
  m_Bmat_rank                                      (0),
  m_Bwp_t__Wy__Bwp                                 (NULL), // to be deleted on destructor
  m_Bop_t__Wy__Bop                                 (NULL), // to be deleted on destructor
  m_Bwp_t__Wy__Bwp__inv                            (NULL), // to be deleted on destructor
  m_Bop_t__Wy__Bop__inv                            (NULL), // to be deleted on destructor
  m_a_y_modifier                                   (0.),
  m_b_y_modifier                                   (0.)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": allOutputsAreScalar = " << allOutputsAreScalar
                            << ": key-debug"
                            << ", some entities just created (not yet populated)"
                            << ", m_Smat_uw.numRowsLocal() = "   << m_Smat_uw.numRowsLocal()
                            << ", m_Smat_uw.numCols() = "        << m_Smat_uw.numCols()
                            << ", m_Smat_uw_t.numRowsLocal() = " << m_Smat_uw_t.numRowsLocal()
                            << ", m_Smat_uw_t.numCols() = "      << m_Smat_uw_t.numCols()
                            << ", m_Zvec_hat_vu.sizeLocal() = "  << m_Zvec_hat_vu.sizeLocal()
                            << ", m_Smat_u.numRowsLocal() = "    << m_Smat_u.numRowsLocal()
                            << ", m_Smat_u.numCols() = "         << m_Smat_u.numCols()
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

  if (allOutputsAreScalar) {
    // Do nothing
    // ppp: set m_Zvec_hat_vu?
  }
  else {
    //********************************************************************************
    // Form 'P_K' matrix
    //********************************************************************************
    D_M PK(m_u_space.zeroVector());
    for (unsigned int i = 0; i < s.m_paper_p_eta; ++i) {
      for (unsigned int j = 0; j < e.m_paper_n; ++j) {
        unsigned int row = j + (e.m_paper_n*i);
        unsigned int col = (j*s.m_paper_p_eta)+i;
        PK(row,col) = 1.;
      }
    }

    if (m_env.checkingLevel() >= 1) {
      // Check transpose operation
      D_M PKt(PK.transpose());
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", tests on PK"
                                << ": PK.numRowsLocal() = "  << PK.numRowsLocal()
                                << ", PK.numCols() = "       << PK.numCols()
                                << ": PKt.numRowsLocal() = " << PKt.numRowsLocal()
                                << ", PKt.numCols() = "      << PKt.numCols()
                                << std::endl;
      }

      D_M matShouldBeI1( PK * PKt );
      D_M matI1        (m_u_space.zeroVector());
      for (unsigned int i = 0; i < matI1.numRowsLocal(); ++i) {
        matI1(i,i) = 1.;
      }
      matShouldBeI1 -= matI1;
      double auxNorm1 = matShouldBeI1.normFrob();

      D_M matShouldBeI2( PKt * PK );
      D_M matI2        (m_u_space.zeroVector());
      for (unsigned int i = 0; i < matI2.numRowsLocal(); ++i) {
        matI2(i,i) = 1.;
      }
      matShouldBeI2 -= matI2;
      double auxNorm2 = matShouldBeI2.normFrob();

      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", tests on PK"
                                << ": matShouldBeI1.numRowsLocal() = "  << matShouldBeI1.numRowsLocal()
                                << ", ||matI1||_2^2 = "                 << matI1.normFrob() * matI1.normFrob()
                                << ", ||matShouldBeI1 - matI1||_2^2 = " << auxNorm1 * auxNorm1
                                << "; matShouldBeI2.numRowsLocal() = "  << matShouldBeI2.numRowsLocal()
                                << ", ||matI2||_2^2 = "                 << matI2.normFrob() * matI2.normFrob()
                                << ", ||matShouldBeI2 - matI2||_2^2 = " << auxNorm2 * auxNorm2
                                << std::endl;
      }
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished forming 'P_K'"
                              << std::endl;
    }

    //********************************************************************************
    // Form 'Kmat_interp_BlockDiag' matrix
    //********************************************************************************
    D_M Kmat_interp_BlockDiag       (m_env,e.m_y_space.map(),m_u_size); // Formed with 'experimentModel.Kmats_interp()' // rr0: use D_M& and copy from 'experimentModel' object
    D_M Kmat_interp_BlockDiag_permut(m_env,e.m_y_space.map(),m_u_size);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": Kmat_interp_BlockDiag.numRowsLocal() = " << Kmat_interp_BlockDiag.numRowsLocal()
                              << ", Kmat_interp_BlockDiag.numCols() = "      << Kmat_interp_BlockDiag.numCols()
                              << std::endl;
    }

    Kmat_interp_BlockDiag.fillWithBlocksDiagonally(0,0,e.m_experimentModel.Kmats_interp(),true,true);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished forming 'Kmat_interp_BlockDiag'"
                              << std::endl;
    }

    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      Kmat_interp_BlockDiag.subWriteContents("Kmat_interp_BlockDiag",
                                             "mat_Kmat_interp_BlockDiag",
                                             "m",
                                             tmpSet);
    }

    //********************************************************************************
    // Compute 'Kmat_interp_BlockDiag_permut' matrix
    //********************************************************************************
    Kmat_interp_BlockDiag_permut = Kmat_interp_BlockDiag * (PK.transpose());

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'Kmat_interp_BlockDiag_permut'"
                              << std::endl;
    }

    //********************************************************************************
    // Form 'B' matrix
    //********************************************************************************
    m_Bmat_with_permut    = new D_M(m_env,e.m_y_space.map(),m_vu_size); // to be deleted on destructor
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", m_Bmat_with_permut just created (not yet populated)"
                              << ", numRowsLocal() = " << m_Bmat_with_permut->numRowsLocal()
                              << ", numCols() = "      << m_Bmat_with_permut->numCols()
                              << std::endl;
    }
    m_Bmat_without_permut = new D_M(m_env,e.m_y_space.map(),m_vu_size); // to be deleted on destructor
    m_Bmat_rank           = std::min(m_Bmat_with_permut->numRowsGlobal(),m_Bmat_with_permut->numCols()); // Might be smaller
    m_Bwp_t__Wy__Bwp      = new D_M(m_vu_space.zeroVector()); // to be deleted on destructor
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", m_Bwp_t__Wy__Bwp just created (not yet populated)"
                              << ", numRowsLocal() = " << m_Bwp_t__Wy__Bwp->numRowsLocal()
                              << ", numCols() = "      << m_Bwp_t__Wy__Bwp->numCols()
                              << std::endl;
    }
    m_Bop_t__Wy__Bop      = new D_M(m_vu_space.zeroVector()); // to be deleted on destructor
    m_Bwp_t__Wy__Bwp__inv = new D_M(m_vu_space.zeroVector()); // to be deleted on destructor
    m_Bop_t__Wy__Bop__inv = new D_M(m_vu_space.zeroVector()); // to be deleted on destructor

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor():"
                              << "\n  m_Bmat_with_permut->numRowsLocal() = "          << m_Bmat_with_permut->numRowsLocal()
                              << "\n  m_Bmat_with_permut->numCols() = "               << m_Bmat_with_permut->numCols()
                              << "\n  m_Dmat_BlockDiag_permut->numRowsLocal() = "     << e.m_Dmat_BlockDiag_permut->numRowsLocal()
                              << "\n  m_Dmat_BlockDiag_permut->numCols() = "          << e.m_Dmat_BlockDiag_permut->numCols()
                              << "\n  m_PD->numRowsLocal() = "                        << e.m_PD->numRowsLocal()
                              << "\n  m_PD->numCols() = "                             << e.m_PD->numCols()
                              << "\n  Kmat_interp_BlockDiag_permut.numRowsLocal() = " << Kmat_interp_BlockDiag_permut.numRowsLocal()
                              << "\n  Kmat_interp_BlockDiag_permut.numCols() = "      << Kmat_interp_BlockDiag_permut.numCols()
                              << "\n  PK.numRowsLocal() = "                           << PK.numRowsLocal()
                              << "\n  PK.numCols() = "                                << PK.numCols()
                              << "\n  gcmOptionsObj.m_ov.m_nuggetValueForBtWyB = "    << gcmOptionsObj.m_ov.m_nuggetValueForBtWyB
                              << "\n  gcmOptionsObj.m_ov.m_nuggetValueForBtWyBInv = " << gcmOptionsObj.m_ov.m_nuggetValueForBtWyBInv
                              << std::endl;
    }

    std::vector<const P_M* > twoMats(2, (P_M*) NULL);

    twoMats[0] = e.m_Dmat_BlockDiag_permut;
    twoMats[1] = &Kmat_interp_BlockDiag_permut;
    m_Bmat_with_permut->fillWithBlocksHorizontally(0,0,twoMats,true,true);

    twoMats[0] = e.m_Dmat_BlockDiag;
    twoMats[1] = &Kmat_interp_BlockDiag;
    m_Bmat_without_permut->fillWithBlocksHorizontally(0,0,twoMats,true,true);

    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bmat_without_permut->subWriteContents("B_op",
                                              "mat_B_op",
                                              "m",
                                              tmpSet);
    }

    //********************************************************************************
    // Form 'Bwp^T' matrix
    //********************************************************************************
    D_M Bwp_t(m_env,m_vu_space.map(),e.m_paper_n_y);
    Bwp_t.fillWithTranspose(0,0,*m_Bmat_with_permut,true,true);

    if (m_env.checkingLevel() >= 1) {
      // Check transpose operation
      D_M Btt(m_env,e.m_y_space.map(),m_vu_size);
      Btt.fillWithTranspose(0,0,Bwp_t,true,true);
      Btt -= *m_Bmat_with_permut;
      double btDiffNorm = Btt.normFrob();
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": ||Btt - B||_2 = " << btDiffNorm
                                << std::endl;
      }
    }

    double bRank14 = 0.;
    if (m_Bmat_with_permut->numRowsGlobal() >= m_Bmat_with_permut->numCols()) {
      m_Bmat_rank = m_Bmat_with_permut->rank(0.,1.e-8 ); // todo: should be an option
      bRank14     = m_Bmat_with_permut->rank(0.,1.e-14);
    }
    else {
      m_Bmat_rank = Bwp_t.rank(0.,1.e-8); // todo: should be an option
      bRank14     = Bwp_t.rank(0.,1.e-14);
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished forming 'm_Bmat'"
                              << ", m_Bmat_with_permut->numRowsLocal() = "  << m_Bmat_with_permut->numRowsLocal()
                              << ", m_Bmat_with_permut->numCols() = "       << m_Bmat_with_permut->numCols()
                              << ", m_Bmat_with_permut->rank(0.,1.e-8) = "  << m_Bmat_rank
                              << ", m_Bmat_with_permut->rank(0.,1.e-14) = " << bRank14
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bmat_with_permut->subWriteContents("B_wp",
                                           "mat_B_wp",
                                           "m",
                                           tmpSet);
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      D_M Bmat_filtered(*m_Bmat_with_permut);
      Bmat_filtered.setPrintHorizontally(false);
      Bmat_filtered.filterSmallValues(1.e-6);
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": Bmat_filtered.numRowsLocal() = " << Bmat_filtered.numRowsLocal()
                              << ", Bmat_filtered.numCols() = "      << Bmat_filtered.numCols()
                              << ", Bmat_filtered contents =\n"      << Bmat_filtered
                              << std::endl;
    }

    //********************************************************************************
    // Compute 'Bwp^T W_y Bwp' matrix
    //********************************************************************************
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor():"
                              << "\n  m_Bwp_t__Wy__Bwp->numRowsLocal() = "     << m_Bwp_t__Wy__Bwp->numRowsLocal()
                              << "\n  m_Bwp_t__Wy__Bwp->numCols() = "          << m_Bwp_t__Wy__Bwp->numCols()
                              << "\n  Bwp_t.numRowsLocal() = "                 << Bwp_t.numRowsLocal()
                              << "\n  Bwp_t.numCols() = "                      << Bwp_t.numCols()
                              << "\n  m_Wy->numRowsLocal() = "                 << e.m_Wy->numRowsLocal()
                              << "\n  m_Wy->numCols() = "                      << e.m_Wy->numCols()
                              << "\n  m_Bmat_with_permut->numRowsLocal() = "   << m_Bmat_with_permut->numRowsLocal()
                              << "\n  m_Bmat_with_permut->numCols() = "        << m_Bmat_with_permut->numCols()
                              << std::endl;
    }
    *m_Bwp_t__Wy__Bwp = Bwp_t * (*e.m_Wy * *m_Bmat_with_permut);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Bwp_t__Wy__Bwp'"
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bwp_t__Wy__Bwp->subWriteContents("Bwp_t__Wy__Bwp",
                                         "mat_Bwp_t__Wy__Bwp",
                                         "m",
                                         tmpSet);
    }

    if (m_env.displayVerbosity() >= 4) {
      double       Bwp_t__Wy__Bwp__LnDeterminant = m_Bwp_t__Wy__Bwp->lnDeterminant();
      unsigned int Bwp_t__Wy__Bwp__Rank          = m_Bwp_t__Wy__Bwp->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bwp_t__Wy__Bwp__Rank14        = m_Bwp_t__Wy__Bwp->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Bwp_t__Wy__Bwp->numRowsLocal() = "  << m_Bwp_t__Wy__Bwp->numRowsLocal()
                                << ", m_Bwp_t__Wy__Bwp->numCols() = "       << m_Bwp_t__Wy__Bwp->numCols()
                                << ", m_Bwp_t__Wy__Bwp->lnDeterminant() = " << Bwp_t__Wy__Bwp__LnDeterminant
                                << ", m_Bwp_t__Wy__Bwp->rank(0.,1.e-8) = "  << Bwp_t__Wy__Bwp__Rank
                                << ", m_Bwp_t__Wy__Bwp->rank(0.,1.e-14) = " << Bwp_t__Wy__Bwp__Rank14
                                << std::endl;
      }
    }

    //********************************************************************************
    // Compute 'Bwp^T W_y Bwp' inverse
    //********************************************************************************
    *m_Bwp_t__Wy__Bwp__inv = m_Bwp_t__Wy__Bwp->inverse(); // inversion savings
    if (m_env.displayVerbosity() >= 4) {
      double Bwp_t__Wy__Bwp__inv__LnDeterminant = m_Bwp_t__Wy__Bwp__inv->lnDeterminant();
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": finished computing 'm_Bwp_t__Wy__Bwp__inv'"
                                << ", m_Bwp_t__Wy__Bwp__inv->lnDeterminant() = " << Bwp_t__Wy__Bwp__inv__LnDeterminant
                                << std::endl;
      }
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bwp_t__Wy__Bwp__inv->subWriteContents("Bwp_t__Wy__Bwp__inv",
                                              "mat_Bwp_t__Wy__Bwp__inv",
                                              "m",
                                              tmpSet);
    }

    if (m_env.displayVerbosity() >= 4) {
      double       Bwp_t__Wy__Bwp__InvLnDeterminant = m_Bwp_t__Wy__Bwp__inv->lnDeterminant();
      unsigned int Bwp_t__Wy__Bwp__InvRank          = m_Bwp_t__Wy__Bwp__inv->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bwp_t__Wy__Bwp__InvRank14        = m_Bwp_t__Wy__Bwp__inv->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", m_Bwp_t__Wy__Bwp__inv->numRowsLocal() = "  << m_Bwp_t__Wy__Bwp__inv->numRowsLocal()
                                << ", m_Bwp_t__Wy__Bwp__inv->numCols() = "       << m_Bwp_t__Wy__Bwp__inv->numCols()
                                << ": m_Bwp_t__Wy__Bwp__inv->lnDeterminant() = " << Bwp_t__Wy__Bwp__InvLnDeterminant
                                << ": m_Bwp_t__Wy__Bwp__inv->rank(0.,1.e-8) = "  << Bwp_t__Wy__Bwp__InvRank
                                << ": m_Bwp_t__Wy__Bwp__inv->rank(0.,1.e-14) = " << Bwp_t__Wy__Bwp__InvRank14
                                << std::endl;
      }
    }

    if (m_env.checkingLevel() >= 1) {
      // Check transpose operation
      D_M matShouldBeI1(*m_Bwp_t__Wy__Bwp      * *m_Bwp_t__Wy__Bwp__inv);
      D_M matShouldBeI2(*m_Bwp_t__Wy__Bwp__inv * *m_Bwp_t__Wy__Bwp     );
      D_M matI         (m_vu_space.zeroVector());
      for (unsigned int i = 0; i < matI.numRowsLocal(); ++i) {
        matI(i,i) = 1.;
      }
      matShouldBeI1 -= matI;
      double auxNorm1 = matShouldBeI1.normFrob();
      matShouldBeI2 -= matI;
      double auxNorm2 = matShouldBeI2.normFrob();
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", tests on m_Bt_Wy_B"
                                << ": matShouldBeI1.numRowsLocal() = " << matShouldBeI1.numRowsLocal()
                              //<< ", ||matShouldBeI1||_2^2 = "        << matShouldBeI1.normFrob() * matShouldBeI1.normFrob()
                              //<< ", ||matShouldBeI2||_2^2 = "        << matShouldBeI2.normFrob() * matShouldBeI2.normFrob()
                                << ", ||matI||_2^2 = "                 << matI.normFrob() * matI.normFrob()
                                << ", ||matShouldBeI1 - matI||_2^2 = " << auxNorm1 * auxNorm1
                                << ", ||matShouldBeI2 - matI||_2^2 = " << auxNorm2 * auxNorm2
                                << std::endl;
      }
    }

    //********************************************************************************
    // Form 'Bop^T' matrix
    //********************************************************************************
    D_M Bop_t(m_env,m_vu_space.map(),e.m_paper_n_y);
    Bop_t.fillWithTranspose(0,0,*m_Bmat_without_permut,true,true);

    //********************************************************************************
    // Compute 'Bop^T W_y Bop' matrix
    //********************************************************************************
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor():"
                              << "\n  m_Bop_t__Wy__Bop->numRowsLocal() = "     << m_Bop_t__Wy__Bop->numRowsLocal()
                              << "\n  m_Bop_t__Wy__Bop->numCols() = "          << m_Bop_t__Wy__Bop->numCols()
                              << "\n  Bop_t.numRowsLocal() = "                 << Bop_t.numRowsLocal()
                              << "\n  Bop_t.numCols() = "                      << Bop_t.numCols()
                              << "\n  m_Wy->numRowsLocal() = "                 << e.m_Wy->numRowsLocal()
                              << "\n  m_Wy->numCols() = "                      << e.m_Wy->numCols()
                              << "\n  m_Bmat_with_permut->numRowsLocal() = "   << m_Bmat_with_permut->numRowsLocal()
                              << "\n  m_Bmat_with_permut->numCols() = "        << m_Bmat_with_permut->numCols()
                              << std::endl;
    }

    *m_Bop_t__Wy__Bop = Bop_t * (*e.m_Wy * *m_Bmat_without_permut);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Bop_t__Wy__Bop before nugget'"
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bop_t__Wy__Bop->subWriteContents("Bop_t__Wy__Bop__beforePermut_beforeNugget",
                                         "mat_Bop_t__Wy__Bop__beforePermut_beforeNugget",
                                         "m",
                                         tmpSet);
    }

    // Just for checking
    if (m_env.displayVerbosity() >= 4) {
      double       Bop_t__Wy__Bop__beforeNugget__LnDeterminant = m_Bop_t__Wy__Bop->lnDeterminant();
      unsigned int Bop_t__Wy__Bop__beforeNugget__Rank          = m_Bop_t__Wy__Bop->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bop_t__Wy__Bop__beforeNugget__Rank14        = m_Bop_t__Wy__Bop->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", m_Bop_t__Wy__Bop__beforeNugget.numRowsLocal() = "  << m_Bop_t__Wy__Bop->numRowsLocal()
                                << ", m_Bop_t__Wy__Bop__beforeNugget.numCols() = "       << m_Bop_t__Wy__Bop->numCols()
                                << ": m_Bop_t__Wy__Bop__beforeNugget.lnDeterminant() = " << Bop_t__Wy__Bop__beforeNugget__LnDeterminant
                                << ", m_Bop_t__Wy__Bop__beforeNugget.rank(0.,1.e-8) = "  << Bop_t__Wy__Bop__beforeNugget__Rank
                                << ", m_Bop_t__Wy__Bop__beforeNugget.rank(0.,1.e-14) = " << Bop_t__Wy__Bop__beforeNugget__Rank14
                                << std::endl;
      }
    }

    //********************************************************************************
    // Apply permutation only now
    //********************************************************************************
    twoMats[0] = e.m_PD;
    twoMats[1] = &PK;
    D_M matP(m_vu_space.zeroVector());
    matP.fillWithBlocksDiagonally(0,0,twoMats,true,true);
    *m_Bop_t__Wy__Bop = matP * (*m_Bop_t__Wy__Bop * (matP.transpose()));

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Bop_t__Wy__Bop after permutation'"
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bop_t__Wy__Bop->subWriteContents("Bop_t__Wy__Bop__afterPermut_beforeNugget",
                                         "mat_Bop_t__Wy__Bop__afterPermut_beforeNugget",
                                         "m",
                                         tmpSet);
    }

    if (m_env.displayVerbosity() >= 4) {
      double       Bop_t__Wy__Bop__LnDeterminant = m_Bop_t__Wy__Bop->lnDeterminant();
      unsigned int Bop_t__Wy__Bop__Rank          = m_Bop_t__Wy__Bop->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bop_t__Wy__Bop__Rank14        = m_Bop_t__Wy__Bop->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Bop_t__Wy__Bop->numRowsLocal() = "  << m_Bop_t__Wy__Bop->numRowsLocal()
                                << ", m_Bop_t__Wy__Bop->numCols() = "       << m_Bop_t__Wy__Bop->numCols()
                                << ", m_Bop_t__Wy__Bop->lnDeterminant() = " << Bop_t__Wy__Bop__LnDeterminant
                                << ", m_Bop_t__Wy__Bop->rank(0.,1.e-8) = "  << Bop_t__Wy__Bop__Rank
                                << ", m_Bop_t__Wy__Bop->rank(0.,1.e-14) = " << Bop_t__Wy__Bop__Rank14
                                << std::endl;
      }
    }

    //********************************************************************************
    // Add nugget only now
    //********************************************************************************
    if (gcmOptionsObj.m_ov.m_nuggetValueForBtWyB != 0.) {
      for (unsigned int i = 0; i < m_Bop_t__Wy__Bop->numRowsLocal(); ++i) {
        (*m_Bop_t__Wy__Bop)(i,i) += gcmOptionsObj.m_ov.m_nuggetValueForBtWyB;
      }
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Bop_t__Wy__Bop after permutation and after nugget'"
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bop_t__Wy__Bop->subWriteContents("Bop_t__Wy__Bop__afterPermut_afterNugget",
                                         "mat_Bop_t__Wy__Bop__afterPermut_afterNugget",
                                         "m",
                                         tmpSet);
    }

    // Just for checking
    if (m_env.displayVerbosity() >= 4) {
      double       Bop_t__Wy__Bop__afterNugget__LnDeterminant = m_Bop_t__Wy__Bop->lnDeterminant();
      unsigned int Bop_t__Wy__Bop__afterNugget__Rank          = m_Bop_t__Wy__Bop->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bop_t__Wy__Bop__afterNugget__Rank14        = m_Bop_t__Wy__Bop->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", m_Bop_t__Wy__Bop__afterNugget.numRowsLocal() = "  << m_Bop_t__Wy__Bop->numRowsLocal()
                                << ", m_Bop_t__Wy__Bop__afterNugget.numCols() = "       << m_Bop_t__Wy__Bop->numCols()
                                << ": m_Bop_t__Wy__Bop__afterNugget.lnDeterminant() = " << Bop_t__Wy__Bop__afterNugget__LnDeterminant
                                << ", m_Bop_t__Wy__Bop__afterNugget.rank(0.,1.e-8) = "  << Bop_t__Wy__Bop__afterNugget__Rank
                                << ", m_Bop_t__Wy__Bop__afterNugget.rank(0.,1.e-14) = " << Bop_t__Wy__Bop__afterNugget__Rank14
                                << std::endl;
      }
    }

    //********************************************************************************
    // Compute 'Bop^T W_y Bop' inverse
    //********************************************************************************
    *m_Bop_t__Wy__Bop__inv = m_Bop_t__Wy__Bop->inverse(); // inversion savings

    if (m_env.displayVerbosity() >= 4) {
      double       Bop_t__Wy__Bop__inv__beforeNugget__LnDeterminant = m_Bop_t__Wy__Bop__inv->lnDeterminant();
      unsigned int Bop_t__Wy__Bop__beforeNugget__InvRank            = m_Bop_t__Wy__Bop__inv->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bop_t__Wy__Bop__beforeNugget__InvRank14          = m_Bop_t__Wy__Bop__inv->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", m_Bop_t__Wy__Bop__inv__beforeNugget.numRowsLocal() = "  << m_Bop_t__Wy__Bop__inv->numRowsLocal()
                                << ", m_Bop_t__Wy__Bop__inv__beforeNugget.numCols() = "       << m_Bop_t__Wy__Bop__inv->numCols()
                                << ": m_Bop_t__Wy__Bop__inv__beforeNugget.lnDeterminant() = " << Bop_t__Wy__Bop__inv__beforeNugget__LnDeterminant
                                << ", m_Bop_t__Wy__Bop__inv__beforeNugget.rank(0.,1.e-8) = "  << Bop_t__Wy__Bop__beforeNugget__InvRank
                                << ", m_Bop_t__Wy__Bop__inv__beforeNugget.rank(0.,1.e-14) = " << Bop_t__Wy__Bop__beforeNugget__InvRank14
                                << std::endl;
      }
    }

    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bop_t__Wy__Bop__inv->subWriteContents("Bop_t__Wy__Bop__inv__beforeNuggetForInv",
                                              "mat_Bop_t__Wy__Bop__inv__beforeNuggetForInv",
                                              "m",
                                              tmpSet);
    }

    // Add nugget
    if (gcmOptionsObj.m_ov.m_nuggetValueForBtWyBInv != 0.) {
      for (unsigned int i = 0; i < m_Bop_t__Wy__Bop__inv->numRowsLocal(); ++i) {
        (*m_Bop_t__Wy__Bop__inv)(i,i) += gcmOptionsObj.m_ov.m_nuggetValueForBtWyBInv;
      }
    }

    if (m_env.displayVerbosity() >= 4) {
      double       Bop_t__Wy__Bop__inv__afterNugget__LnDeterminant = m_Bop_t__Wy__Bop__inv->lnDeterminant();
      unsigned int Bop_t__Wy__Bop__afterNugget__InvRank            = m_Bop_t__Wy__Bop__inv->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int Bop_t__Wy__Bop__afterNugget__InvRank14          = m_Bop_t__Wy__Bop__inv->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", m_Bop_t__Wy__Bop__inv__afterNugget.numRowsLocal() = "  << m_Bop_t__Wy__Bop__inv->numRowsLocal()
                                << ", m_Bop_t__Wy__Bop__inv__afterNugget.numCols() = "       << m_Bop_t__Wy__Bop__inv->numCols()
                                << ": m_Bop_t__Wy__Bop__inv__afterNugget.lnDeterminant() = " << Bop_t__Wy__Bop__inv__afterNugget__LnDeterminant
                                << ", m_Bop_t__Wy__Bop__inv__afterNugget.rank(0.,1.e-8) = "  << Bop_t__Wy__Bop__afterNugget__InvRank
                                << ", m_Bop_t__Wy__Bop__inv__afterNugget.rank(0.,1.e-14) = " << Bop_t__Wy__Bop__afterNugget__InvRank14
                                << std::endl;
      }
    }

    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Bop_t__Wy__Bop__inv->subWriteContents("Bop_t__Wy__Bop__inv__afterNuggetForInv",
                                              "mat_Bop_t__Wy__Bop__inv__afterNuggetForInv",
                                              "m",
                                              tmpSet);
    }

    if (m_env.checkingLevel() >= 1) {
      // Check transpose operation
      D_M matShouldBeI1(*m_Bop_t__Wy__Bop      * *m_Bop_t__Wy__Bop__inv);
      D_M matShouldBeI2(*m_Bop_t__Wy__Bop__inv * *m_Bop_t__Wy__Bop     );
      D_M matI         (m_vu_space.zeroVector());
      for (unsigned int i = 0; i < matI.numRowsLocal(); ++i) {
        matI(i,i) = 1.;
      }
      matShouldBeI1 -= matI;
      double auxNorm1 = matShouldBeI1.normFrob();
      matShouldBeI2 -= matI;
      double auxNorm2 = matShouldBeI2.normFrob();
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", tests on m_Bt_Wy_B"
                                << ": matShouldBeI1.numRowsLocal() = " << matShouldBeI1.numRowsLocal()
                              //<< ", ||matShouldBeI1||_2^2 = "        << matShouldBeI1.normFrob() * matShouldBeI1.normFrob()
                              //<< ", ||matShouldBeI2||_2^2 = "        << matShouldBeI2.normFrob() * matShouldBeI2.normFrob()
                                << ", ||matI||_2^2 = "                 << matI.normFrob() * matI.normFrob()
                                << ", ||matShouldBeI1 - matI||_2^2 = " << auxNorm1 * auxNorm1
                                << ", ||matShouldBeI2 - matI||_2^2 = " << auxNorm2 * auxNorm2
                                << std::endl;
      }
    }

    //********************************************************************************
    // Compute exponent modifiers
    //********************************************************************************
    m_a_y_modifier = ((double) (e.m_paper_n_y - m_Bmat_rank)) / 2.;
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_a_y_modifier = "   << m_a_y_modifier
                              << std::endl;
    }

    D_V yVec_transformed(e.m_experimentStorage.yVec_transformed());
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Zvec_hat_vu.sizeLocal() = "             << m_Zvec_hat_vu.sizeLocal()
                              << ", m_Bop_t__Wy__Bop__inv->numRowsLocal() = " << m_Bop_t__Wy__Bop__inv->numRowsLocal()
                              << ", m_Bop_t__Wy__Bop__inv->numCols() = "      << m_Bop_t__Wy__Bop__inv->numCols()
                              << ", Bop_t.numRowsLocal() = "                  << Bop_t.numRowsLocal()
                              << ", Bop_t.numCols() = "                       << Bop_t.numCols()
                              << ", m_Wy->numRowsLocal() = "                  << e.m_Wy->numRowsLocal()
                              << ", m_Wy->numCols() = "                       << e.m_Wy->numCols()
                              << ", yVec_transformed.sizeLocal() = "          << yVec_transformed.sizeLocal()
                              << std::endl;
    }
    m_Zvec_hat_vu = *m_Bop_t__Wy__Bop__inv * (Bwp_t * (*e.m_Wy * yVec_transformed)); // aqui
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Zvec_hat_vu.subWriteContents("Zvec_hat_vu",
                                     "vec_Zvec_hat_vu",
                                     "m",
                                     tmpSet);
    }
    D_V tmpVec2(yVec_transformed - (*m_Bmat_with_permut * m_Zvec_hat_vu));
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ", yVec_transformed = "  << yVec_transformed.sizeLocal()
                              << ", m_Zvec_hat_vu = "     << m_Zvec_hat_vu
                              << ", tmpVec2 (initial) = " << tmpVec2
                              << std::endl;
    }

    tmpVec2 = *e.m_Wy * tmpVec2;
    m_b_y_modifier = scalarProduct(yVec_transformed,tmpVec2) / 2.;
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": tmpVec2 (final) = " << tmpVec2
                              << ", m_b_y_modifier = "  << m_b_y_modifier
                              << std::endl;
    }
  } // if (allOutputsAreScalar)

  //********************************************************************************
  // Instantiate Smat spaces
  //********************************************************************************
  unsigned int sumDims = 0;
  for (unsigned int i = 0; i < m_Smat_u_is.size(); ++i) {
    m_Rmat_u_is[i] = new D_M(e.m_paper_n_space.zeroVector()); // to be deleted on destructor
    m_Smat_u_is[i] = new D_M(e.m_paper_n_space.zeroVector()); // to be deleted on destructor
    sumDims += e.m_paper_n;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating the m_Smat_u_i spaces"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumDims, m_u_size, "'sumDims' and 'm_u_size' should be equal");

  unsigned int sumNumRows = 0;
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < m_Smat_uw_is.size(); ++i) {
    m_Rmat_uw_is[i] = new D_M(m_env,e.m_paper_n_space.map(),s.m_paper_m); // Yes, 'u' only // to be deleted on destructor
    m_Smat_uw_is[i] = new D_M(m_env,e.m_paper_n_space.map(),s.m_paper_m); // Yes, 'u' only // to be deleted on destructor
    sumNumRows += e.m_paper_n_space.dimLocal();
    sumNumCols += s.m_paper_m;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating the m_Smat_uw_i matrices"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumNumRows, m_u_size, "'sumNumRows' and 'm_u_size' should be equal");
  queso_require_equal_to_msg(sumNumCols, s.m_w_size, "'sumNumCols' and 'm_w_size' should be equal");

  //********************************************************************************
  // Instantiate 'u_hat_u_asterisk' matrices
  //********************************************************************************
  sumNumRows = 0;
  sumNumCols = 0;
  for (unsigned int i = 0; i < m_Smat_u_hat_u_asterisk_is.size(); ++i) {
    m_Rmat_u_hat_u_asterisk_is[i] = new D_M(m_env, e.m_paper_n_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    m_Smat_u_hat_u_asterisk_is[i] = new D_M(m_env, e.m_paper_n_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    sumNumRows += e.m_paper_n_space.dimLocal();
    sumNumCols += 1;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating the m_Smat_u_hat_u_asterisk_i matrices"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumNumRows, m_u_size, "'sumNumRows' and 'm_u_size' should be equal");
  queso_require_equal_to_msg(sumNumCols, s.m_paper_p_eta, "'sumNumCols (1)' and 'm_paper_p_eta' should be equal");

  //********************************************************************************
  // Instantiate 'w_hat_u_asterisk' matrices
  //********************************************************************************
  sumNumRows = 0;
  sumNumCols = 0;
  for (unsigned int i = 0; i < m_Smat_w_hat_u_asterisk_is.size(); ++i) {
    m_Rmat_w_hat_u_asterisk_is[i] = new D_M(m_env, s.m_paper_m_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    m_Smat_w_hat_u_asterisk_is[i] = new D_M(m_env, s.m_paper_m_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    sumNumRows += s.m_paper_m_space.dimLocal();
    sumNumCols += 1;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating the m_Smat_w_hat_u_asterisk_i matrices"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumNumRows, s.m_w_size, "'sumNumRows' and 'm_w_size' should be equal");
  queso_require_equal_to_msg(sumNumCols, s.m_paper_p_eta, "'sumNumCols (2)' and 'm_paper_p_eta' should be equal");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~GcmJointInfo()
{
  for (unsigned int i = 0; i < m_Smat_w_hat_u_asterisk_is.size(); ++i) {
    delete m_Smat_w_hat_u_asterisk_is[i]; // to be deleted on destructor
    m_Smat_w_hat_u_asterisk_is[i] = NULL;
    delete m_Rmat_w_hat_u_asterisk_is[i]; // to be deleted on destructor
    m_Rmat_w_hat_u_asterisk_is[i] = NULL;
  }

  for (unsigned int i = 0; i < m_Smat_u_hat_u_asterisk_is.size(); ++i) {
    delete m_Smat_u_hat_u_asterisk_is[i]; // to be deleted on destructor
    m_Smat_u_hat_u_asterisk_is[i] = NULL;
    delete m_Rmat_u_hat_u_asterisk_is[i]; // to be deleted on destructor
    m_Rmat_u_hat_u_asterisk_is[i] = NULL;
  }

  for (unsigned int i = 0; i < m_Smat_uw_is.size(); ++i) {
    delete m_Smat_uw_is[i]; // to be deleted on destructor
    m_Smat_uw_is[i] = NULL;
    delete m_Rmat_uw_is[i]; // to be deleted on destructor
    m_Rmat_uw_is[i] = NULL;
  }

  for (unsigned int i = 0; i < m_Smat_u_is.size(); ++i) {
    delete m_Smat_u_is[i]; // to be deleted on destructor
    m_Smat_u_is[i] = NULL;
    delete m_Rmat_u_is[i]; // to be deleted on destructor
    m_Rmat_u_is[i] = NULL;
  }

  delete m_Bop_t__Wy__Bop__inv;
  delete m_Bwp_t__Wy__Bwp__inv;
  delete m_Bop_t__Wy__Bop;
  delete m_Bwp_t__Wy__Bwp;
  delete m_Bmat_without_permut;
  delete m_Bmat_with_permut;
}

}  // End namespace QUESO

template class QUESO::GcmJointInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
