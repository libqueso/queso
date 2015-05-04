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

#include <queso/GcmSimulationInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::GcmSimulationInfo(
  const GpmsaComputerModelOptions&                  gcmOptionsObj,
  bool                                                     allOutputsAreScalar,
  const SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationStorage,
  const SimulationModel  <S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationModel)
  :
  m_env                                            (simulationStorage.env()),
  m_simulationStorage                              (simulationStorage),
  m_simulationModel                                (simulationModel),
  m_paper_p_x                                      (simulationStorage.scenarioSpace().dimLocal()),
  m_paper_xs_asterisks_standard                    (simulationModel.xs_asterisks_standard()),
  m_paper_ts_asterisks_standard                    (simulationModel.ts_asterisks_standard()),
  m_paper_p_t                                      (simulationStorage.parameterSpace().dimLocal()),
  m_paper_m                                        (simulationStorage.numSimulations()),
  m_paper_n_eta                                    (simulationStorage.outputSpace().dimLocal()),
  m_paper_p_eta                                    (simulationModel.numBasis()),
  m_paper_m_space                                  (m_env, "paper_m_", m_paper_m, NULL),
  m_1lambdaEtaDim                                  (1), // '1' in paper
  m_1lambdaEtaSpace                                (m_env, "1lambdaEta_", m_1lambdaEtaDim, NULL),
  m_1lambdaEtaMins                                 (m_env,m_1lambdaEtaSpace.map(),0.),
  m_1lambdaEtaMaxs                                 (m_env,m_1lambdaEtaSpace.map(),+INFINITY),
  m_1lambdaEtaDomain                               ("1lambdaEta_",m_1lambdaEtaSpace,m_1lambdaEtaMins,m_1lambdaEtaMaxs),
  m_1lambdaEtaGammaAVec                            (m_env,m_1lambdaEtaSpace.map(),simulationModel.optionsObj().m_ov.m_a_eta),
  m_1lambdaEtaGammaBVec                            (m_env,m_1lambdaEtaSpace.map(),1./simulationModel.optionsObj().m_ov.m_b_eta), // Yes, 1./...
  m_1lambdaEtaPriorRv                              ("1lambdaEta_",m_1lambdaEtaDomain,m_1lambdaEtaGammaAVec,m_1lambdaEtaGammaBVec),
  m_like_previous1                                 (m_1lambdaEtaSpace.zeroVector()),
  m_tmp_1lambdaEtaVec                              (m_1lambdaEtaSpace.zeroVector()),
  m_2lambdaWDim                                    (m_paper_p_eta), // 'p_eta' in paper
  m_2lambdaWSpace                                  (m_env, "2lambdaW_", m_2lambdaWDim, NULL),
  m_2lambdaWMins                                   (m_env,m_2lambdaWSpace.map(),0.),
  m_2lambdaWMaxs                                   (m_env,m_2lambdaWSpace.map(),+INFINITY),
  m_2lambdaWDomain                                 ("2lambdaW_",m_2lambdaWSpace,m_2lambdaWMins,m_2lambdaWMaxs),
  m_2lambdaWGammaAVec                              (m_env,m_2lambdaWSpace.map(),simulationModel.optionsObj().m_ov.m_a_w),
  m_2lambdaWGammaBVec                              (m_env,m_2lambdaWSpace.map(),1./simulationModel.optionsObj().m_ov.m_b_w), // Yes, 1./...
  m_2lambdaWPriorRv                                ("2lambdaW_",m_2lambdaWDomain,m_2lambdaWGammaAVec,m_2lambdaWGammaBVec),
  m_like_previous2                                 (m_2lambdaWSpace.zeroVector()),
  m_tmp_2lambdaWVec                                (m_2lambdaWSpace.zeroVector()),
  m_3rhoWDim                                       (m_paper_p_eta * (m_paper_p_x + m_paper_p_t)), // 'p_eta * (p_x + p_t)' in paper
  m_3rhoWSpace                                     (m_env, "3rhoW_", m_3rhoWDim, NULL),
  m_3rhoWMins                                      (m_env,m_3rhoWSpace.map(),0.),
  m_3rhoWMaxs                                      (m_env,m_3rhoWSpace.map(),1.),
  m_3rhoWDomain                                    ("3rhoW_",m_3rhoWSpace,m_3rhoWMins,m_3rhoWMaxs),
  m_3rhoWBetaAVec                                  (m_env,m_3rhoWSpace.map(),simulationModel.optionsObj().m_ov.m_a_rho_w),
  m_3rhoWBetaBVec                                  (m_env,m_3rhoWSpace.map(),simulationModel.optionsObj().m_ov.m_b_rho_w),
  m_3rhoWPriorRv                                   ("3rhoW_",m_3rhoWDomain,m_3rhoWBetaAVec,m_3rhoWBetaBVec),
  m_like_previous3                                 (m_3rhoWSpace.zeroVector()),
  m_tmp_3rhoWVec                                   (m_3rhoWSpace.zeroVector()),
  m_4lambdaSDim                                    (m_paper_p_eta), // 'p_eta' in paper
  m_4lambdaSSpace                                  (m_env, "4lambdaS_", m_4lambdaSDim, NULL),
  m_4lambdaSMins                                   (m_env,m_4lambdaSSpace.map(),0.),
  m_4lambdaSMaxs                                   (m_env,m_4lambdaSSpace.map(),+INFINITY),
  m_4lambdaSDomain                                 ("4lambdaS_",m_4lambdaSSpace,m_4lambdaSMins,m_4lambdaSMaxs),
  m_4lambdaSGammaAVec                              (m_env,m_4lambdaSSpace.map(),simulationModel.optionsObj().m_ov.m_a_s),
  m_4lambdaSGammaBVec                              (m_env,m_4lambdaSSpace.map(),1./simulationModel.optionsObj().m_ov.m_b_s), // Yes, 1./...
  m_4lambdaSPriorRv                                ("4lambdaS_",m_4lambdaSDomain,m_4lambdaSGammaAVec,m_4lambdaSGammaBVec),
  m_like_previous4                                 (m_4lambdaSSpace.zeroVector()),
  m_tmp_4lambdaSVec                                (m_4lambdaSSpace.zeroVector()),
  m_eta_size                                       (m_paper_m*m_paper_n_eta),
  m_eta_space                                      (m_env, "eta_", m_eta_size, NULL),
  m_w_size                                         (m_paper_m*m_paper_p_eta),
  m_w_space                                        (m_env, "w_", m_w_size, NULL),
  m_unique_w_space                                 (m_env, "unique_w_",  m_paper_p_eta, NULL),
  m_Zvec_hat_w                                     (m_w_space.zeroVector()),
  m_rho_w_space                                    (m_env, "rho_w_", m_paper_p_x+m_paper_p_t, NULL),
  m_tmp_rho_w_vec                                  (m_rho_w_space.zeroVector()),
  m_Rmat_w_is                                      (m_paper_p_eta, (Q_M*) NULL), // to be deleted on destructor
  m_Smat_w_is                                      (m_paper_p_eta, (Q_M*) NULL), // to be deleted on destructor
  m_Smat_w                                         (m_w_space.zeroVector()),
  m_Smat_w_hat                                     (m_w_space.zeroVector()),
  m_Rmat_w_hat_w_asterisk_is                       (m_paper_p_eta, (Q_M*) NULL), // to be deleted on destructor
  m_Smat_w_hat_w_asterisk_is                       (m_paper_p_eta, (Q_M*) NULL), // to be deleted on destructor
  m_Smat_w_hat_w_asterisk                          (m_env,m_w_space.map(), m_paper_p_eta),
  m_Smat_w_hat_w_asterisk_t                        (m_env,m_unique_w_space.map(),m_w_size),
  m_Smat_w_asterisk_w_asterisk                     (m_unique_w_space.zeroVector()),
  m_Kmat                                           (simulationModel.Kmat()),
  m_Kmat_eta                                       (simulationModel.Kmat_eta()),
  m_Kmat_rank                                      (std::min(m_Kmat.numRowsGlobal(),m_Kmat.numCols())), // Might be smaller
  m_Kt_K                                           (NULL), // to be deleted on destructor
  m_Kt_K_inv                                       (NULL), // to be deleted on destructor
  m_a_eta_modifier                                 (0.),
  m_b_eta_modifier                                 (0.),
  m_predW_counter                                  (MiscUintDebugMessage(0,NULL)),
  m_predW_summingRVs_unique_w_meanVec              (m_unique_w_space.zeroVector()),
  m_predW_summingRVs_mean_of_unique_w_covMatrices  (m_unique_w_space.zeroVector()),
  m_predW_summingRVs_covMatrix_of_unique_w_means   (m_unique_w_space.zeroVector()),
  m_predW_summingRVs_corrMatrix_of_unique_w_means  (m_unique_w_space.zeroVector())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": key-debug"
                            << ", some entities just created (not yet populated)"
                            << ", m_Zvec_hat_w.sizeLocal() = " << m_Zvec_hat_w.sizeLocal()
                            << ", m_Smat_w.numRowsLocal() = "  << m_Smat_w.numRowsLocal()
                            << ", m_Smat_w.numCols() = "       << m_Smat_w.numCols()
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

  //********************************************************************************
  // Print information
  //********************************************************************************
  unsigned int m_Kmat_rank = m_Kmat.rank(0.,1.e-8 ); // todo: should be an option
  unsigned int kRank14     = m_Kmat.rank(0.,1.e-14);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": m_Kmat.numRowsLocal() = "  << m_Kmat.numRowsLocal()
                            << ", m_Kmat.numCols() = "       << m_Kmat.numCols()
                            << ", m_Kmat.rank(0.,1.e-8) = "  << m_Kmat_rank
                            << ", m_Kmat.rank(0.,1.e-14) = " << kRank14
                            << std::endl;
  }
  if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
    m_Kmat.subWriteContents("K",
                            "mat_K",
                            "m",
                            tmpSet);
  }

  Q_V etaVec_transformed(simulationModel.etaVec_transformed("Gp.h.004"));
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": m_Zvec_hat_w.sizeLocal() = "       << m_Zvec_hat_w.sizeLocal()
                            << ", etaVec_transformed.sizeLocal() = " << etaVec_transformed.sizeLocal()
                            << std::endl;
  }

  if (allOutputsAreScalar) {
    m_Zvec_hat_w = etaVec_transformed;
  }
  else {
    //********************************************************************************
    // Form 'K^T' matrix
    //********************************************************************************
    Q_M Kt(m_env,m_w_space.map(),m_eta_size);
    Kt.fillWithTranspose(0,0,m_Kmat,true,true);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": Kt.numRowsLocal() = "     << Kt.numRowsLocal()
                              << ", Kt.numCols() = "          << Kt.numCols()
                              << std::endl;
    }

    //********************************************************************************
    // Compute 'K^T K' matrix, and its inverse
    //********************************************************************************
    m_Kt_K     = new Q_M(m_w_space.zeroVector()); // to be deleted on destructor
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", m_Kt_K just created (not yet populated)"
                              << ", numRowsLocal() = " << m_Kt_K->numRowsLocal()
                              << ", numCols() = "      << m_Kt_K->numCols()
                              << std::endl;
    }
    m_Kt_K_inv = new Q_M(m_w_space.zeroVector()); // to be deleted on destructor
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", m_Kt_K_inv just created (not yet populated)"
                              << ", numRowsLocal() = " << m_Kt_K_inv->numRowsLocal()
                              << ", numCols() = "      << m_Kt_K_inv->numCols()
                              << std::endl;
    }

    *m_Kt_K = Kt * m_Kmat;
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      m_Kt_K->setPrintHorizontally(false);
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Kt_K'"
                            //<< "\n m_Kt_K =\n" << m_Kt_K
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Kt_K->subWriteContents("Kt_K",
                               "mat_Kt_K",
                               "m",
                               tmpSet);
    }

    if (m_env.displayVerbosity() >= 4) {
      double       ktKLnDeterminant = m_Kt_K->lnDeterminant();
      unsigned int ktKRank          = m_Kt_K->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int ktKRank14        = m_Kt_K->rank(0.,1.e-14);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Kt_K->numRowsLocal() = "  << m_Kt_K->numRowsLocal()
                                << ", m_Kt_K->numCols() = "       << m_Kt_K->numCols()
                                << ", m_Kt_K->lnDeterminant() = " << ktKLnDeterminant
                                << ", m_Kt_K->rank(0.,1.e-8) = "  << ktKRank
                                << ", m_Kt_K->rank(0.,1.e-14) = " << ktKRank14
                                << std::endl;
      }
    }

    *m_Kt_K_inv = m_Kt_K->inverse(); // inversion savings

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      m_Kt_K_inv->setPrintHorizontally(false);
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished computing 'm_Kt_K_inv'"
                            //<< "\n m_Kt_K_inv =\n" << m_Kt_K_inv
                              << std::endl;
    }
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Kt_K_inv->subWriteContents("Kt_K_inv",
                                   "mat_Kt_K_inv",
                                   "m",
                                   tmpSet);
    }

    if (m_env.displayVerbosity() >= 4) {
      double       ktKInvLnDeterminant = m_Kt_K_inv->lnDeterminant();
      unsigned int ktKInvRank          = m_Kt_K_inv->rank(0.,1.e-8 ); // todo: should be an option
      unsigned int ktKInvRank14        = m_Kt_K_inv->rank(0.,1.e-14);
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": m_Kt_K_inv->numRowsLocal() = "  << m_Kt_K_inv->numRowsLocal()
                                << ", m_Kt_K_inv->numCols() = "       << m_Kt_K_inv->numCols()
                                << ", m_Kt_K_inv->lnDeterminant() = " << ktKInvLnDeterminant
                                << ", m_Kt_K_inv->rank(0.,1.e-8) = "  << ktKInvRank
                                << ", m_Kt_K_inv->rank(0.,1.e-14) = " << ktKInvRank14
                                << std::endl;
      }
    }

    //********************************************************************************
    // Compute exponent modifiers
    //********************************************************************************
    m_a_eta_modifier = ((double) m_paper_m) * ((double) (m_paper_n_eta - m_paper_p_eta)) / 2.;
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_a_eta_modifier = " << m_a_eta_modifier
                              << std::endl;
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_Kt_K_inv->numRowsLocal() = "     << m_Kt_K_inv->numRowsLocal()
                              << ", m_Kt_K_inv->numCols() = "          << m_Kt_K_inv->numCols()
                              << ", Kt.numRowsLocal() = "              << Kt.numRowsLocal()
                              << ", Kt.numCols() = "                   << Kt.numCols()
                              << std::endl;
    }

    m_Zvec_hat_w = (*m_Kt_K_inv) * (Kt * etaVec_transformed);
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Zvec_hat_w.subWriteContents("Zvec_hat_w",
                                    "vec_Zvec_hat_w",
                                    "m",
                                    tmpSet);
    }
    Q_V tmpVec1(etaVec_transformed - (m_Kmat * m_Zvec_hat_w));
    m_b_eta_modifier = scalarProduct(etaVec_transformed,tmpVec1) / 2.;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_b_eta_modifier = " << m_b_eta_modifier
                              << std::endl;
    }
  } // if (allOutputsAreScalar)

  //********************************************************************************
  // Instantiate Smat spaces
  //********************************************************************************
  unsigned int sumDims = 0;
  for (unsigned int i = 0; i < m_Smat_w_is.size(); ++i) {
    m_Rmat_w_is[i] = new Q_M(m_paper_m_space.zeroVector()); // to be deleted on destructor
    m_Smat_w_is[i] = new Q_M(m_paper_m_space.zeroVector()); // to be deleted on destructor
    sumDims += m_paper_m_space.dimLocal();
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating m_Smat_w_is"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumDims, m_w_size, "'sumDims' and 'm_w_size' should be equal");

  unsigned int sumNumRows = 0;
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < m_Smat_w_hat_w_asterisk_is.size(); ++i) {
    m_Rmat_w_hat_w_asterisk_is[i] = new Q_M(m_env, m_paper_m_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    m_Smat_w_hat_w_asterisk_is[i] = new Q_M(m_env, m_paper_m_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    sumNumRows += m_paper_m_space.dimLocal();
    sumNumCols += 1; // m_paper_p_eta;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating the m_Smat_w_hat_w_asterisk_i matrices"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumNumRows, m_w_size, "'sumNumRows' and 'm_w_size' should be equal");
  queso_require_equal_to_msg(sumNumCols, (m_paper_p_eta), "'sumNumCols' and 'm_paper_p_eta*m_paper_p_eta' should be equal");

  //********************************************************************************
  // Display information
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "KEY In GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << "\n KEY  m_paper_p_x = "            << m_paper_p_x
                            << "\n KEY  m_paper_p_t = "            << m_paper_p_t
                            << "\n KEY  m_paper_m = "              << m_paper_m
                            << "\n KEY  m_paper_n_eta = "          << m_paper_n_eta
                            << "\n KEY  m_paper_p_eta = "          << m_paper_p_eta
                            << "\n KEY  m_1lambdaEtaDim = "        << m_1lambdaEtaDim
                            << ", m_1lambdaEtaGammaAVec = "        << m_1lambdaEtaGammaAVec
                            << ", m_1lambdaEtaGammaBVec = "        << m_1lambdaEtaGammaBVec
                            << "\n KEY  m_2lambdaWDim   = "        << m_2lambdaWDim
                            << ", m_2lambdaWGammaAVec = "          << m_2lambdaWGammaAVec
                            << ", m_2lambdaWGammaBVec = "          << m_2lambdaWGammaBVec
                            << "\n KEY  m_3rhoWDim      = "        << m_3rhoWDim
                            << ", m_3rhoWBetaAVec = "              << m_3rhoWBetaAVec
                            << ", m_3rhoWBetaBVec = "              << m_3rhoWBetaBVec
                            << "\n KEY  m_4lambdaSDim   = "        << m_4lambdaSDim
                            << ", m_4lambdaSGammaAVec = "          << m_4lambdaSGammaAVec
                            << ", m_4lambdaSGammaBVec = "          << m_4lambdaSGammaBVec
                            << "\n KEY  full 'eta' vector size = " << m_paper_m * m_paper_n_eta // = simulationModel.etaVec_transformed("Gp.h.001").sizeLocal()
                            << std::endl;
  }

  //********************************************************************************
  // Make checks
  //********************************************************************************
  queso_require_equal_to_msg(simulationModel.etaVec_transformed("Gp.h.002").sizeLocal(), m_eta_size, "incompatible calculations for 'eta' vector size");

  queso_require_equal_to_msg(m_paper_p_x, simulationStorage.scenarioSpace().dimLocal(), "'m_paper_p_x' and 'simulationStorage.scenarioSpace().dimLocal()' should be equal");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }
}

template <class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>::~GcmSimulationInfo()
{
  for (unsigned int i = 0; i < m_Smat_w_hat_w_asterisk_is.size(); ++i) {
    delete m_Smat_w_hat_w_asterisk_is[i]; // to be deleted on destructor
    m_Smat_w_hat_w_asterisk_is[i] = NULL;
    delete m_Rmat_w_hat_w_asterisk_is[i]; // to be deleted on destructor
    m_Rmat_w_hat_w_asterisk_is[i] = NULL;
  }

  for (unsigned int i = 0; i < m_Smat_w_is.size(); ++i) {
    delete m_Smat_w_is[i]; // to be deleted on destructor
    m_Smat_w_is[i] = NULL;
    delete m_Rmat_w_is[i]; // to be deleted on destructor
    m_Rmat_w_is[i] = NULL;
  }

  delete m_Kt_K_inv; // to be deleted on destructor
  delete m_Kt_K;     // to be deleted on destructor
}

}  // End namespace QUESO

template class QUESO::GcmSimulationInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
