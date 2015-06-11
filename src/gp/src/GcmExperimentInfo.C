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

#include <queso/GcmExperimentInfo.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M>
GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::GcmExperimentInfo(
  const GpmsaComputerModelOptions&          gcmOptionsObj,
        bool                                       allOutputsAreScalar,
  const ExperimentStorage<S_V,S_M,D_V,D_M>& experimentStorage,
  const ExperimentModel  <S_V,S_M,D_V,D_M>& experimentModel,
  const BaseVectorRV     <P_V,P_M>&         thetaPriorRv)
  :
  m_env                       (experimentStorage.env()),
  m_experimentStorage         (experimentStorage),
  m_experimentModel           (experimentModel),
  m_paper_p_x                 (experimentStorage.scenarioSpace().dimLocal()),
  m_paper_n                   (experimentStorage.numExperiments()),
  m_paper_xs_standard         (experimentStorage.xs_standard()),
  m_paper_n_ys_transformed    (experimentStorage.n_ys_transformed()),
  m_paper_n_y                 (experimentStorage.n_y()),
  m_paper_p_delta             (experimentModel.numBasis()),
  m_paper_F                   (experimentModel.numBasisGroups()),
  m_paper_Gs                  (experimentModel.Gs()),
  m_paper_n_space             (m_env, "paper_n_", m_paper_n, NULL),
  m_5lambdaYDim               (MiscUintDebugMessage(1,NULL)), // "GcmExperimentInfo::consctructor(), m_5lambdaYDim")), // '1' in paper
  m_5lambdaYSpace             (m_env, "5lambdaY_", m_5lambdaYDim, NULL),
  m_5lambdaYMins              (m_env,m_5lambdaYSpace.map(),0.),
  m_5lambdaYMaxs              (m_env,m_5lambdaYSpace.map(),+INFINITY),
  m_5lambdaYDomain            ("5lambdaY_",m_5lambdaYSpace,m_5lambdaYMins,m_5lambdaYMaxs),
  m_5lambdaYGammaAVec         (m_env,m_5lambdaYSpace.map(),experimentModel.optionsObj().m_ov.m_a_y),
  m_5lambdaYGammaBVec         (m_env,m_5lambdaYSpace.map(),1./experimentModel.optionsObj().m_ov.m_b_y), // Yes, 1./...
  m_5lambdaYPriorRv           ("5lambdaY_",m_5lambdaYDomain,m_5lambdaYGammaAVec,m_5lambdaYGammaBVec),
  m_like_previous5            (m_5lambdaYSpace.zeroVector()),
  m_tmp_5lambdaYVec           (m_5lambdaYSpace.zeroVector()),
  m_6lambdaVDim               (MiscUintDebugMessage(m_paper_F,NULL)), // "GcmExperimentInfo::consctructor(), m_6lambdaVDim")), // 'F' in paper
  m_6lambdaVSpace             (m_env, "6lambdaV_", m_6lambdaVDim, NULL),
  m_6lambdaVMins              (m_env,m_6lambdaVSpace.map(),0.),
  m_6lambdaVMaxs              (m_env,m_6lambdaVSpace.map(),+INFINITY),
  m_6lambdaVDomain            ("6lambdaV_",m_6lambdaVSpace,m_6lambdaVMins,m_6lambdaVMaxs),
  m_6lambdaVGammaAVec         (m_env,m_6lambdaVSpace.map(),experimentModel.optionsObj().m_ov.m_a_v),
  m_6lambdaVGammaBVec         (m_env,m_6lambdaVSpace.map(),1./experimentModel.optionsObj().m_ov.m_b_v), // Yes, 1./...
  m_6lambdaVPriorRv           ("6lambdaV_",m_6lambdaVDomain,m_6lambdaVGammaAVec,m_6lambdaVGammaBVec),
  m_like_previous6            (m_6lambdaVSpace.zeroVector()),
  m_tmp_6lambdaVVec           (m_6lambdaVSpace.zeroVector()),
  m_7rhoVDim                  (MiscUintDebugMessage(m_paper_F * m_paper_p_x,NULL)), // "GcmExperimentInfo::consctructor(), m_7lrhoVDim")), // 'F * p_x' in paper
  m_7rhoVSpace                (m_env, "7rhoV_", m_7rhoVDim, NULL),
  m_7rhoVMins                 (m_env,m_7rhoVSpace.map(),0.),
  m_7rhoVMaxs                 (m_env,m_7rhoVSpace.map(),1.),
  m_7rhoVDomain               ("7rhoV_",m_7rhoVSpace,m_7rhoVMins,m_7rhoVMaxs),
  m_7rhoVBetaAVec             (m_env,m_7rhoVSpace.map(),experimentModel.optionsObj().m_ov.m_a_rho_v),
  m_7rhoVBetaBVec             (m_env,m_7rhoVSpace.map(),experimentModel.optionsObj().m_ov.m_b_rho_v),
  m_7rhoVPriorRv              ("7rhoV_",m_7rhoVDomain,m_7rhoVBetaAVec,m_7rhoVBetaBVec),
  m_like_previous7            (m_7rhoVSpace.zeroVector()),
  m_tmp_7rhoVVec              (m_7rhoVSpace.zeroVector()),
  m_8thetaDim                 (MiscUintDebugMessage(thetaPriorRv.imageSet().vectorSpace().dimLocal(),NULL)), // "GcmExperimentInfo::consctructor(), m_8thetaDim")),
  m_8thetaSpace               (m_env, "8theta_", m_8thetaDim, NULL),
  m_8thetaPriorRv             (thetaPriorRv),
  m_like_previous8            (m_8thetaSpace.zeroVector()),
  m_tmp_8thetaVec             (m_8thetaSpace.zeroVector()),
  m_v_size                    (m_paper_n*m_paper_p_delta),
  m_v_space                   (m_env, "v_", m_v_size, NULL),
  m_unique_v_space            (m_env, "unique_v_",  m_paper_p_delta, NULL),
  m_rho_v_space               (m_env, "rho_v_", m_paper_p_x, NULL),
  m_y_space                   (m_env, "n_y_", m_paper_n_y, NULL),
  m_tmp_rho_v_vec             (m_rho_v_space.zeroVector()),
  m_Imat_v_i_spaces           (m_paper_F, (VectorSpace<D_V,D_M>*) NULL), // to be deleted on destructor
  m_Imat_v_is                 (m_paper_F, (D_M*) NULL),                         // to be deleted on destructor
  m_Rmat_v_i_spaces           (m_paper_F, (VectorSpace<D_V,D_M>*) NULL), // to be deleted on destructor
  m_Rmat_v_is                 (m_paper_F, (D_M*) NULL),                         // to be deleted on destructor
  m_Smat_v_i_spaces           (m_paper_F, (VectorSpace<D_V,D_M>*) NULL), // to be deleted on destructor
  m_Smat_v_is                 (m_paper_F, (D_M*) NULL),                         // to be deleted on destructor
  m_Smat_v                    (m_v_space.zeroVector()),
  m_Rmat_v_hat_v_asterisk_is  (m_paper_p_delta, (D_M*) NULL),                   // to be deleted on destructor
  m_Smat_v_hat_v_asterisk_is  (m_paper_p_delta, (D_M*) NULL),                   // to be deleted on destructor
  m_Smat_v_hat_v_asterisk     (m_env, m_v_space.map(),        m_paper_p_delta),
  m_Smat_v_hat_v_asterisk_t   (m_env, m_unique_v_space.map(), m_v_size),
  m_PD                        (NULL), // to be deleted on destructor
  m_Dmat_BlockDiag            (NULL),
  m_Dmat_BlockDiag_permut     (NULL), // to be deleted on destructor
  m_Wy        (NULL),
  m_Smat_v_asterisk_v_asterisk(m_unique_v_space.zeroVector())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << ": key-debug"
                            << ", some entities just created (not yet populated)"
                            << ", m_Smat_v.numRowsLocal() = " << m_Smat_v.numRowsLocal()
                            << ", m_Smat_v.numCols() = "      << m_Smat_v.numCols()
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

  if (allOutputsAreScalar) {
    // Do nothing
  }
  else {
    m_PD                    = new D_M(m_v_space.zeroVector()); // to be deleted on destructor
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": key-debug"
                              << ", m_PD just created (not yet populated)"
                              << ", numRowsLocal() = " << m_PD->numRowsLocal()
                              << ", numCols() = "      << m_PD->numCols()
                              << std::endl;
    }
    m_Dmat_BlockDiag        = &experimentModel.Dmat_BlockDiag();
    m_Dmat_BlockDiag_permut = new D_M(m_env,m_y_space.map(),m_v_size); // to be deleted on destructor
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": key-debug"
                              << ", m_Dmat_BlockDiag just created (not yet populated)"
                              << ", numRowsLocal() = " << m_Dmat_BlockDiag->numRowsLocal()
                              << ", numCols() = "      << m_Dmat_BlockDiag->numCols()
                              << std::endl;
    }
    m_Wy    = &experimentStorage.Wy();

    //********************************************************************************
    // Print information
    //********************************************************************************
    if (gcmOptionsObj.m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_ov.m_dataOutputAllowedSet.end()) {
      m_Wy->subWriteContents("Wy",
                             "mat_Wy",
                             "m",
                             tmpSet);
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      D_M Wy_filtered(*m_Wy);
      Wy_filtered.setPrintHorizontally(false);
      Wy_filtered.filterSmallValues(1.e-6);
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": Wy_filtered.numRowsLocal() = " << Wy_filtered.numRowsLocal()
                              << ", Wy_filtered.numCols() = "      << Wy_filtered.numCols()
                              << ", Wy_filtered contents =\n"      << Wy_filtered
                              << std::endl;
    }

    //********************************************************************************
    // Form 'P_D' matrix
    //********************************************************************************
    for (unsigned int i = 0; i < m_paper_p_delta; ++i) {
      for (unsigned int j = 0; j < m_paper_n; ++j) {
        unsigned int row = j + (m_paper_n*i);
        unsigned int col = (j*m_paper_p_delta)+i;
        (*m_PD)(row,col) = 1.;
      }
    }

    if (m_env.checkingLevel() >= 1) {
      // Check transpose operation
      D_M PDt(m_PD->transpose());
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", tests on m_PD"
                                << ": m_PD->numRowsLocal() = " << m_PD->numRowsLocal()
                                << ", m_PD->numCols() = "      << m_PD->numCols()
                                << ": PDt.numRowsLocal() = "  << PDt.numRowsLocal()
                                << ", PDt.numCols() = "       << PDt.numCols()
                                << std::endl;
      }

      D_M matShouldBeI1( *m_PD * PDt );
      D_M matI1        (m_v_space.zeroVector());
      for (unsigned int i = 0; i < matI1.numRowsLocal(); ++i) {
        matI1(i,i) = 1.;
      }
      matShouldBeI1 -= matI1;
      double auxNorm1 = matShouldBeI1.normFrob();

      D_M matShouldBeI2( PDt * *m_PD );
      D_M matI2        (m_v_space.zeroVector());
      for (unsigned int i = 0; i < matI2.numRowsLocal(); ++i) {
        matI2(i,i) = 1.;
      }
      matShouldBeI2 -= matI2;
      double auxNorm2 = matShouldBeI2.normFrob();

      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ", tests on m_PD"
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
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": finished forming 'P_D'"
                              << std::endl;
    }

    //********************************************************************************
    // Compute 'Dmat_BlockDiag_permut' matrix
    //********************************************************************************
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": m_Dmat_BlockDiag->numRowsLocal() = " << m_Dmat_BlockDiag->numRowsLocal()
                              << ", m_Dmat_BlockDiag->numCols() = "      << m_Dmat_BlockDiag->numCols()
                              << ", m_PD->numRowsLocal() = "             << m_PD->numRowsLocal()
                              << ", m_PD->numCols() = "                  << m_PD->numCols()
                              << std::endl;
    }

    *m_Dmat_BlockDiag_permut = *m_Dmat_BlockDiag * (m_PD->transpose());

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": finished computing 'm_Dmat_BlockDiag_permut'"
                              << std::endl;
    }
  } // if (allOutputsAreScalar)

  //********************************************************************************
  // Instantiate Smat spaces
  //********************************************************************************
  for (unsigned int i = 0; i < m_Imat_v_i_spaces.size(); ++i) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": before instantiating a m_Imat_v_i space"
                              << ", m_paper_Gs[" << i << "] = " << m_paper_Gs[i]
                              << std::endl;
    }
    m_Imat_v_i_spaces[i] = new VectorSpace<D_V,D_M>(m_env, "Imat_v_i_spaces_", m_paper_Gs[i], NULL); // to be deleted on destructor
    D_V unitVec(m_Imat_v_i_spaces[i]->zeroVector());
    unitVec.cwSet(1.);
    m_Imat_v_is[i] = new D_M(unitVec); // to be deleted on destructor
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << ": before instantiating the m_Rmat_v_i spaces"
                            << ", m_paper_n = " << m_paper_n
                            << std::endl;
  }
  for (unsigned int i = 0; i < m_Rmat_v_i_spaces.size(); ++i) {
    m_Rmat_v_i_spaces[i] = new VectorSpace<D_V,D_M>(m_env, "Smat_v_i_spaces_", m_paper_n, NULL); // to be deleted on destructor
    m_Rmat_v_is[i] = new D_M(m_Rmat_v_i_spaces[i]->zeroVector()); // to be deleted on destructor
  }

  unsigned int sumDims = 0;
  for (unsigned int i = 0; i < m_Smat_v_i_spaces.size(); ++i) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                              << ": before instantiating a m_Smat_v_i space"
                              << ", m_paper_Gs[" << i << "] = " << m_paper_Gs[i]
                              << std::endl;
    }
    m_Smat_v_i_spaces[i] = new VectorSpace<D_V,D_M>(m_env, "Smat_v_i_spaces_", m_paper_n*m_paper_Gs[i], NULL); // to be deleted on destructor
    sumDims += m_paper_n*m_paper_Gs[i];
    m_Smat_v_is[i] = new D_M(m_Smat_v_i_spaces[i]->zeroVector()); // to be deleted on destructor
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << ": finished instantiating the m_Smat_v_i spaces"
                            << ", m_paper_n = "       << m_paper_n
                            << ", m_paper_p_delta = " << m_paper_p_delta
                            << ", sumDims = "         << sumDims
                            << std::endl;
  }
  queso_require_equal_to_msg(sumDims, m_v_size, "'sumDims' and 'm_v_size' should be equal");

  //********************************************************************************
  // Instantiate 'v_hat_v_asterisk' matrices
  //********************************************************************************
  unsigned int sumNumRows = 0;
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < m_Smat_v_hat_v_asterisk_is.size(); ++i) {
    m_Rmat_v_hat_v_asterisk_is[i] = new D_M(m_env, m_paper_n_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    m_Smat_v_hat_v_asterisk_is[i] = new D_M(m_env, m_paper_n_space.map(), (unsigned int) 1); // to be deleted on destructor; Yes, only 1 column
    sumNumRows += m_paper_n_space.dimLocal();
    sumNumCols += 1;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << ": finished instantiating the m_Smat_v_hat_v_asterisk_i matrices"
                            << std::endl;
  }
  queso_require_equal_to_msg(sumNumRows, m_v_size, "'sumNumRows' and 'm_v_size' should be equal");
  queso_require_equal_to_msg(sumNumCols, m_paper_p_delta, "'sumNumCols' and 'm_paper_p_delta' should be equal");

  //********************************************************************************
  // Display information
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "KEY In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << "\n KEY  m_paper_p_x = "                     << m_paper_p_x
                            << "\n KEY  m_paper_n = "                       << m_paper_n
                            << "\n KEY  m_paper_n_ys_transformed.size() = " << m_paper_n_ys_transformed.size()
                            << "\n KEY  m_paper_n_ys_transformed =";
    for (unsigned int i = 0; i < m_paper_n_ys_transformed.size(); ++i) {
      *m_env.subDisplayFile() << " " << m_paper_n_ys_transformed[i];
    }
    *m_env.subDisplayFile() << "\n KEY  m_paper_n_y = "                     << m_paper_n_y
                            << "\n KEY  m_paper_p_delta = "                 << m_paper_p_delta
                            << "\n KEY  m_paper_F = "                       << m_paper_F
                            << "\n KEY  m_paper_Gs.size() = "               << m_paper_Gs.size()
                            << "\n KEY  m_paper_Gs =";
    for (unsigned int i = 0; i < m_paper_Gs.size(); ++i) {
      *m_env.subDisplayFile() << " " << m_paper_Gs[i];
    }
    *m_env.subDisplayFile() << "\n KEY  m_5lambdaYDim   = "                 << m_5lambdaYDim
                            << ", m_5lambdaYGammaAVec = "                   << m_5lambdaYGammaAVec
                            << ", m_5lambdaYGammaBVec = "                   << m_5lambdaYGammaBVec
                            << "\n KEY  m_6lambdaVDim   = "                 << m_6lambdaVDim
                            << ", m_6lambdaVGammaAVec = "                   << m_6lambdaVGammaAVec
                            << ", m_6lambdaVGammaBVec = "                   << m_6lambdaVGammaBVec
                            << "\n KEY  m_7rhoVDim      = "                 << m_7rhoVDim
                            << ", m_7rhoVBetaAVec = "                       << m_7rhoVBetaAVec
                            << ", m_7rhoVBetaBVec = "                       << m_7rhoVBetaBVec
                            << "\n KEY  m_8thetaDim     = "                 << m_8thetaDim
                            << "\n KEY  full 'y'   vector size = "          << m_paper_n_y // = experimentStorage->yVec_transformed().sizeLocal()
                            << std::endl;
  }

  //********************************************************************************
  // Make checks
  //********************************************************************************
  queso_require_equal_to_msg(experimentStorage.yVec_transformed().sizeLocal(), m_paper_n_y, "incompatible calculations for 'y' vector size");

  queso_require_equal_to_msg(m_paper_n, m_paper_n_ys_transformed.size(), "'m_paper_n' and 'm_paper_n_ys_transformed.size()' should be equal");

  queso_require_equal_to_msg(m_paper_F, m_paper_Gs.size(), "'m_paper_F' and 'm_paper_Gs.size()' should be equal");

  unsigned int sumGs = 0;
  for (unsigned int i = 0; i < m_paper_F; ++i) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
       *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                               << ": m_paper_Gs[" << i << "] = " << m_paper_Gs[i]
                               << std::endl;
    }
    sumGs += m_paper_Gs[i];
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << ": sumGs = " << sumGs
                            << std::endl;
  }
  queso_require_equal_to_msg(m_paper_p_delta, sumGs, "'m_paper_p_delta' and 'sumGs' should be equal");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::constructor()"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M>
GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>::~GcmExperimentInfo()
{
  for (unsigned int i = 0; i < m_Smat_v_hat_v_asterisk_is.size(); ++i) {
    delete m_Smat_v_hat_v_asterisk_is[i]; // to be deleted on destructor
    m_Smat_v_hat_v_asterisk_is[i] = NULL;
    delete m_Rmat_v_hat_v_asterisk_is[i]; // to be deleted on destructor
    m_Rmat_v_hat_v_asterisk_is[i] = NULL;
  }

  for (unsigned int i = 0; i < m_Smat_v_i_spaces.size(); ++i) {
    delete m_Smat_v_is[i]; // to be deleted on destructor
    m_Smat_v_is[i] = NULL;
    delete m_Smat_v_i_spaces[i]; // to be deleted on destructor
    m_Smat_v_i_spaces[i] = NULL;

    delete m_Rmat_v_is[i]; // to be deleted on destructor
    m_Rmat_v_is[i] = NULL;
    delete m_Rmat_v_i_spaces[i]; // to be deleted on destructor
    m_Rmat_v_i_spaces[i] = NULL;

    delete m_Imat_v_is[i]; // to be deleted on destructor
    m_Imat_v_is[i] = NULL;
    delete m_Imat_v_i_spaces[i]; // to be deleted on destructor
    m_Imat_v_i_spaces[i] = NULL;
  }

  delete m_Dmat_BlockDiag_permut;
  delete m_PD;
}

}  // End namespace QUESO

template class QUESO::GcmExperimentInfo<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
