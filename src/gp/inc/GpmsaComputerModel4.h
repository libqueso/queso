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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_GCM_4_H__
#define __UQ_GCM_4_H__

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(
  const P_V&         input_1lambdaEtaVec,
  const P_V&         input_2lambdaWVec,
  const P_V&         input_3rhoWVec,
  const P_V&         input_4lambdaSVec,
  const P_V&         input_5lambdaYVec,
  const P_V&         input_6lambdaVVec,
  const P_V&         input_7rhoVVec,
  const P_V&         input_8thetaVec,
        unsigned int outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                            << ", outerCounter = " << outerCounter
                            << ": m_formCMatrix = "                                        << m_formCMatrix
                            << ", m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC = " << m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC
                            << ", m_Cmat = "                                               << m_z->m_Cmat
                            << ", m_cMatIsRankDefficient = "                               << m_cMatIsRankDefficient
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO((m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC == true),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)",
                      "'m_useTildeLogicForRankDefficientC' should be 'false'");

  //********************************************************************************
  // Form '\Sigma_z' matrix
  //********************************************************************************
  // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
  // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
  // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
  // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
  this->formSigma_z(input_2lambdaWVec,
                    input_3rhoWVec,
                    input_4lambdaSVec,
                    input_6lambdaVVec,
                    input_7rhoVVec,
                    input_8thetaVec,
                    outerCounter);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                            << ", outerCounter = "           << outerCounter
                            << ": finished forming 'm_tmp_Smat_z'"
                            << "\n m_tmp_Smat_z contents = " << m_z->m_tmp_Smat_z
                            << std::endl;
  }

  if (m_env.displayVerbosity() >= 4) {
    double       sigmaZLnDeterminant = m_z->m_tmp_Smat_z.lnDeterminant();
    unsigned int sigmaZRank          = m_z->m_tmp_Smat_z.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int sigmaZRank14        = m_z->m_tmp_Smat_z.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                              << ", outerCounter = "                 << outerCounter
                              << ", m_tmp_Smat_z.numRowsLocal() = "  << m_z->m_tmp_Smat_z.numRowsLocal()
                              << ", m_tmp_Smat_z.numCols() = "       << m_z->m_tmp_Smat_z.numCols()
                              << ", m_tmp_Smat_z.lnDeterminant() = " << sigmaZLnDeterminant
                              << ", m_tmp_Smat_z.rank(0.,1.e-8) = "  << sigmaZRank
                              << ", m_tmp_Smat_z.rank(0.,1.e-14) = " << sigmaZRank14
                              << std::endl;
    }
  }


  //********************************************************************************
  // Form '\Sigma_{extra}' matrix
  //********************************************************************************
  if (m_allOutputsAreScalar) {
    // ppp
  }
  else {
    m_z->m_tmp_Smat_extra.cwSet(0.);
    m_z->m_tmp_Smat_extra.cwSet(                                         0,                                    0,(1./input_5lambdaYVec  [0]) * *m_j->m_Bop_t__Wy__Bop__inv);
    m_z->m_tmp_Smat_extra.cwSet(m_j->m_Bop_t__Wy__Bop__inv->numRowsLocal(),m_j->m_Bop_t__Wy__Bop__inv->numCols(),(1./input_1lambdaEtaVec[0]) * *m_s->m_Kt_K_inv           );
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                            << ", outerCounter = "               << outerCounter
                            << ": finished forming 'm_tmp_Smat_extra'"
                            << ", input_5lambdaYVec[0] = "       << input_5lambdaYVec[0]
                            << ", input_1lambdaEtaVec[0] = "     << input_1lambdaEtaVec[0]
                            << "\n m_Bop_t__Wy__Bop__inv = "     << *m_j->m_Bop_t__Wy__Bop__inv
                            << "\n m_Kt_K_inv = "                << *m_s->m_Kt_K_inv
                            << "\n m_tmp_Smat_extra contents = " << m_z->m_tmp_Smat_extra
                            << std::endl;
  }

  if (m_env.displayVerbosity() >= 4) {
    double       extraLnDeterminant = m_z->m_tmp_Smat_extra.lnDeterminant();
    unsigned int extraRank          = m_z->m_tmp_Smat_extra.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int extraRank14        = m_z->m_tmp_Smat_extra.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                              << ", outerCounter = "                     << outerCounter
                              << ", m_tmp_Smat_extra.numRowsLocal() = "  << m_z->m_tmp_Smat_extra.numRowsLocal()
                              << ", m_tmp_Smat_extra.numCols() = "       << m_z->m_tmp_Smat_extra.numCols()
                              << ", m_tmp_Smat_extra.lnDeterminant() = " << extraLnDeterminant
                              << ", m_tmp_Smat_extra.rank(0.,1.e-8) = "  << extraRank
                              << ", m_tmp_Smat_extra.rank(0.,1.e-14) = " << extraRank14
                              << std::endl;
    }
  }

  //********************************************************************************
  // Compute '\Sigma_z_hat' matrix
  //********************************************************************************
  m_z->m_tmp_Smat_z_hat = m_z->m_tmp_Smat_z + m_z->m_tmp_Smat_extra;

  if (m_env.displayVerbosity() >= 4) {
    double       zHatLnDeterminant = m_z->m_tmp_Smat_z_hat.lnDeterminant();
    unsigned int zHatRank          = m_z->m_tmp_Smat_z_hat.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int zHatRank14        = m_z->m_tmp_Smat_z_hat.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                              << ", outerCounter = "                     << outerCounter
                              << ", m_tmp_Smat_z_hat.numRowsLocal() = "  << m_z->m_tmp_Smat_z_hat.numRowsLocal()
                              << ", m_tmp_Smat_z_hat.numCols() = "       << m_z->m_tmp_Smat_z_hat.numCols()
                              << ", m_tmp_Smat_z_hat.lnDeterminant() = " << zHatLnDeterminant
                              << ", m_tmp_Smat_z_hat.rank(0.,1.e-8) = "  << zHatRank
                              << ", m_tmp_Smat_z_hat.rank(0.,1.e-14) = " << zHatRank14
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(1)"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}
									
// Case with no experimental data // checar
template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(
  const P_V&         input_1lambdaEtaVec,
  const P_V&         input_2lambdaWVec,
  const P_V&         input_3rhoWVec,
  const P_V&         input_4lambdaSVec,
        unsigned int outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                            << ", outerCounter = " << outerCounter
                            << ": m_formCMatrix = "                                        << m_formCMatrix
                            << ", m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC = " << m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC
                            << ", m_Cmat = "                                               << m_z->m_Cmat
                            << ", m_cMatIsRankDefficient = "                               << m_cMatIsRankDefficient
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO((m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC == true),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)",
                      "'m_useTildeLogicForRankDefficientC' should be 'false'");

  //********************************************************************************
  // Form '\Sigma_z' matrix
  //********************************************************************************
  // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
  // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
  // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
  // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
  this->formSigma_z(input_2lambdaWVec,
                    input_3rhoWVec,
                    input_4lambdaSVec,
                    outerCounter);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                            << ", outerCounter = "           << outerCounter
                            << ": finished forming 'm_tmp_Smat_z'"
                            << "\n m_tmp_Smat_z contents = " << m_z->m_tmp_Smat_z
                            << std::endl;
  }

  if (m_env.displayVerbosity() >= 4) {
    double       sigmaZLnDeterminant = m_z->m_tmp_Smat_z.lnDeterminant();
    unsigned int sigmaZRank          = m_z->m_tmp_Smat_z.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int sigmaZRank14        = m_z->m_tmp_Smat_z.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                              << ", outerCounter = "                 << outerCounter
                              << ", m_tmp_Smat_z.numRowsLocal() = "  << m_z->m_tmp_Smat_z.numRowsLocal()
                              << ", m_tmp_Smat_z.numCols() = "       << m_z->m_tmp_Smat_z.numCols()
                              << ", m_tmp_Smat_z.lnDeterminant() = " << sigmaZLnDeterminant
                              << ", m_tmp_Smat_z.rank(0.,1.e-8) = "  << sigmaZRank
                              << ", m_tmp_Smat_z.rank(0.,1.e-14) = " << sigmaZRank14
                              << std::endl;
    }
  }

#if 0 // Case with no experimental data // checar
  //********************************************************************************
  // Form '\Sigma_{extra}' matrix
  //********************************************************************************
  if (m_allOutputsAreScalar) {
    // ppp
  }
  else {
    m_z->m_tmp_Smat_extra.cwSet(0.);
    m_z->m_tmp_Smat_extra.cwSet(                                         0,                                    0,(1./input_5lambdaYVec  [0]) * *m_j->m_Bop_t__Wy__Bop__inv);
    m_z->m_tmp_Smat_extra.cwSet(m_j->m_Bop_t__Wy__Bop__inv->numRowsLocal(),m_j->m_Bop_t__Wy__Bop__inv->numCols(),(1./input_1lambdaEtaVec[0]) * *m_s->m_Kt_K_inv           );
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                            << ", outerCounter = "               << outerCounter
                            << ": finished forming 'm_tmp_Smat_extra'"
                            << ", input_5lambdaYVec[0] = "       << input_5lambdaYVec[0]
                            << ", input_1lambdaEtaVec[0] = "     << input_1lambdaEtaVec[0]
                            << "\n m_Bop_t__Wy__Bop__inv = "     << *m_j->m_Bop_t__Wy__Bop__inv
                            << "\n m_Kt_K_inv = "                << *m_s->m_Kt_K_inv
                            << "\n m_tmp_Smat_extra contents = " << m_z->m_tmp_Smat_extra
                            << std::endl;
  }

  if (m_env.displayVerbosity() >= 4) {
    double       extraLnDeterminant = m_z->m_tmp_Smat_extra.lnDeterminant();
    unsigned int extraRank          = m_z->m_tmp_Smat_extra.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int extraRank14        = m_z->m_tmp_Smat_extra.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                              << ", outerCounter = "                     << outerCounter
                              << ", m_tmp_Smat_extra.numRowsLocal() = "  << m_z->m_tmp_Smat_extra.numRowsLocal()
                              << ", m_tmp_Smat_extra.numCols() = "       << m_z->m_tmp_Smat_extra.numCols()
                              << ", m_tmp_Smat_extra.lnDeterminant() = " << extraLnDeterminant
                              << ", m_tmp_Smat_extra.rank(0.,1.e-8) = "  << extraRank
                              << ", m_tmp_Smat_extra.rank(0.,1.e-14) = " << extraRank14
                              << std::endl;
    }
  }

  //********************************************************************************
  // Compute '\Sigma_z_hat' matrix
  //********************************************************************************
  m_z->m_tmp_Smat_z_hat = m_z->m_tmp_Smat_z + m_z->m_tmp_Smat_extra;

  if (m_env.displayVerbosity() >= 4) {
    double       zHatLnDeterminant = m_z->m_tmp_Smat_z_hat.lnDeterminant();
    unsigned int zHatRank          = m_z->m_tmp_Smat_z_hat.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int zHatRank14        = m_z->m_tmp_Smat_z_hat.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                              << ", outerCounter = "                     << outerCounter
                              << ", m_tmp_Smat_z_hat.numRowsLocal() = "  << m_z->m_tmp_Smat_z_hat.numRowsLocal()
                              << ", m_tmp_Smat_z_hat.numCols() = "       << m_z->m_tmp_Smat_z_hat.numCols()
                              << ", m_tmp_Smat_z_hat.lnDeterminant() = " << zHatLnDeterminant
                              << ", m_tmp_Smat_z_hat.rank(0.,1.e-8) = "  << zHatRank
                              << ", m_tmp_Smat_z_hat.rank(0.,1.e-14) = " << zHatRank14
                              << std::endl;
    }
  }
#endif
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_hat(2)"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}
									
template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat(
  const P_V&         input_1lambdaEtaVec,
  const P_V&         input_2lambdaWVec,
  const P_V&         input_3rhoWVec,
  const P_V&         input_4lambdaSVec,
  const P_V&         input_5lambdaYVec,
  const P_V&         input_6lambdaVVec,
  const P_V&         input_7rhoVVec,
  const P_V&         input_8thetaVec,
        unsigned int outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                            << ", outerCounter = " << outerCounter
                            << ": m_formCMatrix = "                                        << m_formCMatrix
                            << ", m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC = " << m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC
                            << ", m_Cmat = "                                               << m_z->m_Cmat
                            << ", m_cMatIsRankDefficient = "                               << m_cMatIsRankDefficient
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_formCMatrix == false,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()",
                      "'m_Cmat' should have been requested");

  UQ_FATAL_TEST_MACRO(m_z->m_Cmat == NULL,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()",
                      "'m_Cmat' should have been formed");

  UQ_FATAL_TEST_MACRO(m_cMatIsRankDefficient == false,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()",
                      "'m_Cmat' should be rank defficient");

  UQ_FATAL_TEST_MACRO((m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC == false),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()",
                      "'m_useTildeLogicForRankDefficientC' should be 'true'");

  //********************************************************************************
  // Form '\Sigma_z' matrix
  //********************************************************************************
  // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
  // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
  // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
  // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
  this->formSigma_z(input_2lambdaWVec,
                    input_3rhoWVec,
                    input_4lambdaSVec,
                    input_6lambdaVVec,
                    input_7rhoVVec,
                    input_8thetaVec,
                    outerCounter);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                            << ", outerCounter = "           << outerCounter
                            << ": finished forming 'm_tmp_Smat_z'"
                            << "\n m_tmp_Smat_z contents = " << m_z->m_tmp_Smat_z
                            << std::endl;
  }

  if (m_env.displayVerbosity() >= 4) {
    double       sigmaZLnDeterminant = m_z->m_tmp_Smat_z.lnDeterminant();
    unsigned int sigmaZRank          = m_z->m_tmp_Smat_z.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int sigmaZRank14        = m_z->m_tmp_Smat_z.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                              << ", outerCounter = "                 << outerCounter
                              << ", m_tmp_Smat_z.numRowsLocal() = "  << m_z->m_tmp_Smat_z.numRowsLocal()
                              << ", m_tmp_Smat_z.numCols() = "       << m_z->m_tmp_Smat_z.numCols()
                              << ", m_tmp_Smat_z.lnDeterminant() = " << sigmaZLnDeterminant
                              << ", m_tmp_Smat_z.rank(0.,1.e-8) = "  << sigmaZRank
                              << ", m_tmp_Smat_z.rank(0.,1.e-14) = " << sigmaZRank14
                              << std::endl;
    }
  }

  //********************************************************************************
  // Form 'L . \Sigma_z . L^T' matrix
  //********************************************************************************
  m_zt->m_tmp_Smat_z_tilde = m_zt->m_Lmat * (m_z->m_tmp_Smat_z * m_zt->m_Lmat_t);

  if (m_env.displayVerbosity() >= 4) {
    double       sigmaZTildeLnDeterminant = m_zt->m_tmp_Smat_z_tilde.lnDeterminant();
    unsigned int sigmaZTildeRank          = m_zt->m_tmp_Smat_z_tilde.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int sigmaZTildeRank14        = m_zt->m_tmp_Smat_z_tilde.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                              << ", outerCounter = "                        << outerCounter
                              << ", m_tmp_Smat_z_tilde->numRowsLocal() = "  << m_zt->m_tmp_Smat_z_tilde.numRowsLocal()
                              << ", m_tmp_Smat_z_tilde->numCols() = "       << m_zt->m_tmp_Smat_z_tilde.numCols()
                              << ", m_tmp_Smat_z_tilde->lnDeterminant() = " << sigmaZTildeLnDeterminant
                              << ", m_tmp_Smat_z_tilde->rank(0.,1.e-8) = "  << sigmaZTildeRank
                              << ", m_tmp_Smat_z_tilde->rank(0.,1.e-14) = " << sigmaZTildeRank14
                              << std::endl;
    }
  }

  //********************************************************************************
  // Form '\Sigma_{extra}_tilde' matrix
  //********************************************************************************
  m_zt->m_tmp_Smat_extra_tilde.cwSet(0.);
  m_zt->m_tmp_Smat_extra_tilde.cwSet(                                           0,                                      0,(1./input_5lambdaYVec  [0]) * (m_jt->m_Btildet_Wy_Btilde_inv));
  m_zt->m_tmp_Smat_extra_tilde.cwSet(m_jt->m_Btildet_Wy_Btilde_inv.numRowsLocal(),m_jt->m_Btildet_Wy_Btilde_inv.numCols(),(1./input_1lambdaEtaVec[0]) * (m_st->m_Ktildet_Ktilde_inv)   );

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                            << ", outerCounter = "                     << outerCounter
                            << ": finished forming 'm_tmp_Smat_extra_tilde'"
                            << ", input_5lambdaYVec[0] = "             << input_5lambdaYVec[0]
                            << ", input_1lambdaEtaVec[0] = "           << input_1lambdaEtaVec[0]
                            << "\n m_Btildet_Wy_Btilde_inv = "         << m_jt->m_Btildet_Wy_Btilde_inv
                            << "\n m_Ktildet_Ktilde_inv = "            << m_st->m_Ktildet_Ktilde_inv
                            << "\n m_tmp_Smat_extra_tilde contents = " << m_zt->m_tmp_Smat_extra_tilde
                            << std::endl;
  }

  if (m_env.displayVerbosity() >= 4) {
    double       extraTildeLnDeterminant = m_zt->m_tmp_Smat_extra_tilde.lnDeterminant();
    unsigned int extraTildeRank          = m_zt->m_tmp_Smat_extra_tilde.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int extraTildeRank14        = m_zt->m_tmp_Smat_extra_tilde.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                              << ", outerCounter = "                           << outerCounter
                              << ", m_tmp_Smat_extra_tilde.numRowsLocal() = "  << m_zt->m_tmp_Smat_extra_tilde.numRowsLocal()
                              << ", m_tmp_Smat_extra_tilde.numCols() = "       << m_zt->m_tmp_Smat_extra_tilde.numCols()
                              << ", m_tmp_Smat_extra_tilde.lnDeterminant() = " << extraTildeLnDeterminant
                              << ", m_tmp_Smat_extra_tilde.rank(0.,1.e-8) = "  << extraTildeRank
                              << ", m_tmp_Smat_extra_tilde.rank(0.,1.e-14) = " << extraTildeRank14
                              << std::endl;
    }
  }

  //********************************************************************************
  // Compute '\Sigma_z_tilde_hat' matrix
  //********************************************************************************
  m_zt->m_tmp_Smat_z_tilde_hat = m_zt->m_tmp_Smat_z_tilde + m_zt->m_tmp_Smat_extra_tilde;

  if (m_env.displayVerbosity() >= 4) {
    double       zTildeHatLnDeterminant = m_zt->m_tmp_Smat_z_tilde_hat.lnDeterminant();
    unsigned int zTildeHatRank          = m_zt->m_tmp_Smat_z_tilde_hat.rank(0.,1.e-8 ); // todo: should be an option
    unsigned int zTildeHatRank14        = m_zt->m_tmp_Smat_z_tilde_hat.rank(0.,1.e-14);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                              << ", outerCounter = "                            << outerCounter
                              << ", m_tmp_Smat_z_tilde_hat->numRowsLocal() = "  << m_zt->m_tmp_Smat_z_tilde_hat.numRowsLocal()
                              << ", m_tmp_Smat_z_tilde_hat->numCols() = "       << m_zt->m_tmp_Smat_z_tilde_hat.numCols()
                              << ", m_tmp_Smat_z_tilde_hat->lnDeterminant() = " << zTildeHatLnDeterminant
                              << ", m_tmp_Smat_z_tilde_hat->rank(0.,1.e-8) = "  << zTildeHatRank
                              << ", m_tmp_Smat_z_tilde_hat->rank(0.,1.e-14) = " << zTildeHatRank14
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z_tilde_hat()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}
									
template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(
  const P_V&         input_2lambdaWVec,
  const P_V&         input_3rhoWVec,
  const P_V&         input_4lambdaSVec,
  const P_V&         input_6lambdaVVec,
  const P_V&         input_7rhoVVec,
  const P_V&         input_8thetaVec,
        unsigned int outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

  this->memoryCheck(90);

  // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
  // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
  // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
  // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
  //********************************************************************************
  // Compute '\Sigma' matrices
  // \Sigma_v:
  // --> Uses page 576-b
  // --> Uses R(all x's;\rho_v_i[size p_x]) and formula (2) to "each pair" of experimental input settings
  // --> \Sigma_v_i = (1/\lambda_v_i).I_|G_i| [X] R(...) is (n.|G_i|) x (n.|G_i|), i = 1,...,F
  // --> \Sigma_v is (n.p_delta) x (n.p_delta) 
  // \Sigma_u:
  // --> Uses page 576-b
  // --> Uses R(all x's,one \theta;\rho_w_i[size p_x+p_t]) and formula (1) to "each pair" of experimental input settings (correlations depend only on x dimensions)
  // --> \Sigma_u_i = (1/\lambda_w_i).R(...) is n x n, i = 1,...,p_eta
  // --> \Sigma_u is (n.p_eta) x (n.p_eta) 
  // \Sigma_w:
  // --> Uses page 575-a
  // --> Uses R(all x^*'s,all t^*'s;\rho_w_i[size p_x+p_t]) and formula (1) to "each pair" of input settings in the design
  // --> \Sigma_w_i = (1/\lambda_w_i).R(...) is m x m, i = 1,...,p_eta
  // --> \Sigma_w is (m.p_eta) x (m.p_eta) 
  // \Sigma_u,w:
  // --> Uses page 577-a
  // --> Uses R(all x's,one \theta,all x^*'s,all t^*'s;\rho_w_i[size p_x+p_t]) and formula (1)
  // --> \Sigma_u,w_i = (1/\lambda_w_i).R(...) is n x m, i = 1,...,p_eta
  // --> \Sigma_u,w is (n.p_eta) x (m.p_eta) 
  //********************************************************************************
  unsigned int initialPos = 0;
  for (unsigned int i = 0; i < m_e->m_Smat_v_i_spaces.size(); ++i) {
    input_7rhoVVec.cwExtract(initialPos,m_e->m_tmp_rho_v_vec);
    initialPos += m_e->m_tmp_rho_v_vec.sizeLocal();
    m_e->m_Rmat_v_is[i]->cwSet(0.);
    this->fillR_formula2_for_Sigma_v(m_e->m_paper_xs_standard,
                                     m_e->m_tmp_rho_v_vec,
                                     *(m_e->m_Rmat_v_is[i]), // IMPORTANT-28
                                     outerCounter);

    m_e->m_Smat_v_is[i]->cwSet(0.);
    m_e->m_Smat_v_is[i]->fillWithTensorProduct(0,0,*(m_e->m_Imat_v_is[i]),*(m_e->m_Rmat_v_is[i]),true,true); // IMPORTANT-28
    *(m_e->m_Smat_v_is[i]) *= (1./input_6lambdaVVec[i]);
  }
  m_e->m_Smat_v.cwSet(0.);
  m_e->m_Smat_v.fillWithBlocksDiagonally(0,0,m_e->m_Smat_v_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_v'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_e->m_Smat_v.subWriteContents("Sigma_v",
                                     "mat_Sigma_v",
                                     "m",
                                     tmpSet);
    }
  }

  this->memoryCheck(91);

  initialPos = 0;
  for (unsigned int i = 0; i < m_j->m_Smat_u_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_j->m_Rmat_u_is[i]->cwSet(0.);
    this->fillR_formula1_for_Sigma_u(m_e->m_paper_xs_standard,
                                     input_8thetaVec,
                                     m_s->m_tmp_rho_w_vec,
                                     *(m_j->m_Rmat_u_is[i]),
                                     outerCounter);

    m_j->m_Smat_u_is[i]->cwSet(0.);
    *(m_j->m_Smat_u_is[i]) = (1./input_2lambdaWVec[i]) * *(m_j->m_Rmat_u_is[i]);
    for (unsigned int j = 0; j < m_j->m_Smat_u_is[i]->numRowsLocal(); ++j) {
      (*(m_j->m_Smat_u_is[i]))(j,j) +=  1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
    }
  }
  m_j->m_Smat_u.cwSet(0.);
  m_j->m_Smat_u.fillWithBlocksDiagonally(0,0,m_j->m_Smat_u_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_u'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_j->m_Smat_u.subWriteContents("Sigma_u",
                                     "mat_Sigma_u",
                                     "m",
                                     tmpSet);
    }
  }

  this->memoryCheck(92);

  initialPos = 0;
  for (unsigned int i = 0; i < m_s->m_Smat_w_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_s->m_Rmat_w_is[i]->cwSet(0.); // This matrix is square: m_paper_m X m_paper_m
    this->fillR_formula1_for_Sigma_w(m_s->m_paper_xs_asterisks_standard, // IMPORTANT
                                     m_s->m_paper_ts_asterisks_standard,
                                     m_s->m_tmp_rho_w_vec,
                                     *(m_s->m_Rmat_w_is[i]),
                                     outerCounter);

    m_s->m_Smat_w_is[i]->cwSet(0.);
    *(m_s->m_Smat_w_is[i]) = (1./input_2lambdaWVec[i]) * *(m_s->m_Rmat_w_is[i]);
    for (unsigned int j = 0; j < m_s->m_Smat_w_is[i]->numRowsLocal(); ++j) {
      (*(m_s->m_Smat_w_is[i]))(j,j) +=  1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
    }
  }
  m_s->m_Smat_w.cwSet(0.);
  m_s->m_Smat_w.fillWithBlocksDiagonally(0,0,m_s->m_Smat_w_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_w'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_s->m_Smat_w.subWriteContents("Sigma_w",
                                     "mat_Sigma_w",
                                     "m",
                                     tmpSet);
    }
  }

  this->memoryCheck(93);

  initialPos = 0;
  for (unsigned int i = 0; i < m_j->m_Smat_uw_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_j->m_Rmat_uw_is[i]->cwSet(0.);
    this->fillR_formula1_for_Sigma_uw(m_e->m_paper_xs_standard,
                                      input_8thetaVec,
                                      m_s->m_paper_xs_asterisks_standard,
                                      m_s->m_paper_ts_asterisks_standard,
                                      m_s->m_tmp_rho_w_vec,*(m_j->m_Rmat_uw_is[i]),
                                      outerCounter);

    m_j->m_Smat_uw_is[i]->cwSet(0.);
    *(m_j->m_Smat_uw_is[i]) = (1./input_2lambdaWVec[i]) * *(m_j->m_Rmat_uw_is[i]);
  }
  m_j->m_Smat_uw.cwSet(0.);
  m_j->m_Smat_uw.fillWithBlocksDiagonally(0,0,m_j->m_Smat_uw_is,true,true);
  m_j->m_Smat_uw_t.fillWithTranspose(0,0,m_j->m_Smat_uw,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_j->m_Smat_uw'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_j->m_Smat_uw.subWriteContents("Sigma_uw",
                                      "mat_Sigma_uw",
                                      "m",
                                      tmpSet);
    }
  }

  this->memoryCheck(94);

  //********************************************************************************
  // Form '\Sigma_z' matrix
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << ": m_tmp_Smat_z.numRowsLocal() = " << m_z->m_tmp_Smat_z.numRowsLocal()
                            << ", m_tmp_Smat_z.numCols() = "      << m_z->m_tmp_Smat_z.numCols()
                            << ", m_Smat_v.numRowsLocal() = "     << m_e->m_Smat_v.numRowsLocal()
                            << ", m_Smat_v.numCols() = "          << m_e->m_Smat_v.numCols()
                            << ", m_Smat_u.numRowsLocal() = "     << m_j->m_Smat_u.numRowsLocal()
                            << ", m_Smat_u.numCols() = "          << m_j->m_Smat_u.numCols()
                            << ", m_Smat_w.numRowsLocal() = "     << m_s->m_Smat_w.numRowsLocal()
                            << ", m_Smat_w.numCols() = "          << m_s->m_Smat_w.numCols()
                            << ", m_Smat_uw.numRowsLocal() = "    << m_j->m_Smat_uw.numRowsLocal()
                            << ", m_Smat_uw.numCols() = "         << m_j->m_Smat_uw.numCols()
                            << ", m_Smat_v_i_spaces.size() = "    << m_e->m_Smat_v_i_spaces.size()
                            << ", m_Smat_u_is.size() = "          << m_j->m_Smat_u_is.size()
                            << ", m_Smat_w_is.size() = "          << m_s->m_Smat_w_is.size()
                            << std::endl;
  }

  this->memoryCheck(95);

  m_z->m_tmp_Smat_z.cwSet(0.);
  if (m_allOutputsAreScalar) {
    // ppp
  }
  else {
    m_z->m_tmp_Smat_z.cwSet(0,0,m_e->m_Smat_v);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal(),                             m_e->m_Smat_v.numCols(),                        m_j->m_Smat_u);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal(),                             m_e->m_Smat_v.numCols()+m_j->m_Smat_u.numCols(),m_j->m_Smat_uw);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal()+m_j->m_Smat_u.numRowsLocal(),m_e->m_Smat_v.numCols(),                        m_j->m_Smat_uw_t);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal()+m_j->m_Smat_u.numRowsLocal(),m_e->m_Smat_v.numCols()+m_j->m_Smat_u.numCols(),m_s->m_Smat_w);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

// Case with no experimental data // checar
template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
  GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(
  const P_V&         input_2lambdaWVec,
  const P_V&         input_3rhoWVec,
  const P_V&         input_4lambdaSVec,
        unsigned int outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());

  this->memoryCheck(90);

  // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
  // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
  // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
  // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
  //********************************************************************************
  // Compute '\Sigma' matrices
  // \Sigma_v:
  // --> Uses page 576-b
  // --> Uses R(all x's;\rho_v_i[size p_x]) and formula (2) to "each pair" of experimental input settings
  // --> \Sigma_v_i = (1/\lambda_v_i).I_|G_i| [X] R(...) is (n.|G_i|) x (n.|G_i|), i = 1,...,F
  // --> \Sigma_v is (n.p_delta) x (n.p_delta) 
  // \Sigma_u:
  // --> Uses page 576-b
  // --> Uses R(all x's,one \theta;\rho_w_i[size p_x+p_t]) and formula (1) to "each pair" of experimental input settings (correlations depend only on x dimensions)
  // --> \Sigma_u_i = (1/\lambda_w_i).R(...) is n x n, i = 1,...,p_eta
  // --> \Sigma_u is (n.p_eta) x (n.p_eta) 
  // \Sigma_w:
  // --> Uses page 575-a
  // --> Uses R(all x^*'s,all t^*'s;\rho_w_i[size p_x+p_t]) and formula (1) to "each pair" of input settings in the design
  // --> \Sigma_w_i = (1/\lambda_w_i).R(...) is m x m, i = 1,...,p_eta
  // --> \Sigma_w is (m.p_eta) x (m.p_eta) 
  // \Sigma_u,w:
  // --> Uses page 577-a
  // --> Uses R(all x's,one \theta,all x^*'s,all t^*'s;\rho_w_i[size p_x+p_t]) and formula (1)
  // --> \Sigma_u,w_i = (1/\lambda_w_i).R(...) is n x m, i = 1,...,p_eta
  // --> \Sigma_u,w is (n.p_eta) x (m.p_eta) 
  //********************************************************************************
#if 0 // Case with no experimental data // checar
  unsigned int initialPos = 0;
  for (unsigned int i = 0; i < m_e->m_Smat_v_i_spaces.size(); ++i) {
    input_7rhoVVec.cwExtract(initialPos,m_e->m_tmp_rho_v_vec);
    initialPos += m_e->m_tmp_rho_v_vec.sizeLocal();
    m_e->m_Rmat_v_is[i]->cwSet(0.);
    this->fillR_formula2_for_Sigma_v(m_e->m_paper_xs_standard,
                                     m_e->m_tmp_rho_v_vec,
                                     *(m_e->m_Rmat_v_is[i]),
                                     outerCounter);

    m_e->m_Smat_v_is[i]->cwSet(0.);
    m_e->m_Smat_v_is[i]->fillWithTensorProduct(0,0,*(m_e->m_Imat_v_is[i]),*(m_e->m_Rmat_v_is[i]),true,true);
    *(m_e->m_Smat_v_is[i]) *= (1./input_6lambdaVVec[i]);
  }
  m_e->m_Smat_v.cwSet(0.);
  m_e->m_Smat_v.fillWithBlocksDiagonally(0,0,m_e->m_Smat_v_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_v'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_e->m_Smat_v.subWriteContents("Sigma_v",
                                     "mat_Sigma_v",
                                     "m",
                                     tmpSet);
    }
  }
  this->memoryCheck(91);

  initialPos = 0;
  for (unsigned int i = 0; i < m_j->m_Smat_u_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_j->m_Rmat_u_is[i]->cwSet(0.);
    this->fillR_formula1_for_Sigma_u(m_e->m_paper_xs_standard,
                                     input_8thetaVec,
                                     m_s->m_tmp_rho_w_vec,
                                     *(m_j->m_Rmat_u_is[i]),
                                     outerCounter);

    m_j->m_Smat_u_is[i]->cwSet(0.);
    *(m_j->m_Smat_u_is[i]) = (1./input_2lambdaWVec[i]) * *(m_j->m_Rmat_u_is[i]);
    for (unsigned int j = 0; j < m_j->m_Smat_u_is[i]->numRowsLocal(); ++j) {
      (*(m_j->m_Smat_u_is[i]))(j,j) +=  1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
    }
  }
  m_j->m_Smat_u.cwSet(0.);
  m_j->m_Smat_u.fillWithBlocksDiagonally(0,0,m_j->m_Smat_u_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_u'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_j->m_Smat_u.subWriteContents("Sigma_u",
                                     "mat_Sigma_u",
                                     "m",
                                     tmpSet);
    }
  }

  this->memoryCheck(92);

  initialPos = 0;
  for (unsigned int i = 0; i < m_s->m_Smat_w_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_s->m_Rmat_w_is[i]->cwSet(0.); // This matrix is square: m_paper_m X m_paper_m
    this->fillR_formula1_for_Sigma_w(m_s->m_paper_xs_asterisks_standard,
                                     m_s->m_paper_ts_asterisks_standard,
                                     m_s->m_tmp_rho_w_vec,
                                     *(m_s->m_Rmat_w_is[i]),
                                     outerCounter);

    m_s->m_Smat_w_is[i]->cwSet(0.);
    *(m_s->m_Smat_w_is[i]) = (1./input_2lambdaWVec[i]) * *(m_s->m_Rmat_w_is[i]);
    for (unsigned int j = 0; j < m_s->m_Smat_w_is[i]->numRowsLocal(); ++j) {
      (*(m_s->m_Smat_w_is[i]))(j,j) +=  1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
    }
  }
  m_s->m_Smat_w.cwSet(0.);
  m_s->m_Smat_w.fillWithBlocksDiagonally(0,0,m_s->m_Smat_w_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_w'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_s->m_Smat_w.subWriteContents("Sigma_w",
                                     "mat_Sigma_w",
                                     "m",
                                     tmpSet);
    }
  }

  this->memoryCheck(93);

  initialPos = 0;
  for (unsigned int i = 0; i < m_j->m_Smat_uw_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_j->m_Rmat_uw_is[i]->cwSet(0.);
    this->fillR_formula1_for_Sigma_uw(m_e->m_paper_xs_standard,
                                      input_8thetaVec,
                                      m_s->m_paper_xs_asterisks_standard,
                                      m_s->m_paper_ts_asterisks_standard,
                                      m_s->m_tmp_rho_w_vec,*(m_j->m_Rmat_uw_is[i]),
                                      outerCounter);

    m_j->m_Smat_uw_is[i]->cwSet(0.);
    *(m_j->m_Smat_uw_is[i]) = (1./input_2lambdaWVec[i]) * *(m_j->m_Rmat_uw_is[i]);
  }
  m_j->m_Smat_uw.cwSet(0.);
  m_j->m_Smat_uw.fillWithBlocksDiagonally(0,0,m_j->m_Smat_uw_is,true,true);
  m_j->m_Smat_uw_t.fillWithTranspose(0,0,m_j->m_Smat_uw,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_j->m_Smat_uw'"
                            << std::endl;
  }

  if (outerCounter == 1) {
    if (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end()) {
      m_j->m_Smat_uw.subWriteContents("Sigma_uw",
                                      "mat_Sigma_uw",
                                      "m",
                                      tmpSet);
    }
  }
#endif
  this->memoryCheck(94);

  //********************************************************************************
  // Form '\Sigma_z' matrix
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << ": m_tmp_Smat_z.numRowsLocal() = " << m_z->m_tmp_Smat_z.numRowsLocal()
                            << ", m_tmp_Smat_z.numCols() = "      << m_z->m_tmp_Smat_z.numCols()
                            << ", m_Smat_v.numRowsLocal() = "     << m_e->m_Smat_v.numRowsLocal()
                            << ", m_Smat_v.numCols() = "          << m_e->m_Smat_v.numCols()
                            << ", m_Smat_u.numRowsLocal() = "     << m_j->m_Smat_u.numRowsLocal()
                            << ", m_Smat_u.numCols() = "          << m_j->m_Smat_u.numCols()
                            << ", m_Smat_w.numRowsLocal() = "     << m_s->m_Smat_w.numRowsLocal()
                            << ", m_Smat_w.numCols() = "          << m_s->m_Smat_w.numCols()
                            << ", m_Smat_uw.numRowsLocal() = "    << m_j->m_Smat_uw.numRowsLocal()
                            << ", m_Smat_uw.numCols() = "         << m_j->m_Smat_uw.numCols()
                            << ", m_Smat_v_i_spaces.size() = "    << m_e->m_Smat_v_i_spaces.size()
                            << ", m_Smat_u_is.size() = "          << m_j->m_Smat_u_is.size()
                            << ", m_Smat_w_is.size() = "          << m_s->m_Smat_w_is.size()
                            << std::endl;
  }

  this->memoryCheck(95);

  m_z->m_tmp_Smat_z.cwSet(0.);
  if (m_allOutputsAreScalar) {
    // ppp
  }
  else {
    m_z->m_tmp_Smat_z.cwSet(0,0,m_e->m_Smat_v);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal(),                             m_e->m_Smat_v.numCols(),                        m_j->m_Smat_u);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal(),                             m_e->m_Smat_v.numCols()+m_j->m_Smat_u.numCols(),m_j->m_Smat_uw);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal()+m_j->m_Smat_u.numRowsLocal(),m_e->m_Smat_v.numCols(),                        m_j->m_Smat_uw_t);
    m_z->m_tmp_Smat_z.cwSet(m_e->m_Smat_v.numRowsLocal()+m_j->m_Smat_u.numRowsLocal(),m_e->m_Smat_v.numCols()+m_j->m_Smat_u.numCols(),m_s->m_Smat_w);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(2)"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_w_hat(
  const P_V&         input_1lambdaEtaVec,
  const P_V&         input_2lambdaWVec,
  const P_V&         input_3rhoWVec,
  const P_V&         input_4lambdaSVec,
  const P_V&         input_8thetaVec,
        unsigned int outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_w_hat()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  unsigned int initialPos = 0;
  for (unsigned int i = 0; i < m_s->m_Smat_w_is.size(); ++i) {
    input_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
    initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
    m_s->m_Rmat_w_is[i]->cwSet(0.); // This matrix is square: m_paper_m X m_paper_m
    this->fillR_formula1_for_Sigma_w(m_s->m_paper_xs_asterisks_standard,
                                     m_s->m_paper_ts_asterisks_standard,
                                     m_s->m_tmp_rho_w_vec,
                                     *(m_s->m_Rmat_w_is[i]),
                                     outerCounter);

    m_s->m_Smat_w_is[i]->cwSet(0.);
    *(m_s->m_Smat_w_is[i]) = (1./input_2lambdaWVec[i]) * *(m_s->m_Rmat_w_is[i]);
    for (unsigned int j = 0; j < m_s->m_Smat_w_is[i]->numRowsLocal(); ++j) {
      (*(m_s->m_Smat_w_is[i]))(j,j) +=  1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
    }
  }
  m_s->m_Smat_w.cwSet(0.);
  m_s->m_Smat_w.fillWithBlocksDiagonally(0,0,m_s->m_Smat_w_is,true,true);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_w_hat()"
                            << ", outerCounter = " << outerCounter
                            << ": finished instantiating 'm_Smat_w'"
                            << std::endl;
  }

  if (m_allOutputsAreScalar) {
    // ppp
  }
  else {
    m_s->m_Smat_w_hat = m_s->m_Smat_w + (1./input_1lambdaEtaVec[0]) * (*m_s->m_Kt_K_inv);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_w_hat()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

}  // End namespace QUESO

#endif // __UQ_GCM_4_H__
