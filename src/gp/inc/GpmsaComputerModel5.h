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

#ifndef __UQ_GCM_5_H__
#define __UQ_GCM_5_H__

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v(
  const std::vector<const S_V* >& xVecs,
  const P_V&                      rho_v_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(xVecs.size() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()",
                      "xVecs.size() is wrong");
  for (unsigned int i = 0; i < xVecs.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()",
                        "xVecs[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO((rho_v_vec.sizeLocal() == 0) || (rho_v_vec.sizeLocal() > m_s->m_paper_p_x), // Should be equal to m_Gs[i], for some 'i'
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()",
                      "rho_v_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()",
                      "Rmat.numCols() is wrong");

  S_V vecI(*(xVecs[0]));  
  S_V vecJ(*(xVecs[0]));  
  for (unsigned int i = 0; i < m_e->m_paper_n; ++i) {
    vecI = *(xVecs[i]);
    for (unsigned int j = 0; j < m_e->m_paper_n; ++j) {
      vecJ = *(xVecs[j]);
      Rmat(i,j) = 1.;
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = vecI[k] - vecJ[k];
        Rmat(i,j) *= std::pow(rho_v_vec[k],4.*diffTerm*diffTerm);
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
          *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()"
                                  << ", outerCounter = " << outerCounter
                                  << ": i = "            << i
                                  << ", j = "            << j
                                  << ", k = "            << k
                                  << ", vecI[k] = "      << vecI[k]
                                  << ", vecJ[k] = "      << vecJ[k]
                                  << ", diffTerm = "     << diffTerm
                                  << ", rho_v_vec[k] = " << rho_v_vec[k]
                                  << std::endl;
        }
      }
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()"
                                << ", outerCounter = " << outerCounter
                                << ": i = "            << i
                                << ", j = "            << j
                                << ", Rmat(i,j) = "    << Rmat(i,j)
                                << std::endl;
      }
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u(
  const std::vector<const S_V* >& xVecs,
  const P_V&                      tVec,
  const P_V&                      rho_w_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(xVecs.size() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                      "xVecs.size() is wrong");
  for (unsigned int i = 0; i < xVecs.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                        "xVecs[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVec.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                      "tVec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(rho_w_vec.sizeLocal() != (m_s->m_paper_p_x+m_s->m_paper_p_t),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                      "rho_w_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                      "Rmat.numCols() is wrong");

  S_V vecI(*(xVecs[0]));
  S_V vecJ(*(xVecs[0]));
  for (unsigned int i = 0; i < m_e->m_paper_n; ++i) {
    vecI = *(xVecs[i]);
    for (unsigned int j = 0; j < m_e->m_paper_n; ++j) {
      vecJ = *(xVecs[j]);
      Rmat(i,j) = 1.;
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) { // Yes, just 'p_x', instead of 'p_x + p_t', since 't' is the same for all pairs
        double diffTerm = vecI[k] - vecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[k],4.*diffTerm*diffTerm);
      }
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w(
  const std::vector<const S_V* >& xVecs,
  const std::vector<const P_V* >& tVecs,
  const P_V&                      rho_w_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(xVecs.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                      "xVecs.size() is wrong");
  for (unsigned int i = 0; i < xVecs.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                        "xVecs[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVecs.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                      "tVecs.size() is wrong");
  for (unsigned int i = 0; i < tVecs.size(); ++i) {
    UQ_FATAL_TEST_MACRO(tVecs[i]->sizeLocal() != m_s->m_paper_p_t,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                        "tVecs[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(rho_w_vec.sizeLocal() != (m_s->m_paper_p_x+m_s->m_paper_p_t),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                      "rho_w_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()",
                      "Rmat.numCols() is wrong");

  S_V xVecI(*(xVecs[0]));
  S_V xVecJ(*(xVecs[0]));
  P_V tVecI(*(tVecs[0]));
  P_V tVecJ(*(tVecs[0]));
  for (unsigned int i = 0; i < m_s->m_paper_m; ++i) {
    xVecI = *(xVecs[i]);
    tVecI = *(tVecs[i]);
    for (unsigned int j = 0; j < m_s->m_paper_m; ++j) {
      xVecJ = *(xVecs[j]);
      tVecJ = *(tVecs[j]);
      Rmat(i,j) = 1.;
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = xVecI[k] - xVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[k],4.*diffTerm*diffTerm);
      }
      for (unsigned int k = 0; k < m_s->m_paper_p_t; ++k) {
        double diffTerm = tVecI[k] - tVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[m_s->m_paper_p_x+k],4.*diffTerm*diffTerm);
      }
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw(
  const std::vector<const S_V* >& xVecs1,
  const P_V&                      tVec1,
  const std::vector<const S_V* >& xVecs2,
  const std::vector<const P_V* >& tVecs2,
  const P_V&                      rho_w_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(xVecs1.size() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                      "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs1[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                        "xVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVec1.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u()",
                      "tVec1.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(xVecs2.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                      "xVecs2.size() is wrong");
  for (unsigned int i = 0; i < xVecs2.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs2[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                        "xVecs2[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVecs2.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                      "tVecs2.size() is wrong");
  for (unsigned int i = 0; i < tVecs2.size(); ++i) {
    UQ_FATAL_TEST_MACRO(tVecs2[i]->sizeLocal() != m_s->m_paper_p_t,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                        "tVecs2[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(rho_w_vec.sizeLocal() != (m_s->m_paper_p_x+m_s->m_paper_p_t),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                      "rho_w_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()",
                      "Rmat.numCols() is wrong");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  S_V xVecI(*(xVecs1[0]));
  S_V xVecJ(*(xVecs2[0]));
  P_V tVecI(tVec1);
  P_V tVecJ(*(tVecs2[0]));
  for (unsigned int i = 0; i < m_e->m_paper_n; ++i) {
    xVecI = *(xVecs1[i]);
    tVecI = tVec1;
    for (unsigned int j = 0; j < m_s->m_paper_m; ++j) {
      xVecJ = *(xVecs2[j]);
      tVecJ = *(tVecs2[j]);
      Rmat(i,j) = 1.;
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = xVecI[k] - xVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[k],4.*diffTerm*diffTerm);
      }
      for (unsigned int k = 0; k < m_s->m_paper_p_t; ++k) {
        double diffTerm = tVecI[k] - tVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[m_s->m_paper_p_x+k],4.*diffTerm*diffTerm);
      }
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_uw()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk(
  const std::vector<const S_V* >& xVecs1,
  const std::vector<const P_V* >& tVecs1,
  const S_V&                      xVec2,
  const P_V&                      tVec2,
  const P_V&                      rho_v_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(xVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs1[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                        "xVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(tVecs1[i]->sizeLocal() != m_s->m_paper_p_t,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                        "tVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(xVec2.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "xVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(tVec2.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "tVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO((rho_v_vec.sizeLocal() == 0) || (rho_v_vec.sizeLocal() > m_s->m_paper_p_x), // Should be equal to m_Gs[i], for some 'i'
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "rho_v_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_e->m_paper_n,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != 1,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()",
                      "Rmat.numCols() is wrong");

  S_V xVecI(*(xVecs1[0]));
  S_V xVecJ(xVec2);
  P_V tVecI(*(tVecs1[0]));
  P_V tVecJ(tVec2);
  for (unsigned int i = 0; i < m_e->m_paper_n; ++i) {
    xVecI = *(xVecs1[i]);
    tVecI = *(tVecs1[i]);
    for (unsigned int j = 0; j < 1; ++j) { // Yes, only '1'
      xVecJ = xVec2;
      tVecJ = tVec2;
      Rmat(i,j) = 1.;
      // checar
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = xVecI[k] - xVecJ[k];
        Rmat(i,j) *= std::pow(rho_v_vec[k],4.*diffTerm*diffTerm);
      }
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula2_for_Sigma_v_hat_v_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk(
  const std::vector<const S_V* >& xVecs1,
  const std::vector<const P_V* >& tVecs1,
  const S_V&                      xVec2,
  const P_V&                      tVec2,
  const P_V&                      rho_w_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }
  UQ_FATAL_TEST_MACRO(xVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs1[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                        "xVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(tVecs1[i]->sizeLocal() != m_s->m_paper_p_t,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                        "tVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(xVec2.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "xVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(tVec2.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "tVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(rho_w_vec.sizeLocal() != (m_s->m_paper_p_x+m_s->m_paper_p_t),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "rho_w_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != 1,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()",
                      "Rmat.numCols() is wrong");

  S_V xVecI(*(xVecs1[0]));
  S_V xVecJ(xVec2);
  P_V tVecI(*(tVecs1[0]));
  P_V tVecJ(tVec2);
  for (unsigned int i = 0; i < m_s->m_paper_m; ++i) {
    xVecI = *(xVecs1[i]);
    tVecI = *(tVecs1[i]);
    for (unsigned int j = 0; j < 1; ++j) { // Yes, only '1'
      xVecJ = xVec2;
      tVecJ = tVec2;
      Rmat(i,j) = 1.;
      // checar
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = xVecI[k] - xVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[k],4.*diffTerm*diffTerm);
      }
      for (unsigned int k = 0; k < m_s->m_paper_p_t; ++k) {
        double diffTerm = tVecI[k] - tVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[m_s->m_paper_p_x+k],4.*diffTerm*diffTerm);
      }
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_u_hat_u_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk(
  const std::vector<const S_V* >& xVecs1,
  const std::vector<const P_V* >& tVecs1,
  const S_V&                      xVec2,
  const P_V&                      tVec2,
  const P_V&                      rho_w_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }
  UQ_FATAL_TEST_MACRO(xVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs1[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                        "xVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(tVecs1[i]->sizeLocal() != m_s->m_paper_p_t,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                        "tVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(xVec2.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "xVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(tVec2.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "tVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(rho_w_vec.sizeLocal() != (m_s->m_paper_p_x+m_s->m_paper_p_t),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "rho_w_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != 1,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()",
                      "Rmat.numCols() is wrong");

  S_V xVecI(*(xVecs1[0]));
  S_V xVecJ(xVec2);
  P_V tVecI(*(tVecs1[0]));
  P_V tVecJ(tVec2);
  for (unsigned int i = 0; i < m_s->m_paper_m; ++i) {
    xVecI = *(xVecs1[i]);
    tVecI = *(tVecs1[i]);
    for (unsigned int j = 0; j < 1; ++j) { // Yes, only '1'
      xVecJ = xVec2;
      tVecJ = tVec2;
      Rmat(i,j) = 1.;
      // checar
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = xVecI[k] - xVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[k],4.*diffTerm*diffTerm);
      }
      for (unsigned int k = 0; k < m_s->m_paper_p_t; ++k) {
        double diffTerm = tVecI[k] - tVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[m_s->m_paper_p_x+k],4.*diffTerm*diffTerm);
      }
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_u_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk(
  const std::vector<const S_V* >& xVecs1,
  const std::vector<const P_V* >& tVecs1,
  const S_V&                      xVec2,
  const P_V&                      tVec2,
  const P_V&                      rho_w_vec,
        D_M&                      Rmat,
        unsigned int              outerCounter)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(xVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(xVecs1[i]->sizeLocal() != m_s->m_paper_p_x,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                        "xVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(tVecs1.size() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    UQ_FATAL_TEST_MACRO(tVecs1[i]->sizeLocal() != m_s->m_paper_p_t,
                        m_env.worldRank(),
                        "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                        "tVecs1[i]->sizeLocal() is wrong");
  }
  UQ_FATAL_TEST_MACRO(xVec2.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "xVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(tVec2.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "tVec2.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(rho_w_vec.sizeLocal() != (m_s->m_paper_p_x+m_s->m_paper_p_t),
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "rho_w_vec.sizeLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numRowsLocal() != m_s->m_paper_m,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "Rmat.numRowsLocal() is wrong");
  UQ_FATAL_TEST_MACRO(Rmat.numCols() != 1,
                      m_env.worldRank(),
                      "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()",
                      "Rmat.numCols() is wrong");

  S_V xVecI(*(xVecs1[0]));
  S_V xVecJ(xVec2);
  P_V tVecI(*(tVecs1[0]));
  P_V tVecJ(tVec2);
  for (unsigned int i = 0; i < m_s->m_paper_m; ++i) {
    xVecI = *(xVecs1[i]);
    tVecI = *(tVecs1[i]);
    for (unsigned int j = 0; j < 1; ++j) { // Yes, only '1'
      xVecJ = xVec2;
      tVecJ = tVec2;
      Rmat(i,j) = 1.;
      for (unsigned int k = 0; k < m_s->m_paper_p_x; ++k) {
        double diffTerm = xVecI[k] - xVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[k],4.*diffTerm*diffTerm);
      }
      for (unsigned int k = 0; k < m_s->m_paper_p_t; ++k) {
        double diffTerm = tVecI[k] - tVecJ[k];
        Rmat(i,j) *= std::pow(rho_w_vec[m_s->m_paper_p_x+k],4.*diffTerm*diffTerm);
      }
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::fillR_formula1_for_Sigma_w_hat_w_asterisk()"
                            << ", outerCounter = " << outerCounter
                            << std::endl;
  }

  return;
}

}  // End namespace QUESO

#endif // __UQ_GCM_5_H__
