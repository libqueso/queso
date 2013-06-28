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

#ifndef __UQ_GCM_3_H__
#define __UQ_GCM_3_H__

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSpaceClass<P_V,P_M>&
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalSpace() const
{
  UQ_FATAL_TEST_MACRO(m_t == NULL,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalSpace()",
                      "m_t is NULL");
  return m_t->m_totalSpace;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSpaceClass<P_V,P_M>&
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::unique_vu_space() const
{
  UQ_FATAL_TEST_MACRO(m_j == NULL,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::unique_vu_space()",
                      "m_j is NULL");
  return m_j->m_unique_vu_space;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const uqBaseVectorRVClass<P_V,P_M>& 
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalPriorRv() const
{
  UQ_FATAL_TEST_MACRO(m_t == NULL,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalPriorRv()",
                      "m_t is NULL");
  return m_t->m_totalPriorRv;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const uqGenericVectorRVClass<P_V,P_M>& 
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalPostRv() const
{
  UQ_FATAL_TEST_MACRO(m_t == NULL,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalPostRv()",
                      "m_t is NULL");
  return m_t->m_totalPostRv;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::memoryCheck(unsigned int codePositionId)
{
#if 0
  std::cout << "Entering memoryCheck(), m_like_counter = " << m_like_counter << ", codePositionId = " << codePositionId << std::endl;

  double sumZ = 0.;
  for (unsigned int i = 0; i < m_z->m_tmp_Smat_z.numRowsLocal(); ++i) {
    //std::cout << "i = " << i << std::endl;
    for (unsigned int j = 0; j < m_z->m_tmp_Smat_z.numCols(); ++j) {
      //std::cout << "j = " << j << std::endl;
      sumZ += m_z->m_tmp_Smat_z(i,j);
    }
  }
  //std::cout << "Aqui 000-000, sumZ = " << sumZ << std::endl;

  double sumV = 0.;
  for (unsigned int i = 0; i < m_e->m_Smat_v.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < m_e->m_Smat_v.numCols(); ++j) {
      sumV += m_e->m_Smat_v(i,j);
    }
  }

  double sumU = 0.;
  for (unsigned int i = 0; i < m_j->m_Smat_u.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < m_j->m_Smat_u.numCols(); ++j) {
      sumU += m_j->m_Smat_u(i,j);
    }
  }

  double sumW = 0.;
  for (unsigned int i = 0; i < m_s->m_Smat_w.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < m_s->m_Smat_w.numCols(); ++j) {
      sumW += m_s->m_Smat_w(i,j);
    }
  }

  double sumUW = 0.;
  for (unsigned int i = 0; i < m_j->m_Smat_uw.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < m_j->m_Smat_uw.numCols(); ++j) {
      sumUW += m_j->m_Smat_uw(i,j);
    }
  }

  double sumUW_T = 0.;
  for (unsigned int i = 0; i < m_j->m_Smat_uw_t.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < m_j->m_Smat_uw_t.numCols(); ++j) {
      sumUW_T += m_j->m_Smat_uw_t(i,j);
    }
  }

  std::cout << "Leaving memoryCheck(), m_like_counter = " << m_like_counter << ", codePositionId = " << codePositionId << std::endl;
#endif
  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::generatePriorSeq()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::generatePriorSeq()..."
                            << ": m_optionsObj->m_prefix.c_str() = "          << m_optionsObj->m_prefix.c_str()
                            << ", m_optionsObj->m_ov.m_priorSeqNumSamples = " << m_optionsObj->m_ov.m_priorSeqNumSamples
                            << std::endl;
  }

  uqSequenceOfVectorsClass<P_V,P_M> priorSeq(m_t->m_totalSpace,m_optionsObj->m_ov.m_priorSeqNumSamples,m_optionsObj->m_prefix+"priorSeq");
  P_V totalSample(m_t->m_totalSpace.zeroVector());
  for (unsigned int sampleId = 0; sampleId < m_optionsObj->m_ov.m_priorSeqNumSamples; ++sampleId) {
    m_t->m_totalPriorRv.realizer().realization(totalSample);
    priorSeq.setPositionValues(sampleId,totalSample);
  }
  priorSeq.unifiedWriteContents(m_optionsObj->m_ov.m_priorSeqDataOutputFileName,
                                m_optionsObj->m_ov.m_priorSeqDataOutputFileType);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::generatePriorSeq()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
double
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::staticLikelihoodRoutine(
  const P_V&  totalValues,
  const P_V*  totalDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  return ((uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>*) functionDataPtr)->likelihoodRoutine(totalValues,
                                                                                                            totalDirection,
                                                                                                            functionDataPtr,
                                                                                                            gradVector,
                                                                                                            hessianMatrix,
                                                                                                            hessianEffect);
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
double
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(
  const P_V&  totalValues,
  const P_V*  totalDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  m_like_counter++;
  //std::cout << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(), m_like_counter = " << m_like_counter << std::endl;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()..."
                            << ": m_like_counter = "            << m_like_counter
                            << ", totalValues = "               << totalValues
                            << ", m_env.subComm().NumProc() = " << m_env.subComm().NumProc()
                            << ", my subRank = "                << m_env.subRank()
                            << ", m_formCMatrix = "             << m_formCMatrix
                            << ", m_cMatIsRankDefficient = "    << m_cMatIsRankDefficient
                            << std::endl;
  }

  double lnLikelihoodValue = 0.;

  if (totalDirection  &&
      functionDataPtr &&
      gradVector      &&
      hessianMatrix   &&
      hessianEffect) {
    // Just to eliminate INTEL compiler warnings
  }

  //********************************************************************************
  // Total values = (\lambda_eta[1],\lambda_w[p_eta],\rho_w[(p_x+p_t).p_eta],\lambda_y[1],\lambda_v[F],\rho_v[F.p_x],\theta)
  //********************************************************************************
  UQ_FATAL_TEST_MACRO((m_s->m_1lambdaEtaSpace.dimLocal() != 1                                                       ) ||
                      (m_s->m_2lambdaWSpace.dimLocal()   != m_s->m_paper_p_eta                                      ) ||
                      (m_s->m_3rhoWSpace.dimLocal()      != (m_s->m_paper_p_eta*(m_s->m_paper_p_x+m_s->m_paper_p_t))) ||
                      (m_s->m_4lambdaSSpace.dimLocal()   != m_s->m_paper_p_eta                                      ),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()",
                      "inconsistent 'm_s' space dimensions");
  if (m_thereIsExperimentalData) {
    UQ_FATAL_TEST_MACRO((m_e->m_5lambdaYSpace.dimLocal() != 1                                ) ||
                        (m_e->m_6lambdaVSpace.dimLocal() != m_e->m_paper_F                   ) ||
                        (m_e->m_7rhoVSpace.dimLocal()    != (m_e->m_paper_F*m_s->m_paper_p_x)),
                        m_env.worldRank(),
                        "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()",
                        "inconsistent 'm_e' space dimensions");
  }

  this->memoryCheck(50);

  unsigned int currPosition = 0;
  totalValues.cwExtract(currPosition,m_s->m_tmp_1lambdaEtaVec); // Total of '1' in paper
  currPosition += m_s->m_tmp_1lambdaEtaVec.sizeLocal();
  totalValues.cwExtract(currPosition,m_s->m_tmp_2lambdaWVec);   // Total of 'p_eta' in paper
  currPosition += m_s->m_tmp_2lambdaWVec.sizeLocal();
  totalValues.cwExtract(currPosition,m_s->m_tmp_3rhoWVec);      // Total of 'p_eta*(p_x+p_t)' in paper
  currPosition += m_s->m_tmp_3rhoWVec.sizeLocal();
  totalValues.cwExtract(currPosition,m_s->m_tmp_4lambdaSVec);   // Total of 'p_eta' in matlab code
  currPosition += m_s->m_tmp_4lambdaSVec.sizeLocal();

  if (m_thereIsExperimentalData) {
    totalValues.cwExtract(currPosition,m_e->m_tmp_5lambdaYVec);   // Total of '1' in paper
    currPosition += m_e->m_tmp_5lambdaYVec.sizeLocal();
    totalValues.cwExtract(currPosition,m_e->m_tmp_6lambdaVVec);   // Total of 'F' in paper
    currPosition += m_e->m_tmp_6lambdaVVec.sizeLocal();
    totalValues.cwExtract(currPosition,m_e->m_tmp_7rhoVVec);      // Total of 'F*p_x' in paper
    currPosition += m_e->m_tmp_7rhoVVec.sizeLocal();
    totalValues.cwExtract(currPosition,m_e->m_tmp_8thetaVec);     // Application specific
    currPosition += m_e->m_tmp_8thetaVec.sizeLocal();
  }
  UQ_FATAL_TEST_MACRO(currPosition != totalValues.sizeLocal(),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()",
                      "'currPosition' and 'totalValues.sizeLocal()' should be equal");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
                            << ", m_like_counter = " << m_like_counter
                            << ": finished extracting components from 'totalValues'"
                            << ", m_tmp_1lambdaEtaVec = " << m_s->m_tmp_1lambdaEtaVec
                            << ", m_tmp_2lambdaWVec = "   << m_s->m_tmp_2lambdaWVec
                            << ", m_tmp_3rhoWVec = "      << m_s->m_tmp_3rhoWVec
                            << ", m_tmp_4lambdaSVec = "   << m_s->m_tmp_4lambdaSVec;
    if (m_thereIsExperimentalData) {
      *m_env.subDisplayFile() << ", m_tmp_5lambdaYVec = " << m_e->m_tmp_5lambdaYVec
                              << ", m_tmp_6lambdaVVec = " << m_e->m_tmp_6lambdaVVec
                              << ", m_tmp_7rhoVVec = "    << m_e->m_tmp_7rhoVVec
                              << ", m_tmp_8thetaVec = "   << m_e->m_tmp_8thetaVec;
    }
    *m_env.subDisplayFile() << std::endl;
  }

  //********************************************************************************
  // Check if current 'totalValues' has any common values with previous 'totalValues' (todo)
  //********************************************************************************
  bool here_1_repeats = false;
  bool here_2_repeats = false;
  bool here_3_repeats = false;
  bool here_4_repeats = false;
  bool here_5_repeats = false;
  bool here_6_repeats = false;
  bool here_7_repeats = false;
  bool here_8_repeats = false;
  if ((m_optionsObj->m_ov.m_checkAgainstPreviousSample) &&
      (m_like_counter == 1                            )) {
    here_1_repeats = (m_s->m_like_previous1 == m_s->m_tmp_1lambdaEtaVec);
    here_2_repeats = (m_s->m_like_previous2 == m_s->m_tmp_2lambdaWVec);
    here_3_repeats = (m_s->m_like_previous3 == m_s->m_tmp_3rhoWVec);
    here_4_repeats = (m_s->m_like_previous2 == m_s->m_tmp_4lambdaSVec);
    if (m_thereIsExperimentalData) {
      here_5_repeats = (m_e->m_like_previous5 == m_e->m_tmp_5lambdaYVec);
      here_6_repeats = (m_e->m_like_previous6 == m_e->m_tmp_6lambdaVVec);
      here_7_repeats = (m_e->m_like_previous7 == m_e->m_tmp_7rhoVVec);
      here_8_repeats = (m_e->m_like_previous8 == m_e->m_tmp_8thetaVec);
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
                              << ", m_like_counter = "     << m_like_counter
                              << "\n  m_like_previous1 = " << m_s->m_like_previous1
                              << "\n  m_like_previous2 = " << m_s->m_like_previous2
                              << "\n  m_like_previous3 = " << m_s->m_like_previous3;
      if (m_thereIsExperimentalData) {
        *m_env.subDisplayFile() << "\n  m_like_previous5 = " << m_e->m_like_previous5
                                << "\n  m_like_previous6 = " << m_e->m_like_previous6
                                << "\n  m_like_previous7 = " << m_e->m_like_previous7
                                << "\n  m_like_previous8 = " << m_e->m_like_previous8;
      }
      *m_env.subDisplayFile() << std::endl;
    }
    if (here_1_repeats ||
        here_2_repeats ||
        here_3_repeats ||
        here_4_repeats ||
        here_5_repeats ||
        here_6_repeats ||
        here_7_repeats ||
	here_8_repeats) {}; // just to remove compiler warning
  }

  lnLikelihoodValue = 0.;
  if ((m_formCMatrix         ) &&
      (m_cMatIsRankDefficient)) {
    //********************************************************************************
    // 'm_Cmat' is rank defficient
    //********************************************************************************
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                              << ", m_like_counter = " << m_like_counter
                              << ": going through true 'm_cMatIsRankDefficient' case"
                              << std::endl;
    }

    this->memoryCheck(57);

    if (m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC) {
      //********************************************************************************
      // Compute '\Sigma_z_tilde_hat' matrix
      //********************************************************************************
      // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
      // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
      // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
      // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
      // Then fill m_tmp_Smat_z
      // Fill m_Rmat_extra
      // Then fill m_tmp_Smat_z_tilde_hat
      this->formSigma_z_tilde_hat(m_s->m_tmp_1lambdaEtaVec,
                                  m_s->m_tmp_2lambdaWVec,
                                  m_s->m_tmp_3rhoWVec,
                                  m_s->m_tmp_4lambdaSVec,
                                  m_e->m_tmp_5lambdaYVec,
                                  m_e->m_tmp_6lambdaVVec,
                                  m_e->m_tmp_7rhoVVec,
                                  m_e->m_tmp_8thetaVec,
                                  m_like_counter);

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                               << m_like_counter
                                << ": finished computing 'm_tmp_Smat_z_tilde_hat' =\n" << m_zt->m_tmp_Smat_z_tilde_hat
                                << std::endl;
      }

      this->memoryCheck(58);

      //********************************************************************************
      // Compute the determinant of '\Sigma_z_tilde_hat' matrix
      //********************************************************************************
      double Smat_z_tilde_hat_lnDeterminant = m_zt->m_tmp_Smat_z_tilde_hat.lnDeterminant();

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                                               << m_like_counter
                                << ": finished computing 'm_tmp_Smat_z_tilde_hat->lnDeterminant()' = " << Smat_z_tilde_hat_lnDeterminant
                                << std::endl;
      }
      lnLikelihoodValue += -0.5*Smat_z_tilde_hat_lnDeterminant;

      this->memoryCheck(59);

      //********************************************************************************
      // Compute Gaussian contribution
      //********************************************************************************
      double tmpValue1 = scalarProduct(m_zt->m_Zvec_tilde_hat,m_zt->m_tmp_Smat_z_tilde_hat.invertMultiply(m_zt->m_Zvec_tilde_hat)); // inversion savings
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue1 = " << tmpValue1
                                << std::endl;
      }
      lnLikelihoodValue += -0.5*tmpValue1;

      this->memoryCheck(60);

      //********************************************************************************
      // Include effect of exponent modifiers
      //********************************************************************************
      double tmpValue2 = m_st->m_a_eta_modifier_tilde*std::log(m_s->m_tmp_1lambdaEtaVec[0]) - m_st->m_b_eta_modifier_tilde*m_s->m_tmp_1lambdaEtaVec[0];
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue2 = " << tmpValue2
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue2;

      this->memoryCheck(61);

      double tmpValue3 = m_jt->m_a_y_modifier_tilde*std::log(m_e->m_tmp_5lambdaYVec[0]) - m_jt->m_b_y_modifier_tilde*m_e->m_tmp_5lambdaYVec[0];
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue3 = " << tmpValue3
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue3;

      this->memoryCheck(62);
    }
    else { // if (m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC)
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)",
                          "incomplete code for situation 'm_useTildeLogicForRankDefficientC == false'");
    }
  }
  else { // if (m_formCMatrix) && (m_cMatIsRankDefficient)
    //********************************************************************************
    // 'm_Cmat' (i) does not exist or (ii) is full rank
    //********************************************************************************
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                              << ", m_like_counter = " << m_like_counter
                              << ": going through case where C matrix (i) does not exist or (ii) is full rank"
                              << ", m_thereIsExperimentalData = " << m_thereIsExperimentalData
                              << std::endl;
    }

    this->memoryCheck(51);

    //********************************************************************************
    // Compute '\Sigma_z_hat' matrix
    //********************************************************************************
    // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
    // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
    // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
    // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
    // Then fill m_tmp_Smat_z
    // Fill m_Rmat_extra
    // Then fill m_tmp_Smat_z_hat
    
    if (m_thereIsExperimentalData) {
      this->formSigma_z_hat(m_s->m_tmp_1lambdaEtaVec,
                            m_s->m_tmp_2lambdaWVec,
                            m_s->m_tmp_3rhoWVec,
                            m_s->m_tmp_4lambdaSVec,
                            m_e->m_tmp_5lambdaYVec,
                            m_e->m_tmp_6lambdaVVec,
                            m_e->m_tmp_7rhoVVec,
                            m_e->m_tmp_8thetaVec,
                            m_like_counter);
    }
    else {
      UQ_FATAL_TEST_MACRO(true, // (m_thereIsExperimentalData == false)
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)",
                          "incomplete code for situation 'm_thereIsExperimentalData == false'");

      this->formSigma_z_hat(m_s->m_tmp_1lambdaEtaVec,
                            m_s->m_tmp_2lambdaWVec,
                            m_s->m_tmp_3rhoWVec,
                            m_s->m_tmp_4lambdaSVec,
                            m_like_counter);
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                              << ", m_like_counter = "                     << m_like_counter
                              << ": finished computing 'm_tmp_Smat_z_hat' =\n" << m_z->m_tmp_Smat_z_hat
                              << std::endl;
    }

    this->memoryCheck(52);

    //********************************************************************************
    // Compute the determinant of '\Sigma_z_hat' matrix
    //********************************************************************************
    double Smat_z_hat_lnDeterminant = m_z->m_tmp_Smat_z_hat.lnDeterminant();

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                              << ", m_like_counter = "                                        << m_like_counter
                              << ": finished computing 'm_tmp_Smat_z_hat.lnDeterminant()' = " << Smat_z_hat_lnDeterminant
                              << std::endl;
    }
    lnLikelihoodValue += -0.5*Smat_z_hat_lnDeterminant;

    this->memoryCheck(53);

    //********************************************************************************
    // Compute Gaussian contribution
    //********************************************************************************
    double tmpValue1 = scalarProduct(m_z->m_Zvec_hat,m_z->m_tmp_Smat_z_hat.invertMultiply(m_z->m_Zvec_hat)); // inversion savings
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                              << ", m_like_counter = "                << m_like_counter
                              << ": finished computing 'tmpValue1 = " << tmpValue1
                              << std::endl;
    }
    lnLikelihoodValue += -0.5*tmpValue1;

    this->memoryCheck(54);

    if (m_allOutputsAreScalar) {
      // Do nothing
    }
    else {
      //********************************************************************************
      // Include effect of exponent modifiers
      //********************************************************************************
      double tmpValue2 = m_s->m_a_eta_modifier*std::log(m_s->m_tmp_1lambdaEtaVec[0]) - m_s->m_b_eta_modifier*m_s->m_tmp_1lambdaEtaVec[0];
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue2 = " << tmpValue2
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue2;

      this->memoryCheck(55);

      double tmpValue3 = m_j->m_a_y_modifier*std::log(m_e->m_tmp_5lambdaYVec[0]) - m_j->m_b_y_modifier*m_e->m_tmp_5lambdaYVec[0];
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue3 = " << tmpValue3
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue3;

      this->memoryCheck(56);

      //for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
      //  for (unsigned int j = 0; j < (m_paper_p_x + m_paper_p_t); ++j) {
      //  }
      //}
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
                            << ", m_like_counter = " << m_like_counter
                            << ": finished computing ln(likelihood)"
                            << ", lnLikelihoodValue = " << lnLikelihoodValue
                            << std::endl;
  }

  this->memoryCheck(63);

  //******************************************************************************
  // Prepare to return
  //******************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
                            << ", m_like_counter = " << m_like_counter
                            << ": starting saving current samples as previous"
                            << std::endl;
  }

  m_t->m_like_previousTotal = totalValues;
  m_s->m_like_previous1 = m_s->m_tmp_1lambdaEtaVec;
  m_s->m_like_previous2 = m_s->m_tmp_2lambdaWVec;
  m_s->m_like_previous3 = m_s->m_tmp_3rhoWVec;
  m_s->m_like_previous2 = m_s->m_tmp_4lambdaSVec;
  m_e->m_like_previous5 = m_e->m_tmp_5lambdaYVec;
  m_e->m_like_previous6 = m_e->m_tmp_6lambdaVVec;
  m_e->m_like_previous7 = m_e->m_tmp_7rhoVVec;
  m_e->m_like_previous8 = m_e->m_tmp_8thetaVec;

  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);

  //std::cout << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(), m_like_counter = " << m_like_counter << std::endl;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
                            << ": m_like_counter = "    << m_like_counter
                            << ", totalValues = "       << totalValues
                            << ", lnLikelihoodValue = " << lnLikelihoodValue
                            << " after "                << totalTime
                            << " seconds"
                            << std::endl;
  }

  if (m_env.subRank() == 0) {
#if 0
    std::cout << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
              << ", m_like_counter = "      << m_like_counter
              << ": totalValues = "         << totalValues
              << ", lnLikelihoodValue = "   << lnLikelihoodValue
              << " after "                  << totalTime
              << " seconds"
              << std::endl;
#else
    //sprintf(syncMsg,"In likelihoodRoutine(), likeCount = %u, total = %11.4e, lnL = %11.4e, time = %11.4e",
    //                m_like_counter,
    //                totalValues[0],
    //                lnLikelihoodValue,
    //                totalTime);
    //m_env.inter0Comm().syncPrintDebugMsg(syncMsg, 0, 1000);
#endif
  }

  m_env.subComm().Barrier();

  if (gradVector) {
    if (gradVector->sizeLocal() >= 4) {
      (*gradVector)[0] = lnLikelihoodValue;
      (*gradVector)[1] = 0.;
      (*gradVector)[2] = 0.;
      (*gradVector)[3] = 0.;
    }
  }

  if (m_like_counter == 0) {
    std::cout << "Exiting in likelihoodRoutine(), on purpose..." << std::endl;
    sleep(1);
    exit(1);
  }

  return lnLikelihoodValue;
}

#endif // __UQ_GCM_3_H__
