//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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

#ifndef __UQ_GCM_2_H__
#define __UQ_GCM_2_H__

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::uqGpmsaComputerModelClass(
  const char*                                               prefix,
  const uqGcmOptionsValuesClass*                            alternativeOptionsValues, // dakota
  const uqSimulationStorageClass <S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationStorage,
  const uqSimulationModelClass   <S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationModel,
  const uqExperimentStorageClass <S_V,S_M,D_V,D_M>*         experimentStorage,
  const uqExperimentModelClass   <S_V,S_M,D_V,D_M>*         experimentModel,
  const uqBaseVectorRVClass      <P_V,P_M>*                 thetaPriorRv)
  // Handle case of no experiments, that is, experiment pointers == NULL? (todo)
  :
  m_env                     (simulationStorage.env()),
  m_alternativeOptionsValues(),
  m_optionsObj              (NULL),
  m_s                       (NULL),
  m_e                       (NULL),
  m_j                       (NULL),
  m_z                       (NULL),
  m_t                       (NULL),
  m_st                      (NULL),
  m_jt                      (NULL),
  m_zt                      (NULL),
  m_thereIsExperimentalData ((experimentStorage != NULL) && (experimentModel != NULL) && (thetaPriorRv != NULL)),
  m_allOutputsAreScalar     (simulationStorage.outputSpace().dimLocal() == 1), // it might become 'false' if there are experiments
  m_formCMatrix             (true), // it will be updated
  m_cMatIsRankDefficient    (false),
  m_likelihoodFunction      (NULL),
  m_like_counter            (uqMiscUintDebugMessage(0,NULL))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = "                       << prefix
                            << ", alternativeOptionsValues = "     << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << ", m_thereIsExperimentalData = "    << m_thereIsExperimentalData
                            << std::endl;
  }

  if ((experimentStorage == NULL) &&
      (experimentModel   == NULL) &&
      (thetaPriorRv      == NULL)) {
    // Ok
    if (m_allOutputsAreScalar) {
      UQ_FATAL_TEST_MACRO((simulationModel.numBasis() != 1),
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                          "inconsistent numBasis (1)");
    }
  }
  else if ((experimentStorage != NULL) &&
           (experimentModel   != NULL) &&
           (thetaPriorRv      != NULL)) {
    // Ok
    for (unsigned int i = 0; m_allOutputsAreScalar && (i < experimentStorage->numExperiments()); ++i) {
      m_allOutputsAreScalar = m_allOutputsAreScalar && (experimentStorage->n_ys_transformed()[i] == 1);
    }
    if (m_allOutputsAreScalar) {
      UQ_FATAL_TEST_MACRO((simulationModel.numBasis() != 1) || (experimentModel->numBasis() != 1),
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                          "inconsistent numBasis (2)");
    }
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.worldRank(),
                        "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                        "inconsistent experimental information");
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = "                << prefix
                            << ", m_allOutputsAreScalar = " << m_allOutputsAreScalar
                            << std::endl;
  }

  m_formCMatrix = m_formCMatrix && (m_allOutputsAreScalar == false);

  //********************************************************************************
  // Handle options
  //********************************************************************************
  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new uqGpmsaComputerModelOptionsClass(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    //std::cout << "In uqGpmsaComputerModel constructor: scanning options from file..." << std::endl;
    m_optionsObj = new uqGpmsaComputerModelOptionsClass(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }

  m_formCMatrix = m_formCMatrix && m_optionsObj->m_ov.m_formCMatrix;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = "        << prefix
                            << ", 'final' m_formCMatrix = " << m_formCMatrix
                            << std::endl;
  }

  //********************************************************************************
  // Open output file // todo_r: is this necessary???
  //********************************************************************************
  if ((m_optionsObj->m_ov.m_dataOutputFileName                       != "."                                            ) &&
      (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end())) {
    m_env.openOutputFile(m_optionsObj->m_ov.m_dataOutputFileName,
                         UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                         m_optionsObj->m_ov.m_dataOutputAllowedSet,
                         false,
                         m_dataOutputFilePtrSet);
    UQ_FATAL_TEST_MACRO(m_dataOutputFilePtrSet.ofsVar == NULL,
                        m_env.worldRank(),
                        "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                        "data output file could not be created");
  }

  //********************************************************************************
  // Instantiate:
  // --> m_s
  // --> m_e
  // --> m_j
  // --> m_z
  // --> m_t
  // --> m_st
  // --> m_jt
  // --> m_zt
  //********************************************************************************
  //********************************************************************************
  // Instantiate Smat spaces
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
  m_s = new uqGcmSimulationInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>(*m_optionsObj,
                                                              m_allOutputsAreScalar,
                                                              simulationStorage,
                                                              simulationModel);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating m_s"
                            << ", experimentStorage = " << experimentStorage
                            << ", experimentModel = "   << experimentModel
                            << ", thetaPriorRv = "      << thetaPriorRv
                            << std::endl;
  }

  if (m_thereIsExperimentalData) {
    UQ_FATAL_TEST_MACRO(simulationStorage.scenarioSpace().dimLocal() != experimentStorage->scenarioSpace().dimLocal(),
                        m_env.worldRank(),
                        "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                        "inconsistent dimension of scenario space");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": about to instantiate m_e"
                              << std::endl;
    }

    m_e = new uqGcmExperimentInfoClass<S_V,S_M,D_V,D_M,P_V,P_M>(*m_optionsObj,
                                                                m_allOutputsAreScalar,
                                                                *experimentStorage,
                                                                *experimentModel,
                                                                *thetaPriorRv);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_e"
                              << std::endl;
    }

    m_j = new uqGcmJointInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_optionsObj,
                                                                   m_allOutputsAreScalar,
                                                                   *m_s,
                                                                   *m_e);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_j"
                              << std::endl;
    }

    if (m_allOutputsAreScalar) {
      m_z = new uqGcmZInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(m_allOutputsAreScalar,
                                                                 *m_s,
                                                                 *m_e);
    }
    else {
      m_z = new uqGcmZInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(m_formCMatrix,
                                                                 m_allOutputsAreScalar,
                                                                 *m_s,
                                                                 *m_e,
                                                                 *m_j);
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_z"
                              << std::endl;
    }

    m_t = new uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_s,
                                                                   *m_e);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_t"
                              << std::endl;
    }
  }
  else {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": about to instantiate m_z (no experiments)"
                              << std::endl;
    }

    m_z = new uqGcmZInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(m_formCMatrix,
                                                               m_allOutputsAreScalar,
                                                               *m_s);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_z (no experiments)"
                              << std::endl;
    }

    m_t = new uqGcmTotalInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_s);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_t (no experiments)"
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating non-tilde auxiliary structures"
                            << std::endl;
  }

  this->memoryCheck(0);

  if (m_formCMatrix) {
    //********************************************************************************
    // 'm_Cmat' has been formed: check if it is rank defficient
    //********************************************************************************
    if (m_z->m_Cmat_rank < m_z->m_Cmat->numCols()) m_cMatIsRankDefficient = true;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": m_z->m_Cmat_rank = "       << m_z->m_Cmat_rank 
                              << ", m_z->m_Cmat = "            << m_z->m_Cmat
                              << ", m_z->m_Cmat->numCols() = " << m_z->m_Cmat->numCols()
                              << ", m_cMatIsRankDefficient = " << m_cMatIsRankDefficient
                              << std::endl;
    }

    if (m_cMatIsRankDefficient == true) {
      //**********************************************************************************
      // 'm_Cmat' is rank difficient
      //**********************************************************************************
      if (m_optionsObj->m_ov.m_useTildeLogicForRankDefficientC) {
        //********************************************************************************
        // Use tilde logic
        //********************************************************************************
        if (m_thereIsExperimentalData) {
          //******************************************************************************
          // Tilde situation: form 'm_Bmat_tilde'
          // Tilde situation: form 'm_vu_tilde_space'
          // Tilde situation: form 'm_Lbmat'
          //******************************************************************************
          m_jt = new uqGcmJointTildeInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_optionsObj,*m_e,*m_j);

          //******************************************************************************
          // Tilde situation: form 'm_Kmat_tilde'
          // Tilde situation: form 'm_w_tilde_space'
          // Tilde situation: form 'm_Lkmat'
          //******************************************************************************
           m_st = new uqGcmSimulationTildeInfoClass<S_V,S_M,P_V,P_M,Q_V,Q_M>(*m_optionsObj,*m_s);

          //******************************************************************************
          // Tilde situation: form 'm_Cmat_tilde'
          //******************************************************************************
          m_zt = new uqGcmZTildeInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_optionsObj,*m_j,*m_z,*m_st,*m_jt);
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.worldRank(),
                              "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                              "incomplete code for the situation 'm_useTildeLogicForRankDefficientC == true' and 'm_thereIsExperimentalData == false'");
        }
      } // if (m_useTildeLogicForRankDefficientC)
      else {
        //********************************************************************************
        // Do not use tilde logic
        //********************************************************************************
        if (m_thereIsExperimentalData) {
          //******************************************************************************
          // Naive formation of 'm_Cmat_tilde'
          //******************************************************************************
          m_zt = new uqGcmZTildeInfoClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_optionsObj,*m_j,*m_z);
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.worldRank(),
                              "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
                              "incomplete code for the situation 'm_useTildeLogicForRankDefficientC == false' and 'm_thereIsExperimentalData == false'");
        }
      }
    } // if (m_cMatIsRankDefficient)
    else {
      //**********************************************************************************
      // 'm_Cmat' is full rank
      //**********************************************************************************
      // Ok. There is nothing extra to be done
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating tilde auxiliary structures"
                              << std::endl;
    }
  } // if (m_formCMatrix)
  else {
    //********************************************************************************
    // 'm_Cmat' will not be formed
    //********************************************************************************
    if (m_thereIsExperimentalData) {
      // Ok. There is nothing extra to be done
    }
    else {
      // Ok. No experimental data. There is nothing extra to be done // checar
      //UQ_FATAL_TEST_MACRO(true,
      //                    m_env.worldRank(),
      //                    "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
      //                    "incomplete code for the situation 'm_thereIsExperimentalData == false'");
    }
  }

  this->memoryCheck(1);

  //********************************************************************************
  // Generate prior sequence
  //********************************************************************************
  if (m_optionsObj->m_ov.m_priorSeqNumSamples > 0) {
    this->generatePriorSeq();
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished generating prior sequence"
                            << std::endl;
  }

  this->memoryCheck(2);

  //********************************************************************************
  // Instantiate likelihood function object
  //********************************************************************************
  m_likelihoodFunction = new uqGenericScalarFunctionClass<P_V,P_M>
                           ("like_",
                            m_t->m_totalDomain,
                            &staticLikelihoodRoutine,
                            (void *) this,
                            true); // routine computes [ln(function)]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating likelihood Function"
                            << std::endl;
  }

  this->memoryCheck(3);

  //********************************************************************************
  // Leave constructor
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~uqGpmsaComputerModelClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::destructor()..."
                            << std::endl;
  }

  delete m_zt;
  delete m_jt;
  delete m_st;
  delete m_t;
  delete m_z;
  delete m_j;
  delete m_e;
  delete m_s;

  delete m_dataOutputFilePtrSet.ofsVar;
  if (m_optionsObj) delete m_optionsObj;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::destructor()"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings(
  const uqMhOptionsValuesClass* alternativeOptionsValues, // dakota
  const P_V&                    totalInitialValues,
  const P_M*                    totalInitialProposalCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << ", totalInitialValues = "             << totalInitialValues
                            << std::endl;
  }

  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",1,3000000);

  UQ_FATAL_TEST_MACRO(m_t->m_totalPriorRv.imageSet().vectorSpace().dimLocal() != totalInitialValues.sizeLocal(),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",
                      "'m_totalPriorRv' and 'totalInitialValues' should have equal dimensions");

  if (totalInitialProposalCovMatrix) {
    UQ_FATAL_TEST_MACRO(m_t->m_totalPriorRv.imageSet().vectorSpace().dimLocal() != totalInitialProposalCovMatrix->numRowsLocal(),
                        m_env.worldRank(),
                        "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",
                        "'m_totalPriorRv' and 'totalInitialProposalCovMatrix' should have equal dimensions");
    UQ_FATAL_TEST_MACRO(totalInitialProposalCovMatrix->numCols() != totalInitialProposalCovMatrix->numRowsGlobal(),
                        m_env.worldRank(),
                        "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",
                        "'totalInitialProposalCovMatrix' should be a square matrix");
  }

  //std::cout << "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 002"
  //          << std::endl;

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_t->m_solutionDomain = uqInstantiateIntersection(m_t->m_totalPriorRv.pdf().domainSet(),m_likelihoodFunction->domainSet());

  //std::cout << "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 003"
  //          << std::endl;

  m_t->m_solutionPdf = new uqBayesianJointPdfClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                            m_t->m_totalPriorRv.pdf(),
                                                            *m_likelihoodFunction,
                                                            1.,
                                                            *(m_t->m_solutionDomain));

  //std::cout << "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 004"
  //          << std::endl;

  m_t->m_totalPostRv.setPdf(*(m_t->m_solutionPdf));

  //std::cout << "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 005"
  //          << std::endl;

  // Compute output realizer: Metropolis-Hastings approach
  m_t->m_chain = new uqSequenceOfVectorsClass<P_V,P_M>(m_t->m_totalPostRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");

  //std::cout << "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 006"
  //          << std::endl;

  m_t->m_mhSeqGenerator = new uqMetropolisHastingsSGClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(), // dakota
                                                                   alternativeOptionsValues,
                                                                   m_t->m_totalPostRv,
                                                                   totalInitialValues,
                                                                   totalInitialProposalCovMatrix);

  //std::cout << "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 007"
  //          << std::endl;

  m_t->m_mhSeqGenerator->generateSequence(*(m_t->m_chain),NULL,NULL);

  // todo_rr
  // m_totalPostMean
  // m_totalPostMedian
  // m_totalPostMode
  // m_totalPostMaxLnValue
  // m_totalMLE
  // m_totalLikeMaxLnValue

  m_t->m_solutionRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                         *(m_t->m_chain));

  m_t->m_totalPostRv.setRealizer(*(m_t->m_solutionRealizer));

  //m_env.fullComm().syncPrintDebugMsg("In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings(), code place 1",3,3000000);

  //m_env.fullComm().syncPrintDebugMsg("Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",1,3000000);
  m_env.fullComm().Barrier();

  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << ", after "                            << totalTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc(
  const uqMhOptionsValuesClass* alternativeOptionsValues, // dakota
  const P_V&                    totalInitialValues,
  const P_M*                    totalInitialProposalCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << ", totalInitialValues = "             << totalInitialValues
                            << std::endl;
  }

  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()",1,3000000);

  // ppp

  m_env.fullComm().syncPrintDebugMsg("Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()",1,3000000);
  m_env.fullComm().Barrier();

  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << ", after "                            << totalTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(
  const S_V& newScenarioVec,
  const P_V& newParameterVec,
        P_V& vuMeanVec,
        P_M& vuCovMatrix,
        P_V& vMeanVec,
        P_M& vCovMatrix,
        P_V& uMeanVec,
        P_M& uCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                            << ", m_predVU_counter = " << m_j->m_predVU_counter
                            << ", newScenarioVec = "   << newScenarioVec
                            << ", newParameterVec = "  << newParameterVec
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(newScenarioVec.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'newScenarioVec'");

  UQ_FATAL_TEST_MACRO(newParameterVec.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'newParameterVec'");

  UQ_FATAL_TEST_MACRO(vuMeanVec.sizeLocal() != (m_e->m_paper_p_delta+m_s->m_paper_p_eta),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'vuMeanVec'");

  UQ_FATAL_TEST_MACRO(vuCovMatrix.numRowsLocal() != (m_e->m_paper_p_delta + m_s->m_paper_p_eta),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'vuCovMatrix.numRowsLocal()'");

  UQ_FATAL_TEST_MACRO(vuCovMatrix.numCols() != (m_e->m_paper_p_delta + m_s->m_paper_p_eta),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'vuCovMatrix.numCols()'");

  UQ_FATAL_TEST_MACRO(vMeanVec.sizeLocal() != m_e->m_paper_p_delta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'vMeanVec'");

  UQ_FATAL_TEST_MACRO(vCovMatrix.numRowsLocal() != m_e->m_paper_p_delta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'vCovMatrix.numRowsLocal()'");

  UQ_FATAL_TEST_MACRO(vCovMatrix.numCols() != m_e->m_paper_p_delta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'vCovMatrix.numCols()'");

  UQ_FATAL_TEST_MACRO(uMeanVec.sizeLocal() != m_s->m_paper_p_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'uMeanVec'");

  UQ_FATAL_TEST_MACRO(uCovMatrix.numRowsLocal() != m_s->m_paper_p_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'uCovMatrix.numRowsLocal()'");

  UQ_FATAL_TEST_MACRO(uCovMatrix.numCols() != m_s->m_paper_p_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                      "invalid 'uCovMatrix.numCols()'");

  if (m_optionsObj->m_ov.m_predVUsBySamplingRVs) {
  }

  if (m_optionsObj->m_ov.m_predVUsBySummingRVs) {
    unsigned int numSamples = (unsigned int) ((double) m_t->m_totalPostRv.realizer().subPeriod())/((double) m_optionsObj->m_ov.m_predLag);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                              << ": m_t->m_totalPostRv.realizer().subPeriod() = " << m_t->m_totalPostRv.realizer().subPeriod()
                              << ", m_optionsObj->m_ov.m_predLag = "              << m_optionsObj->m_ov.m_predLag
                              << std::endl;
    }

    uqSequenceOfVectorsClass<P_V,P_M> unique_vu_means(m_j->m_unique_vu_space,numSamples,m_optionsObj->m_prefix+"vu_means");
    P_M mean_of_unique_vu_covMatrices(m_j->m_unique_vu_space.zeroVector());

    P_V totalSample(m_t->m_totalSpace.zeroVector());
    P_V muVec1     (m_z->m_z_space.zeroVector());
    P_V muVec2     (m_j->m_vu_space.zeroVector());
    P_M sigmaMat11 (m_z->m_z_space.zeroVector());
    P_M sigmaMat12 (m_env,muVec1.map(),muVec2.sizeGlobal());
    P_M sigmaMat21 (m_env,muVec2.map(),muVec1.sizeGlobal());
    P_M sigmaMat22 (m_j->m_vu_space.zeroVector());

    P_M here_Smat_z_hat_v_asterisk  (m_env, m_z->m_z_space.map(),        m_e->m_paper_p_delta);
    P_M here_Smat_z_hat_v_asterisk_t(m_env, m_e->m_unique_v_space.map(), m_z->m_z_size       );
    P_M here_Smat_z_hat_u_asterisk  (m_env, m_z->m_z_space.map(),        m_s->m_paper_p_eta  );
    P_M here_Smat_z_hat_u_asterisk_t(m_env, m_j->m_unique_u_space.map(), m_z->m_z_size       );

    std::vector<const P_M*> twoMats_uw(2,NULL);
    std::vector<const P_M*> twoMats_12(2,NULL);
    std::vector<const P_M*> twoMats_21(2,NULL);
    std::vector<const P_M*> twoMats_22(2,NULL);
    for (unsigned int sampleId = 0; sampleId < numSamples; ++sampleId) {
      m_j->m_predVU_counter++;

      if (sampleId > 0) {
        for (unsigned int i = 1; i < m_optionsObj->m_ov.m_predLag; ++i) { // Yes, '1'
          m_t->m_totalPostRv.realizer().realization(totalSample);
        }
      }
      m_t->m_totalPostRv.realizer().realization(totalSample);

      unsigned int currPosition = 0;
      totalSample.cwExtract(currPosition,m_s->m_tmp_1lambdaEtaVec); // Total of '1' in paper
      currPosition += m_s->m_tmp_1lambdaEtaVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_s->m_tmp_2lambdaWVec);   // Total of 'p_eta' in paper
      currPosition += m_s->m_tmp_2lambdaWVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_s->m_tmp_3rhoWVec);      // Total of 'p_eta*(p_x+p_t)' in paper
      currPosition += m_s->m_tmp_3rhoWVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_s->m_tmp_4lambdaSVec);   // Total of 'p_eta' in matlab code
      currPosition += m_s->m_tmp_4lambdaSVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_5lambdaYVec);   // Total of '1' in paper
      currPosition += m_e->m_tmp_5lambdaYVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_6lambdaVVec);   // Total of 'F' in paper
      currPosition += m_e->m_tmp_6lambdaVVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_7rhoVVec);      // Total of 'F*p_x' in paper
      currPosition += m_e->m_tmp_7rhoVVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_8thetaVec);     // Application specific
      currPosition += m_e->m_tmp_8thetaVec.sizeLocal();
      UQ_FATAL_TEST_MACRO(currPosition != totalSample.sizeLocal(),
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                          "'currPosition' and 'totalSample.sizeLocal()' should be equal");

      //********************************************************************************
      // Submatrix (1,1): Compute '\Sigma_z_hat' matrix
      //********************************************************************************
      // Fill m_Rmat_v_is,  m_Smat_v_is,  m_Smat_v
      // Fill m_Rmat_u_is,  m_Smat_u_is,  m_Smat_u
      // Fill m_Rmat_w_is,  m_Smat_w_is,  m_Smat_w
      // Fill m_Rmat_uw_is, m_Smat_uw_is, m_Smat_uw
      // Then fill m_tmp_Smat_z
      // Fill m_Rmat_extra
      // Then fill m_tmp_Smat_z_hat
      this->formSigma_z_hat(m_s->m_tmp_1lambdaEtaVec, // todo_rr0: necessary to recompute for every posterior sample?
                            m_s->m_tmp_2lambdaWVec,
                            m_s->m_tmp_3rhoWVec,
                            m_s->m_tmp_4lambdaSVec,
                            m_e->m_tmp_5lambdaYVec,
                            m_e->m_tmp_6lambdaVVec,
                            m_e->m_tmp_7rhoVVec,
                            newParameterVec, //m_e->m_tmp_8thetaVec,
                            m_j->m_predVU_counter);

      //********************************************************************************
      // Submatrix (1,2): Compute '\Sigma_z_hat_v_asterisk' matrix
      // Submatrix (2,1): Compute '\Sigma_z_hat_v_asterisk' transpose matrix
      //********************************************************************************
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                                << ", m_predVU_counter = " << m_j->m_predVU_counter
                                << ": about to populate 'm_Smat_v_hat_v_asterisk'"
                                << ", m_e->m_Smat_v_hat_v_asterisk_is.size() = " << m_e->m_Smat_v_hat_v_asterisk_is.size() // 13
                                << ", m_e->m_tmp_rho_v_vec.sizeLocal() = "       << m_e->m_tmp_rho_v_vec.sizeLocal()       //  1
                                << ", m_e->m_tmp_7rhoVVec.sizeLocal() = "        << m_e->m_tmp_7rhoVVec.sizeLocal()        //  1
                                << std::endl;
      }
      UQ_FATAL_TEST_MACRO((m_e->m_Smat_v_hat_v_asterisk_is.size() * m_e->m_tmp_rho_v_vec.sizeLocal()) != m_e->m_tmp_7rhoVVec.sizeLocal(),
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                          "invalid size for 'v' variables");
      unsigned int initialPos = 0;
      for (unsigned int i = 0; i < m_e->m_Smat_v_hat_v_asterisk_is.size(); ++i) {
        m_e->m_tmp_7rhoVVec.cwExtract(initialPos,m_e->m_tmp_rho_v_vec);
        initialPos += m_e->m_tmp_rho_v_vec.sizeLocal();
        m_e->m_Rmat_v_hat_v_asterisk_is[i]->cwSet(0.);
        this->fillR_formula2_for_Sigma_v_hat_v_asterisk(m_s->m_paper_xs_asterisks_standard,
                                                        m_s->m_paper_ts_asterisks_standard,
                                                        newScenarioVec,
                                                        newParameterVec, //m_e->m_tmp_8thetaVec,
                                                        m_e->m_tmp_rho_v_vec,
                                                        *(m_e->m_Rmat_v_hat_v_asterisk_is[i]), // IMPORTANT-28
                                                        m_j->m_predVU_counter);
        m_e->m_Smat_v_hat_v_asterisk_is[i]->cwSet(0.);
        // IMPORTANT-28
	*(m_e->m_Smat_v_hat_v_asterisk_is[i]) = (1./m_e->m_tmp_6lambdaVVec[i]) * *(m_e->m_Rmat_v_hat_v_asterisk_is[i]);
      }
      m_e->m_Smat_v_hat_v_asterisk.cwSet(0.);
      m_e->m_Smat_v_hat_v_asterisk.fillWithBlocksDiagonally(0,0,m_e->m_Smat_v_hat_v_asterisk_is,true,true);
      m_e->m_Smat_v_hat_v_asterisk_t.fillWithTranspose(0,0,m_e->m_Smat_v_hat_v_asterisk,true,true);

      here_Smat_z_hat_v_asterisk.cwSet(0.); 
      here_Smat_z_hat_v_asterisk.cwSet(0,0,m_e->m_Smat_v_hat_v_asterisk); // checar
      here_Smat_z_hat_v_asterisk_t.fillWithTranspose(0,0,here_Smat_z_hat_v_asterisk,true,true);

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                                << ", m_predVU_counter = " << m_j->m_predVU_counter
                                << ": finished instantiating 'm_Smat_v_hat_v_asterisk'"
                                << std::endl;
      }

      //********************************************************************************
      // Submatrix (1,3): Compute '\Sigma_z_hat_u_asterisk' matrix
      // Submatrix (3,1): Compute '\Sigma_z_hat_u_asterisk' transpose matrix
      //********************************************************************************
      initialPos = 0;
      for (unsigned int i = 0; i < m_j->m_Smat_u_hat_u_asterisk_is.size(); ++i) {
	m_s->m_tmp_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
        initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();

        m_j->m_Rmat_u_hat_u_asterisk_is[i]->cwSet(0.);
        this->fillR_formula1_for_Sigma_u_hat_u_asterisk(m_s->m_paper_xs_asterisks_standard,
                                                        m_s->m_paper_ts_asterisks_standard,
                                                        newScenarioVec,
                                                        newParameterVec, //m_e->m_tmp_8thetaVec,
                                                        m_s->m_tmp_rho_w_vec,
                                                        *(m_j->m_Rmat_u_hat_u_asterisk_is[i]),
                                                        m_j->m_predVU_counter);
        m_j->m_Smat_u_hat_u_asterisk_is[i]->cwSet(0.);
        *(m_j->m_Smat_u_hat_u_asterisk_is[i]) = (1./m_s->m_tmp_2lambdaWVec[i]) * *(m_j->m_Rmat_u_hat_u_asterisk_is[i]);

        m_j->m_Rmat_w_hat_u_asterisk_is[i]->cwSet(0.);
        this->fillR_formula1_for_Sigma_w_hat_u_asterisk(m_s->m_paper_xs_asterisks_standard,
                                                        m_s->m_paper_ts_asterisks_standard,
                                                        newScenarioVec,
                                                        newParameterVec, //m_e->m_tmp_8thetaVec,
                                                        m_s->m_tmp_rho_w_vec,
                                                        *(m_j->m_Rmat_w_hat_u_asterisk_is[i]),
                                                        m_j->m_predVU_counter);
        m_j->m_Smat_w_hat_u_asterisk_is[i]->cwSet(0.);
        *(m_j->m_Smat_w_hat_u_asterisk_is[i]) = (1./m_s->m_tmp_2lambdaWVec[i]) * *(m_j->m_Rmat_w_hat_u_asterisk_is[i]);
      }

      m_j->m_Smat_u_hat_u_asterisk.cwSet(0.);
      m_j->m_Smat_u_hat_u_asterisk.fillWithBlocksDiagonally(0,0,m_j->m_Smat_u_hat_u_asterisk_is,true,true);
      m_j->m_Smat_u_hat_u_asterisk_t.fillWithTranspose(0,0,m_j->m_Smat_u_hat_u_asterisk,true,true);

      m_j->m_Smat_w_hat_u_asterisk.cwSet(0.);
      m_j->m_Smat_w_hat_u_asterisk.fillWithBlocksDiagonally(0,0,m_j->m_Smat_w_hat_u_asterisk_is,true,true);
      m_j->m_Smat_w_hat_u_asterisk_t.fillWithTranspose(0,0,m_j->m_Smat_w_hat_u_asterisk,true,true);

      twoMats_uw[0] = &m_j->m_Smat_u_hat_u_asterisk;
      twoMats_uw[1] = &m_j->m_Smat_w_hat_u_asterisk;
      here_Smat_z_hat_u_asterisk.cwSet(0.);
      here_Smat_z_hat_u_asterisk.fillWithBlocksVertically(m_e->m_v_size,0,twoMats_uw,true,true); // checar
      here_Smat_z_hat_u_asterisk_t.fillWithTranspose(0,0,here_Smat_z_hat_u_asterisk,true,true);

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                                << ", m_predVU_counter = " << m_j->m_predVU_counter
                                << ": finished instantiating '>m_Smat_z_hat_u_asterisk'"
                                << std::endl;
      }

      //********************************************************************************
      // Submatrix (2,2): Compute '\Sigma_v_asterisk_v_asterisk' matrix
      //********************************************************************************
      UQ_FATAL_TEST_MACRO(m_e->m_Smat_v_asterisk_v_asterisk.numRowsLocal() != m_e->m_paper_p_delta,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                          "invalid 'm_Smat_v_asterisk_v_asterisk.numRowsLocal()'");
      UQ_FATAL_TEST_MACRO(m_e->m_tmp_6lambdaVVec.sizeLocal() != m_e->m_paper_p_delta,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                          "invalid 'm_tmp_6lambdaVVec.sizeLocal()'");

      m_e->m_Smat_v_asterisk_v_asterisk.cwSet(0.);
      for (unsigned int i = 0; i < m_e->m_paper_p_delta; ++i) {
        m_e->m_Smat_v_asterisk_v_asterisk(i,i) = 1./m_e->m_tmp_6lambdaVVec[i];
      }

      //********************************************************************************
      // Submatrix (3,3): Compute '\Sigma_u_asterisk_u_asterisk' matrix
      //********************************************************************************
      UQ_FATAL_TEST_MACRO(m_j->m_Smat_u_asterisk_u_asterisk.numRowsLocal() != m_s->m_paper_p_eta,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                          "invalid 'm_Smat_u_asterisk_u_asterisk.numRowsLocal()'");
      UQ_FATAL_TEST_MACRO(m_s->m_tmp_2lambdaWVec.sizeLocal() != m_s->m_paper_p_eta,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)",
                          "invalid 'm_tmp_2lambdaWVec.sizeLocal()'");

      m_j->m_Smat_u_asterisk_u_asterisk.cwSet(0.);
      for (unsigned int i = 0; i < m_s->m_paper_p_eta; ++i) {
        m_j->m_Smat_u_asterisk_u_asterisk(i,i) = 1./m_s->m_tmp_2lambdaWVec[i] + 1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
      }

      //********************************************************************************
      // For the current parameter sample, compute the mean vector and the covariance
      // matrix of the corresponding Gassian RVs 'v' and 'u' conditioned on 'z_hat'
      // Use the already computed (in constructor) 'm_Zvec_hat'
      //********************************************************************************
      P_V unique_vu_vec(m_j->m_unique_vu_space.zeroVector());
      P_M unique_vu_mat(m_j->m_unique_vu_space.zeroVector());

      twoMats_12[0] = &here_Smat_z_hat_v_asterisk_t;
      twoMats_12[1] = &here_Smat_z_hat_u_asterisk_t;
      twoMats_21[0] = &here_Smat_z_hat_v_asterisk;
      twoMats_21[1] = &here_Smat_z_hat_u_asterisk;
      twoMats_22[0] = &m_e->m_Smat_v_asterisk_v_asterisk;
      twoMats_22[1] = &m_j->m_Smat_u_asterisk_u_asterisk;

      sigmaMat11 = m_z->m_tmp_Smat_z_hat;
      sigmaMat12.fillWithBlocksHorizontally(0,0,twoMats_12,true,true);
      sigmaMat21.fillWithBlocksVertically  (0,0,twoMats_21,true,true);
      sigmaMat22.fillWithBlocksDiagonally  (0,0,twoMats_22,true,true);
      uqComputeConditionalGaussianVectorRV(muVec1,muVec2,sigmaMat11,sigmaMat12,sigmaMat21,sigmaMat22,m_z->m_Zvec_hat,unique_vu_vec,unique_vu_mat);

      unique_vu_means.setPositionValues(sampleId,unique_vu_vec);
      m_j->m_predVU_summingRVs_mean_of_unique_vu_covMatrices += unique_vu_mat;
    }

    //********************************************************************************
    // Final calculations
    //********************************************************************************
    m_j->m_predVU_summingRVs_mean_of_unique_vu_covMatrices *= (1./(double) numSamples);

    m_j->m_predVU_summingRVs_unique_vu_meanVec = unique_vu_means.unifiedMeanPlain();

    m_j->m_predVU_summingRVs_covMatrix_of_unique_vu_means.cwSet(0.);
    m_j->m_predVU_summingRVs_corrMatrix_of_unique_vu_means.cwSet(0.);
    uqComputeCovCorrMatricesBetweenVectorSequences(unique_vu_means,
                                                   unique_vu_means,
                                                   unique_vu_means.subSequenceSize(),
                                                   m_j->m_predVU_summingRVs_covMatrix_of_unique_vu_means,
                                                   m_j->m_predVU_summingRVs_corrMatrix_of_unique_vu_means);

    vuMeanVec = m_j->m_predVU_summingRVs_unique_vu_meanVec;
    vuMeanVec.cwExtract(0,vMeanVec);
    vuMeanVec.cwExtract(m_e->m_paper_p_delta,uMeanVec);

    P_M vuCovMatrix(m_j->m_unique_vu_space.zeroVector());
    vuCovMatrix = m_j->m_predVU_summingRVs_mean_of_unique_vu_covMatrices + m_j->m_predVU_summingRVs_covMatrix_of_unique_vu_means;
    vuCovMatrix.cwExtract(0,0,vCovMatrix);                                       // checar
    vuCovMatrix.cwExtract(m_e->m_paper_p_delta,m_e->m_paper_p_delta,uCovMatrix); // checar

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                              << ", m_predVU_counter = " << m_j->m_predVU_counter
                              << ": finished computing all means and covariances"
                              << "\n  vuMeanVec = "                                        << vuMeanVec
                              << "\n  m_predW_summingRVs_covMatrix_of_unique_vu_means = "  << m_j->m_predVU_summingRVs_covMatrix_of_unique_vu_means
                              << "\n  m_predW_summingRVs_mean_of_unique_vu_covMatrices = " << m_j->m_predVU_summingRVs_mean_of_unique_vu_covMatrices
                              << "\n  vuCovMatrix = "                                      << vuCovMatrix
                              << std::endl;

    }

    mean_of_unique_vu_covMatrices *= (1./(double) numSamples);
  }

  if (m_optionsObj->m_ov.m_predVUsAtKeyPoints) {
  }

  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                            << ", m_predVU_counter = " << m_j->m_predVU_counter
                            << ", newScenarioVec = "   << newScenarioVec
                            << ", after "              << totalTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint(
  const S_V& newScenarioVec,
  const P_V& newParameterVec,
  const P_V* forcingSampleVecForDebug, // Usually NULL
        P_V& wMeanVec,
        P_M& wCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                            << ", m_predW_counter = " << m_s->m_predW_counter
                            << ", newScenarioVec = "  << newScenarioVec
                            << ", newParameterVec = " << newParameterVec
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(newScenarioVec.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                      "invalid 'newScenarioVec'");

  UQ_FATAL_TEST_MACRO(newParameterVec.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                      "invalid 'newParameterVec'");

  UQ_FATAL_TEST_MACRO(wMeanVec.sizeLocal() != m_s->m_paper_p_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                      "invalid 'wMeanVec'");

  UQ_FATAL_TEST_MACRO(wCovMatrix.numRowsLocal() != m_s->m_paper_p_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                      "invalid 'wCovMatrix.numRowsLocal()'");

  UQ_FATAL_TEST_MACRO(wCovMatrix.numCols() != m_s->m_paper_p_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                      "invalid 'wCovMatrix.numCols()'");

  if (m_optionsObj->m_ov.m_predWsBySamplingRVs) {
  }

  if (m_optionsObj->m_ov.m_predWsBySummingRVs) {
    unsigned int numSamples = (unsigned int) ((double) m_t->m_totalPostRv.realizer().subPeriod())/((double) m_optionsObj->m_ov.m_predLag);
    if (forcingSampleVecForDebug) {
      numSamples = 1;
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                              << ": m_t->m_totalPostRv.realizer().subPeriod() = " << m_t->m_totalPostRv.realizer().subPeriod()
                              << ", m_optionsObj->m_ov.m_predLag = "              << m_optionsObj->m_ov.m_predLag
                              << ", numSamples = "                                << numSamples
                              << std::endl;
    }

    uqSequenceOfVectorsClass<P_V,P_M> unique_w_means(m_s->m_unique_w_space,numSamples,m_optionsObj->m_prefix+"w_means");

    P_V totalSample(m_t->m_totalSpace.zeroVector());
    P_V muVec1     (m_s->m_unique_w_space.zeroVector());
    P_V muVec2     (m_s->m_w_space.zeroVector());
    P_M sigmaMat11 (m_s->m_unique_w_space.zeroVector());
    P_M sigmaMat12 (m_env,muVec1.map(),muVec2.sizeGlobal());
    P_M sigmaMat21 (m_env,muVec2.map(),muVec1.sizeGlobal());
    P_M sigmaMat22 (m_s->m_w_space.zeroVector());
    for (unsigned int sampleId = 0; sampleId < numSamples; ++sampleId) {
      m_s->m_predW_counter++;

      if (sampleId > 0) {
        for (unsigned int i = 1; i < m_optionsObj->m_ov.m_predLag; ++i) { // Yes, '1'
          m_t->m_totalPostRv.realizer().realization(totalSample);
        }
      }
      m_t->m_totalPostRv.realizer().realization(totalSample);
      if (forcingSampleVecForDebug) {
        totalSample = *forcingSampleVecForDebug;
      }

      unsigned int currPosition = 0;
      totalSample.cwExtract(currPosition,m_s->m_tmp_1lambdaEtaVec); // Total of '1' in paper
      currPosition += m_s->m_tmp_1lambdaEtaVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_s->m_tmp_2lambdaWVec);   // Total of 'p_eta' in paper
      currPosition += m_s->m_tmp_2lambdaWVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_s->m_tmp_3rhoWVec);      // Total of 'p_eta*(p_x+p_t)' in paper
      currPosition += m_s->m_tmp_3rhoWVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_s->m_tmp_4lambdaSVec);   // Total of 'p_eta' in matlab code
      currPosition += m_s->m_tmp_4lambdaSVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_5lambdaYVec);   // Total of '1' in paper
      currPosition += m_e->m_tmp_5lambdaYVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_6lambdaVVec);   // Total of 'F' in paper
      currPosition += m_e->m_tmp_6lambdaVVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_7rhoVVec);      // Total of 'F*p_x' in paper
      currPosition += m_e->m_tmp_7rhoVVec.sizeLocal();
      totalSample.cwExtract(currPosition,m_e->m_tmp_8thetaVec);     // Application specific
      currPosition += m_e->m_tmp_8thetaVec.sizeLocal();
      UQ_FATAL_TEST_MACRO(currPosition != totalSample.sizeLocal(),
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                          "'currPosition' and 'totalSample.sizeLocal()' should be equal");

      //********************************************************************************
      // Submatrix (1,1): Compute '\Sigma_w_hat' matrix
      //********************************************************************************
      // Fill m_Rmat_w_is, m_Smat_w_is, m_Smat_w
      // Then add to m_Smat_w in order to obtain m_Smat_w_hat
      this->formSigma_w_hat(m_s->m_tmp_1lambdaEtaVec,
                            m_s->m_tmp_2lambdaWVec,
                            m_s->m_tmp_3rhoWVec,
                            m_s->m_tmp_4lambdaSVec,
                            newParameterVec, // todo_rr0
                            m_s->m_predW_counter);

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                  << ": m_s->m_Smat_w_hat = " << m_s->m_Smat_w_hat
                                  << std::endl;
        }
      }

      //********************************************************************************
      // Submatrix (1,2): Compute '\Sigma_w_hat_w_asterisk' matrix
      // Submatrix (2,1): Compute '\Sigma_w_hat_w_asterisk' transpose matrix
      //********************************************************************************
      unsigned int initialPos = 0;
      for (unsigned int i = 0; i < m_s->m_Smat_w_hat_w_asterisk_is.size(); ++i) {
        m_s->m_tmp_3rhoWVec.cwExtract(initialPos,m_s->m_tmp_rho_w_vec);
        initialPos += m_s->m_tmp_rho_w_vec.sizeLocal();
        m_s->m_Rmat_w_hat_w_asterisk_is[i]->cwSet(0.); // This matrix is rectangular: m_paper_m X 1
        this->fillR_formula1_for_Sigma_w_hat_w_asterisk(m_s->m_paper_xs_asterisks_standard, // IMPORTANT
                                                        m_s->m_paper_ts_asterisks_standard,
                                                        newScenarioVec,
                                                        newParameterVec,
                                                        m_s->m_tmp_rho_w_vec,
                                                        *(m_s->m_Rmat_w_hat_w_asterisk_is[i]),
                                                        m_s->m_predW_counter);
        m_s->m_Smat_w_hat_w_asterisk_is[i]->cwSet(0.);
        *(m_s->m_Smat_w_hat_w_asterisk_is[i]) = (1./m_s->m_tmp_2lambdaWVec[i]) * *(m_s->m_Rmat_w_hat_w_asterisk_is[i]);
      }
      m_s->m_Smat_w_hat_w_asterisk.cwSet(0.);
      m_s->m_Smat_w_hat_w_asterisk.fillWithBlocksDiagonally(0,0,m_s->m_Smat_w_hat_w_asterisk_is,true,true);
      m_s->m_Smat_w_hat_w_asterisk_t.fillWithTranspose(0,0,m_s->m_Smat_w_hat_w_asterisk,true,true);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
        *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                << ", m_predW_counter = " << m_s->m_predW_counter
                                << ": finished instantiating 'm_Smat_w_hat_w_asterisk'"
                                << std::endl;
      }

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                  << ": m_s->m_Smat_w_hat_w_asterisk = " << m_s->m_Smat_w_hat_w_asterisk
                                  << std::endl;
        }
      }

      //********************************************************************************
      // Submatrix (2,2): Compute '\Sigma_w_asterisk_w_asterisk' matrix
      //********************************************************************************
      UQ_FATAL_TEST_MACRO(m_s->m_Smat_w_asterisk_w_asterisk.numRowsLocal() != m_s->m_paper_p_eta,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                          "invalid 'm_Smat_w_asterisk_w_asterisk.numRowsLocal()'");
      UQ_FATAL_TEST_MACRO(m_s->m_tmp_2lambdaWVec.sizeLocal() != m_s->m_paper_p_eta,
                          m_env.worldRank(),
                          "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()",
                          "invalid 'm_tmp_2lambdaWVec.sizeLocal()'");

      m_s->m_Smat_w_asterisk_w_asterisk.cwSet(0.);
      for (unsigned int i = 0; i < m_s->m_paper_p_eta; ++i) {
        m_s->m_Smat_w_asterisk_w_asterisk(i,i) = 1./m_s->m_tmp_2lambdaWVec[i] + 1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
      }

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                  << ": m_s->m_Smat_w_asterisk_w_asterisk = " << m_s->m_Smat_w_asterisk_w_asterisk
                                  << std::endl;
        }
      }

      //********************************************************************************
      // For the current parameter sample, compute the mean vector and the covariance
      // matrix of the corresponding Gassian RV 'w' conditioned on 'w_hat'
      // Use the already computed (in constructor) 'm_Zvec_hat_w'
      //********************************************************************************
      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                  << ": muVec1 = "            << muVec1
                                  << ", muVec2 = "            << muVec2
                                  << ", m_s->m_Zvec_hat_w = " << m_s->m_Zvec_hat_w
                                  << std::endl;
        }
      }

      P_V unique_w_vec(m_s->m_unique_w_space.zeroVector());
      P_M unique_w_mat(m_s->m_unique_w_space.zeroVector());

      sigmaMat11 = m_s->m_Smat_w_asterisk_w_asterisk;
      sigmaMat12 = m_s->m_Smat_w_hat_w_asterisk_t;
      sigmaMat21 = m_s->m_Smat_w_hat_w_asterisk;
      sigmaMat22 = m_s->m_Smat_w_hat;
      uqComputeConditionalGaussianVectorRV(muVec1,muVec2,sigmaMat11,sigmaMat12,sigmaMat21,sigmaMat22,m_s->m_Zvec_hat_w,unique_w_vec,unique_w_mat);

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                  << ", unique_w_vec = " << unique_w_vec
                                  << ", unique_w_mat = " << unique_w_mat
                                  << std::endl;
        }
      }

      unique_w_means.setPositionValues(sampleId,unique_w_vec);
      m_s->m_predW_summingRVs_mean_of_unique_w_covMatrices += unique_w_mat;

      // aqui: display periodic message
    }

    //********************************************************************************
    // Final calculations
    //********************************************************************************
    m_s->m_predW_summingRVs_mean_of_unique_w_covMatrices *= (1./(double) numSamples);

    m_s->m_predW_summingRVs_unique_w_meanVec = unique_w_means.unifiedMeanPlain();

    m_s->m_predW_summingRVs_covMatrix_of_unique_w_means.cwSet(0.);
    m_s->m_predW_summingRVs_corrMatrix_of_unique_w_means.cwSet(0.);
    uqComputeCovCorrMatricesBetweenVectorSequences(unique_w_means,
                                                   unique_w_means,
                                                   unique_w_means.subSequenceSize(),
                                                   m_s->m_predW_summingRVs_covMatrix_of_unique_w_means,
                                                   m_s->m_predW_summingRVs_corrMatrix_of_unique_w_means);

    wMeanVec   = m_s->m_predW_summingRVs_unique_w_meanVec;
    wCovMatrix = m_s->m_predW_summingRVs_mean_of_unique_w_covMatrices + m_s->m_predW_summingRVs_covMatrix_of_unique_w_means;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
      *m_env.subDisplayFile() << "In uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                              << ", m_predW_counter = " << m_s->m_predW_counter
                              << ": finished computing all means and covariances"
                              << "\n  wMeanVec = "                                        << wMeanVec
                              << "\n  m_predW_summingRVs_covMatrix_of_unique_w_means = "  << m_s->m_predW_summingRVs_covMatrix_of_unique_w_means
                              << "\n  m_predW_summingRVs_mean_of_unique_w_covMatrices = " << m_s->m_predW_summingRVs_mean_of_unique_w_covMatrices
                              << "\n  wCovMatrix = "                                      << wCovMatrix
                              << std::endl;

    }
  }

  if (m_optionsObj->m_ov.m_predWsAtKeyPoints) {
  }

  double totalTime = uqMiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                            << ", m_predW_counter = " << m_s->m_predW_counter
                            << ", newScenarioVec = "  << newScenarioVec
                            << ", newParameterVec = " << newParameterVec
                            << ", after "             << totalTime
                            << " seconds"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults(
  const S_V& newScenarioVec,
  const D_M& newKmat_interp,
  const D_M& newDmat,
        D_V& simulationOutputMeanVec, // todo_rr: pass as pointer
        D_V& discrepancyMeanVec)      // todo_rr: pass as pointer
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()"
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(newScenarioVec.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()",
                      "invalid 'newScenarioVec'");

  UQ_FATAL_TEST_MACRO((newKmat_interp.numRowsLocal() != m_s->m_paper_n_eta) || (newKmat_interp.numCols() != m_s->m_paper_p_eta),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()",
                      "invalid 'newKmat_interp'");

  UQ_FATAL_TEST_MACRO((newDmat.numRowsLocal() != m_s->m_paper_n_eta) || (newDmat.numCols() != m_e->m_paper_p_delta),
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()",
                      "invalid 'newDmat'");

  UQ_FATAL_TEST_MACRO(simulationOutputMeanVec.sizeLocal() != m_s->m_paper_n_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()",
                      "invalid 'simulationOutputMeanVec'");

  UQ_FATAL_TEST_MACRO(discrepancyMeanVec.sizeLocal() != m_s->m_paper_n_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()",
                      "invalid 'discrepancyMeanVec'");

  P_V vMeanVec  (m_e->m_unique_v_space.zeroVector());
  P_M vCovMatrix(m_e->m_unique_v_space.zeroVector());
  P_V uMeanVec  (m_s->m_unique_u_space.zeroVector());
  P_M uCovMatrix(m_s->m_unique_u_space.zeroVector());
#if 0
  this->predictVUsAtGridPoint(newScenarioVec,
                              vMeanVec,
                              vCovMatrix,
                              uMeanVec,
                              uCovMatrix);
#endif
  simulationOutputMeanVec = newKmat_interp * uMeanVec;
  discrepancyMeanVec = newDmat * vMeanVec;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs(
  const S_V& newScenarioVec,
  const P_V& newParameterVec,
        Q_V& simulationOutputMeanVec)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()"
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(newScenarioVec.sizeLocal() != m_s->m_paper_p_x,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()",
                      "invalid 'newScenarioVec'");

  UQ_FATAL_TEST_MACRO(newParameterVec.sizeLocal() != m_s->m_paper_p_t,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()",
                      "invalid 'newParameterVec'");

  UQ_FATAL_TEST_MACRO(simulationOutputMeanVec.sizeLocal() != m_s->m_paper_n_eta,
                      m_env.worldRank(),
                      "uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()",
                      "invalid 'simulationOutputMeanVec'");

  P_V wMeanVec  (m_s->m_unique_w_space.zeroVector());
  P_M wCovMatrix(m_s->m_unique_w_space.zeroVector());
  this->predictWsAtGridPoint(newScenarioVec,
                             newParameterVec,
                             wMeanVec,
                             wCovMatrix);

  // todo_rr (Should one denormalize qoi here???

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving uqGpmsaComputerModelClass<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()"
                            << std::endl;
  }

  return;
}

#endif // __UQ_GCM_2_H__
