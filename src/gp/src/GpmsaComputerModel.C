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

#include <queso/GpmsaComputerModel.h>
#include <queso/GenericScalarFunction.h>
#include <queso/SequentialVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::GpmsaComputerModel(
  const char*                                               prefix,
  const GcmOptionsValues*                            alternativeOptionsValues, // dakota
  const SimulationStorage <S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationStorage,
  const SimulationModel   <S_V,S_M,P_V,P_M,Q_V,Q_M>& simulationModel,
  const ExperimentStorage <S_V,S_M,D_V,D_M>*         experimentStorage,
  const ExperimentModel   <S_V,S_M,D_V,D_M>*         experimentModel,
  const BaseVectorRV      <P_V,P_M>*                 thetaPriorRv)
  // Handle case of no experiments, that is, experiment pointers == NULL? (todo)
  :
  m_env                     (simulationStorage.env()),
  m_optionsObj              (alternativeOptionsValues),
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
  m_like_counter            (MiscUintDebugMessage(0,NULL))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
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
      queso_require_equal_to_msg(simulationModel.numBasis(), 1, "inconsistent numBasis (1)");
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
      queso_require_msg(!((simulationModel.numBasis() != 1) || (experimentModel->numBasis() != 1)), "inconsistent numBasis (2)");
    }
  }
  else {
    queso_error_msg("inconsistent experimental information");
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = "                << prefix
                            << ", m_allOutputsAreScalar = " << m_allOutputsAreScalar
                            << std::endl;
  }

  m_formCMatrix = m_formCMatrix && (m_allOutputsAreScalar == false);

  //********************************************************************************
  // Handle options
  //********************************************************************************
  // If NULL, we create one
  if (m_optionsObj == NULL) {
    GcmOptionsValues * tempOptions = new GcmOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;
  }

  m_formCMatrix = m_formCMatrix && m_optionsObj->m_formCMatrix;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = "        << prefix
                            << ", 'final' m_formCMatrix = " << m_formCMatrix
                            << std::endl;
  }

  //********************************************************************************
  // Open output file // todo_r: is this necessary???
  //********************************************************************************
  if ((m_optionsObj->m_dataOutputFileName                       != "."                                            ) &&
      (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end())) {
    m_env.openOutputFile(m_optionsObj->m_dataOutputFileName,
                         UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                         m_optionsObj->m_dataOutputAllowedSet,
                         false,
                         m_dataOutputFilePtrSet);
    queso_require_msg(m_dataOutputFilePtrSet.ofsVar, "data output file could not be created");
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

  // These old options are deprecated.  We do this to preserve backwards
  // compatibility.
  GpmsaComputerModelOptions * gpmsaComputerModelOptions =
    new GpmsaComputerModelOptions(m_env, prefix, *m_optionsObj);
  m_s = new GcmSimulationInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>(*gpmsaComputerModelOptions,
                                                              m_allOutputsAreScalar, // csri (new GcmSimulationInfo)
                                                              simulationStorage,
                                                              simulationModel);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating m_s"
                            << ", experimentStorage = " << experimentStorage
                            << ", experimentModel = "   << experimentModel
                            << ", thetaPriorRv = "      << thetaPriorRv
                            << std::endl;
  }

  if (m_thereIsExperimentalData) {
    queso_require_equal_to_msg(simulationStorage.scenarioSpace().dimLocal(), experimentStorage->scenarioSpace().dimLocal(), "inconsistent dimension of scenario space");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": about to instantiate m_e"
                              << std::endl;
    }

    m_e = new GcmExperimentInfo<S_V,S_M,D_V,D_M,P_V,P_M>(*gpmsaComputerModelOptions,
                                                                m_allOutputsAreScalar, // csri (new GcmExperimentInfo)
                                                                *experimentStorage,
                                                                *experimentModel,
                                                                *thetaPriorRv);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_e"
                              << std::endl;
    }

    m_j = new GcmJointInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*gpmsaComputerModelOptions,
                                                                   m_allOutputsAreScalar, // csri
                                                                   *m_s,
                                                                   *m_e);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_j"
                              << std::endl;
    }

    if (m_allOutputsAreScalar) {
      m_z = new GcmZInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(m_allOutputsAreScalar, // csri ????
                                                                 *m_s,
                                                                 *m_e);
    }
    else {
      m_z = new GcmZInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(m_formCMatrix,
                                                                 m_allOutputsAreScalar,
                                                                 *m_s,
                                                                 *m_e,
                                                                 *m_j);
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_z"
                              << std::endl;
    }

    m_t = new GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_s,
                                                                   *m_e);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_t"
                              << std::endl;
    }
  }
  else {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": about to instantiate m_z (no experiments)"
                              << std::endl;
    }

    m_z = new GcmZInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(m_formCMatrix,
                                                               m_allOutputsAreScalar,
                                                               *m_s);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_z (no experiments)"
                              << std::endl;
    }

    m_t = new GcmTotalInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*m_s);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished instantiating m_t (no experiments)"
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
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
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
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
      if (m_optionsObj->m_useTildeLogicForRankDefficientC) {
        //********************************************************************************
        // Use tilde logic
        //********************************************************************************
        if (m_thereIsExperimentalData) {
          //******************************************************************************
          // Tilde situation: form 'm_Bmat_tilde'
          // Tilde situation: form 'm_vu_tilde_space'
          // Tilde situation: form 'm_Lbmat'
          //******************************************************************************
          m_jt = new GcmJointTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*gpmsaComputerModelOptions,*m_e,*m_j);

          //******************************************************************************
          // Tilde situation: form 'm_Kmat_tilde'
          // Tilde situation: form 'm_w_tilde_space'
          // Tilde situation: form 'm_Lkmat'
          //******************************************************************************
           m_st = new GcmSimulationTildeInfo<S_V,S_M,P_V,P_M,Q_V,Q_M>(*gpmsaComputerModelOptions,*m_s);

          //******************************************************************************
          // Tilde situation: form 'm_Cmat_tilde'
          //******************************************************************************
          m_zt = new GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*gpmsaComputerModelOptions,*m_j,*m_z,*m_st,*m_jt);
        }
        else {
          queso_error_msg("incomplete code for the situation 'm_useTildeLogicForRankDefficientC == true' and 'm_thereIsExperimentalData == false'");
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
          m_zt = new GcmZTildeInfo<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>(*gpmsaComputerModelOptions,*m_j,*m_z);
        }
        else {
          queso_error_msg("incomplete code for the situation 'm_useTildeLogicForRankDefficientC == false' and 'm_thereIsExperimentalData == false'");
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
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
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

      //                    m_env.worldRank(),
      //                    "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()",
      //                    "incomplete code for the situation 'm_thereIsExperimentalData == false'");
    }
  }

  this->memoryCheck(1);

  //********************************************************************************
  // Generate prior sequence
  //********************************************************************************
  if (m_optionsObj->m_priorSeqNumSamples > 0) {
    this->generatePriorSeq();
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished generating prior sequence"
                            << std::endl;
  }

  this->memoryCheck(2);

  //********************************************************************************
  // Instantiate likelihood function object
  //********************************************************************************
  m_likelihoodFunction = new GenericScalarFunction<P_V,P_M>
                           ("like_",
                            m_t->m_totalDomain,
                            &staticLikelihoodRoutine,
                            (void *) this,
                            true); // routine computes [ln(function)]

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished instantiating likelihood Function"
                            << std::endl;
  }

  this->memoryCheck(3);

  // Done with the old options now, so deallocate
  delete gpmsaComputerModelOptions;

  //********************************************************************************
  // Leave constructor
  //********************************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::~GpmsaComputerModel()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::destructor()..."
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
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::destructor()"
                            << std::endl;
  }
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings(
  const MhOptionsValues* alternativeOptionsValues, // dakota
  const P_V&                    totalInitialValues,
  const P_M*                    totalInitialProposalCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << ", totalInitialValues = "             << totalInitialValues
                            << std::endl;
  }

  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",1,3000000);

  queso_require_equal_to_msg(m_t->m_totalPriorRv.imageSet().vectorSpace().dimLocal(), totalInitialValues.sizeLocal(), "'m_totalPriorRv' and 'totalInitialValues' should have equal dimensions");

  if (totalInitialProposalCovMatrix) {
    queso_require_equal_to_msg(m_t->m_totalPriorRv.imageSet().vectorSpace().dimLocal(), totalInitialProposalCovMatrix->numRowsLocal(), "'m_totalPriorRv' and 'totalInitialProposalCovMatrix' should have equal dimensions");
    queso_require_equal_to_msg(totalInitialProposalCovMatrix->numCols(), totalInitialProposalCovMatrix->numRowsGlobal(), "'totalInitialProposalCovMatrix' should be a square matrix");
  }

#if 0
  P_V currPosition (totalInitialValues);
  currPosition.cwSet(0.);
  P_V epsilonVector             (currPosition);
  P_V plusVectorOfLnLikelihoods (currPosition);
  P_V minusVectorOfLnLikelihoods(currPosition);
  P_V deltaVectorOfLnLikelihoods(currPosition);
  P_V vectorOfLnAbsGrads        (currPosition);

  currPosition = totalInitialValues;
  double referenceValue = staticLikelihoodRoutine(currPosition,
                                                  NULL,
                                                  (void *) this,
                                                  NULL,
                                                  NULL,
                                                  NULL);

  for (unsigned int paramId = 0; paramId < totalInitialValues.sizeLocal(); ++paramId) {
    currPosition = totalInitialValues;
    epsilonVector[paramId] = 1.e-8 * totalInitialValues[paramId];
    if (epsilonVector[paramId] == 0.) epsilonVector[paramId] = 1.e-8;

    currPosition[paramId] = totalInitialValues[paramId] + epsilonVector[paramId];
    plusVectorOfLnLikelihoods[paramId] = staticLikelihoodRoutine(currPosition,
                                                                 NULL,
                                                                 (void *) this,
                                                                 NULL,
                                                                 NULL,
                                                                 NULL);

    currPosition[paramId] = totalInitialValues[paramId] - epsilonVector[paramId];
    minusVectorOfLnLikelihoods[paramId] = staticLikelihoodRoutine(currPosition,
                                                                  NULL,
                                                                  (void *) this,
                                                                  NULL,
                                                                  NULL,
                                                                  NULL);

    deltaVectorOfLnLikelihoods[paramId] = plusVectorOfLnLikelihoods[paramId] - minusVectorOfLnLikelihoods[paramId];
    if (deltaVectorOfLnLikelihoods[paramId] > 0.) {
      vectorOfLnAbsGrads[paramId] =  minusVectorOfLnLikelihoods[paramId] + std::log( std::exp( deltaVectorOfLnLikelihoods[paramId]) - 1. ) - std::log(2.*epsilonVector[paramId]);
    }
    else if (deltaVectorOfLnLikelihoods[paramId] == 0.) {
      vectorOfLnAbsGrads[paramId] = -INFINITY;
    }
    else {
      vectorOfLnAbsGrads[paramId] =  plusVectorOfLnLikelihoods [paramId] + std::log( std::exp(-deltaVectorOfLnLikelihoods[paramId]) - 1. ) - std::log(2.*epsilonVector[paramId]);
    }
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
                            << ": referenceValue = "              << referenceValue
                            << "\n epsilonVector              = " << epsilonVector
                            << "\n plusVectorOfLnLikelihoods  = " << plusVectorOfLnLikelihoods
                            << "\n minusVectorOfLnLikelihoods = " << minusVectorOfLnLikelihoods
                            << "\n deltaVectorOfLnLikelihoods = " << deltaVectorOfLnLikelihoods
                            << "\n vectorOfLnAbsGrads         = " << vectorOfLnAbsGrads
                            << std::endl;
  }
#endif

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 002"
  //          << std::endl;

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_t->m_solutionDomain = InstantiateIntersection(m_t->m_totalPriorRv.pdf().domainSet(),m_likelihoodFunction->domainSet());

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 003"
  //          << std::endl;

  m_t->m_solutionPdf = new BayesianJointPdf<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                            m_t->m_totalPriorRv.pdf(),
                                                            *m_likelihoodFunction,
                                                            1.,
                                                            *(m_t->m_solutionDomain));

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 004"
  //          << std::endl;

  m_t->m_totalPostRv.setPdf(*(m_t->m_solutionPdf));

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 005"
  //          << std::endl;

  // Compute output realizer: Metropolis-Hastings approach
  m_t->m_chain = new SequenceOfVectors<P_V,P_M>(m_t->m_totalPostRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 006"
  //          << std::endl;

  m_t->m_mhSeqGenerator = new MetropolisHastingsSG<P_V,P_M>(m_optionsObj->m_prefix.c_str(), // dakota
                                                                   alternativeOptionsValues,
                                                                   m_t->m_totalPostRv,
                                                                   totalInitialValues,
                                                                   totalInitialProposalCovMatrix);

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
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

  m_t->m_solutionRealizer = new SequentialVectorRealizer<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                         *(m_t->m_chain));

  m_t->m_totalPostRv.setRealizer(*(m_t->m_solutionRealizer));

  //m_env.fullComm().syncPrintDebugMsg("In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings(), code place 1",3,3000000);

  //m_env.fullComm().syncPrintDebugMsg("Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()",1,3000000);
  m_env.fullComm().Barrier();

  double totalTime = MiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()..."
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
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc(
  const MhOptionsValues* alternativeOptionsValues, // dakota
  const P_V&                    totalInitialValues,
  const P_M*                    totalInitialProposalCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << ", totalInitialValues = "             << totalInitialValues
                            << std::endl;
  }

  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()",1,3000000);

  // ppp

  m_env.fullComm().syncPrintDebugMsg("Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()",1,3000000);
  m_env.fullComm().Barrier();

  double totalTime = MiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithLanlMcmc()..."
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
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMLSampling()
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMLSampling()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << ", m_env.subComm().NumProc() = "      << m_env.subComm().NumProc()
                            << ", my subRank = "                     << m_env.subRank()
                            << std::endl;
  }

  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMLSampling()",1,3000000);

  // ppp
  m_t->m_solutionDomain = InstantiateIntersection(m_t->m_totalPriorRv.pdf().domainSet(),m_likelihoodFunction->domainSet());

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 003"
  //          << std::endl;

  m_t->m_solutionPdf = new BayesianJointPdf<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                            m_t->m_totalPriorRv.pdf(),
                                                            *m_likelihoodFunction,
                                                            1.,
                                                            *(m_t->m_solutionDomain));

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 004"
  //          << std::endl;

  m_t->m_totalPostRv.setPdf(*(m_t->m_solutionPdf));

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 005"
  //          << std::endl;

  // Compute output realizer: Metropolis-Hastings approach
  m_t->m_chain = new SequenceOfVectors<P_V,P_M>(m_t->m_totalPostRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 006"
  //          << std::endl;

  m_t->m_mlSampler = new MLSampling<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                    m_t->m_totalPriorRv,
                                                    *m_likelihoodFunction);

  //std::cout << "GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMetropolisHastings()"
  //          << ": passing at point 007"
  //          << std::endl;

  m_t->m_mlSampler->generateSequence(*(m_t->m_chain),NULL,NULL);

  // todo_rr
  // m_totalPostMean
  // m_totalPostMedian
  // m_totalPostMode
  // m_totalPostMaxLnValue
  // m_totalMLE
  // m_totalLikeMaxLnValue

  m_t->m_solutionRealizer = new SequentialVectorRealizer<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                         *(m_t->m_chain));

  m_t->m_totalPostRv.setRealizer(*(m_t->m_solutionRealizer));

  m_env.fullComm().syncPrintDebugMsg("Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMLSampling()",1,3000000);
  m_env.fullComm().Barrier();

  double totalTime = MiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::calibrateWithBayesMLSampling()..."
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
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(
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
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                            << ", m_predVU_counter = " << m_j->m_predVU_counter
                            << ", newScenarioVec = "   << newScenarioVec
                            << ", newParameterVec = "  << newParameterVec
                            << std::endl;
  }

  queso_require_equal_to_msg(newScenarioVec.sizeLocal(), m_s->m_paper_p_x, "invalid 'newScenarioVec'");

  queso_require_equal_to_msg(newParameterVec.sizeLocal(), m_s->m_paper_p_t, "invalid 'newParameterVec'");

  queso_require_equal_to_msg(vuMeanVec.sizeLocal(), (m_e->m_paper_p_delta+m_s->m_paper_p_eta), "invalid 'vuMeanVec'");

  queso_require_equal_to_msg(vuCovMatrix.numRowsLocal(), (m_e->m_paper_p_delta + m_s->m_paper_p_eta), "invalid 'vuCovMatrix.numRowsLocal()'");

  queso_require_equal_to_msg(vuCovMatrix.numCols(), (m_e->m_paper_p_delta + m_s->m_paper_p_eta), "invalid 'vuCovMatrix.numCols()'");

  queso_require_equal_to_msg(vMeanVec.sizeLocal(), m_e->m_paper_p_delta, "invalid 'vMeanVec'");

  queso_require_equal_to_msg(vCovMatrix.numRowsLocal(), m_e->m_paper_p_delta, "invalid 'vCovMatrix.numRowsLocal()'");

  queso_require_equal_to_msg(vCovMatrix.numCols(), m_e->m_paper_p_delta, "invalid 'vCovMatrix.numCols()'");

  queso_require_equal_to_msg(uMeanVec.sizeLocal(), m_s->m_paper_p_eta, "invalid 'uMeanVec'");

  queso_require_equal_to_msg(uCovMatrix.numRowsLocal(), m_s->m_paper_p_eta, "invalid 'uCovMatrix.numRowsLocal()'");

  queso_require_equal_to_msg(uCovMatrix.numCols(), m_s->m_paper_p_eta, "invalid 'uCovMatrix.numCols()'");

  if (m_optionsObj->m_predVUsBySamplingRVs) {
  }

  if (m_optionsObj->m_predVUsBySummingRVs) {
    unsigned int numSamples = (unsigned int) ((double) m_t->m_totalPostRv.realizer().subPeriod())/((double) m_optionsObj->m_predLag);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                              << ": m_t->m_totalPostRv.realizer().subPeriod() = " << m_t->m_totalPostRv.realizer().subPeriod()
                              << ", m_optionsObj->m_predLag = "              << m_optionsObj->m_predLag
                              << std::endl;
    }

    SequenceOfVectors<P_V,P_M> unique_vu_means(m_j->m_unique_vu_space,numSamples,m_optionsObj->m_prefix+"vu_means");
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
        for (unsigned int i = 1; i < m_optionsObj->m_predLag; ++i) { // Yes, '1'
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
      queso_require_equal_to_msg(currPosition, totalSample.sizeLocal(), "'currPosition' and 'totalSample.sizeLocal()' should be equal");

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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                                << ", m_predVU_counter = " << m_j->m_predVU_counter
                                << ": about to populate 'm_Smat_v_hat_v_asterisk'"
                                << ", m_e->m_Smat_v_hat_v_asterisk_is.size() = " << m_e->m_Smat_v_hat_v_asterisk_is.size() // 13
                                << ", m_e->m_tmp_rho_v_vec.sizeLocal() = "       << m_e->m_tmp_rho_v_vec.sizeLocal()       //  1
                                << ", m_e->m_tmp_7rhoVVec.sizeLocal() = "        << m_e->m_tmp_7rhoVVec.sizeLocal()        //  1
                                << std::endl;
      }
      queso_require_equal_to_msg((m_e->m_Smat_v_hat_v_asterisk_is.size() * m_e->m_tmp_rho_v_vec.sizeLocal()), m_e->m_tmp_7rhoVVec.sizeLocal(), "invalid size for 'v' variables");
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
                                << ", m_predVU_counter = " << m_j->m_predVU_counter
                                << ": finished instantiating '>m_Smat_z_hat_u_asterisk'"
                                << std::endl;
      }

      //********************************************************************************
      // Submatrix (2,2): Compute '\Sigma_v_asterisk_v_asterisk' matrix
      //********************************************************************************
      queso_require_equal_to_msg(m_e->m_Smat_v_asterisk_v_asterisk.numRowsLocal(), m_e->m_paper_p_delta, "invalid 'm_Smat_v_asterisk_v_asterisk.numRowsLocal()'");
      queso_require_equal_to_msg(m_e->m_tmp_6lambdaVVec.sizeLocal(), m_e->m_paper_p_delta, "invalid 'm_tmp_6lambdaVVec.sizeLocal()'");

      m_e->m_Smat_v_asterisk_v_asterisk.cwSet(0.);
      for (unsigned int i = 0; i < m_e->m_paper_p_delta; ++i) {
        m_e->m_Smat_v_asterisk_v_asterisk(i,i) = 1./m_e->m_tmp_6lambdaVVec[i];
      }

      //********************************************************************************
      // Submatrix (3,3): Compute '\Sigma_u_asterisk_u_asterisk' matrix
      //********************************************************************************
      queso_require_equal_to_msg(m_j->m_Smat_u_asterisk_u_asterisk.numRowsLocal(), m_s->m_paper_p_eta, "invalid 'm_Smat_u_asterisk_u_asterisk.numRowsLocal()'");
      queso_require_equal_to_msg(m_s->m_tmp_2lambdaWVec.sizeLocal(), m_s->m_paper_p_eta, "invalid 'm_tmp_2lambdaWVec.sizeLocal()'");

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
      ComputeConditionalGaussianVectorRV(muVec1,muVec2,sigmaMat11,sigmaMat12,sigmaMat21,sigmaMat22,m_z->m_Zvec_hat,unique_vu_vec,unique_vu_mat);

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
    ComputeCovCorrMatricesBetweenVectorSequences(unique_vu_means,
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
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
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

  if (m_optionsObj->m_predVUsAtKeyPoints) {
  }

  double totalTime = MiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictVUsAtGridPoint(1)"
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
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint(
  const S_V& newScenarioVec,
  const P_V& newParameterVec,
  const P_V* forcingSampleVecForDebug, // Usually NULL
        P_V& wMeanVec,
        P_M& wCovMatrix)
{
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                            << ", m_predW_counter = " << m_s->m_predW_counter
                            << ", newScenarioVec = "  << newScenarioVec
                            << ", newParameterVec = " << newParameterVec
                            << std::endl;
  }

  queso_require_equal_to_msg(newScenarioVec.sizeLocal(), m_s->m_paper_p_x, "invalid 'newScenarioVec'");

  queso_require_equal_to_msg(newParameterVec.sizeLocal(), m_s->m_paper_p_t, "invalid 'newParameterVec'");

  queso_require_equal_to_msg(wMeanVec.sizeLocal(), m_s->m_paper_p_eta, "invalid 'wMeanVec'");

  queso_require_equal_to_msg(wCovMatrix.numRowsLocal(), m_s->m_paper_p_eta, "invalid 'wCovMatrix.numRowsLocal()'");

  queso_require_equal_to_msg(wCovMatrix.numCols(), m_s->m_paper_p_eta, "invalid 'wCovMatrix.numCols()'");

  if (m_optionsObj->m_predWsBySamplingRVs) {
  }

  if (m_optionsObj->m_predWsBySummingRVs) {
    unsigned int numSamples = (unsigned int) ((double) m_t->m_totalPostRv.realizer().subPeriod())/((double) m_optionsObj->m_predLag);
    if (forcingSampleVecForDebug) {
      numSamples = 1;
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                              << ": m_t->m_totalPostRv.realizer().subPeriod() = " << m_t->m_totalPostRv.realizer().subPeriod()
                              << ", m_optionsObj->m_predLag = "              << m_optionsObj->m_predLag
                              << ", numSamples = "                                << numSamples
                              << std::endl;
    }

    SequenceOfVectors<P_V,P_M> unique_w_means(m_s->m_unique_w_space,numSamples,m_optionsObj->m_prefix+"w_means");

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
        for (unsigned int i = 1; i < m_optionsObj->m_predLag; ++i) { // Yes, '1'
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
      queso_require_equal_to_msg(currPosition, totalSample.sizeLocal(), "'currPosition' and 'totalSample.sizeLocal()' should be equal");

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
          *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                << ", m_predW_counter = " << m_s->m_predW_counter
                                << ": finished instantiating 'm_Smat_w_hat_w_asterisk'"
                                << std::endl;
      }

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                                  << ": m_s->m_Smat_w_hat_w_asterisk = " << m_s->m_Smat_w_hat_w_asterisk
                                  << std::endl;
        }
      }

      //********************************************************************************
      // Submatrix (2,2): Compute '\Sigma_w_asterisk_w_asterisk' matrix
      //********************************************************************************
      queso_require_equal_to_msg(m_s->m_Smat_w_asterisk_w_asterisk.numRowsLocal(), m_s->m_paper_p_eta, "invalid 'm_Smat_w_asterisk_w_asterisk.numRowsLocal()'");
      queso_require_equal_to_msg(m_s->m_tmp_2lambdaWVec.sizeLocal(), m_s->m_paper_p_eta, "invalid 'm_tmp_2lambdaWVec.sizeLocal()'");

      m_s->m_Smat_w_asterisk_w_asterisk.cwSet(0.);
      for (unsigned int i = 0; i < m_s->m_paper_p_eta; ++i) {
        m_s->m_Smat_w_asterisk_w_asterisk(i,i) = 1./m_s->m_tmp_2lambdaWVec[i] + 1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
      }

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
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
          *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
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
      ComputeConditionalGaussianVectorRV(muVec1,muVec2,sigmaMat11,sigmaMat12,sigmaMat21,sigmaMat22,m_s->m_Zvec_hat_w,unique_w_vec,unique_w_mat);

      if (forcingSampleVecForDebug) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
          *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
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
    ComputeCovCorrMatricesBetweenVectorSequences(unique_w_means,
                                                   unique_w_means,
                                                   unique_w_means.subSequenceSize(),
                                                   m_s->m_predW_summingRVs_covMatrix_of_unique_w_means,
                                                   m_s->m_predW_summingRVs_corrMatrix_of_unique_w_means);

    wMeanVec   = m_s->m_predW_summingRVs_unique_w_meanVec;
    wCovMatrix = m_s->m_predW_summingRVs_mean_of_unique_w_covMatrices + m_s->m_predW_summingRVs_covMatrix_of_unique_w_means;

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 4)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
                              << ", m_predW_counter = " << m_s->m_predW_counter
                              << ": finished computing all means and covariances"
                              << "\n  wMeanVec = "                                        << wMeanVec
                              << "\n  m_predW_summingRVs_covMatrix_of_unique_w_means = "  << m_s->m_predW_summingRVs_covMatrix_of_unique_w_means
                              << "\n  m_predW_summingRVs_mean_of_unique_w_covMatrices = " << m_s->m_predW_summingRVs_mean_of_unique_w_covMatrices
                              << "\n  wCovMatrix = "                                      << wCovMatrix
                              << std::endl;

    }
  }

  if (m_optionsObj->m_predWsAtKeyPoints) {
  }

  double totalTime = MiscGetEllapsedSeconds(&timevalBegin);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictWsAtGridPoint()"
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
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults(
  const S_V& newScenarioVec,
  const D_M& newKmat_interp,
  const D_M& newDmat,
        D_V& simulationOutputMeanVec, // todo_rr: pass as pointer
        D_V& discrepancyMeanVec)      // todo_rr: pass as pointer
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()"
                            << std::endl;
  }

  queso_require_equal_to_msg(newScenarioVec.sizeLocal(), m_s->m_paper_p_x, "invalid 'newScenarioVec'");

  queso_require_msg(!((newKmat_interp.numRowsLocal() != m_s->m_paper_n_eta) || (newKmat_interp.numCols() != m_s->m_paper_p_eta)), "invalid 'newKmat_interp'");

  queso_require_msg(!((newDmat.numRowsLocal() != m_s->m_paper_n_eta) || (newDmat.numCols() != m_e->m_paper_p_delta)), "invalid 'newDmat'");

  queso_require_equal_to_msg(simulationOutputMeanVec.sizeLocal(), m_s->m_paper_n_eta, "invalid 'simulationOutputMeanVec'");

  queso_require_equal_to_msg(discrepancyMeanVec.sizeLocal(), m_s->m_paper_n_eta, "invalid 'discrepancyMeanVec'");

  P_V vMeanVec  (m_e->m_unique_v_space.zeroVector());
  P_M vCovMatrix(m_e->m_unique_v_space.zeroVector());

  // Damon: Had to change m_unique_u_space to m_unique_w_space below.  I hope
  // it's right.  m_s doesn't have a m_unique_u_space
  P_V uMeanVec  (m_s->m_unique_w_space.zeroVector());
  P_M uCovMatrix(m_s->m_unique_w_space.zeroVector());
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
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictExperimentResults()"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs(
  const S_V& newScenarioVec,
  const P_V& newParameterVec,
        Q_V& simulationOutputMeanVec)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()"
                            << std::endl;
  }

  queso_require_equal_to_msg(newScenarioVec.sizeLocal(), m_s->m_paper_p_x, "invalid 'newScenarioVec'");

  queso_require_equal_to_msg(newParameterVec.sizeLocal(), m_s->m_paper_p_t, "invalid 'newParameterVec'");

  queso_require_equal_to_msg(simulationOutputMeanVec.sizeLocal(), m_s->m_paper_n_eta, "invalid 'simulationOutputMeanVec'");

  P_V wMeanVec  (m_s->m_unique_w_space.zeroVector());
  P_M wCovMatrix(m_s->m_unique_w_space.zeroVector());
  this->predictWsAtGridPoint(newScenarioVec,
                             newParameterVec,
                             NULL,  // Damon: Adding NULL here; Ernesto left this out
                             wMeanVec,
                             wCovMatrix);

  // todo_rr (Should one denormalize qoi here???

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::predictSimulationOutputs()"
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const VectorSpace<P_V,P_M>&
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalSpace() const
{
  queso_require_msg(m_t, "m_t is NULL");
  return m_t->m_totalSpace;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const VectorSpace<P_V,P_M>&
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::unique_vu_space() const
{
  queso_require_msg(m_j, "m_j is NULL");
  return m_j->m_unique_vu_space;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const BaseVectorRV<P_V,P_M>&
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalPriorRv() const
{
  queso_require_msg(m_t, "m_t is NULL");
  return m_t->m_totalPriorRv;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
const GenericVectorRV<P_V,P_M>&
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::totalPostRv() const
{
  queso_require_msg(m_t, "m_t is NULL");
  return m_t->m_totalPostRv;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
void
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::memoryCheck(unsigned int codePositionId)
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
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::generatePriorSeq()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::generatePriorSeq()..."
                            << ": m_optionsObj->m_prefix.c_str() = "          << m_optionsObj->m_prefix.c_str()
                            << ", m_optionsObj->m_priorSeqNumSamples = " << m_optionsObj->m_priorSeqNumSamples
                            << std::endl;
  }

  SequenceOfVectors<P_V,P_M> priorSeq(m_t->m_totalSpace,m_optionsObj->m_priorSeqNumSamples,m_optionsObj->m_prefix+"priorSeq");
  P_V totalSample(m_t->m_totalSpace.zeroVector());
  for (unsigned int sampleId = 0; sampleId < m_optionsObj->m_priorSeqNumSamples; ++sampleId) {
    m_t->m_totalPriorRv.realizer().realization(totalSample);
    priorSeq.setPositionValues(sampleId,totalSample);
  }
  priorSeq.unifiedWriteContents(m_optionsObj->m_priorSeqDataOutputFileName,
                                m_optionsObj->m_priorSeqDataOutputFileType);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::generatePriorSeq()..."
                            << ": m_optionsObj->m_prefix.c_str() = " << m_optionsObj->m_prefix.c_str()
                            << std::endl;
  }

  return;
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
double
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::staticLikelihoodRoutine(
  const P_V&  totalValues,
  const P_V*  totalDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  return ((GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>*) functionDataPtr)->likelihoodRoutine(totalValues,
                                                                                                            totalDirection,
                                                                                                            functionDataPtr,
                                                                                                            gradVector,
                                                                                                            hessianMatrix,
                                                                                                            hessianEffect);
}

template <class S_V,class S_M,class D_V,class D_M,class P_V,class P_M,class Q_V,class Q_M>
double
GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(
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
  //std::cout << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(), m_like_counter = " << m_like_counter << std::endl;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()..."
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
  queso_require_equal_to(m_s->m_1lambdaEtaSpace.dimLocal(), 1);

  queso_require_equal_to(m_s->m_2lambdaWSpace.dimLocal(),
                         m_s->m_paper_p_eta);

  queso_require_equal_to(m_s->m_3rhoWSpace.dimLocal(),
                         (m_s->m_paper_p_eta*(m_s->m_paper_p_x+m_s->m_paper_p_t)));
  queso_require_equal_to(m_s->m_4lambdaSSpace.dimLocal(), m_s->m_paper_p_eta);

  if (m_thereIsExperimentalData) {
    queso_require_equal_to(m_e->m_5lambdaYSpace.dimLocal(), 1);

    queso_require_equal_to(m_e->m_6lambdaVSpace.dimLocal(),
                           m_e->m_paper_F);

    queso_require_equal_to(m_e->m_7rhoVSpace.dimLocal(),
                           m_e->m_paper_F*m_s->m_paper_p_x);
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
  queso_require_equal_to_msg(currPosition, totalValues.sizeLocal(), "'currPosition' and 'totalValues.sizeLocal()' should be equal");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
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
  if ((m_optionsObj->m_checkAgainstPreviousSample) &&
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
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
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
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                              << ", m_like_counter = " << m_like_counter
                              << ": going through true 'm_cMatIsRankDefficient' case"
                              << std::endl;
    }

    this->memoryCheck(57);

    if (m_optionsObj->m_useTildeLogicForRankDefficientC) {
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
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
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue2 = " << tmpValue2
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue2;

      this->memoryCheck(61);

      double tmpValue3 = m_jt->m_a_y_modifier_tilde*std::log(m_e->m_tmp_5lambdaYVec[0]) - m_jt->m_b_y_modifier_tilde*m_e->m_tmp_5lambdaYVec[0];
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3/*99*/)) {
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue3 = " << tmpValue3
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue3;

      this->memoryCheck(62);
    }
    else { // if (m_optionsObj->m_useTildeLogicForRankDefficientC)
      queso_error_msg("incomplete code for situation 'm_useTildeLogicForRankDefficientC == false'");
    }
  }
  else { // if (m_formCMatrix) && (m_cMatIsRankDefficient)
    //********************************************************************************
    // 'm_Cmat' (i) does not exist or (ii) is full rank
    //********************************************************************************
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
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
      queso_error_msg("incomplete code for situation 'm_thereIsExperimentalData == false'");

      this->formSigma_z_hat(m_s->m_tmp_1lambdaEtaVec,
                            m_s->m_tmp_2lambdaWVec,
                            m_s->m_tmp_3rhoWVec,
                            m_s->m_tmp_4lambdaSVec,
                            m_like_counter);
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                              << ", m_like_counter = "                         << m_like_counter
                              << ": finished computing 'm_tmp_Smat_z_hat' =\n" << m_z->m_tmp_Smat_z_hat
                              << std::endl;
    }

    this->memoryCheck(52);

    //********************************************************************************
    // Compute the determinant of '\Sigma_z_hat' matrix
    //********************************************************************************
    double Smat_z_hat_lnDeterminant = m_z->m_tmp_Smat_z_hat.lnDeterminant();

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
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
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
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
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
                                << ", m_like_counter = "                << m_like_counter
                                << ": finished computing 'tmpValue2 = " << tmpValue2
                                << std::endl;
      }
      lnLikelihoodValue += tmpValue2;

      this->memoryCheck(55);

      double tmpValue3 = m_j->m_a_y_modifier*std::log(m_e->m_tmp_5lambdaYVec[0]) - m_j->m_b_y_modifier*m_e->m_tmp_5lambdaYVec[0];
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(non-tilde)"
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
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
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
    *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
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

  double totalTime = MiscGetEllapsedSeconds(&timevalBegin);

  //std::cout << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine(), m_like_counter = " << m_like_counter << std::endl;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
                            << ": m_like_counter = "    << m_like_counter
                            << ", totalValues = "       << totalValues
                            << ", lnLikelihoodValue = " << lnLikelihoodValue
                            << " after "                << totalTime
                            << " seconds"
                            << std::endl;
  }

  if (m_env.subRank() == 0) {
#if 0
    std::cout << "Leaving GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::likelihoodRoutine()"
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
    sleep(1);
    queso_error_msg("Exiting in likelihoodRoutine(), on purpose...");
  }

  return lnLikelihoodValue;
}

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
                            << ", m_optionsObj->m_useTildeLogicForRankDefficientC = " << m_optionsObj->m_useTildeLogicForRankDefficientC
                            << ", m_Cmat = "                                               << m_z->m_Cmat
                            << ", m_cMatIsRankDefficient = "                               << m_cMatIsRankDefficient
                            << std::endl;
  }

  queso_require_msg(!m_optionsObj->m_useTildeLogicForRankDefficientC,
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

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());
  if (outerCounter == 1) {
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
      m_z->m_tmp_Smat_z.subWriteContents("Sigma_z",
          "mat_Sigma_z",
          "m",
          tmpSet);
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

  if (outerCounter == 1) {
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
      m_z->m_tmp_Smat_extra.subWriteContents("Sigma_extra",
          "mat_Sigma_extra",
          "m",
          tmpSet);
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

  if (outerCounter == 1) {
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
      m_z->m_tmp_Smat_z_hat.subWriteContents("Sigma_z_hat",
          "mat_Sigma_z_hat",
          "m",
          tmpSet);
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
                            << ", m_optionsObj->m_useTildeLogicForRankDefficientC = " << m_optionsObj->m_useTildeLogicForRankDefficientC
                            << ", m_Cmat = "                                               << m_z->m_Cmat
                            << ", m_cMatIsRankDefficient = "                               << m_cMatIsRankDefficient
                            << std::endl;
  }

  queso_require_msg(!m_optionsObj->m_useTildeLogicForRankDefficientC,
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
                            << ", m_optionsObj->m_useTildeLogicForRankDefficientC = " << m_optionsObj->m_useTildeLogicForRankDefficientC
                            << ", m_Cmat = "                                               << m_z->m_Cmat
                            << ", m_cMatIsRankDefficient = "                               << m_cMatIsRankDefficient
                            << std::endl;
  }

  queso_require_msg(m_formCMatrix, "'m_Cmat' should have been requested");

  queso_require_msg(m_z->m_Cmat, "'m_Cmat' should have been formed");

  queso_require_msg(m_cMatIsRankDefficient, "'m_Cmat' should be rank defficient");

  queso_require_msg(m_optionsObj->m_useTildeLogicForRankDefficientC, "'m_useTildeLogicForRankDefficientC' should be 'true'");

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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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

    if ((outerCounter == 1) &&
        (i            == 0)) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                                << ", outerCounter = " << outerCounter
                                << ": during the calculation of 'm_Smat_w'"
                                << "\n  m_s->m_paper_xs_asterisks_standard.size() = " << m_s->m_paper_xs_asterisks_standard.size();
        for (unsigned int tmpId = 0; tmpId < m_s->m_paper_xs_asterisks_standard.size(); ++tmpId) {
          *m_env.subDisplayFile() << "\n  m_s->m_paper_xs_asterisks_standard[" << tmpId << "] = "
                                  << *(m_s->m_paper_xs_asterisks_standard[tmpId])
                                  << std::endl;
        }
        *m_env.subDisplayFile() << "\n  m_s->m_paper_ts_asterisks_standard.size() = " << m_s->m_paper_ts_asterisks_standard.size();
        for (unsigned int tmpId = 0; tmpId < m_s->m_paper_ts_asterisks_standard.size(); ++tmpId) {
          *m_env.subDisplayFile() << "\n  m_s->m_paper_ts_asterisks_standard[" << tmpId << "] = "
                                  << *(m_s->m_paper_ts_asterisks_standard[tmpId])
                                  << std::endl;
        }
        *m_env.subDisplayFile() << "\n  m_s->m_tmp_rho_w_vec                      = " << m_s->m_tmp_rho_w_vec
                                << "\n  (*(m_s->m_Rmat_w_is[i=0]))(0,0)           = " << (*(m_s->m_Rmat_w_is[i]))(0,0)
                                << "\n  input_2lambdaWVec[i=0]                    = " << input_2lambdaWVec[i]
                                << "\n  (*(m_s->m_Smat_w_is[i=0]))(0,0)           = " << (*(m_s->m_Smat_w_is[i]))(0,0)
                                << "\n  m_s->m_tmp_4lambdaSVec[i=0]               = " << m_s->m_tmp_4lambdaSVec[i]
                                << std::endl;
      }
    }

    for (unsigned int j = 0; j < m_s->m_Smat_w_is[i]->numRowsLocal(); ++j) {
      (*(m_s->m_Smat_w_is[i]))(j,j) +=  1/m_s->m_tmp_4lambdaSVec[i]; // lambda_s
    }

    if ((outerCounter == 1) &&
        (i            == 0)) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
        *m_env.subDisplayFile() << "In GpmsaComputerModel<S_V,S_M,D_V,D_M,P_V,P_M,Q_V,Q_M>::formSigma_z(1)"
                                << ", outerCounter = " << outerCounter
                                << ": during the calculation of 'm_Smat_w'"
                                << "\n  (*(m_s->m_Smat_w_is[i=0]))(0,0) = " << (*(m_s->m_Smat_w_is[i]))(0,0)
                                << std::endl;
      }
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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
    if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
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

  queso_require_equal_to_msg(xVecs.size(), m_e->m_paper_n, "xVecs.size() is wrong");
  for (unsigned int i = 0; i < xVecs.size(); ++i) {
    queso_require_equal_to_msg(xVecs[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs[i]->sizeLocal() is wrong");
  }
  queso_require_msg(!((rho_v_vec.sizeLocal() == 0) || (rho_v_vec.sizeLocal() > m_s->m_paper_p_x)), "rho_v_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_e->m_paper_n, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), m_e->m_paper_n, "Rmat.numCols() is wrong");

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

  queso_require_equal_to_msg(xVecs.size(), m_e->m_paper_n, "xVecs.size() is wrong");
  for (unsigned int i = 0; i < xVecs.size(); ++i) {
    queso_require_equal_to_msg(xVecs[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVec.sizeLocal(), m_s->m_paper_p_t, "tVec.sizeLocal() is wrong");
  queso_require_equal_to_msg(rho_w_vec.sizeLocal(), (m_s->m_paper_p_x+m_s->m_paper_p_t), "rho_w_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_e->m_paper_n, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), m_e->m_paper_n, "Rmat.numCols() is wrong");

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

  queso_require_equal_to_msg(xVecs.size(), m_s->m_paper_m, "xVecs.size() is wrong");
  for (unsigned int i = 0; i < xVecs.size(); ++i) {
    queso_require_equal_to_msg(xVecs[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVecs.size(), m_s->m_paper_m, "tVecs.size() is wrong");
  for (unsigned int i = 0; i < tVecs.size(); ++i) {
    queso_require_equal_to_msg(tVecs[i]->sizeLocal(), m_s->m_paper_p_t, "tVecs[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(rho_w_vec.sizeLocal(), (m_s->m_paper_p_x+m_s->m_paper_p_t), "rho_w_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_s->m_paper_m, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), m_s->m_paper_m, "Rmat.numCols() is wrong");

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

  queso_require_equal_to_msg(xVecs1.size(), m_e->m_paper_n, "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    queso_require_equal_to_msg(xVecs1[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVec1.sizeLocal(), m_s->m_paper_p_t, "tVec1.sizeLocal() is wrong");
  queso_require_equal_to_msg(xVecs2.size(), m_s->m_paper_m, "xVecs2.size() is wrong");
  for (unsigned int i = 0; i < xVecs2.size(); ++i) {
    queso_require_equal_to_msg(xVecs2[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs2[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVecs2.size(), m_s->m_paper_m, "tVecs2.size() is wrong");
  for (unsigned int i = 0; i < tVecs2.size(); ++i) {
    queso_require_equal_to_msg(tVecs2[i]->sizeLocal(), m_s->m_paper_p_t, "tVecs2[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(rho_w_vec.sizeLocal(), (m_s->m_paper_p_x+m_s->m_paper_p_t), "rho_w_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_e->m_paper_n, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), m_s->m_paper_m, "Rmat.numCols() is wrong");

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

  queso_require_equal_to_msg(xVecs1.size(), m_s->m_paper_m, "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    queso_require_equal_to_msg(xVecs1[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVecs1.size(), m_s->m_paper_m, "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    queso_require_equal_to_msg(tVecs1[i]->sizeLocal(), m_s->m_paper_p_t, "tVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(xVec2.sizeLocal(), m_s->m_paper_p_x, "xVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(tVec2.sizeLocal(), m_s->m_paper_p_t, "tVec2.sizeLocal() is wrong");
  queso_require_msg(!((rho_v_vec.sizeLocal() == 0) || (rho_v_vec.sizeLocal() > m_s->m_paper_p_x)), "rho_v_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_e->m_paper_n, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), 1, "Rmat.numCols() is wrong");

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
  queso_require_equal_to_msg(xVecs1.size(), m_s->m_paper_m, "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    queso_require_equal_to_msg(xVecs1[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVecs1.size(), m_s->m_paper_m, "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    queso_require_equal_to_msg(tVecs1[i]->sizeLocal(), m_s->m_paper_p_t, "tVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(xVec2.sizeLocal(), m_s->m_paper_p_x, "xVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(tVec2.sizeLocal(), m_s->m_paper_p_t, "tVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(rho_w_vec.sizeLocal(), (m_s->m_paper_p_x+m_s->m_paper_p_t), "rho_w_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_s->m_paper_m, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), 1, "Rmat.numCols() is wrong");

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
  queso_require_equal_to_msg(xVecs1.size(), m_s->m_paper_m, "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    queso_require_equal_to_msg(xVecs1[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVecs1.size(), m_s->m_paper_m, "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    queso_require_equal_to_msg(tVecs1[i]->sizeLocal(), m_s->m_paper_p_t, "tVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(xVec2.sizeLocal(), m_s->m_paper_p_x, "xVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(tVec2.sizeLocal(), m_s->m_paper_p_t, "tVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(rho_w_vec.sizeLocal(), (m_s->m_paper_p_x+m_s->m_paper_p_t), "rho_w_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_s->m_paper_m, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), 1, "Rmat.numCols() is wrong");

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

  queso_require_equal_to_msg(xVecs1.size(), m_s->m_paper_m, "xVecs1.size() is wrong");
  for (unsigned int i = 0; i < xVecs1.size(); ++i) {
    queso_require_equal_to_msg(xVecs1[i]->sizeLocal(), m_s->m_paper_p_x, "xVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(tVecs1.size(), m_s->m_paper_m, "tVecs1.size() is wrong");
  for (unsigned int i = 0; i < tVecs1.size(); ++i) {
    queso_require_equal_to_msg(tVecs1[i]->sizeLocal(), m_s->m_paper_p_t, "tVecs1[i]->sizeLocal() is wrong");
  }
  queso_require_equal_to_msg(xVec2.sizeLocal(), m_s->m_paper_p_x, "xVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(tVec2.sizeLocal(), m_s->m_paper_p_t, "tVec2.sizeLocal() is wrong");
  queso_require_equal_to_msg(rho_w_vec.sizeLocal(), (m_s->m_paper_p_x+m_s->m_paper_p_t), "rho_w_vec.sizeLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numRowsLocal(), m_s->m_paper_m, "Rmat.numRowsLocal() is wrong");
  queso_require_equal_to_msg(Rmat.numCols(), 1, "Rmat.numCols() is wrong");

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

template class QUESO::GpmsaComputerModel<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
