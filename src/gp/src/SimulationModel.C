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

#include <queso/SimulationModel.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::SimulationModel(
  const char*                                              prefix,
  const SmOptionsValues*                            alternativeOptionsValues, // dakota
  const SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>& simulStorage)
  :
  m_env                     (simulStorage.env()),
  m_optionsObj              (alternativeOptionsValues),
  m_paper_p_x               (simulStorage.scenarioSpace().dimLocal()),
  m_paper_p_t               (simulStorage.parameterSpace().dimLocal()),
  m_paper_m                 (simulStorage.numSimulations()),
  m_paper_n_eta             (simulStorage.outputSpace().dimLocal()),
  m_p_x_space               (m_env, "p_x_", m_paper_p_x, NULL),
  m_xSeq_original           (m_p_x_space, m_paper_m, "simul_xSeq_original"),
  m_xSeq_original_mins      (m_p_x_space.zeroVector()),
  m_xSeq_original_maxs      (m_p_x_space.zeroVector()),
  m_xSeq_original_ranges    (m_p_x_space.zeroVector()),
  m_xSeq_standard           (m_p_x_space, m_paper_m, "simul_xSeq_standard"),
  m_xSeq_standard_mins      (m_p_x_space.zeroVector()),
  m_xSeq_standard_maxs      (m_p_x_space.zeroVector()),
  m_xSeq_standard_ranges    (m_p_x_space.zeroVector()),
  m_xs_asterisks_standard   (simulStorage.xs_asterisks_original().size(),(const S_V*) NULL), // to be deleted on destructor
  m_p_t_space               (m_env, "p_t_", m_paper_p_t, NULL),
  m_tSeq_original           (m_p_t_space, m_paper_m, "simul_tSeq_original"),
  m_tSeq_mins               (m_p_t_space.zeroVector()),
  m_tSeq_maxs               (m_p_t_space.zeroVector()),
  m_tSeq_ranges             (m_p_t_space.zeroVector()),
  m_tSeq_standard           (m_p_t_space, m_paper_m, "simul_tSeq_standard"),
  m_ts_asterisks_standard   (simulStorage.ts_asterisks_original().size(),(const P_V*) NULL), // to be deleted on destructor
  m_n_eta_space             (m_env, "n_eta_", m_paper_n_eta, NULL),
  m_etaSeq_original         (m_n_eta_space, m_paper_m, "simul_etaSeq_original"),
  m_etaSeq_original_mean    (m_n_eta_space.zeroVector()),
  m_etaSeq_original_std     (m_n_eta_space.zeroVector()),
#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
  m_etaSeq_chunkMeans       (simulStorage.chunkSizes().size(),0.),
  m_etaSeq_chunkStds        (simulStorage.chunkSizes().size(),0.),
#else
  m_etaSeq_allMean          (0.),
  m_etaSeq_allStd           (0.),
#endif
  m_etaSeq_transformed      (m_n_eta_space, m_paper_m, "simul_etaSeq_transformed"),
  m_etaSeq_transformed_mean (m_n_eta_space.zeroVector()),
  m_etaSeq_transformed_std  (m_n_eta_space.zeroVector()),
  m_eta_space               (m_env, "m_eta_simul_model", m_paper_m*m_paper_n_eta, NULL),
  m_etaVec_transformed      (m_eta_space.zeroVector()),
  m_etaMat_transformed      (m_env, m_n_eta_space.map(), m_paper_m),
  m_m_space                 (m_env, "m_", m_paper_m, NULL),
  m_m_unitVec               (m_env, m_m_space.map(), 1.),
  m_m_Imat                  (m_m_unitVec),
  m_paper_p_eta             (0),
  m_p_eta_space             (NULL),                       // to be deleted on destructor
  m_Kmat_eta                (NULL),                       // to be deleted on destructor
  m_kvec_is                 (m_paper_p_eta, (Q_V*) NULL), // to be deleted on destructor
  m_Kmat_is                 (m_paper_p_eta, (Q_M*) NULL), // to be deleted on destructor
  m_Kmat                    (NULL)                        // to be deleted on destructor
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = "     << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << "\n  m_paper_p_x = "                << m_paper_p_x
                            << "\n  m_paper_p_t = "                << m_paper_p_t
                            << "\n  m_paper_m = "                  << m_paper_m
                            << "\n  m_paper_n_eta ="               << m_paper_n_eta
                            << "\n  m_paper_p_eta ="               << m_paper_p_eta
                            << "\n  m_eta_space.dimLocal() = "     << m_eta_space.dimLocal()
                            << std::endl;
  }

  // If NULL, we create one
  if (m_optionsObj == NULL) {
    SmOptionsValues * tempOptions = new SmOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;
  }

  // We'll need to remove this later
  m_simulationModelOptions = new SimulationModelOptions(m_env, prefix);


  queso_require_less_equal_msg(m_paper_p_eta, m_paper_m, "'m_paper_p_eta' cannot be bigger than 'm_paper_m'");

  queso_require_less_equal_msg(m_paper_p_eta, m_paper_n_eta, "'m_paper_p_eta' cannot be bigger than 'm_paper_n_eta'");

  queso_require_equal_to_msg(m_paper_m, simulStorage.xs_asterisks_original().size(), "invalid m_paper_m (1)");

  queso_require_equal_to_msg(m_paper_m, simulStorage.ts_asterisks_original().size(), "invalid m_paper_m (2)");

  //***********************************************************************
  // Standardize 'xs_asterisks_original'
  //***********************************************************************
  // Form 'm_xSeq_original', compute 'm_xSeq_original_mins' and 'm_xSeq_original_ranges', form 'm_xSeq_standard'
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_xSeq_original.setPositionValues(i,*(simulStorage.xs_asterisks_original()[i]));
  }
  m_xSeq_original_mins = m_xSeq_original.subMinPlain();
  m_xSeq_original_maxs = m_xSeq_original.subMaxPlain();
  m_xSeq_original_ranges = m_xSeq_original_maxs - m_xSeq_original_mins;
  S_V tmpXVec(m_p_x_space.zeroVector());
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_xSeq_original.getPositionValues(i,tmpXVec);
    tmpXVec -= m_xSeq_original_mins;
    for (unsigned int j = 0; j < m_paper_p_x; ++j) {
      if (m_xSeq_original_ranges[j] == 0.) { // CSRI, 2013-Aug-06, with Laura
        tmpXVec[j] = 0.5;
      }
      else {
        tmpXVec[j] /= m_xSeq_original_ranges[j];
      }
    }
    m_xSeq_standard.setPositionValues(i,tmpXVec);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished setting 'm_xSeq_original'"
                            << ", computing mins of 'm_xSeq_original'"
                            << ", computing ranges of 'm_xSeq_original'"
                            << ", computing 'm_xSeq_standard'"
                            << std::endl;
  }

  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_xSeq_standard.getPositionValues(i,tmpXVec);
    m_xs_asterisks_standard[i] = new S_V(tmpXVec); // to be deleted on destructor
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()";
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_xSeq_original.getPositionValues(i,tmpXVec);
      *m_env.subDisplayFile() << "\n  m_xSeq_original[" << i << "] = " << tmpXVec;
    }
    *m_env.subDisplayFile() << "\n  m_xSeq_original_mins = "   << m_xSeq_original_mins
                            << "\n  m_xSeq_original_maxs = "   << m_xSeq_original_maxs
                            << "\n  m_xSeq_original_ranges = " << m_xSeq_original_ranges
                            << std::endl;
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_xSeq_standard.getPositionValues(i,tmpXVec);
      *m_env.subDisplayFile() << "\n  m_xSeq_standard[" << i << "] = " << tmpXVec;
    }
    m_xSeq_standard_mins = m_xSeq_standard.subMinPlain();
    m_xSeq_standard_maxs = m_xSeq_standard.subMaxPlain();
    m_xSeq_standard_ranges = m_xSeq_standard_maxs - m_xSeq_standard_mins;
    *m_env.subDisplayFile() << "\n  m_xSeq_standard_mins = "   << m_xSeq_standard_mins
                            << "\n  m_xSeq_standard_maxs = "   << m_xSeq_standard_maxs
                            << "\n  m_xSeq_standard_ranges = " << m_xSeq_standard_ranges
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished forming 'm_xs_asterisks_standard'"
                            << std::endl;
  }

  //***********************************************************************
  // Standardize 'ts_asterisks_original'
  //***********************************************************************
  // Form 'm_tSeq_original', compute 'm_tSeq_mins' and 'm_tSeq_ranges', form 'm_tSeq_standard'
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_tSeq_original.setPositionValues(i,*(simulStorage.ts_asterisks_original()[i]));
  }
  m_tSeq_mins = m_tSeq_original.subMinPlain();
  m_tSeq_maxs = m_tSeq_original.subMaxPlain();
  m_tSeq_ranges = m_tSeq_maxs - m_tSeq_mins;
  S_V tmpTVec(m_p_t_space.zeroVector());
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_tSeq_original.getPositionValues(i,tmpTVec);
    tmpTVec -= m_tSeq_mins;
    for (unsigned int j = 0; j < m_paper_p_t; ++j) {
      if (m_tSeq_ranges[j] == 0.) { // CSRI, 2013-Aug-06, with Laura
        tmpTVec[j] = 0.5;
      }
      else {
        tmpTVec[j] /= m_tSeq_ranges[j];
      }
    }
    m_tSeq_standard.setPositionValues(i,tmpTVec);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished setting 'm_tSeq_original'"
                            << ", computing mins of 'm_tSeq_original'"
                            << ", computing ranges of 'm_tSeq_original'"
                            << ", computing 'm_tSeq_standard'"
                            << std::endl;
  }

  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_tSeq_standard.getPositionValues(i,tmpTVec);
    m_ts_asterisks_standard[i] = new S_V(tmpTVec); // to be deleted on destructor
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor():";
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_tSeq_original.getPositionValues(i,tmpTVec);
      *m_env.subDisplayFile() << "\n  m_tSeq_original[" << i << "] = " << tmpTVec;
    }
    *m_env.subDisplayFile() << "\n  m_tSeq_original_mins = "   << m_tSeq_mins
                            << "\n  m_tSeq_original_maxs = "   << m_tSeq_maxs
                            << "\n  m_tSeq_original_ranges = " << m_tSeq_ranges
                            << std::endl;
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_tSeq_standard.getPositionValues(i,tmpTVec);
      *m_env.subDisplayFile() << "\n  m_tSeq_standard[" << i << "] = " << tmpTVec;
    }
    m_tSeq_mins = m_tSeq_standard.subMinPlain();
    m_tSeq_maxs = m_tSeq_standard.subMaxPlain();
    m_tSeq_ranges = m_tSeq_maxs - m_tSeq_mins;
    *m_env.subDisplayFile() << "\n  m_tSeq_standard_mins = "   << m_tSeq_mins
                            << "\n  m_tSeq_standard_maxs = "   << m_tSeq_maxs
                            << "\n  m_tSeq_standard_ranges = " << m_tSeq_ranges
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished forming 'm_ts_asterisks_standard'"
                            << std::endl;
  }

  //***********************************************************************
  // Form 'etaVec_transformed' and 'etaMat_transformed' matrix
  //***********************************************************************
  // Form 'm_etaSeq_original', compute 'm_etaSeq_original_mean' and 'm_etaSeq_original_std', form 'm_etaSeq_transformed'
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_original.setPositionValues(i,simulStorage.outputVec_original(i));
  }
  m_etaSeq_original_mean = m_etaSeq_original.subMeanPlain();
  m_etaSeq_original.subSampleStd(0,m_paper_m,m_etaSeq_original_mean,m_etaSeq_original_std);
  Q_V tmpEtaVec(m_n_eta_space.zeroVector());
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_original.getPositionValues(i,tmpEtaVec);
    tmpEtaVec -= m_etaSeq_original_mean;
    //for (unsigned int j = 0; j < m_paper_n_eta; ++j) {
    //  tmpEtaVec[j] /= m_etaSeq_original_std[j];
    //}
    m_etaSeq_transformed.setPositionValues(i,tmpEtaVec);
  }

#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
  unsigned int cumulativeSize = 0;
  for (unsigned int chunkId = 0; chunkId < m_etaSeq_chunkStds.size(); ++chunkId) {
    m_etaSeq_chunkMeans[chunkId] = 0.;
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
      for (unsigned int j = cumulativeSize; j < (cumulativeSize + simulStorage.chunkSizes()[chunkId]); ++j) {
        m_etaSeq_chunkMeans[chunkId] += tmpEtaVec[j];
      }
    }
    m_etaSeq_chunkMeans[chunkId] /= ((double) (m_paper_m * simulStorage.chunkSizes()[chunkId]));
    cumulativeSize += simulStorage.chunkSizes()[chunkId];
  }
  queso_require_equal_to_msg(cumulativeSize, m_paper_n_eta, "inconsistency in 'cumulativeSize' (1)");

  cumulativeSize = 0;
  for (unsigned int chunkId = 0; chunkId < m_etaSeq_chunkStds.size(); ++chunkId) {
    m_etaSeq_chunkStds[chunkId] = 0.;
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
      for (unsigned int j = cumulativeSize; j < (cumulativeSize + simulStorage.chunkSizes()[chunkId]); ++j) {
        double diff = tmpEtaVec[j] - m_etaSeq_chunkMeans[chunkId];
        m_etaSeq_chunkStds[chunkId] += diff * diff;
      }
    }
    m_etaSeq_chunkStds[chunkId] = sqrt( m_etaSeq_chunkStds[chunkId] / ((double) (m_paper_m * simulStorage.chunkSizes()[chunkId] - 1)) );
    cumulativeSize += simulStorage.chunkSizes()[chunkId];
  }
  queso_require_equal_to_msg(cumulativeSize, m_paper_n_eta, "inconsistency in 'cumulativeSize' (2)");

  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
    cumulativeSize = 0;
    for (unsigned int chunkId = 0; chunkId < m_etaSeq_chunkStds.size(); ++chunkId) {
      for (unsigned int j = cumulativeSize; j < (cumulativeSize + simulStorage.chunkSizes()[chunkId]); ++j) {
        tmpEtaVec[j] /= m_etaSeq_chunkStds[chunkId];
      }
      cumulativeSize += simulStorage.chunkSizes()[chunkId];
    }
    queso_require_equal_to_msg(cumulativeSize, m_paper_n_eta, "inconsistency in 'cumulativeSize' (3)");
    m_etaSeq_transformed.setPositionValues(i,tmpEtaVec);
  }
#else
  m_etaSeq_allMean = 0.;
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
    for (unsigned int j = 0; j < m_paper_n_eta; ++j) {
      m_etaSeq_allMean += tmpEtaVec[j];
    }
  }
  m_etaSeq_allMean /= ((double) (m_paper_m * m_paper_n_eta));

  m_etaSeq_allStd = 0.;
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
    for (unsigned int j = 0; j < m_paper_n_eta; ++j) {
      double diff = tmpEtaVec[j] - m_etaSeq_allMean;
      m_etaSeq_allStd += diff * diff;
    }
  }
  m_etaSeq_allStd = sqrt( m_etaSeq_allStd / ((double) (m_paper_m * m_paper_n_eta - 1)) );

  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
    tmpEtaVec /= m_etaSeq_allStd;
    m_etaSeq_transformed.setPositionValues(i,tmpEtaVec);
  }
#endif

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished setting 'm_etaSeq_original'"
                            << ", computing means of 'm_etaSeq_original'"
                            << ", computing sample stds of 'm_etaSeq_original'"
#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
                            << ", computing 'm_etaSeq_chunkMeans'"
                            << ", computing 'm_etaSeq_chunkStds'"
#else
                            << ", computing 'm_etaSeq_allmean' = " << m_etaSeq_allMean
                            << ", computing 'm_etaSeq_allStd' = "  << m_etaSeq_allStd
#endif
                            << ", computing 'm_etaSeq_transformed'"
                            << std::endl;
  }

#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": 'm_etaSeq_chunkMeans' =";
    for (unsigned int chunkId = 0; chunkId < m_etaSeq_chunkStds.size(); ++chunkId) {
      *m_env.subDisplayFile() << " " << m_etaSeq_chunkMeans[chunkId];
    }
    *m_env.subDisplayFile() << std::endl
                            << ": 'm_etaSeq_chunkStds' =";
    for (unsigned int chunkId = 0; chunkId < m_etaSeq_chunkStds.size(); ++chunkId) {
      *m_env.subDisplayFile() << " " << m_etaSeq_chunkStds[chunkId];
    }
    *m_env.subDisplayFile() << std::endl;
  }
#endif

  // Form 'etaVec_transformed' and 'etaMat_transformed' matrix
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": before forming 'm_etaVec_transformed' and 'm_etaMat_transformed"
                            << "\n  m_paper_p_eta = "                          << m_paper_p_eta
                            << "\n  m_paper_m = "                              << m_paper_m
                            << "\n  m_etaSeq_original.subSequenceSize() = "    << m_etaSeq_original.subSequenceSize()
                            << "\n  m_etaSeq_transformed.subSequenceSize() = " << m_etaSeq_transformed.subSequenceSize()
                            << "\n  m_etaMat_transformed.numRowsLocal() = "    << m_etaMat_transformed.numRowsLocal()
                            << "\n  m_etaMat_transformed.numCols() = "         << m_etaMat_transformed.numCols()
                            << std::endl;
  }

  unsigned int cumulativeEtaVecPosition = 0;
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
    m_etaVec_transformed.cwSet(cumulativeEtaVecPosition,tmpEtaVec);
    cumulativeEtaVecPosition += tmpEtaVec.sizeLocal();
    m_etaMat_transformed.setColumn(i,tmpEtaVec);
  }
  queso_require_equal_to_msg(cumulativeEtaVecPosition, m_etaVec_transformed.sizeLocal(), "inconsistency in the assemble of 'm_etaVec_transformed'");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()";
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_etaSeq_original.getPositionValues(i,tmpEtaVec);
      *m_env.subDisplayFile() << "\n  m_etaSeq_original[" << i << "] = " << tmpEtaVec;
    }
    *m_env.subDisplayFile() << "\n  m_etaSeq_original_mean = " << m_etaSeq_original_mean
                            << "\n  m_etaSeq_original_std = "  << m_etaSeq_original_std
                            << std::endl;
    for (unsigned int i = 0; i < m_paper_m; ++i) {
      m_etaSeq_transformed.getPositionValues(i,tmpEtaVec);
      *m_env.subDisplayFile() << "\n  m_etaSeq_transformed[" << i << "] = " << tmpEtaVec;
    }
    m_etaSeq_transformed_mean = m_etaSeq_transformed.subMeanPlain();
    m_etaSeq_transformed.subSampleStd(0,m_paper_m,m_etaSeq_transformed_mean,m_etaSeq_transformed_std);
    *m_env.subDisplayFile() << "\n  m_etaSeq_transformed_mean = " << m_etaSeq_transformed_mean
                            << "\n  m_etaSeq_transformed_std = "  << m_etaSeq_transformed_std
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    m_etaMat_transformed.setPrintHorizontally(false);
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished forming 'm_etaVec_transformed' and 'm_etaMat_transformed'"
                            << "\n (key-debug) m_etaMat_transformed =\n" << m_etaMat_transformed
                            << std::endl;
  }

  //***********************************************************************
  // Compute SVD of 'etaMat_transformed' matrix
  //***********************************************************************
  if (m_etaMat_transformed.numRowsGlobal() >= m_etaMat_transformed.numCols()) {
    //***********************************************************************
    // Case 'etaMat.numRows >= etaMat.numCols': perform svd on 'etaMat_transformed'
    //***********************************************************************
    Q_M svdU_mat (m_env, m_n_eta_space.map(), m_paper_m);
    Q_V svdS_vec (m_m_space.zeroVector());
    Q_M svdVt_mat(m_m_space.zeroVector());
    //std::cout << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
    //          << ": before m_etaMat_transformed.svd()"
    //          << std::endl;
    int iRC = m_etaMat_transformed.svd(svdU_mat, svdS_vec, svdVt_mat);
    //std::cout << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
    //          << ": after m_etaMat_transformed.svd(), iRC = " << iRC
    //          << std::endl;
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << "\n iRC from m_etaMat_transformed.svd(1) = " << iRC
                              << std::endl;
    }
    queso_require_equal_to_msg(iRC, 0, "svd(1) failed");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished performing svd(1) on 'm_etaMat_transformed'"
                              << "\n svdU_mat  =\n" << svdU_mat
                              << "\n svdS_vec  =\n" << svdS_vec
                              << "\n svdVt_mat =\n" << svdVt_mat
                              << std::endl;
    }

    //***********************************************************************
    // Case 'etaMat.numRows >= etaMat.numCols': determine 'm_paper_p_eta'
    // [U,S,V]=svd(ysimStd,0);
    // lam=diag(S).^2/sum(diag(S).^2);
    // lam=cumsum(lam);
    // pu=sum(lam<pcpct)+1;
    // Ksim=U(:,1:pu)*S(1:pu,1:pu)./sqrt(m);
    //***********************************************************************
    m_paper_p_eta = computePEta(svdS_vec);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ", m_optionsObj->m_cdfThresholdForPEta = "       << m_optionsObj->m_cdfThresholdForPEta
                              << ", m_optionsObj->m_zeroRelativeSingularValue = " << m_optionsObj->m_zeroRelativeSingularValue
                              << ": m_paper_p_eta(1) = "                               << m_paper_p_eta
                              << std::endl;
    }

    m_p_eta_space = new VectorSpace<Q_V,Q_M>(m_env, "p_eta_", m_paper_p_eta, NULL);      // to be deleted on destructor
    m_Kmat_eta    = new Q_M                        (m_env, m_n_eta_space.map(), m_paper_p_eta); // to be deleted on destructor

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", m_Kmat_eta just created (situation 1; not yet populated)"
                              << ", numRowsLocal = " << m_Kmat_eta->numRowsLocal()
                              << ", numCols = "      << m_Kmat_eta->numCols()
                              << std::endl;
    }

    //***********************************************************************
    // Case 'etaMat.numRows >= etaMat.numCols': form 'Kmat_eta' matrix
    //***********************************************************************
    Q_M matU_tmp(m_env, m_n_eta_space.map(), m_paper_p_eta);
    Q_V vecU_tmp(m_n_eta_space.zeroVector());
    for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
      svdU_mat.getColumn(i, vecU_tmp);
      matU_tmp.setColumn(i, vecU_tmp);
    }

    Q_M matS_tmp(m_p_eta_space->zeroVector());
    for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
      matS_tmp(i,i) = svdS_vec[i];
    }

    (*m_Kmat_eta) = matU_tmp * matS_tmp;

     if (m_env.identifyingString() == "CASLexample") {
       // IMPORTANT: temporary, just to match with CASL example of the Matlab version of GPMSA
       if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
         *m_env.subDisplayFile() << "IMPORTANT In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                 << ": multiplying first and third columns of 'm_Kmat_eta' by '-1' in order to match the results of Matlab GPMSA"
                                 << std::endl;
       }
       for (unsigned int i = 0; i < m_Kmat_eta->numRowsLocal(); ++i) {
         (*m_Kmat_eta)(i,0) *= -1.;
       }
       for (unsigned int i = 0; i < m_Kmat_eta->numRowsLocal(); ++i) {
         (*m_Kmat_eta)(i,2) *= -1.;
       }
     }
  }
  else {
    //***********************************************************************
    // Case 'etaMat.numRows < etaMat.numCols': perform svd on 'etaMat_transformed'
    //***********************************************************************
    Q_M svdU_mat (m_env, m_m_space.map(), m_paper_n_eta);
    Q_V svdS_vec (m_n_eta_space.zeroVector());
    Q_M svdVt_mat(m_n_eta_space.zeroVector());
    Q_M etaMat_transfored_transpose(svdU_mat);
    etaMat_transfored_transpose.fillWithTranspose(0,0,m_etaMat_transformed,true,true);
    int iRC = etaMat_transfored_transpose.svd(svdU_mat, svdS_vec, svdVt_mat);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << "\n iRC from m_etaMat_transformed.svd(2) = " << iRC
                              << std::endl;
    }
    queso_require_equal_to_msg(iRC, 0, "svd(2) failed");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": finished performing svd(2) on 'm_etaMat_transformed'"
                              << "\n svdU_mat  =\n" << svdU_mat
                              << "\n svdS_vec  =\n" << svdS_vec
                              << "\n svdVt_mat =\n" << svdVt_mat
                              << std::endl;
    }

    //***********************************************************************
    // Case 'etaMat.numRows < etaMat.numCols': determine 'm_paper_p_eta'
    // [U,S,V]=svd(ysimStd,0);
    // lam=diag(S).^2/sum(diag(S).^2);
    // lam=cumsum(lam);
    // pu=sum(lam<pcpct)+1;
    // Ksim=U(:,1:pu)*S(1:pu,1:pu)./sqrt(m);
    //***********************************************************************
    m_paper_p_eta = computePEta(svdS_vec);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ", m_optionsObj->m_cdfThresholdForPEta = "       << m_optionsObj->m_cdfThresholdForPEta
                              << ", m_optionsObj->m_zeroRelativeSingularValue = " << m_optionsObj->m_zeroRelativeSingularValue
                              << ": m_paper_p_eta(2) = "                               << m_paper_p_eta
                              << std::endl;
    }

    m_p_eta_space = new VectorSpace<Q_V,Q_M>(m_env, "p_eta_", m_paper_p_eta, NULL);      // to be deleted on destructor
    m_Kmat_eta    = new Q_M                        (m_env, m_n_eta_space.map(), m_paper_p_eta); // to be deleted on destructor

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", m_Kmat_eta just created (situation 2; not yet populated)"
                              << ", numRowsLocal = " << m_Kmat_eta->numRowsLocal()
                              << ", numCols = "      << m_Kmat_eta->numCols()
                              << std::endl;
    }

    //***********************************************************************
    // Case 'etaMat.numRows < etaMat.numCols': form 'Kmat_eta' matrix
    //***********************************************************************
    Q_M svdV_mat(svdVt_mat.transpose());
    Q_M matV_tmp(m_env, m_n_eta_space.map(), m_paper_p_eta);
    Q_V vecV_tmp(m_n_eta_space.zeroVector());
    for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
      svdV_mat.getColumn(i, vecV_tmp);
      matV_tmp.setColumn(i, vecV_tmp);
    }

    Q_M matS_tmp(m_p_eta_space->zeroVector());
    for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
      matS_tmp(i,i) = svdS_vec[i];
    }

    (*m_Kmat_eta) = matV_tmp * matS_tmp;

    if (m_env.identifyingString() == "towerExampleToMatchGPMSA") {
      // IMPORTANT: temporary, just to match with tower example of the Matlab version of GPMSA
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
        *m_env.subDisplayFile() << "IMPORTANT In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                                << ": multiplying first column of 'm_Kmat_eta' by '-1' in order to match the results of Matlab GPMSA"
                                << std::endl;
      }
      for (unsigned int i = 0; i < m_Kmat_eta->numRowsLocal(); ++i) {
        (*m_Kmat_eta)(i,0) *= -1.;
      }
    }
  }
  (*m_Kmat_eta) /= sqrt(m_paper_m); // Scaling (check)

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    m_Kmat_eta->setPrintHorizontally(false);
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished forming 'm_Kmat_eta'"
                            << "\n (key-debug) *m_Kmat_eta =\n" << *m_Kmat_eta
                            << std::endl;
  }

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());
  if (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end()) {
    m_Kmat_eta->subWriteContents("Kmat_eta",
                                 "mat_K_eta",
                                 "m",
                                 tmpSet);
  }


  m_kvec_is.resize(m_paper_p_eta, (Q_V*) NULL);                        // to be deleted on destructor
  m_Kmat_is.resize(m_paper_p_eta, (Q_M*) NULL);                        // to be deleted on destructor
  m_Kmat = new Q_M(m_env, m_eta_space.map(), m_paper_m*m_paper_p_eta); // to be deleted on destructor
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": key-debug"
                            << ", m_Kmat just created (not yet populated)"
                            << ", numRowsLocal = " << m_Kmat->numRowsLocal()
                            << ", numCols = "      << m_Kmat->numCols()
                            << std::endl;
  }

  //***********************************************************************
  // Form 'kvec_is'
  //***********************************************************************
  for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
    m_kvec_is[i] = new Q_V(m_n_eta_space.zeroVector());
    m_Kmat_eta->getColumn(i,*(m_kvec_is[i]));
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished computing 'm_kvec_is'"
                            << std::endl;
  }

  //***********************************************************************
  // Compute 'Kmat_is' (via tensor product)
  //***********************************************************************
  for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
    m_Kmat_is[i] = new Q_M(m_env, m_eta_space.map(), m_paper_m);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                              << ": key-debug"
                              << ", before 'm_Kmat_is[i]->fillWithTensorProduct()"
                              << ", m_Kmat_is[" << i << "]->numRowsLocal() = " << m_Kmat_is[i]->numRowsLocal()
                              << ", m_Kmat_is[" << i << "]->numCols() = "      << m_Kmat_is[i]->numCols()
                              << ", m_m_Imat.numRowsLocal() = "                << m_m_Imat.numRowsLocal()
                              << ", m_m_Imat.numCols() = "                     << m_m_Imat.numCols()
                              << ", m_kvec_is[" << i << "]->sizeLocal() = "    << m_kvec_is[i]->sizeLocal()
                              << std::endl;
    }
    m_Kmat_is[i]->fillWithTensorProduct(0,0,m_m_Imat,*(m_kvec_is[i]),true,true);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": finished computing 'm_Kmat_is'"
                            << std::endl;
  }

  //***********************************************************************
  // Form 'Kmat' matrix
  //***********************************************************************
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": before 'm_Kmat->fillWithBlocksDiagonal()"
                            << "\n  m_Kmat->numRowsLocal() = " << m_Kmat->numRowsLocal()
                            << "\n  m_Kmat->numCols() = "      << m_Kmat->numCols()
                            << std::endl;
  }
  m_Kmat->fillWithBlocksHorizontally(0,0,m_Kmat_is,true,true);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::~SimulationModel()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()..."
                            << std::endl;
  }

  delete m_Kmat; // to be deleted on destructor
  for (unsigned int i = 0; i < m_paper_p_eta; ++i) {
    delete m_Kmat_is[i]; // to be deleted on destructor
    delete m_kvec_is[i]; // to be deleted on destructor
  }
  delete m_Kmat_eta;    // to be deleted on destructor
  delete m_p_eta_space; // to be deleted on destructor
  for (unsigned int i = 0; i < m_paper_m; ++i) {
    delete m_xs_asterisks_standard[i]; // to be deleted on destructor
    delete m_ts_asterisks_standard[i]; // to be deleted on destructor
  }
  if (m_optionsObj) delete m_optionsObj;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()"
                            << std::endl;
  }
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
unsigned int
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::computePEta(const Q_V& svdS_vec)
{
  unsigned int result = 0;
  if (m_optionsObj->m_cdfThresholdForPEta <= 0.) {
    result = m_optionsObj->m_p_eta;
  }
  else if (m_optionsObj->m_cdfThresholdForPEta >= 1.) {
    for (unsigned int i = 0; i < svdS_vec.sizeLocal(); ++i) {
      if ((svdS_vec[i]/svdS_vec[0]) > m_optionsObj->m_zeroRelativeSingularValue) {
        result++;
      }
    }
  }
  else {
    Q_V auxVec(svdS_vec);
    double auxSum = 0.;
    for (unsigned int i = 0; i < auxVec.sizeLocal(); ++i) {
      double auxTerm = auxVec[i];
      auxVec[i] = auxTerm * auxTerm;
      auxSum += auxVec[i];
    }
    for (unsigned int i = 0; i < auxVec.sizeLocal(); ++i) {
      auxVec[i] /= auxSum;
    }
    Q_V auxCumSum(m_m_space.zeroVector());
    auxCumSum[0] = auxVec[0];
    for (unsigned int i = 1; i < auxVec.sizeLocal(); ++i) {
      auxCumSum[i] = auxCumSum[i-1] + auxVec[i];
    }
    for (unsigned int i = 0; i < auxVec.sizeLocal(); ++i) {
      if (auxCumSum[i] < m_optionsObj->m_cdfThresholdForPEta) {
        result++;
      }
    }
    result += 1;
  }

  return result;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
unsigned int
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::numBasis() const
{
  return m_paper_p_eta;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const std::vector<const S_V* >&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::xs_asterisks_standard() const
{
  return m_xs_asterisks_standard;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const S_V&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::xSeq_original_mins() const
{
  return m_xSeq_original_mins;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const S_V&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::xSeq_original_ranges() const
{
  return m_xSeq_original_ranges;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const std::vector<const P_V* >&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::ts_asterisks_standard() const
{
  return m_ts_asterisks_standard;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_V&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaSeq_original_mean() const
{
  return m_etaSeq_original_mean;
}

#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
double
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaSeq_chunkStd(unsigned int chunkId) const
{
  queso_require_less_equal_msg(chunkId, m_etaSeq_chunkStds.size(), "'chunkId' is too big");

  return m_etaSeq_chunkStds[chunkId];
}
#else
template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
double
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaSeq_allStd() const
{
  return m_etaSeq_allStd;
}
#endif

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_V&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaVec_transformed(const std::string& debugString) const
{
  //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
  //  *m_env.subDisplayFile() << "Entering SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaVec_transformed()"
  //                          << ": debugString = " << debugString
  //                          << ", m_etaVec_transformed = " << m_etaVec_transformed
  //                          << std::endl;
  //}

  return m_etaVec_transformed;
}
template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_V&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::basisVec(unsigned int basisId) const
{
  return *(m_kvec_is[basisId]);
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_M&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::Kmat_eta() const
{
  return *m_Kmat_eta;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_M&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::Kmat() const
{
  return *m_Kmat;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const SimulationModelOptions&
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::optionsObj() const
{
  // This method returns the old SimulationModelOptions object.  We're using
  // the new SmOptionsValues so this method is deprecated.
  queso_deprecated();
  return *m_simulationModelOptions;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
void
SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO

template class QUESO::SimulationModel<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
