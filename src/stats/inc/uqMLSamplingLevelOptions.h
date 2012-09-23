//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#ifndef __UQ_MULTI_LEVEL_SAMPLING_LEVEL_OPTIONS_H__
#define __UQ_MULTI_LEVEL_SAMPLING_LEVEL_OPTIONS_H__

#define LEVEL_REF_ID 0

#include <uqEnvironment.h>
#include <uqSequenceStatisticalOptions.h>
#define UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
#define UQ_ML_SAMPLING_L_CHECKPOINT_OUTPUT_FILE_NAME_ODV                      UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#endif
#define UQ_ML_SAMPLING_L_STOP_AT_END_ODV                                      0
#define UQ_ML_SAMPLING_L_DATA_OUTPUT_FILE_NAME_ODV                            UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_DATA_OUTPUT_ALLOW_ODV                                ""
#define UQ_ML_SAMPLING_L_LOAD_BALANCE_ALGORITHM_ID_ODV                        2
#define UQ_ML_SAMPLING_L_LOAD_BALANCE_TRESHOLD_ODV                            1.
#define UQ_ML_SAMPLING_L_MIN_EFFECTIVE_SIZE_RATIO_ODV                         0.85
#define UQ_ML_SAMPLING_L_MAX_EFFECTIVE_SIZE_RATIO_ODV                         0.91
#define UQ_ML_SAMPLING_L_SCALE_COV_MATRIX_ODV                                 1
#define UQ_ML_SAMPLING_L_MIN_REJECTION_RATE_ODV                               0.50
#define UQ_ML_SAMPLING_L_MAX_REJECTION_RATE_ODV                               0.75
#define UQ_ML_SAMPLING_L_COV_REJECTION_RATE_ODV                               0.25
#define UQ_ML_SAMPLING_L_TOTALLY_MUTE_ODV                                     1
#define UQ_ML_SAMPLING_L_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV            UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV            UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV                   UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV                   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_RAW_CHAIN_SIZE_ODV                                   100
#define UQ_ML_SAMPLING_L_RAW_CHAIN_GENERATE_EXTRA_ODV                         0
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DISPLAY_PERIOD_ODV                         500
#define UQ_ML_SAMPLING_L_RAW_CHAIN_MEASURE_RUN_TIMES_ODV                      1
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV                     0
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV                  UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV                  UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV                ""
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#define UQ_ML_SAMPLING_L_RAW_CHAIN_COMPUTE_STATS_ODV                          0
#endif
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_GENERATE_ODV                          0
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DISCARDED_PORTION_ODV                 0.
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_LAG_ODV                               1
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV             UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV           ""
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_COMPUTE_STATS_ODV                     0
#endif
#define UQ_ML_SAMPLING_L_DISPLAY_CANDIDATES_ODV                               0
#define UQ_ML_SAMPLING_L_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV                       1
#define UQ_ML_SAMPLING_L_TK_USE_LOCAL_HESSIAN_ODV                             0
#define UQ_ML_SAMPLING_L_TK_USE_NEWTON_COMPONENT_ODV                          1
#define UQ_ML_SAMPLING_L_DR_MAX_NUM_EXTRA_STAGES_ODV                          0
#define UQ_ML_SAMPLING_L_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV               "1."
#define UQ_ML_SAMPLING_L_DR_DURING_AM_NON_ADAPTIVE_INT_ODV                    1
#define UQ_ML_SAMPLING_L_AM_KEEP_INITIAL_MATRIX_ODV                           0
#define UQ_ML_SAMPLING_L_AM_INIT_NON_ADAPT_INT_ODV                            0
#define UQ_ML_SAMPLING_L_AM_ADAPT_INTERVAL_ODV                                0
#define UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV           0
#define UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV        UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV        UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_AM_ETA_ODV                                           1.
#define UQ_ML_SAMPLING_L_AM_EPSILON_ODV                                       1.e-5

class uqMLSamplingLevelOptionsClass
{
public:
  uqMLSamplingLevelOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
//uqMLSamplingLevelOptionsClass(const uqMLSamplingLevelOptionsClass& inputOptions);
 ~uqMLSamplingLevelOptionsClass();

  const uqBaseEnvironmentClass& env() const;
//void changePrefix     (const char* prefix);
  void scanOptionsValues(const uqMLSamplingLevelOptionsClass* defaultOptions);
  void print            (std::ostream& os) const;

  std::string                        m_prefix;

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  std::string                        m_checkpointOutputFileName;
#endif
  bool                               m_stopAtEnd;
  std::string                        m_dataOutputFileName;
  std::set<unsigned int>             m_dataOutputAllowedSet;
  std::string                        m_str1;
  unsigned int                       m_loadBalanceAlgorithmId;
  double                             m_loadBalanceTreshold;
  double                             m_minEffectiveSizeRatio;
  double                             m_maxEffectiveSizeRatio;
  bool                               m_scaleCovMatrix;
  double                             m_minRejectionRate;
  double                             m_maxRejectionRate;
  double                             m_covRejectionRate;
  bool                               m_totallyMute;
  std::string                        m_initialPositionDataInputFileName;
  std::string                        m_initialPositionDataInputFileType;
  std::string                        m_initialProposalCovMatrixDataInputFileName;
  std::string                        m_initialProposalCovMatrixDataInputFileType;
  std::string                        m_rawChainDataInputFileName;
  std::string                        m_rawChainDataInputFileType;
  unsigned int                       m_rawChainSize;
  bool                               m_rawChainGenerateExtra;
  unsigned int                       m_rawChainDisplayPeriod;
  bool                               m_rawChainMeasureRunTimes;
  unsigned int                       m_rawChainDataOutputPeriod;
  std::string                        m_rawChainDataOutputFileName;
  std::string                        m_rawChainDataOutputFileType;
  std::set<unsigned int>             m_rawChainDataOutputAllowedSet;
  std::string                        m_str2;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_rawChainComputeStats;
  uqSequenceStatisticalOptionsClass* m_rawChainStatisticalOptionsObj;
  bool                               m_rawChainStatOptsInstantiated;
#endif

  bool                               m_filteredChainGenerate;
  double                             m_filteredChainDiscardedPortion; // input or set during run time
  unsigned int                       m_filteredChainLag;              // input or set during run time
  std::string                        m_filteredChainDataOutputFileName;
  std::string                        m_filteredChainDataOutputFileType;
  std::set<unsigned int>             m_filteredChainDataOutputAllowedSet;
  std::string                        m_str3;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_filteredChainComputeStats;
  uqSequenceStatisticalOptionsClass* m_filteredChainStatisticalOptionsObj;
  bool                               m_filteredChainStatOptsInstantiated;
#endif

  bool                               m_displayCandidates;
  bool                               m_putOutOfBoundsInChain;
  bool                               m_tkUseLocalHessian;
  bool                               m_tkUseNewtonComponent;
  unsigned int                       m_drMaxNumExtraStages;
  std::vector<double>                m_drScalesForExtraStages;
  std::string                        m_str4;
  bool                               m_drDuringAmNonAdaptiveInt;
  bool                               m_amKeepInitialMatrix;
  unsigned int                       m_amInitialNonAdaptInterval;
  unsigned int                       m_amAdaptInterval;
  unsigned int                       m_amAdaptedMatricesDataOutputPeriod;
  std::string                        m_amAdaptedMatricesDataOutputFileName;
  std::string                        m_amAdaptedMatricesDataOutputFileType;
  double                             m_amEta;
  double                             m_amEpsilon;

private:
  void   copyOptionsValues(const uqMLSamplingLevelOptionsClass& srcOptions);
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  std::string                   m_option_checkpointOutputFileName;
#endif
  std::string                   m_option_stopAtEnd;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
  std::string                   m_option_loadBalanceAlgorithmId;
  std::string                   m_option_loadBalanceTreshold;
  std::string                   m_option_minEffectiveSizeRatio;
  std::string                   m_option_maxEffectiveSizeRatio;
  std::string                   m_option_scaleCovMatrix;
  std::string                   m_option_minRejectionRate;
  std::string                   m_option_maxRejectionRate;
  std::string                   m_option_covRejectionRate;
  std::string                   m_option_totallyMute;
  std::string                   m_option_initialPosition_dataInputFileName;
  std::string                   m_option_initialPosition_dataInputFileType;
  std::string                   m_option_initialProposalCovMatrix_dataInputFileName;
  std::string                   m_option_initialProposalCovMatrix_dataInputFileType;
  std::string                   m_option_rawChain_dataInputFileName;
  std::string                   m_option_rawChain_dataInputFileType;
  std::string                   m_option_rawChain_size;
  std::string                   m_option_rawChain_generateExtra;
  std::string                   m_option_rawChain_displayPeriod;
  std::string                   m_option_rawChain_measureRunTimes;
  std::string                   m_option_rawChain_dataOutputPeriod;
  std::string                   m_option_rawChain_dataOutputFileName;
  std::string                   m_option_rawChain_dataOutputFileType;
  std::string                   m_option_rawChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_rawChain_computeStats;
#endif

  std::string                   m_option_filteredChain_generate;
  std::string                   m_option_filteredChain_discardedPortion;
  std::string                   m_option_filteredChain_lag;
  std::string                   m_option_filteredChain_dataOutputFileName;
  std::string                   m_option_filteredChain_dataOutputFileType;
  std::string                   m_option_filteredChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_filteredChain_computeStats;
#endif

  std::string                   m_option_displayCandidates;
  std::string                   m_option_putOutOfBoundsInChain;
  std::string                   m_option_tk_useLocalHessian;
  std::string                   m_option_tk_useNewtonComponent;
  std::string                   m_option_dr_maxNumExtraStages;
  std::string                   m_option_dr_listOfScalesForExtraStages;
  std::string                   m_option_dr_duringAmNonAdaptiveInt;
  std::string                   m_option_am_keepInitialMatrix;
  std::string                   m_option_am_initialNonAdaptInterval;
  std::string                   m_option_am_adaptInterval;
  std::string                   m_option_am_adaptedMatrices_dataOutputPeriod;
  std::string                   m_option_am_adaptedMatrices_dataOutputFileName;
  std::string                   m_option_am_adaptedMatrices_dataOutputFileType;
  std::string                   m_option_am_eta;
  std::string                   m_option_am_epsilon;
};

std::ostream& operator<<(std::ostream& os, const uqMLSamplingLevelOptionsClass& obj);
#endif // __UQ_MULTI_LEVEL_SAMPLING_LEVEL_OPTIONS_H__
