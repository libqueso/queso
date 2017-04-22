//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef UQ_MULTI_LEVEL_SAMPLING_LEVEL_OPTIONS_H
#define UQ_MULTI_LEVEL_SAMPLING_LEVEL_OPTIONS_H

#define LEVEL_REF_ID 0

#include <queso/Environment.h>
#include <queso/SequenceStatisticalOptions.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#else
#include <queso/getpot.h>
#endif

#define UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_ML_SAMPLING_L_HELP                                                 ""
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
#define UQ_ML_SAMPLING_L_CHECKPOINT_OUTPUT_FILE_NAME_ODV                      UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#endif
#define UQ_ML_SAMPLING_L_STOP_AT_END_ODV                                      0
#define UQ_ML_SAMPLING_L_DATA_OUTPUT_FILE_NAME_ODV                            UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_DATA_OUTPUT_ALLOW_ALL_ODV                            0
#define UQ_ML_SAMPLING_L_DATA_OUTPUT_ALLOWED_SET_ODV                          ""
#define UQ_ML_SAMPLING_L_LOAD_BALANCE_ALGORITHM_ID_ODV                        2
#define UQ_ML_SAMPLING_L_LOAD_BALANCE_TRESHOLD_ODV                            1.
#define UQ_ML_SAMPLING_L_MIN_EFFECTIVE_SIZE_RATIO_ODV                         0.85
#define UQ_ML_SAMPLING_L_MAX_EFFECTIVE_SIZE_RATIO_ODV                         0.91
#define UQ_ML_SAMPLING_L_SCALE_COV_MATRIX_ODV                                 1
#define UQ_ML_SAMPLING_L_MIN_REJECTION_RATE_ODV                               0.50
#define UQ_ML_SAMPLING_L_MAX_REJECTION_RATE_ODV                               0.75
#define UQ_ML_SAMPLING_L_COV_REJECTION_RATE_ODV                               0.25
#define UQ_ML_SAMPLING_L_MIN_ACCEPTABLE_ETA_ODV                               0.
#define UQ_ML_SAMPLING_L_TOTALLY_MUTE_ODV                                     1
#define UQ_ML_SAMPLING_L_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV            UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV            UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_INITIAL_POSITION_USE_PREVIOUS_LEVEL_LIKELIHOOD_ODV   0
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV                   UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV                   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_LIST_OF_DISABLED_PARAMETERS_ODV                      ""
#define UQ_ML_SAMPLING_L_INITIAL_VALUES_OF_DISABLED_PARAMETERS_ODV            ""
#define UQ_ML_SAMPLING_L_RAW_CHAIN_SIZE_ODV                                   100
#define UQ_ML_SAMPLING_L_RAW_CHAIN_GENERATE_EXTRA_ODV                         0
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DISPLAY_PERIOD_ODV                         500
#define UQ_ML_SAMPLING_L_RAW_CHAIN_MEASURE_RUN_TIMES_ODV                      1
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV                     0
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV                  UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV                  UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV                  0
#define UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV                ""
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#define UQ_ML_SAMPLING_L_RAW_CHAIN_COMPUTE_STATS_ODV                          0
#endif
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_GENERATE_ODV                          0
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DISCARDED_PORTION_ODV                 0.
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_LAG_ODV                               1
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV             UQ_ML_SAMPLING_L_FILENAME_FOR_NO_FILE
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV             0
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
#define UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV        0
#define UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV      ""
#define UQ_ML_SAMPLING_L_AM_ETA_ODV                                           1.
#define UQ_ML_SAMPLING_L_AM_EPSILON_ODV                                       1.e-5
#define UQ_ML_SAMPLING_L_DO_LOGIT_TRANSFORM                                   0
#define UQ_ML_SAMPLING_L_ALGORITHM                                            "random_walk"
#define UQ_ML_SAMPLING_L_TK                                                   "random_walk"
#define UQ_ML_SAMPLING_L_UPDATE_INTERVAL                                      1

namespace QUESO {

/*!\file MLSamplingLevelOptions.h
   \brief Classes to allow options to be passed to the Multilevel algorithm, per level.*/

/*!\class MLSamplingLevelOptions
 * \brief This class provides options for each level of the Multilevel sequence generator if no input file is available.
 *
 * Multilevel sequence generator expects options for its methods to operate in each level
 * of the sequece generated. This class provides default values for such options if no
 * input file is available. */

class MLSamplingLevelOptions
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor: reads options from the input file.
  MLSamplingLevelOptions(const BaseEnvironment& env, const char* prefix);

  //MLSamplingLevelOptions(const MLSamplingLevelOptions& inputOptions);

  //! Destructor
  virtual ~MLSamplingLevelOptions();
  //@}

 //! @name Misc method
  //@{
  //! Access to the environment.
  const BaseEnvironment& env() const;

  //void changePrefix     (const char* prefix);
  //@}

  //! @name I/O methods
  //@{
  //!  It prints the option values.
  //! It scans the option values from the options input file.
  void scanOptionsValues(const MLSamplingLevelOptions* defaultOptions);

  void print            (std::ostream& os) const;
  //@}
  std::string                        m_prefix;

  //! If non-empty string, print options and values to output file
  std::string m_help;

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  //! Name of checkpoint output file
  std::string                        m_checkpointOutputFileName;
#endif
  //! Stop at end of such level.
  bool                               m_stopAtEnd;

  //! Name of generic output file.
  std::string                        m_dataOutputFileName;

  //! subEnvs that will write to generic output file.
  bool                               m_dataOutputAllowAll;

  //! subEnvs that will write to generic output file.
  std::set<unsigned int>             m_dataOutputAllowedSet;

  //! subEnvs that will write to generic output file.
  std::string                        m_str1;

  //! Perform load balancing with chosen algorithm (0 = no balancing).
  unsigned int                       m_loadBalanceAlgorithmId;

  //! Perform load balancing if load unbalancing ratio > threshold.
  double                             m_loadBalanceTreshold;

  //! Minimum allowed effective size ratio wrt previous level.
  double                             m_minEffectiveSizeRatio;

  //! Maximum allowed effective size ratio wrt previous level.
  double                             m_maxEffectiveSizeRatio;

  //! Whether or not scale proposal covariance matrix.
  bool                               m_scaleCovMatrix;

  //! minimum allowed attempted rejection rate at current level
  double                             m_minRejectionRate;

  //! maximum allowed attempted rejection rate at current level.
  double                             m_maxRejectionRate;

  //! c.o.v. for judging attempted rejection rate at current level.
  double                             m_covRejectionRate;

  //! minimum acceptable eta,
  /*! Used in the GPMSA code.*/
  double                             m_minAcceptableEta; // gpmsa1

  //! Whether or not to be totally mute (no printout message).
  bool                               m_totallyMute;

  //! Name of input file for initial position.
  std::string                        m_initialPositionDataInputFileName;

  //! Type of input file for initial position.
  std::string                        m_initialPositionDataInputFileType;

  //! Name of input file for initial proposal covariance matrix.
  std::string                        m_initialProposalCovMatrixDataInputFileName;

  //! Type of input file for initial proposal covariance matrix.
  std::string                        m_initialProposalCovMatrixDataInputFileType;

  //! Use previous level likelihood for initial chain position instead of re-computing it from target pdf
  bool                               m_initialPositionUsePreviousLevelLikelihood;  // ml_likelihood_caching


  std::set<unsigned int>             m_parameterDisabledSet; // gpmsa2
  std::string                        m_str2; // gpmsa2
  std::vector<double>                m_initialValuesOfDisabledParameters; // gpmsa2
  std::string                        m_str3; // gpmsa2

  //! Name of input file for raw chain.
  std::string                        m_rawChainDataInputFileName;


  //! Type of input file for raw chain
  std::string                        m_rawChainDataInputFileType;

  //! Size of raw chain
  unsigned int                       m_rawChainSize;

  //! Generate extra information about raw chain.
  bool                               m_rawChainGenerateExtra;

  //! Period of message display during raw chain generation.
  unsigned int                       m_rawChainDisplayPeriod;

  //!  Whether or not to measure run times.
  bool                               m_rawChainMeasureRunTimes;

  //! Period of message display during raw chain generation.
  unsigned int                       m_rawChainDataOutputPeriod;

  //! Name of output file for raw chain.
  std::string                        m_rawChainDataOutputFileName;

  //! Type of output file for raw chain.
  /*!
   * See MhOptionsValues::m_rawChainDataOutputFileType
   */
  std::string                        m_rawChainDataOutputFileType;

  //! Whether or not subEnvs will write to output file for raw chain.
  bool                               m_rawChainDataOutputAllowAll;

  //! subEnvs that will write to output file for raw chain.
  std::set<unsigned int>             m_rawChainDataOutputAllowedSet;

  std::string                        m_str4;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  //! Compute statistics on raw chain.
  bool                               m_rawChainComputeStats;
  SequenceStatisticalOptions*        m_rawChainStatisticalOptionsObj;
  bool                               m_rawChainStatOptsInstantiated;
#endif
  //! Whether or not to generate filtered chain.
  bool                               m_filteredChainGenerate;

  //! Initial discarded portion for chain filtering.
  double                             m_filteredChainDiscardedPortion; // input or set during run time

  //! Spacing for chain filtering.
  unsigned int                       m_filteredChainLag;              // input or set during run time

  //! Name of output file for filtered chain.
  std::string                        m_filteredChainDataOutputFileName;

  //! Type of output file for filtered chain.
  std::string                        m_filteredChainDataOutputFileType;

  //! Whether or not subEnvs will write to output file for filtered chain.
  bool                               m_filteredChainDataOutputAllowAll;

  //! subEnvs that will write to output file for filtered chain.
  std::set<unsigned int>             m_filteredChainDataOutputAllowedSet;
  std::string                        m_str5;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  //! Compute statistics on filtered chain.
  bool                               m_filteredChainComputeStats;
  SequenceStatisticalOptions*        m_filteredChainStatisticalOptionsObj;
  bool                               m_filteredChainStatOptsInstantiated;
#endif
  //! Display candidates generated in the core MH algorithm.
  bool                               m_displayCandidates;

  //! Put 'out of bound' candidates in chain as well.
  bool                               m_putOutOfBoundsInChain;

  //! Whether or not 'proposal' uses local Hessian.
  bool                               m_tkUseLocalHessian;

  //! Whether or not 'proposal' uses Newton component.
  bool                               m_tkUseNewtonComponent;

  //! 'dr' maximum number of extra stages.
  unsigned int                       m_drMaxNumExtraStages;

  //! 'dr' list of scales for proposal covariance matrices from 2nd stage on.
  std::vector<double>                m_drScalesForExtraStages;
  std::string                        m_str6;

  //! Whether or not 'dr' is used during 'am' non adaptive interval.
  bool                               m_drDuringAmNonAdaptiveInt;

  //! Whether or not 'am' will keep initial (given) matrix.
  bool                               m_amKeepInitialMatrix;

  //! 'am' initial non adaptation interval
  unsigned int                       m_amInitialNonAdaptInterval;

  //! 'am' adaptation interval.
  unsigned int                       m_amAdaptInterval;

  //! Period for outputing 'am' adapted matrices.
  unsigned int                       m_amAdaptedMatricesDataOutputPeriod;

  //! Name of output file for 'am' adapted matrices.
  std::string                        m_amAdaptedMatricesDataOutputFileName;

  //! Type of output file for 'am' adapted matrices.
  std::string                        m_amAdaptedMatricesDataOutputFileType;

  //! Whether or not subEnvs will write to output file for 'am' adapted matrices.
  bool                               m_amAdaptedMatricesDataOutputAllowAll;

  //! subEnvs that will write to output file for 'am' adapted matrices.
  std::set<unsigned int>             m_amAdaptedMatricesDataOutputAllowedSet;


  std::string                        m_str7;

  //! 'am' eta.
  double                             m_amEta;

  //! 'am' epsilon.
  double                             m_amEpsilon;

  //! Whether or not a logit transform will be done for bounded domains
  bool m_doLogitTransform;

  //! Which algorithm to use for sampling
  std::string m_algorithm;

  //! Which transition kernel to use for sampling
  std::string m_tk;

  //! How often to call the TK's updateTK method
  unsigned int m_updateInterval;

private:
  //! Copies the option values from \c srcOptions to \c this.
  void   copyOptionsValues(const MLSamplingLevelOptions& srcOptions);
  void getAllOptions();
  void defineAllOptions();

  const BaseEnvironment&        m_env;

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  BoostInputOptionsParser * m_parser;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  std::string                   m_option_help;

#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  std::string                   m_option_checkpointOutputFileName;
#endif
  std::string                   m_option_stopAtEnd;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowAll;
  std::string                   m_option_dataOutputAllowedSet;
  std::string                   m_option_loadBalanceAlgorithmId;
  std::string                   m_option_loadBalanceTreshold;
  std::string                   m_option_minEffectiveSizeRatio;
  std::string                   m_option_maxEffectiveSizeRatio;
  std::string                   m_option_scaleCovMatrix;
  std::string                   m_option_minRejectionRate;
  std::string                   m_option_maxRejectionRate;
  std::string                   m_option_covRejectionRate;
  std::string                   m_option_minAcceptableEta;  // gpmsa1
  std::string                   m_option_totallyMute;
  std::string                   m_option_initialPosition_dataInputFileName;
  std::string                   m_option_initialPosition_dataInputFileType;
  std::string                   m_option_initialProposalCovMatrix_dataInputFileName;
  std::string                   m_option_initialProposalCovMatrix_dataInputFileType;
  std::string                   m_option_initialPositionUsePreviousLevelLikelihood;  // ml_likelihood_caching
  std::string                   m_option_listOfDisabledParameters;  // gpmsa2
  std::string                   m_option_initialValuesOfDisabledParameters;  // gpmsa2
  std::string                   m_option_rawChain_dataInputFileName;
  std::string                   m_option_rawChain_dataInputFileType;
  std::string                   m_option_rawChain_size;
  std::string                   m_option_rawChain_generateExtra;
  std::string                   m_option_rawChain_displayPeriod;
  std::string                   m_option_rawChain_measureRunTimes;
  std::string                   m_option_rawChain_dataOutputPeriod;
  std::string                   m_option_rawChain_dataOutputFileName;

  //! Option name for MLSamplingLevelOptions::m_rawChainDataOutputFileType.  Option name is m_prefix + "ml_rawChain_dataOutputFileType"
  std::string                   m_option_rawChain_dataOutputFileType;

  std::string                   m_option_rawChain_dataOutputAllowAll;
  std::string                   m_option_rawChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_rawChain_computeStats;
#endif

  std::string                   m_option_filteredChain_generate;
  std::string                   m_option_filteredChain_discardedPortion;
  std::string                   m_option_filteredChain_lag;
  std::string                   m_option_filteredChain_dataOutputFileName;
  std::string                   m_option_filteredChain_dataOutputFileType;
  std::string                   m_option_filteredChain_dataOutputAllowAll;
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
  std::string                   m_option_am_adaptedMatrices_dataOutputAllowAll;
  std::string                   m_option_am_adaptedMatrices_dataOutputAllowedSet;
  std::string                   m_option_am_eta;
  std::string                   m_option_am_epsilon;
  std::string                   m_option_doLogitTransform;
  std::string                   m_option_algorithm;
  std::string                   m_option_tk;
  std::string                   m_option_updateInterval;

  void checkOptions(const BaseEnvironment * env);

  friend std::ostream & operator<<(std::ostream & os,
      const MLSamplingLevelOptions & obj);
};


}  // End namespace QUESO

#endif // UQ_MULTI_LEVEL_SAMPLING_LEVEL_OPTIONS_H
