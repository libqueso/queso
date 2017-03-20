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

#include <queso/MLSamplingLevelOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

MLSamplingLevelOptions::MLSamplingLevelOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
    m_prefix                                   ((std::string)(prefix) + ""),
    m_help                                     (UQ_ML_SAMPLING_L_HELP),
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
    m_checkpointOutputFileName                 (UQ_ML_SAMPLING_L_CHECKPOINT_OUTPUT_FILE_NAME_ODV),
#endif
    m_stopAtEnd                                (UQ_ML_SAMPLING_L_STOP_AT_END_ODV),
    m_dataOutputFileName                       (UQ_ML_SAMPLING_L_DATA_OUTPUT_FILE_NAME_ODV),
    m_dataOutputAllowAll                       (UQ_ML_SAMPLING_L_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_dataOutputAllowedSet                     (),
    m_str1                                     (""),
    m_loadBalanceAlgorithmId                   (UQ_ML_SAMPLING_L_LOAD_BALANCE_ALGORITHM_ID_ODV),
    m_loadBalanceTreshold                      (UQ_ML_SAMPLING_L_LOAD_BALANCE_TRESHOLD_ODV),
    m_minEffectiveSizeRatio                    (UQ_ML_SAMPLING_L_MIN_EFFECTIVE_SIZE_RATIO_ODV),
    m_maxEffectiveSizeRatio                    (UQ_ML_SAMPLING_L_MAX_EFFECTIVE_SIZE_RATIO_ODV),
    m_scaleCovMatrix                           (UQ_ML_SAMPLING_L_SCALE_COV_MATRIX_ODV),
    m_minRejectionRate                         (UQ_ML_SAMPLING_L_MIN_REJECTION_RATE_ODV),
    m_maxRejectionRate                         (UQ_ML_SAMPLING_L_MAX_REJECTION_RATE_ODV),
    m_covRejectionRate                         (UQ_ML_SAMPLING_L_COV_REJECTION_RATE_ODV),
    m_minAcceptableEta                         (UQ_ML_SAMPLING_L_MIN_ACCEPTABLE_ETA_ODV), // gpmsa1
    m_totallyMute                              (UQ_ML_SAMPLING_L_TOTALLY_MUTE_ODV),
    m_initialPositionDataInputFileName         (UQ_ML_SAMPLING_L_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV),
    m_initialPositionDataInputFileType         (UQ_ML_SAMPLING_L_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV),
    m_initialProposalCovMatrixDataInputFileName(UQ_ML_SAMPLING_L_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV),
    m_initialProposalCovMatrixDataInputFileType(UQ_ML_SAMPLING_L_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV),
    m_initialPositionUsePreviousLevelLikelihood(UQ_ML_SAMPLING_L_INITIAL_POSITION_USE_PREVIOUS_LEVEL_LIKELIHOOD_ODV),  // ml_likelihood_caching
  //m_parameterDisabledSet                     (), // gpmsa2
    m_str2                                     (""),
    m_initialValuesOfDisabledParameters        (0),
    m_str3                                     (""),
    m_rawChainDataInputFileName                (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV),
    m_rawChainDataInputFileType                (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV),
    m_rawChainSize                             (UQ_ML_SAMPLING_L_RAW_CHAIN_SIZE_ODV),
    m_rawChainGenerateExtra                    (UQ_ML_SAMPLING_L_RAW_CHAIN_GENERATE_EXTRA_ODV),
    m_rawChainDisplayPeriod                    (UQ_ML_SAMPLING_L_RAW_CHAIN_DISPLAY_PERIOD_ODV),
    m_rawChainMeasureRunTimes                  (UQ_ML_SAMPLING_L_RAW_CHAIN_MEASURE_RUN_TIMES_ODV),
    m_rawChainDataOutputPeriod                 (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV),
    m_rawChainDataOutputFileName               (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
    m_rawChainDataOutputFileType               (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
    m_rawChainDataOutputAllowAll               (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_rawChainDataOutputAllowedSet             (),
    m_str4                                     (""),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_rawChainComputeStats                     (UQ_ML_SAMPLING_L_RAW_CHAIN_COMPUTE_STATS_ODV),
    m_rawChainStatisticalOptionsObj            (NULL),
    m_rawChainStatOptsInstantiated             (false),
#endif
    m_filteredChainGenerate                    (UQ_ML_SAMPLING_L_FILTERED_CHAIN_GENERATE_ODV),
    m_filteredChainDiscardedPortion            (UQ_ML_SAMPLING_L_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
    m_filteredChainLag                         (UQ_ML_SAMPLING_L_FILTERED_CHAIN_LAG_ODV),
    m_filteredChainDataOutputFileName          (UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
    m_filteredChainDataOutputFileType          (UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
    m_filteredChainDataOutputAllowAll          (UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_filteredChainDataOutputAllowedSet        (),
    m_str5                                     (""),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_filteredChainComputeStats                (UQ_ML_SAMPLING_L_FILTERED_CHAIN_COMPUTE_STATS_ODV),
    m_filteredChainStatisticalOptionsObj       (NULL),
    m_filteredChainStatOptsInstantiated        (false),
#endif
    m_displayCandidates                        (UQ_ML_SAMPLING_L_DISPLAY_CANDIDATES_ODV),
    m_putOutOfBoundsInChain                    (UQ_ML_SAMPLING_L_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV),
    m_tkUseLocalHessian                        (UQ_ML_SAMPLING_L_TK_USE_LOCAL_HESSIAN_ODV),
    m_tkUseNewtonComponent                     (UQ_ML_SAMPLING_L_TK_USE_NEWTON_COMPONENT_ODV),
    m_drMaxNumExtraStages                      (UQ_ML_SAMPLING_L_DR_MAX_NUM_EXTRA_STAGES_ODV),
    m_drScalesForExtraStages                   (0),
    m_str6                                     ("1. "),
    m_drDuringAmNonAdaptiveInt                 (UQ_ML_SAMPLING_L_DR_DURING_AM_NON_ADAPTIVE_INT_ODV),
    m_amKeepInitialMatrix                      (UQ_ML_SAMPLING_L_AM_KEEP_INITIAL_MATRIX_ODV),
    m_amInitialNonAdaptInterval                (UQ_ML_SAMPLING_L_AM_INIT_NON_ADAPT_INT_ODV),
    m_amAdaptInterval                          (UQ_ML_SAMPLING_L_AM_ADAPT_INTERVAL_ODV),
    m_amAdaptedMatricesDataOutputPeriod        (UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV),
    m_amAdaptedMatricesDataOutputFileName      (UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV),
    m_amAdaptedMatricesDataOutputFileType      (UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV),
    m_amAdaptedMatricesDataOutputAllowAll      (UQ_ML_SAMPLING_L_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_amAdaptedMatricesDataOutputAllowedSet    (),
    m_str7                                     (""),
    m_amEta                                    (UQ_ML_SAMPLING_L_AM_ETA_ODV),
    m_amEpsilon                                (UQ_ML_SAMPLING_L_AM_EPSILON_ODV),
    m_env                                      (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(env.optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help                                      (m_prefix + "help"                                      ),
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
    m_option_checkpointOutputFileName                  (m_prefix + "checkpointOutputFileName"                  ),
#endif
    m_option_stopAtEnd                                 (m_prefix + "stopAtEnd"                                 ),
    m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
    m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
    m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
    m_option_loadBalanceAlgorithmId                    (m_prefix + "loadBalanceAlgorithmId"                    ),
    m_option_loadBalanceTreshold                       (m_prefix + "loadBalanceTreshold"                       ),
    m_option_minEffectiveSizeRatio                     (m_prefix + "minEffectiveSizeRatio"                     ),
    m_option_maxEffectiveSizeRatio                     (m_prefix + "maxEffectiveSizeRatio"                     ),
    m_option_scaleCovMatrix                            (m_prefix + "scaleCovMatrix"                            ),
    m_option_minRejectionRate                          (m_prefix + "minRejectionRate"                          ),
    m_option_maxRejectionRate                          (m_prefix + "maxRejectionRate"                          ),
    m_option_covRejectionRate                          (m_prefix + "covRejectionRate"                          ),
    m_option_minAcceptableEta                          (m_prefix + "minAcceptableEta"                          ), // gpmsa1
    m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
    m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
    m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
    m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
    m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
    m_option_initialPositionUsePreviousLevelLikelihood (m_prefix + "initialPositionUsePreviousLevelLikelihood" ), // ml_likelihood_caching
    m_option_listOfDisabledParameters                  (m_prefix + "listOfDisabledParameters"                  ), // gpmsa2
    m_option_initialValuesOfDisabledParameters         (m_prefix + "initialValuesOfDisabledParameters"         ), // gpmsa2
    m_option_rawChain_dataInputFileName                (m_prefix + "rawChain_dataInputFileName"                ),
    m_option_rawChain_dataInputFileType                (m_prefix + "rawChain_dataInputFileType"                ),
    m_option_rawChain_size                             (m_prefix + "rawChain_size"                             ),
    m_option_rawChain_generateExtra                    (m_prefix + "rawChain_generateExtra"                    ),
    m_option_rawChain_displayPeriod                    (m_prefix + "rawChain_displayPeriod"                    ),
    m_option_rawChain_measureRunTimes                  (m_prefix + "rawChain_measureRunTimes"                  ),
    m_option_rawChain_dataOutputPeriod                 (m_prefix + "rawChain_dataOutputPeriod"                 ),
    m_option_rawChain_dataOutputFileName               (m_prefix + "rawChain_dataOutputFileName"               ),
    m_option_rawChain_dataOutputFileType               (m_prefix + "rawChain_dataOutputFileType"               ),
    m_option_rawChain_dataOutputAllowAll               (m_prefix + "rawChain_dataOutputAllowAll"               ),
    m_option_rawChain_dataOutputAllowedSet             (m_prefix + "rawChain_dataOutputAllowedSet"             ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_option_rawChain_computeStats                     (m_prefix + "rawChain_computeStats"                     ),
#endif
    m_option_filteredChain_generate                    (m_prefix + "filteredChain_generate"                    ),
    m_option_filteredChain_discardedPortion            (m_prefix + "filteredChain_discardedPortion"            ),
    m_option_filteredChain_lag                         (m_prefix + "filteredChain_lag"                         ),
    m_option_filteredChain_dataOutputFileName          (m_prefix + "filteredChain_dataOutputFileName"          ),
    m_option_filteredChain_dataOutputFileType          (m_prefix + "filteredChain_dataOutputFileType"          ),
    m_option_filteredChain_dataOutputAllowAll          (m_prefix + "filteredChain_dataOutputAllowAll"          ),
    m_option_filteredChain_dataOutputAllowedSet        (m_prefix + "filteredChain_dataOutputAllowedSet"        ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_option_filteredChain_computeStats                (m_prefix + "filteredChain_computeStats"                ),
#endif
    m_option_displayCandidates                         (m_prefix + "displayCandidates"                         ),
    m_option_putOutOfBoundsInChain                     (m_prefix + "putOutOfBoundsInChain"                     ),
    m_option_tk_useLocalHessian                        (m_prefix + "tk_useLocalHessian"                        ),
    m_option_tk_useNewtonComponent                     (m_prefix + "tk_useNewtonComponent"                     ),
    m_option_dr_maxNumExtraStages                      (m_prefix + "dr_maxNumExtraStages"                      ),
    m_option_dr_listOfScalesForExtraStages             (m_prefix + "dr_listOfScalesForExtraStages"             ),
    m_option_dr_duringAmNonAdaptiveInt                 (m_prefix + "dr_duringAmNonAdaptiveInt"                 ),
    m_option_am_keepInitialMatrix                      (m_prefix + "am_keepInitialMatrix"                      ),
    m_option_am_initialNonAdaptInterval                (m_prefix + "am_initialNonAdaptInterval"                ),
    m_option_am_adaptInterval                          (m_prefix + "am_adaptInterval"                          ),
    m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "amAdaptedMatrices_dataOutputPeriod"        ),
    m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "amAdaptedMatrices_dataOutputFileName"      ),
    m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "amAdaptedMatrices_dataOutputFileType"      ),
    m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "amAdaptedMatrices_dataOutputAllowAll"      ),
    m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "amAdaptedMatrices_dataOutputAllowedSet"    ),
    m_option_am_eta                                    (m_prefix + "am_eta"                                    ),
    m_option_am_epsilon                                (m_prefix + "am_epsilon"                                ),
    m_option_doLogitTransform                          (m_prefix + "doLogitTransform"                          ),
    m_option_algorithm                                 (m_prefix + "algorithm"                                 ),
    m_option_tk                                        (m_prefix + "tk"                                        ),
    m_option_updateInterval                            (m_prefix + "updateInterval"                            )
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  this->defineAllOptions();
  m_parser->scanInputFile();
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  this->getAllOptions();

  checkOptions(&env);
}

void
MLSamplingLevelOptions::defineAllOptions()
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string >(m_option_help,                                       UQ_ML_SAMPLING_L_HELP                      , "produce help message for Bayesian Markov chain distr. calculator");
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  m_parser->registerOption<std::string >(m_option_checkpointOutputFileName,                   m_checkpointOutputFileName                 , "name of checpoint output file"                                   );
#endif
  m_parser->registerOption<bool        >(m_option_stopAtEnd,                                  m_stopAtEnd                                , "stop at end of such level"                                       );
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,                         m_dataOutputFileName                       , "name of generic output file"                                     );
  m_parser->registerOption<bool        >(m_option_dataOutputAllowAll,                         m_dataOutputAllowAll                       , "subEnvs that will write to generic output file"                  );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,                       m_str1                                     , "subEnvs that will write to generic output file"                  );
  m_parser->registerOption<unsigned int>(m_option_loadBalanceAlgorithmId,                     m_loadBalanceAlgorithmId                   , "Perform load balancing with chosen algorithm (0 = no balancing)" );
  m_parser->registerOption<double      >(m_option_loadBalanceTreshold,                        m_loadBalanceTreshold                      , "Perform load balancing if load unbalancing ratio > treshold"     );
  m_parser->registerOption<double      >(m_option_minEffectiveSizeRatio,                      m_minEffectiveSizeRatio                    , "minimum allowed effective size ratio wrt previous level"         );
  m_parser->registerOption<double      >(m_option_maxEffectiveSizeRatio,                      m_maxEffectiveSizeRatio                    , "maximum allowed effective size ratio wrt previous level"         );
  m_parser->registerOption<bool        >(m_option_scaleCovMatrix,                             m_scaleCovMatrix                           , "scale proposal covariance matrix"                                );
  m_parser->registerOption<double      >(m_option_minRejectionRate,                           m_minRejectionRate                         , "minimum allowed attempted rejection rate at current level"       );
  m_parser->registerOption<double      >(m_option_maxRejectionRate,                           m_maxRejectionRate                         , "maximum allowed attempted rejection rate at current level"       );
  m_parser->registerOption<double      >(m_option_covRejectionRate,                           m_covRejectionRate                         , "c.o.v. for judging attempted rejection rate at current level"    );
  m_parser->registerOption<double      >(m_option_minAcceptableEta,                           m_minAcceptableEta                         , "min acceptable eta"                                              );
  m_parser->registerOption<bool        >(m_option_totallyMute,                                m_totallyMute                              , "totally mute (no printout message)"                              );
  m_parser->registerOption<std::string >(m_option_initialPosition_dataInputFileName,          m_initialPositionDataInputFileName         , "name of input file for initial position"                         );
  m_parser->registerOption<std::string >(m_option_initialPosition_dataInputFileType,          m_initialPositionDataInputFileType         , "type of input file for initial position"                         );
  m_parser->registerOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName, "name of input file for initial proposal covariance matrix"       );
  m_parser->registerOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType, "type of input file for initial proposal covariance matrix"       );
  m_parser->registerOption<bool        >(m_option_initialPositionUsePreviousLevelLikelihood,  m_initialPositionUsePreviousLevelLikelihood, "use previous level likelihood for initial chain position instead of re-computing from target pdf");
  m_parser->registerOption<std::string >(m_option_listOfDisabledParameters,                   m_str2                                     , "list of disabled parameters"                                     );  // gpmsa2
  m_parser->registerOption<std::string >(m_option_initialValuesOfDisabledParameters,          m_str3                                     , "initial values of disabled parameters"                           );  // gpmsa2
  m_parser->registerOption<std::string >(m_option_rawChain_dataInputFileName,                 m_rawChainDataInputFileName                , "name of input file for raw chain "                               );
  m_parser->registerOption<std::string >(m_option_rawChain_dataInputFileType,                 m_rawChainDataInputFileType                , "type of input file for raw chain "                               );
  m_parser->registerOption<unsigned int>(m_option_rawChain_size,                              m_rawChainSize                             , "size of raw chain"                                               );
  m_parser->registerOption<bool        >(m_option_rawChain_generateExtra,                     m_rawChainGenerateExtra                    , "generate extra information about raw chain"                      );
  m_parser->registerOption<unsigned int>(m_option_rawChain_displayPeriod,                     m_rawChainDisplayPeriod                    , "period of message display during raw chain generation"           );
  m_parser->registerOption<bool        >(m_option_rawChain_measureRunTimes,                   m_rawChainMeasureRunTimes                  , "measure run times"                                               );
  m_parser->registerOption<unsigned int>(m_option_rawChain_dataOutputPeriod,                  m_rawChainDataOutputPeriod                 , "period of message display during raw chain generation"           );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputFileName,                m_rawChainDataOutputFileName               , "name of output file for raw chain "                              );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputFileType,                m_rawChainDataOutputFileType               , "type of output file for raw chain "                              );
  m_parser->registerOption<bool        >(m_option_rawChain_dataOutputAllowAll,                m_rawChainDataOutputAllowAll               , "subEnvs that will write to output file for raw chain"            );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputAllowedSet,              m_str4                                     , "subEnvs that will write to output file for raw chain"            );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_rawChain_computeStats,                      m_rawChainComputeStats                     , "compute statistics on raw chain"                                 );
#endif
  m_parser->registerOption<bool        >(m_option_filteredChain_generate,                     m_filteredChainGenerate                    , "generate filtered chain"                                         );
  m_parser->registerOption<double      >(m_option_filteredChain_discardedPortion,             m_filteredChainDiscardedPortion            , "initial discarded portion for chain filtering"                   );
  m_parser->registerOption<unsigned int>(m_option_filteredChain_lag,                          m_filteredChainLag                         , "spacing for chain filtering"                                     );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputFileName,           m_filteredChainDataOutputFileName          , "name of output file for filtered chain"                          );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputFileType,           m_filteredChainDataOutputFileType          , "type of output file for filtered chain"                          );
  m_parser->registerOption<bool        >(m_option_filteredChain_dataOutputAllowAll,           m_filteredChainDataOutputAllowAll          , "subEnvs that will write to output file for filtered chain"       );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputAllowedSet,         m_str5                                     , "subEnvs that will write to output file for filtered chain"       );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_filteredChain_computeStats,                 m_filteredChainComputeStats                , "compute statistics on filtered chain"                            );
#endif
  m_parser->registerOption<bool        >(m_option_displayCandidates,                          m_displayCandidates                        , "display candidates generated in the core MH algorithm"           );
  m_parser->registerOption<bool        >(m_option_putOutOfBoundsInChain,                      m_putOutOfBoundsInChain                    , "put 'out of bound' candidates in chain as well"                  );
  m_parser->registerOption<bool        >(m_option_tk_useLocalHessian,                         m_tkUseLocalHessian                        , "'proposal' use local Hessian"                                    );
  m_parser->registerOption<bool        >(m_option_tk_useNewtonComponent,                      m_tkUseNewtonComponent                     , "'proposal' use Newton component"                                 );
  m_parser->registerOption<unsigned int>(m_option_dr_maxNumExtraStages,                       m_drMaxNumExtraStages                      , "'dr' maximum number of extra stages"                             );
  m_parser->registerOption<std::string >(m_option_dr_listOfScalesForExtraStages,              m_str6                                     , "'dr' list of scales for proposal cov matrices from 2nd stage on" );
  m_parser->registerOption<bool        >(m_option_dr_duringAmNonAdaptiveInt,                  m_drDuringAmNonAdaptiveInt                 , "'dr' used during 'am' non adaptive interval"                     );
  m_parser->registerOption<bool        >(m_option_am_keepInitialMatrix,                       m_amKeepInitialMatrix                      , "'am' keep initial (given) matrix"                                );
  m_parser->registerOption<unsigned int>(m_option_am_initialNonAdaptInterval,                 m_amInitialNonAdaptInterval                , "'am' initial non adaptation interval"                            );
  m_parser->registerOption<unsigned int>(m_option_am_adaptInterval,                           m_amAdaptInterval                          , "'am' adaptation interval"                                        );
  m_parser->registerOption<unsigned int>(m_option_am_adaptedMatrices_dataOutputPeriod,        m_amAdaptedMatricesDataOutputPeriod        , "period for outputing 'am' adapted matrices"                      );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileName,      m_amAdaptedMatricesDataOutputFileName      , "name of output file for 'am' adapted matrices"                   );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileType,      m_amAdaptedMatricesDataOutputFileType      , "type of output file for 'am' adapted matrices"                   );
  m_parser->registerOption<bool        >(m_option_am_adaptedMatrices_dataOutputAllowAll,      m_amAdaptedMatricesDataOutputAllowAll      , "type of output file for 'am' adapted matrices"                   );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputAllowedSet,    m_str7                                     , "type of output file for 'am' adapted matrices"                   );
  m_parser->registerOption<double      >(m_option_am_eta,                                     m_amEta                                    , "'am' eta"                                                        );
  m_parser->registerOption<double      >(m_option_am_epsilon,                                 m_amEpsilon                                , "'am' epsilon"                                                    );
  m_parser->registerOption<bool        >(m_option_doLogitTransform,                           UQ_ML_SAMPLING_L_DO_LOGIT_TRANSFORM        , "flag for doing logit transform for bounded domains"              );
  m_parser->registerOption<std::string >(m_option_algorithm,                                  UQ_ML_SAMPLING_L_ALGORITHM                 , "which algorithm to use for sampling"                             );
  m_parser->registerOption<std::string >(m_option_tk,                                         UQ_ML_SAMPLING_L_TK                        , "which transition kernel to use for sampling"                     );
  m_parser->registerOption<unsigned int>(m_option_updateInterval,                             UQ_ML_SAMPLING_L_UPDATE_INTERVAL           , "how often to call TK's updateTK method"                          );
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
}

void
MLSamplingLevelOptions::getAllOptions()
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->getOption<std::string >(m_option_help,                                       m_help);
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  m_parser->getOption<std::string >(m_option_checkpointOutputFileName,                   m_checkpointOutputFileName                 );
#endif
  m_parser->getOption<bool        >(m_option_stopAtEnd,                                  m_stopAtEnd                                );
  m_parser->getOption<std::string >(m_option_dataOutputFileName,                         m_dataOutputFileName                       );
  m_parser->getOption<bool        >(m_option_dataOutputAllowAll,                         m_dataOutputAllowAll                       );
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,                       m_dataOutputAllowedSet);
  m_parser->getOption<unsigned int>(m_option_loadBalanceAlgorithmId,                     m_loadBalanceAlgorithmId                   );
  m_parser->getOption<double      >(m_option_loadBalanceTreshold,                        m_loadBalanceTreshold                      );
  m_parser->getOption<double      >(m_option_minEffectiveSizeRatio,                      m_minEffectiveSizeRatio                    );
  m_parser->getOption<double      >(m_option_maxEffectiveSizeRatio,                      m_maxEffectiveSizeRatio                    );
  m_parser->getOption<bool        >(m_option_scaleCovMatrix,                             m_scaleCovMatrix                           );
  m_parser->getOption<double      >(m_option_minRejectionRate,                           m_minRejectionRate                         );
  m_parser->getOption<double      >(m_option_maxRejectionRate,                           m_maxRejectionRate                         );
  m_parser->getOption<double      >(m_option_covRejectionRate,                           m_covRejectionRate                         );
  m_parser->getOption<double      >(m_option_minAcceptableEta,                           m_minAcceptableEta                         );
  m_parser->getOption<bool        >(m_option_totallyMute,                                m_totallyMute                              );
  m_parser->getOption<std::string >(m_option_initialPosition_dataInputFileName,          m_initialPositionDataInputFileName         );
  m_parser->getOption<std::string >(m_option_initialPosition_dataInputFileType,          m_initialPositionDataInputFileType         );
  m_parser->getOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName);
  m_parser->getOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType);
  m_parser->getOption<bool        >(m_option_initialPositionUsePreviousLevelLikelihood,  m_initialPositionUsePreviousLevelLikelihood);
  m_parser->getOption<std::set<unsigned int> >(m_option_listOfDisabledParameters,                   m_parameterDisabledSet);  // gpmsa2
  m_parser->getOption<std::vector<double> >(m_option_initialValuesOfDisabledParameters,          m_initialValuesOfDisabledParameters);  // gpmsa2
  m_parser->getOption<std::string >(m_option_rawChain_dataInputFileName,                 m_rawChainDataInputFileName                );
  m_parser->getOption<std::string >(m_option_rawChain_dataInputFileType,                 m_rawChainDataInputFileType                );
  m_parser->getOption<unsigned int>(m_option_rawChain_size,                              m_rawChainSize                             );
  m_parser->getOption<bool        >(m_option_rawChain_generateExtra,                     m_rawChainGenerateExtra                    );
  m_parser->getOption<unsigned int>(m_option_rawChain_displayPeriod,                     m_rawChainDisplayPeriod                    );
  m_parser->getOption<bool        >(m_option_rawChain_measureRunTimes,                   m_rawChainMeasureRunTimes                  );
  m_parser->getOption<unsigned int>(m_option_rawChain_dataOutputPeriod,                  m_rawChainDataOutputPeriod                 );
  m_parser->getOption<std::string >(m_option_rawChain_dataOutputFileName,                m_rawChainDataOutputFileName               );
  m_parser->getOption<std::string >(m_option_rawChain_dataOutputFileType,                m_rawChainDataOutputFileType               );
  m_parser->getOption<bool        >(m_option_rawChain_dataOutputAllowAll,                m_rawChainDataOutputAllowAll               );
  m_parser->getOption<std::set<unsigned int> >(m_option_rawChain_dataOutputAllowedSet,              m_rawChainDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool        >(m_option_rawChain_computeStats,                      m_rawChainComputeStats                     );
#endif
  m_parser->getOption<bool        >(m_option_filteredChain_generate,                     m_filteredChainGenerate                    );
  m_parser->getOption<double      >(m_option_filteredChain_discardedPortion,             m_filteredChainDiscardedPortion            );
  m_parser->getOption<unsigned int>(m_option_filteredChain_lag,                          m_filteredChainLag                         );
  m_parser->getOption<std::string >(m_option_filteredChain_dataOutputFileName,           m_filteredChainDataOutputFileName          );
  m_parser->getOption<std::string >(m_option_filteredChain_dataOutputFileType,           m_filteredChainDataOutputFileType          );
  m_parser->getOption<bool        >(m_option_filteredChain_dataOutputAllowAll,           m_filteredChainDataOutputAllowAll          );
  m_parser->getOption<std::set<unsigned int> >(m_option_filteredChain_dataOutputAllowedSet,         m_filteredChainDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool        >(m_option_filteredChain_computeStats,                 m_filteredChainComputeStats                );
#endif
  m_parser->getOption<bool        >(m_option_displayCandidates,                          m_displayCandidates                        );
  m_parser->getOption<bool        >(m_option_putOutOfBoundsInChain,                      m_putOutOfBoundsInChain                    );
  m_parser->getOption<bool        >(m_option_tk_useLocalHessian,                         m_tkUseLocalHessian                        );
  m_parser->getOption<bool        >(m_option_tk_useNewtonComponent,                      m_tkUseNewtonComponent                     );
  m_parser->getOption<unsigned int>(m_option_dr_maxNumExtraStages,                       m_drMaxNumExtraStages                      );
  m_parser->getOption<std::vector<double> >(m_option_dr_listOfScalesForExtraStages,              m_drScalesForExtraStages);
  m_parser->getOption<bool        >(m_option_dr_duringAmNonAdaptiveInt,                  m_drDuringAmNonAdaptiveInt                 );
  m_parser->getOption<bool        >(m_option_am_keepInitialMatrix,                       m_amKeepInitialMatrix                      );
  m_parser->getOption<unsigned int>(m_option_am_initialNonAdaptInterval,                 m_amInitialNonAdaptInterval                );
  m_parser->getOption<unsigned int>(m_option_am_adaptInterval,                           m_amAdaptInterval                          );
  m_parser->getOption<unsigned int>(m_option_am_adaptedMatrices_dataOutputPeriod,        m_amAdaptedMatricesDataOutputPeriod        );
  m_parser->getOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileName,      m_amAdaptedMatricesDataOutputFileName      );
  m_parser->getOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileType,      m_amAdaptedMatricesDataOutputFileType      );
  m_parser->getOption<bool        >(m_option_am_adaptedMatrices_dataOutputAllowAll,      m_amAdaptedMatricesDataOutputAllowAll      );
  m_parser->getOption<std::set<unsigned int> >(m_option_am_adaptedMatrices_dataOutputAllowedSet,    m_amAdaptedMatricesDataOutputAllowedSet);
  m_parser->getOption<double      >(m_option_am_eta,                                     m_amEta                                    );
  m_parser->getOption<double      >(m_option_am_epsilon,                                 m_amEpsilon                                );
  m_parser->getOption<bool        >(m_option_doLogitTransform,                           m_doLogitTransform);
  m_parser->getOption<std::string >(m_option_algorithm,                                  m_algorithm);
  m_parser->getOption<std::string >(m_option_tk,                                         m_tk);
  m_parser->getOption<unsigned int>(m_option_updateInterval,                             m_updateInterval);
#else
  m_help = m_env.input()(m_option_help,                                       UQ_ML_SAMPLING_L_HELP                      );
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  m_checkpointOutputFileName                  = m_env.input()(m_option_checkpointOutputFileName,                   m_checkpointOutputFileName                 );
#endif
  m_stopAtEnd                                 = m_env.input()(m_option_stopAtEnd,                                  m_stopAtEnd                                );
  m_dataOutputFileName                        = m_env.input()(m_option_dataOutputFileName,                         m_dataOutputFileName                       );
  m_dataOutputAllowAll                        = m_env.input()(m_option_dataOutputAllowAll,                         m_dataOutputAllowAll                       );

  // Is the empty set (string) by default
  unsigned int size = m_env.input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env.input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

  m_loadBalanceAlgorithmId                    = m_env.input()(m_option_loadBalanceAlgorithmId,                     m_loadBalanceAlgorithmId                   );
  m_loadBalanceTreshold                       = m_env.input()(m_option_loadBalanceTreshold,                        m_loadBalanceTreshold                      );
  m_minEffectiveSizeRatio                     = m_env.input()(m_option_minEffectiveSizeRatio,                      m_minEffectiveSizeRatio                    );
  m_maxEffectiveSizeRatio                     = m_env.input()(m_option_maxEffectiveSizeRatio,                      m_maxEffectiveSizeRatio                    );
  m_scaleCovMatrix                            = m_env.input()(m_option_scaleCovMatrix,                             m_scaleCovMatrix                           );
  m_minRejectionRate                          = m_env.input()(m_option_minRejectionRate,                           m_minRejectionRate                         );
  m_maxRejectionRate                          = m_env.input()(m_option_maxRejectionRate,                           m_maxRejectionRate                         );
  m_covRejectionRate                          = m_env.input()(m_option_covRejectionRate,                           m_covRejectionRate                         );
  m_minAcceptableEta                          = m_env.input()(m_option_minAcceptableEta,                           m_minAcceptableEta                         );
  m_totallyMute                               = m_env.input()(m_option_totallyMute,                                m_totallyMute                              );
  m_initialPositionDataInputFileName          = m_env.input()(m_option_initialPosition_dataInputFileName,          m_initialPositionDataInputFileName         );
  m_initialPositionDataInputFileType          = m_env.input()(m_option_initialPosition_dataInputFileType,          m_initialPositionDataInputFileType         );
  m_initialProposalCovMatrixDataInputFileName = m_env.input()(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName);
  m_initialProposalCovMatrixDataInputFileType = m_env.input()(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType);
  m_initialPositionUsePreviousLevelLikelihood = m_env.input()(m_option_initialPositionUsePreviousLevelLikelihood,  m_initialPositionUsePreviousLevelLikelihood);

  // Is the empty set (string) by default
  size = m_env.input().vector_variable_size(m_option_listOfDisabledParameters);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env.input()(m_option_listOfDisabledParameters, i, i);
    m_parameterDisabledSet.insert(allowed);
  }

  // Is the empty set (string) by default
  size = m_env.input().vector_variable_size(m_option_initialValuesOfDisabledParameters);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int value = m_env.input()(m_option_initialValuesOfDisabledParameters, i, i);
    m_initialValuesOfDisabledParameters.push_back(value);
  }

  m_rawChainDataInputFileName                 = m_env.input()(m_option_rawChain_dataInputFileName,                 m_rawChainDataInputFileName                );
  m_rawChainDataInputFileType                 = m_env.input()(m_option_rawChain_dataInputFileType,                 m_rawChainDataInputFileType                );
  m_rawChainSize                              = m_env.input()(m_option_rawChain_size,                              m_rawChainSize                             );
  m_rawChainGenerateExtra                     = m_env.input()(m_option_rawChain_generateExtra,                     m_rawChainGenerateExtra                    );
  m_rawChainDisplayPeriod                     = m_env.input()(m_option_rawChain_displayPeriod,                     m_rawChainDisplayPeriod                    );
  m_rawChainMeasureRunTimes                   = m_env.input()(m_option_rawChain_measureRunTimes,                   m_rawChainMeasureRunTimes                  );
  m_rawChainDataOutputPeriod                  = m_env.input()(m_option_rawChain_dataOutputPeriod,                  m_rawChainDataOutputPeriod                 );
  m_rawChainDataOutputFileName                = m_env.input()(m_option_rawChain_dataOutputFileName,                m_rawChainDataOutputFileName               );
  m_rawChainDataOutputFileType                = m_env.input()(m_option_rawChain_dataOutputFileType,                m_rawChainDataOutputFileType               );
  m_rawChainDataOutputAllowAll                = m_env.input()(m_option_rawChain_dataOutputAllowAll,                m_rawChainDataOutputAllowAll               );

  // Is the empty set (string) by default
  size = m_env.input().vector_variable_size(m_option_rawChain_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env.input()(m_option_rawChain_dataOutputAllowedSet, i, i);
    m_rawChainDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChainComputeStats                      = m_env.input()(m_option_rawChain_computeStats,                      m_rawChainComputeStats                     );
#endif
  m_filteredChainGenerate                     = m_env.input()(m_option_filteredChain_generate,                     m_filteredChainGenerate                    );
  m_filteredChainDiscardedPortion             = m_env.input()(m_option_filteredChain_discardedPortion,             m_filteredChainDiscardedPortion            );
  m_filteredChainLag                          = m_env.input()(m_option_filteredChain_lag,                          m_filteredChainLag                         );
  m_filteredChainDataOutputFileName           = m_env.input()(m_option_filteredChain_dataOutputFileName,           m_filteredChainDataOutputFileName          );
  m_filteredChainDataOutputFileType           = m_env.input()(m_option_filteredChain_dataOutputFileType,           m_filteredChainDataOutputFileType          );
  m_filteredChainDataOutputAllowAll           = m_env.input()(m_option_filteredChain_dataOutputAllowAll,           m_filteredChainDataOutputAllowAll          );

  // Is the empty set (string) by default
  size = m_env.input().vector_variable_size(m_option_filteredChain_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env.input()(m_option_filteredChain_dataOutputAllowedSet, i, i);
    m_filteredChainDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_filteredChainComputeStats                 = m_env.input()(m_option_filteredChain_computeStats,                 m_filteredChainComputeStats                );
#endif
  m_displayCandidates                         = m_env.input()(m_option_displayCandidates,                          m_displayCandidates                        );
  m_putOutOfBoundsInChain                     = m_env.input()(m_option_putOutOfBoundsInChain,                      m_putOutOfBoundsInChain                    );
  m_tkUseLocalHessian                         = m_env.input()(m_option_tk_useLocalHessian,                         m_tkUseLocalHessian                        );
  m_tkUseNewtonComponent                      = m_env.input()(m_option_tk_useNewtonComponent,                      m_tkUseNewtonComponent                     );
  m_drMaxNumExtraStages                       = m_env.input()(m_option_dr_maxNumExtraStages,                       m_drMaxNumExtraStages                      );

  // Is the empty set (string) by default
  size = m_env.input().vector_variable_size(m_option_dr_listOfScalesForExtraStages);

  m_drScalesForExtraStages.clear();
  for (unsigned int i = 0; i < size; i++) {
    unsigned int value = m_env.input()(m_option_dr_listOfScalesForExtraStages, i, i);
    m_drScalesForExtraStages.push_back(value);
  }

  m_drDuringAmNonAdaptiveInt                  = m_env.input()(m_option_dr_duringAmNonAdaptiveInt,                  m_drDuringAmNonAdaptiveInt                 );
  m_amKeepInitialMatrix                       = m_env.input()(m_option_am_keepInitialMatrix,                       m_amKeepInitialMatrix                      );
  m_amInitialNonAdaptInterval                 = m_env.input()(m_option_am_initialNonAdaptInterval,                 m_amInitialNonAdaptInterval                );
  m_amAdaptInterval                           = m_env.input()(m_option_am_adaptInterval,                           m_amAdaptInterval                          );
  m_amAdaptedMatricesDataOutputPeriod         = m_env.input()(m_option_am_adaptedMatrices_dataOutputPeriod,        m_amAdaptedMatricesDataOutputPeriod        );
  m_amAdaptedMatricesDataOutputFileName       = m_env.input()(m_option_am_adaptedMatrices_dataOutputFileName,      m_amAdaptedMatricesDataOutputFileName      );
  m_amAdaptedMatricesDataOutputFileType       = m_env.input()(m_option_am_adaptedMatrices_dataOutputFileType,      m_amAdaptedMatricesDataOutputFileType      );
  m_amAdaptedMatricesDataOutputAllowAll       = m_env.input()(m_option_am_adaptedMatrices_dataOutputAllowAll,      m_amAdaptedMatricesDataOutputAllowAll      );

  size = m_env.input().vector_variable_size(m_option_am_adaptedMatrices_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = m_env.input()(m_option_am_adaptedMatrices_dataOutputAllowedSet, i, i);
    m_amAdaptedMatricesDataOutputAllowedSet.insert(allowed);
  }

  m_amEta                                     = m_env.input()(m_option_am_eta,                                     m_amEta                                    );
  m_amEpsilon                                 = m_env.input()(m_option_am_epsilon,                                 m_amEpsilon                                );
  m_doLogitTransform                          = m_env.input()(m_option_doLogitTransform,                           UQ_ML_SAMPLING_L_DO_LOGIT_TRANSFORM        );
  m_algorithm                                 = m_env.input()(m_option_algorithm,                                  UQ_ML_SAMPLING_L_ALGORITHM                 );
  m_tk                                        = m_env.input()(m_option_tk,                                         UQ_ML_SAMPLING_L_TK                        );
  m_updateInterval                            = m_env.input()(m_option_updateInterval,                             UQ_ML_SAMPLING_L_UPDATE_INTERVAL           );
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
}

void
MLSamplingLevelOptions::copyOptionsValues(const MLSamplingLevelOptions& srcOptions)
{
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
  m_checkpointOutputFileName                  = srcOptions.m_checkpointOutputFileName;
#endif
  m_stopAtEnd                                 = srcOptions.m_stopAtEnd;
  m_dataOutputFileName                        = srcOptions.m_dataOutputFileName;
  m_dataOutputAllowAll                        = srcOptions.m_dataOutputAllowAll;
  m_dataOutputAllowedSet                      = srcOptions.m_dataOutputAllowedSet;
  m_str1                                      = srcOptions.m_str1;
  m_loadBalanceAlgorithmId                    = srcOptions.m_loadBalanceAlgorithmId;
  m_loadBalanceTreshold                       = srcOptions.m_loadBalanceTreshold;
  m_minEffectiveSizeRatio                     = srcOptions.m_minEffectiveSizeRatio;
  m_maxEffectiveSizeRatio                     = srcOptions.m_maxEffectiveSizeRatio;
  m_scaleCovMatrix                            = srcOptions.m_scaleCovMatrix;
  m_minRejectionRate                          = srcOptions.m_minRejectionRate;
  m_maxRejectionRate                          = srcOptions.m_maxRejectionRate;
  m_covRejectionRate                          = srcOptions.m_covRejectionRate;
  m_minAcceptableEta                          = srcOptions.m_minAcceptableEta; // gpmsa1
  m_totallyMute                               = srcOptions.m_totallyMute;
  m_initialPositionDataInputFileName          = srcOptions.m_initialPositionDataInputFileName;
  m_initialPositionDataInputFileType          = srcOptions.m_initialPositionDataInputFileType;
  m_initialProposalCovMatrixDataInputFileName = srcOptions.m_initialProposalCovMatrixDataInputFileName;
  m_initialProposalCovMatrixDataInputFileType = srcOptions.m_initialProposalCovMatrixDataInputFileType;
  m_initialPositionUsePreviousLevelLikelihood = srcOptions.m_initialPositionUsePreviousLevelLikelihood;
  m_parameterDisabledSet                      = srcOptions.m_parameterDisabledSet; // gpmsa2
  m_str2                                      = srcOptions.m_str2; // gpmsa2
  m_initialValuesOfDisabledParameters         = srcOptions.m_initialValuesOfDisabledParameters;
  m_str3                                      = srcOptions.m_str3;
  m_rawChainDataInputFileName                 = srcOptions.m_rawChainDataInputFileName;
  m_rawChainDataInputFileType                 = srcOptions.m_rawChainDataInputFileType;
  m_rawChainSize                              = srcOptions.m_rawChainSize;
//std::cout << "In copy(), rawChainSize = " << m_rawChainSize << std::endl;
  m_rawChainGenerateExtra                     = srcOptions.m_rawChainGenerateExtra;
  m_rawChainDisplayPeriod                     = srcOptions.m_rawChainDisplayPeriod;
  m_rawChainMeasureRunTimes                   = srcOptions.m_rawChainMeasureRunTimes;
  m_rawChainDataOutputPeriod                  = srcOptions.m_rawChainDataOutputPeriod;
  m_rawChainDataOutputFileName                = srcOptions.m_rawChainDataOutputFileName;
  m_rawChainDataOutputFileType                = srcOptions.m_rawChainDataOutputFileType;
  m_rawChainDataOutputAllowAll                = srcOptions.m_rawChainDataOutputAllowAll;
  m_rawChainDataOutputAllowedSet              = srcOptions.m_rawChainDataOutputAllowedSet;
  m_str4                                      = srcOptions.m_str4;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChainComputeStats                      = srcOptions.m_rawChainComputeStats;
  m_rawChainStatisticalOptionsObj             = NULL; // Yes, 'NULL'
  m_rawChainStatOptsInstantiated              = false;
#endif
  m_filteredChainGenerate                     = srcOptions.m_filteredChainGenerate;
  m_filteredChainDiscardedPortion             = srcOptions.m_filteredChainDiscardedPortion;
  m_filteredChainLag                          = srcOptions.m_filteredChainLag;
  m_filteredChainDataOutputFileName           = srcOptions.m_filteredChainDataOutputFileName;
  m_filteredChainDataOutputFileType           = srcOptions.m_filteredChainDataOutputFileType;
  m_filteredChainDataOutputAllowAll           = srcOptions.m_filteredChainDataOutputAllowAll;
  m_filteredChainDataOutputAllowedSet         = srcOptions.m_filteredChainDataOutputAllowedSet;
  m_str5                                      = srcOptions.m_str5;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_filteredChainComputeStats                 = srcOptions.m_filteredChainComputeStats;
  m_filteredChainStatisticalOptionsObj        = NULL; // Yes, 'NULL'
  m_filteredChainStatOptsInstantiated         = false;
#endif
  m_displayCandidates                         = srcOptions.m_displayCandidates;
  m_putOutOfBoundsInChain                     = srcOptions.m_putOutOfBoundsInChain;
  m_tkUseLocalHessian                         = srcOptions.m_tkUseLocalHessian;
  m_tkUseNewtonComponent                      = srcOptions.m_tkUseNewtonComponent;
  m_drMaxNumExtraStages                       = srcOptions.m_drMaxNumExtraStages;
  m_drScalesForExtraStages                    = srcOptions.m_drScalesForExtraStages;
  m_str6                                      = srcOptions.m_str6;
  m_drDuringAmNonAdaptiveInt                  = srcOptions.m_drDuringAmNonAdaptiveInt;
  m_amKeepInitialMatrix                       = srcOptions.m_amKeepInitialMatrix;
  m_amInitialNonAdaptInterval                 = srcOptions.m_amInitialNonAdaptInterval;
  m_amAdaptInterval                           = srcOptions.m_amAdaptInterval;
  m_amAdaptedMatricesDataOutputPeriod         = srcOptions.m_amAdaptedMatricesDataOutputPeriod;
  m_amAdaptedMatricesDataOutputFileName       = srcOptions.m_amAdaptedMatricesDataOutputFileName;
  m_amAdaptedMatricesDataOutputFileType       = srcOptions.m_amAdaptedMatricesDataOutputFileType;
  m_amAdaptedMatricesDataOutputAllowAll       = srcOptions.m_amAdaptedMatricesDataOutputAllowAll;
  m_amAdaptedMatricesDataOutputAllowedSet     = srcOptions.m_amAdaptedMatricesDataOutputAllowedSet;
  m_str7                                      = srcOptions.m_str7;
  m_amEta                                     = srcOptions.m_amEta;
  m_amEpsilon                                 = srcOptions.m_amEpsilon;
  m_doLogitTransform                          = srcOptions.m_doLogitTransform;
  m_algorithm                                 = srcOptions.m_algorithm;
  m_tk                                        = srcOptions.m_tk;
  m_updateInterval                            = srcOptions.m_updateInterval;

  return;
}

MLSamplingLevelOptions::~MLSamplingLevelOptions()
{
  //std::cout << "In MLSamplingLevelOptions::destructor()"
  //          << ": m_filteredChainStatOptsInstantiated = " << m_filteredChainStatOptsInstantiated
  //          << ", m_rawChainStatOptsInstantiated = "      << m_rawChainStatOptsInstantiated
  //          << std::endl;
  //sleep(1);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_filteredChainStatOptsInstantiated) delete m_filteredChainStatisticalOptionsObj;
  if (m_rawChainStatOptsInstantiated     ) delete m_rawChainStatisticalOptionsObj;
#endif
}

void
MLSamplingLevelOptions::scanOptionsValues(const MLSamplingLevelOptions* defaultOptions)
{
  queso_deprecated();

  // FIXME
  if (defaultOptions) this->copyOptionsValues(*defaultOptions);

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  // Replace the parser since default values changed
  if (m_parser) {
    delete m_parser;
    m_parser = new BoostInputOptionsParser(m_env.optionsInputFileName());
  }

  // FIXME: Does this work with GetPot?
  this->defineAllOptions();
  m_parser->scanInputFile();
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  this->getAllOptions();

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_rawChainComputeStats) {
    m_rawChainStatisticalOptionsObj = new SequenceStatisticalOptions(m_env,m_prefix + "rawChain_");
    m_rawChainStatOptsInstantiated  = true;
  }
  if (m_filteredChainComputeStats) {
    m_filteredChainStatisticalOptionsObj = new SequenceStatisticalOptions(m_env,m_prefix + "filteredChain_");
    m_filteredChainStatOptsInstantiated  = true;
  }
#endif
}

void
MLSamplingLevelOptions::checkOptions(const BaseEnvironment * env)
{
  char tmpStr[64];

  // DM: Print here because I don't know where the class is instantiated
  if (m_help != "") {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*this) << std::endl;
    }
  }

  if (m_dataOutputAllowAll) {
    m_dataOutputAllowedSet.clear();
    m_dataOutputAllowedSet.insert(m_env.subId());
  }

  // DM: Not sure what this is for
  m_str1.clear();
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",(int)(*setIt));
    m_str1 += tmpStr;
    m_str1 += " ";
  }

  queso_require_less_msg(m_minEffectiveSizeRatio, 1.0, "option `" << m_option_minEffectiveSizeRatio << "` must be less than 1.0");
  queso_require_less_msg(m_maxEffectiveSizeRatio, 1.0, "option `" << m_option_maxEffectiveSizeRatio << "` must be less than 1.0");
  queso_require_less_msg(m_minRejectionRate, 1.0, "option `" << m_option_minRejectionRate << "` must be less than 1.0");
  queso_require_less_msg(m_maxRejectionRate, 1.0, "option `" << m_option_maxRejectionRate << "` must be less than 1.0");
  queso_require_less_msg(m_covRejectionRate, 1.0, "option `" << m_option_covRejectionRate << "` must be less than 1.0");

  // DM: Not sure what this is for
  m_str2.clear();
  for (std::set<unsigned int>::iterator setIt = m_parameterDisabledSet.begin(); setIt != m_parameterDisabledSet.end(); ++setIt) { // gpmsa2
    sprintf(tmpStr,"%d",(int)(*setIt));
    m_str2 += tmpStr;
    m_str2 += " ";
  }

  // DM: Not sure what this is for
  m_str3.clear();
  for (unsigned int i = 0; i < m_initialValuesOfDisabledParameters.size(); ++i) {
    sprintf(tmpStr,"%e",m_initialValuesOfDisabledParameters[i]);
    m_str3 += tmpStr;
    m_str3 += " ";
  }

  if (m_rawChainDataOutputAllowAll) {
    m_rawChainDataOutputAllowedSet.clear();
    m_rawChainDataOutputAllowedSet.insert(m_env.subId());
  }

  // DM: Not sure what this is for
  m_str4.clear();
  for (std::set<unsigned int>::iterator setIt = m_rawChainDataOutputAllowedSet.begin(); setIt != m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",(int)(*setIt));
    m_str4 += tmpStr;
    m_str4 += " ";
  }

  if (m_filteredChainGenerate == true) {
    queso_require_greater_equal_msg(m_filteredChainLag, 2, "option `" << m_option_filteredChain_lag << "` must be at least 2");
  }

  if (m_filteredChainDataOutputAllowAll) {
    m_filteredChainDataOutputAllowedSet.clear();
    m_filteredChainDataOutputAllowedSet.insert(m_env.subId());
  }

  // DM: Not sure what this is for
  m_str5.clear();
  for (std::set<unsigned int>::iterator setIt = m_filteredChainDataOutputAllowedSet.begin(); setIt != m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",(int)(*setIt));
    m_str5 += tmpStr;
    m_str5 += " ";
  }

  if (m_drMaxNumExtraStages > 0) {
    unsigned int size = m_drScalesForExtraStages.size();

    // If we asked more stages than scales provided, pad with ones
    if (m_drMaxNumExtraStages > size) {
      for (unsigned int i = size; i < m_drMaxNumExtraStages; i++) {
        double scale = 1.0;
        m_drScalesForExtraStages.push_back(scale);
      }
    }
  }

  // DM: Not sure what this is for
  m_str6.clear();
  for (unsigned int i = 0; i < m_drScalesForExtraStages.size(); ++i) {
    sprintf(tmpStr,"%e",m_drScalesForExtraStages[i]);
    m_str6 += tmpStr;
    m_str6 += " ";
  }

  if (m_amAdaptedMatricesDataOutputAllowAll) {
    m_amAdaptedMatricesDataOutputAllowedSet.clear();
    m_amAdaptedMatricesDataOutputAllowedSet.insert(m_env.subId());
  }

  // DM: Not sure what this is for
  m_str7.clear();
  for (std::set<unsigned int>::iterator setIt = m_amAdaptedMatricesDataOutputAllowedSet.begin(); setIt != m_amAdaptedMatricesDataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",(int)(*setIt));
    m_str7 += tmpStr;
    m_str7 += " ";
  }
}

void
MLSamplingLevelOptions::print(std::ostream& os) const
{
  os <<        "m_prefix"                         << " = " << m_prefix
#ifdef ML_CODE_HAS_NEW_RESTART_CAPABILITY
#else
     << "\n" << m_option_checkpointOutputFileName << " = " << m_checkpointOutputFileName
#endif
     << "\n" << m_option_stopAtEnd                << " = " << m_stopAtEnd
     << "\n" << m_option_dataOutputFileName       << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowAll       << " = " << m_dataOutputAllowAll
     << "\n" << m_option_dataOutputAllowedSet     << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_loadBalanceAlgorithmId                     << " = " << m_loadBalanceAlgorithmId
     << "\n" << m_option_loadBalanceTreshold                        << " = " << m_loadBalanceTreshold
     << "\n" << m_option_minEffectiveSizeRatio                      << " = " << m_minEffectiveSizeRatio
     << "\n" << m_option_maxEffectiveSizeRatio                      << " = " << m_maxEffectiveSizeRatio
     << "\n" << m_option_scaleCovMatrix                             << " = " << m_scaleCovMatrix
     << "\n" << m_option_minRejectionRate                           << " = " << m_minRejectionRate
     << "\n" << m_option_maxRejectionRate                           << " = " << m_maxRejectionRate
     << "\n" << m_option_covRejectionRate                           << " = " << m_covRejectionRate
     << "\n" << m_option_minAcceptableEta                           << " = " << m_minAcceptableEta // gpmsa1
     << "\n" << m_option_totallyMute                                << " = " << m_totallyMute
     << "\n" << m_option_initialPosition_dataInputFileName          << " = " << m_initialPositionDataInputFileName
     << "\n" << m_option_initialPosition_dataInputFileType          << " = " << m_initialPositionDataInputFileType
     << "\n" << m_option_initialProposalCovMatrix_dataInputFileName << " = " << m_initialProposalCovMatrixDataInputFileName
     << "\n" << m_option_initialProposalCovMatrix_dataInputFileType << " = " << m_initialProposalCovMatrixDataInputFileType
     << "\n" << m_option_initialPositionUsePreviousLevelLikelihood  << " = " << m_initialPositionUsePreviousLevelLikelihood
     << "\n" << m_option_listOfDisabledParameters                   << " = "; // gpmsa2
  for (std::set<unsigned int>::iterator setIt = m_parameterDisabledSet.begin(); setIt != m_parameterDisabledSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_initialValuesOfDisabledParameters          << " = "; // gpmsa2
  for (unsigned int i = 0; i < m_initialValuesOfDisabledParameters.size(); ++i) {
    os << m_initialValuesOfDisabledParameters[i] << " ";
  }
  os << "\n" << m_option_rawChain_dataInputFileName                 << " = " << m_rawChainDataInputFileName
     << "\n" << m_option_rawChain_dataInputFileType                 << " = " << m_rawChainDataInputFileType
     << "\n" << m_option_rawChain_size                              << " = " << m_rawChainSize
     << "\n" << m_option_rawChain_generateExtra                     << " = " << m_rawChainGenerateExtra
     << "\n" << m_option_rawChain_displayPeriod                     << " = " << m_rawChainDisplayPeriod
     << "\n" << m_option_rawChain_measureRunTimes                   << " = " << m_rawChainMeasureRunTimes
     << "\n" << m_option_rawChain_dataOutputPeriod                  << " = " << m_rawChainDataOutputPeriod
     << "\n" << m_option_rawChain_dataOutputFileName                << " = " << m_rawChainDataOutputFileName
     << "\n" << m_option_rawChain_dataOutputFileType                << " = " << m_rawChainDataOutputFileType
     << "\n" << m_option_rawChain_dataOutputAllowAll                << " = " << m_rawChainDataOutputAllowAll
     << "\n" << m_option_rawChain_dataOutputAllowedSet              << " = ";
  for (std::set<unsigned int>::iterator setIt = m_rawChainDataOutputAllowedSet.begin(); setIt != m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << m_option_rawChain_computeStats                      << " = " << m_rawChainComputeStats
#endif
     << "\n" << m_option_filteredChain_generate                     << " = " << m_filteredChainGenerate
     << "\n" << m_option_filteredChain_discardedPortion             << " = " << m_filteredChainDiscardedPortion
     << "\n" << m_option_filteredChain_lag                          << " = " << m_filteredChainLag
     << "\n" << m_option_filteredChain_dataOutputFileName           << " = " << m_filteredChainDataOutputFileName
     << "\n" << m_option_filteredChain_dataOutputFileType           << " = " << m_filteredChainDataOutputFileType
     << "\n" << m_option_filteredChain_dataOutputAllowAll           << " = " << m_filteredChainDataOutputAllowAll
     << "\n" << m_option_filteredChain_dataOutputAllowedSet         << " = ";
  for (std::set<unsigned int>::iterator setIt = m_filteredChainDataOutputAllowedSet.begin(); setIt != m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << m_option_filteredChain_computeStats                 << " = " << m_filteredChainComputeStats
#endif
     << "\n" << m_option_displayCandidates                          << " = " << m_displayCandidates
     << "\n" << m_option_putOutOfBoundsInChain                      << " = " << m_putOutOfBoundsInChain
     << "\n" << m_option_tk_useLocalHessian                         << " = " << m_tkUseLocalHessian
     << "\n" << m_option_tk_useNewtonComponent                      << " = " << m_tkUseNewtonComponent
     << "\n" << m_option_dr_maxNumExtraStages                       << " = " << m_drMaxNumExtraStages
     << "\n" << m_option_dr_listOfScalesForExtraStages              << " = ";
  for (unsigned int i = 0; i < m_drScalesForExtraStages.size(); ++i) {
    os << m_drScalesForExtraStages[i] << " ";
  }
  os << "\n" << m_option_dr_duringAmNonAdaptiveInt                  << " = " << m_drDuringAmNonAdaptiveInt
     << "\n" << m_option_am_keepInitialMatrix                       << " = " << m_amKeepInitialMatrix
     << "\n" << m_option_am_initialNonAdaptInterval                 << " = " << m_amInitialNonAdaptInterval
     << "\n" << m_option_am_adaptInterval                           << " = " << m_amAdaptInterval
     << "\n" << m_option_am_adaptedMatrices_dataOutputPeriod        << " = " << m_amAdaptedMatricesDataOutputPeriod
     << "\n" << m_option_am_adaptedMatrices_dataOutputFileName      << " = " << m_amAdaptedMatricesDataOutputFileName
     << "\n" << m_option_am_adaptedMatrices_dataOutputFileType      << " = " << m_amAdaptedMatricesDataOutputFileType
     << "\n" << m_option_am_adaptedMatrices_dataOutputAllowAll      << " = " << m_amAdaptedMatricesDataOutputAllowAll
     << "\n" << m_option_am_adaptedMatrices_dataOutputAllowedSet    << " = ";
  for (std::set<unsigned int>::iterator setIt = m_amAdaptedMatricesDataOutputAllowedSet.begin(); setIt != m_amAdaptedMatricesDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_am_eta                                     << " = " << m_amEta
     << "\n" << m_option_am_epsilon                                 << " = " << m_amEpsilon
     << "\n" << m_option_doLogitTransform                           << " = " << m_doLogitTransform
     << "\n" << m_option_algorithm                                  << " = " << m_algorithm
     << "\n" << m_option_tk                                         << " = " << m_tk
     << "\n" << m_option_updateInterval                             << " = " << m_updateInterval
     << std::endl;

  return;
}

const BaseEnvironment&
MLSamplingLevelOptions::env() const
{
  return m_env;
}

std::ostream& operator<<(std::ostream& os, const MLSamplingLevelOptions& obj)
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  obj.print(os);
  return os;
}

}  // End namespace QUESO
