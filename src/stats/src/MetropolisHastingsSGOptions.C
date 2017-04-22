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

#include <queso/Environment.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/Miscellaneous.h>

// -------------------------------------------------
// MhOptionsValues --------------------------
// -------------------------------------------------

namespace QUESO {

// Default constructor -----------------------------
MhOptionsValues::MhOptionsValues(
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  const SsOptionsValues* alternativeRawSsOptionsValues,
  const SsOptionsValues* alternativeFilteredSsOptionsValues
#endif
  )
  :
    m_prefix                                           ("mh_"),
    m_help(UQ_MH_SG_HELP),
    m_dataOutputFileName                       (UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV),
    m_dataOutputAllowAll                       (UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_dataOutputAllowedSet                     (),
    m_totallyMute                              (UQ_MH_SG_TOTALLY_MUTE_ODV),
    m_initialPositionDataInputFileName         (UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV),
    m_initialPositionDataInputFileType         (UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV),
    m_initialProposalCovMatrixDataInputFileName(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV),
    m_initialProposalCovMatrixDataInputFileType(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV),
  //m_parameterDisabledSet                     (),
    m_rawChainDataInputFileName                (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV),
    m_rawChainDataInputFileType                (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV),
    m_rawChainSize                             (UQ_MH_SG_RAW_CHAIN_SIZE_ODV),
    m_rawChainGenerateExtra                    (UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV),
    m_rawChainDisplayPeriod                    (UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV),
    m_rawChainMeasureRunTimes                  (UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV),
    m_rawChainDataOutputPeriod                 (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV),
    m_rawChainDataOutputFileName               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
    m_rawChainDataOutputFileType               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
    m_rawChainDataOutputAllowAll               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_rawChainDataOutputAllowedSet             (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_rawChainComputeStats                     (UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV),
#endif
    m_filteredChainGenerate                    (UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV),
    m_filteredChainDiscardedPortion            (UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
    m_filteredChainLag                         (UQ_MH_SG_FILTERED_CHAIN_LAG_ODV),
    m_filteredChainDataOutputFileName          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
    m_filteredChainDataOutputFileType          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
    m_filteredChainDataOutputAllowAll          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_filteredChainDataOutputAllowedSet        (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_filteredChainComputeStats                (UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV),
#endif
    m_displayCandidates                        (UQ_MH_SG_DISPLAY_CANDIDATES_ODV),
    m_putOutOfBoundsInChain                    (UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV),
    m_tkUseLocalHessian                        (UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV),
    m_tkUseNewtonComponent                     (UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV),
    m_drMaxNumExtraStages                      (UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV),
    m_drScalesForExtraStages                   (0),
    m_drDuringAmNonAdaptiveInt                 (UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV),
    m_amKeepInitialMatrix                      (UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV),
    m_amInitialNonAdaptInterval                (UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV),
    m_amAdaptInterval                          (UQ_MH_SG_AM_ADAPT_INTERVAL_ODV),
    m_amAdaptedMatricesDataOutputPeriod        (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV),
    m_amAdaptedMatricesDataOutputFileName      (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV),
    m_amAdaptedMatricesDataOutputFileType      (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV),
    m_amAdaptedMatricesDataOutputAllowAll      (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_amAdaptedMatricesDataOutputAllowedSet    (),
    m_amEta                                    (UQ_MH_SG_AM_ETA_ODV),
    m_amEpsilon                                (UQ_MH_SG_AM_EPSILON_ODV),
    m_enableBrooksGelmanConvMonitor            (UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR),
    m_BrooksGelmanLag                          (UQ_MH_SG_BROOKS_GELMAN_LAG),
    m_outputLogLikelihood                      (UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD),
    m_outputLogTarget                          (UQ_MH_SG_OUTPUT_LOG_TARGET),
    m_doLogitTransform                         (UQ_MH_SG_DO_LOGIT_TRANSFORM),
    m_algorithm                                (UQ_MH_SG_ALGORITHM),
    m_tk                                       (UQ_MH_SG_TK),
    m_updateInterval                           (UQ_MH_SG_UPDATE_INTERVAL),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_alternativeRawSsOptionsValues            (),
    m_alternativeFilteredSsOptionsValues       (),
#endif
    m_env(NULL),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help                                      (m_prefix + "help"                                      ),
    m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
    m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
    m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
    m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
    m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
    m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
    m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
    m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
    m_option_listOfDisabledParameters                  (m_prefix + "listOfDisabledParameters"                  ),
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
    m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "am_adaptedMatrices_dataOutputPeriod"       ),
    m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "am_adaptedMatrices_dataOutputFileName"     ),
    m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "am_adaptedMatrices_dataOutputFileType"     ),
    m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "am_adaptedMatrices_dataOutputAllowAll"     ),
    m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "am_adaptedMatrices_dataOutputAllowedSet"   ),
    m_option_am_eta                                    (m_prefix + "am_eta"                                    ),
    m_option_am_epsilon                                (m_prefix + "am_epsilon"                                ),
    m_option_enableBrooksGelmanConvMonitor             (m_prefix + "enableBrooksGelmanConvMonitor"             ),
    m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           ),
    m_option_outputLogLikelihood                       (m_prefix + "outputLogLikelihood"                       ),
    m_option_outputLogTarget                           (m_prefix + "outputLogTarget"                           ),
    m_option_doLogitTransform                          (m_prefix + "doLogitTransform"                          ),
    m_option_algorithm                                 (m_prefix + "algorithm"                                 ),
    m_option_tk                                        (m_prefix + "tk"                                        ),
    m_option_updateInterval                            (m_prefix + "updateInterval"                            )
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativeRawSsOptionsValues     ) m_alternativeRawSsOptionsValues      = *alternativeRawSsOptionsValues;
  if (alternativeFilteredSsOptionsValues) m_alternativeFilteredSsOptionsValues = *alternativeFilteredSsOptionsValues;
#endif
}

MhOptionsValues::MhOptionsValues(
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  const SsOptionsValues* alternativeRawSsOptionsValues,
  const SsOptionsValues* alternativeFilteredSsOptionsValues,
#endif
  const BaseEnvironment * env,
  const char * prefix
  )
  :
    m_prefix                                           ((std::string)(prefix) + "mh_"),
    m_help(UQ_MH_SG_HELP),
    m_dataOutputFileName                       (UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV),
    m_dataOutputAllowAll                       (UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_dataOutputAllowedSet                     (),
    m_totallyMute                              (UQ_MH_SG_TOTALLY_MUTE_ODV),
    m_initialPositionDataInputFileName         (UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV),
    m_initialPositionDataInputFileType         (UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV),
    m_initialProposalCovMatrixDataInputFileName(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV),
    m_initialProposalCovMatrixDataInputFileType(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV),
  //m_parameterDisabledSet                     (),
    m_rawChainDataInputFileName                (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV),
    m_rawChainDataInputFileType                (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV),
    m_rawChainSize                             (UQ_MH_SG_RAW_CHAIN_SIZE_ODV),
    m_rawChainGenerateExtra                    (UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV),
    m_rawChainDisplayPeriod                    (UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV),
    m_rawChainMeasureRunTimes                  (UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV),
    m_rawChainDataOutputPeriod                 (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV),
    m_rawChainDataOutputFileName               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
    m_rawChainDataOutputFileType               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
    m_rawChainDataOutputAllowAll               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_rawChainDataOutputAllowedSet             (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_rawChainComputeStats                     (UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV),
#endif
    m_filteredChainGenerate                    (UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV),
    m_filteredChainDiscardedPortion            (UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
    m_filteredChainLag                         (UQ_MH_SG_FILTERED_CHAIN_LAG_ODV),
    m_filteredChainDataOutputFileName          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
    m_filteredChainDataOutputFileType          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
    m_filteredChainDataOutputAllowAll          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_filteredChainDataOutputAllowedSet        (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_filteredChainComputeStats                (UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV),
#endif
    m_displayCandidates                        (UQ_MH_SG_DISPLAY_CANDIDATES_ODV),
    m_putOutOfBoundsInChain                    (UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV),
    m_tkUseLocalHessian                        (UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV),
    m_tkUseNewtonComponent                     (UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV),
    m_drMaxNumExtraStages                      (UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV),
    m_drScalesForExtraStages                   (0),
    m_drDuringAmNonAdaptiveInt                 (UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV),
    m_amKeepInitialMatrix                      (UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV),
    m_amInitialNonAdaptInterval                (UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV),
    m_amAdaptInterval                          (UQ_MH_SG_AM_ADAPT_INTERVAL_ODV),
    m_amAdaptedMatricesDataOutputPeriod        (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV),
    m_amAdaptedMatricesDataOutputFileName      (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV),
    m_amAdaptedMatricesDataOutputFileType      (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV),
    m_amAdaptedMatricesDataOutputAllowAll      (UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV),
  //m_amAdaptedMatricesDataOutputAllowedSet    (),
    m_amEta                                    (UQ_MH_SG_AM_ETA_ODV),
    m_amEpsilon                                (UQ_MH_SG_AM_EPSILON_ODV),
    m_enableBrooksGelmanConvMonitor            (UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR),
    m_BrooksGelmanLag                          (UQ_MH_SG_BROOKS_GELMAN_LAG),
    m_outputLogLikelihood                      (UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD),
    m_outputLogTarget                          (UQ_MH_SG_OUTPUT_LOG_TARGET),
    m_doLogitTransform                         (UQ_MH_SG_DO_LOGIT_TRANSFORM),
    m_algorithm                                (UQ_MH_SG_ALGORITHM),
    m_tk                                       (UQ_MH_SG_TK),
    m_updateInterval                           (UQ_MH_SG_UPDATE_INTERVAL),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_alternativeRawSsOptionsValues            (),
    m_alternativeFilteredSsOptionsValues       (),
#endif
    m_env(env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help                                      (m_prefix + "help"                                      ),
    m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
    m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
    m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
    m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
    m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
    m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
    m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
    m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
    m_option_listOfDisabledParameters                  (m_prefix + "listOfDisabledParameters"                  ),
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
    m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "am_adaptedMatrices_dataOutputPeriod"       ),
    m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "am_adaptedMatrices_dataOutputFileName"     ),
    m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "am_adaptedMatrices_dataOutputFileType"     ),
    m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "am_adaptedMatrices_dataOutputAllowAll"     ),
    m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "am_adaptedMatrices_dataOutputAllowedSet"   ),
    m_option_am_eta                                    (m_prefix + "am_eta"                                    ),
    m_option_am_epsilon                                (m_prefix + "am_epsilon"                                ),
    m_option_enableBrooksGelmanConvMonitor             (m_prefix + "enableBrooksGelmanConvMonitor"             ),
    m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           ),
    m_option_outputLogLikelihood                       (m_prefix + "outputLogLikelihood"                       ),
    m_option_outputLogTarget                           (m_prefix + "outputLogTarget"                           ),
    m_option_doLogitTransform                          (m_prefix + "doLogitTransform"                          ),
    m_option_algorithm                                 (m_prefix + "algorithm"                                 ),
    m_option_tk                                        (m_prefix + "tk"                                        ),
    m_option_updateInterval                            (m_prefix + "updateInterval"                            )
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativeRawSsOptionsValues     ) m_alternativeRawSsOptionsValues      = *alternativeRawSsOptionsValues;
  if (alternativeFilteredSsOptionsValues) m_alternativeFilteredSsOptionsValues = *alternativeFilteredSsOptionsValues;
#endif

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string >(m_option_help,                                       UQ_MH_SG_HELP,                                                 "produce help msg for Bayesian Metropolis-Hastings"          );
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,                         UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                           , "name of generic output file"                                );
  m_parser->registerOption<bool        >(m_option_dataOutputAllowAll,                         UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV                           , "allow all subEnvs write to a generic output file"           );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,                       UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV                         , "subEnvs that will write to generic output file"             );
  m_parser->registerOption<bool        >(m_option_totallyMute,                                UQ_MH_SG_TOTALLY_MUTE_ODV                                    , "totally mute (no printout msg)"                             );
  m_parser->registerOption<std::string >(m_option_initialPosition_dataInputFileName,          UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV           , "name of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_initialPosition_dataInputFileType,          UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV           , "type of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileName, UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV, "name of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileType, UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV, "type of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_listOfDisabledParameters,                   UQ_MH_SG_LIST_OF_DISABLED_PARAMETERS_ODV                     , "list of disabled parameters"                                );
  m_parser->registerOption<std::string >(m_option_rawChain_dataInputFileName,                 UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV                  , "name of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_rawChain_dataInputFileType,                 UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV                  , "type of input file for raw chain "                          );
  m_parser->registerOption<unsigned int>(m_option_rawChain_size,                              UQ_MH_SG_RAW_CHAIN_SIZE_ODV                                  , "size of raw chain"                                          );
  m_parser->registerOption<bool        >(m_option_rawChain_generateExtra,                     UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV                        , "generate extra information about raw chain"                 );
  m_parser->registerOption<unsigned int>(m_option_rawChain_displayPeriod,                     UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV                        , "period of msg display during raw chain generation"          );
  m_parser->registerOption<bool        >(m_option_rawChain_measureRunTimes,                   UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV                     , "measure run times"                                          );
  m_parser->registerOption<unsigned int>(m_option_rawChain_dataOutputPeriod,                  UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV                    , "period of msg display during raw chain generation"          );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputFileName,                UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV                 , "name of output file for raw chain "                         );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputFileType,                UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV                 , "type of output file for raw chain "                         );
  m_parser->registerOption<bool        >(m_option_rawChain_dataOutputAllowAll,                UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV                 , "allow all subEnvs to write raw chain to an output file"     );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputAllowedSet,              UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV               , "subEnvs that will write raw chain to output file"           );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_rawChain_computeStats,                      UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV                         , "compute statistics on raw chain"                            );
#endif
  m_parser->registerOption<bool        >(m_option_filteredChain_generate,                     UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV                         , "generate filtered chain"                                    );
  m_parser->registerOption<double      >(m_option_filteredChain_discardedPortion,             UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV                , "initial discarded portion for chain filtering"              );
  m_parser->registerOption<unsigned int>(m_option_filteredChain_lag,                          UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                              , "spacing for chain filtering"                                );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputFileName,           UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV            , "name of output file for filtered chain"                     );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputFileType,           UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV            , "type of output file for filtered chain"                     );
  m_parser->registerOption<bool        >(m_option_filteredChain_dataOutputAllowAll,           UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV            , "allow all subEnvs to write filt chain to an output file"    );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputAllowedSet,         UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV          , "subEnvs that will write filt chain to output file"          );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_filteredChain_computeStats,                 UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV                    , "compute statistics on filtered chain"                       );
#endif
  m_parser->registerOption<bool        >(m_option_displayCandidates,                          UQ_MH_SG_DISPLAY_CANDIDATES_ODV                              , "display candidates in the core MH algorithm"                );
  m_parser->registerOption<bool        >(m_option_putOutOfBoundsInChain,                      UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV                      , "put 'out of bound' candidates in chain as well"             );
  m_parser->registerOption<bool        >(m_option_tk_useLocalHessian,                         UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                            , "'proposal' use local Hessian"                               );
  m_parser->registerOption<bool        >(m_option_tk_useNewtonComponent,                      UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV                         , "'proposal' use Newton component"                            );
  m_parser->registerOption<unsigned int>(m_option_dr_maxNumExtraStages,                       UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV                         , "'dr' maximum number of extra stages"                        );
  m_parser->registerOption<std::string >(m_option_dr_listOfScalesForExtraStages,              UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV              , "'dr' scales for prop cov matrices from 2nd stage on"        );
  m_parser->registerOption<bool        >(m_option_dr_duringAmNonAdaptiveInt,                  UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV                   , "'dr' used during 'am' non adaptive interval"                );
  m_parser->registerOption<bool        >(m_option_am_keepInitialMatrix,                       UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV                          , "'am' keep initial (given) matrix"                           );
  m_parser->registerOption<unsigned int>(m_option_am_initialNonAdaptInterval,                 UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                           , "'am' initial non adaptation interval"                       );
  m_parser->registerOption<unsigned int>(m_option_am_adaptInterval,                           UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                               , "'am' adaptation interval"                                   );
  m_parser->registerOption<unsigned int>(m_option_am_adaptedMatrices_dataOutputPeriod,        UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV          , "period for outputting 'am' adapted matrices"                );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileName,      UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV       , "name of output file for 'am' adapted matrices"              );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileType,      UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV       , "type of output file for 'am' adapted matrices"              );
  m_parser->registerOption<bool        >(m_option_am_adaptedMatrices_dataOutputAllowAll,      UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV       , "type of output file for 'am' adapted matrices"              );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputAllowedSet,    UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV     , "type of output file for 'am' adapted matrices"              );
  m_parser->registerOption<double      >(m_option_am_eta,                                     UQ_MH_SG_AM_ETA_ODV                                          , "'am' eta"                                                   );
  m_parser->registerOption<double      >(m_option_am_epsilon,                                 UQ_MH_SG_AM_EPSILON_ODV                                      , "'am' epsilon"                                               );
  m_parser->registerOption<unsigned int>(m_option_enableBrooksGelmanConvMonitor,              UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR                   , "assess convergence using Brooks-Gelman metric"              );
  m_parser->registerOption<unsigned int>(m_option_BrooksGelmanLag,                            UQ_MH_SG_BROOKS_GELMAN_LAG                                   , "number of chain positions before starting to compute metric");
  m_parser->registerOption<bool        >(m_option_outputLogLikelihood,                        UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD                               , "flag to toggle output of log likelihood values"             );
  m_parser->registerOption<bool        >(m_option_outputLogTarget,                            UQ_MH_SG_OUTPUT_LOG_TARGET                                   , "flag to toggle output of log target values"                 );
  m_parser->registerOption<bool        >(m_option_doLogitTransform,                           UQ_MH_SG_DO_LOGIT_TRANSFORM                                  , "flag to toggle logit transform for bounded domains"         );
  m_parser->registerOption<std::string >(m_option_algorithm,                                  UQ_MH_SG_ALGORITHM                                           , "which MCMC algorithm to use"                                );
  m_parser->registerOption<std::string >(m_option_tk,                                         UQ_MH_SG_TK                                                  , "which MCMC transition kernel to use"                        );
  m_parser->registerOption<unsigned int>(m_option_updateInterval,                             UQ_MH_SG_UPDATE_INTERVAL                                     , "how often to call updateTK method"                          );

  m_parser->scanInputFile();

  m_parser->getOption<std::string >(m_option_help,                                       m_help);
  m_parser->getOption<std::string >(m_option_dataOutputFileName,                         m_dataOutputFileName);
  m_parser->getOption<bool        >(m_option_dataOutputAllowAll,                         m_dataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,                       m_dataOutputAllowedSet);
  m_parser->getOption<bool        >(m_option_totallyMute,                                m_totallyMute);
  m_parser->getOption<std::string >(m_option_initialPosition_dataInputFileName,          m_initialPositionDataInputFileName);
  m_parser->getOption<std::string >(m_option_initialPosition_dataInputFileType,          m_initialPositionDataInputFileType);
  m_parser->getOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName);
  m_parser->getOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType);
  m_parser->getOption<std::set<unsigned int> >(m_option_listOfDisabledParameters,                   m_parameterDisabledSet);
  m_parser->getOption<std::string >(m_option_rawChain_dataInputFileName,                 m_rawChainDataInputFileName);
  m_parser->getOption<std::string >(m_option_rawChain_dataInputFileType,                 m_rawChainDataInputFileType);
  m_parser->getOption<unsigned int>(m_option_rawChain_size,                              m_rawChainSize);
  m_parser->getOption<bool        >(m_option_rawChain_generateExtra,                     m_rawChainGenerateExtra);
  m_parser->getOption<unsigned int>(m_option_rawChain_displayPeriod,                     m_rawChainDisplayPeriod);
  m_parser->getOption<bool        >(m_option_rawChain_measureRunTimes,                   m_rawChainMeasureRunTimes);
  m_parser->getOption<unsigned int>(m_option_rawChain_dataOutputPeriod,                  m_rawChainDataOutputPeriod);
  m_parser->getOption<std::string >(m_option_rawChain_dataOutputFileName,                m_rawChainDataOutputFileName);
  m_parser->getOption<std::string >(m_option_rawChain_dataOutputFileType,                m_rawChainDataOutputFileType);
  m_parser->getOption<bool        >(m_option_rawChain_dataOutputAllowAll,                m_rawChainDataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_rawChain_dataOutputAllowedSet,              m_rawChainDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool        >(m_option_rawChain_computeStats,                      m_rawChain_computeStats);
#endif
  m_parser->getOption<bool        >(m_option_filteredChain_generate,                     m_filteredChainGenerate);
  m_parser->getOption<double      >(m_option_filteredChain_discardedPortion,             m_filteredChainDiscardedPortion);
  m_parser->getOption<unsigned int>(m_option_filteredChain_lag,                          m_filteredChainLag);
  m_parser->getOption<std::string >(m_option_filteredChain_dataOutputFileName,           m_filteredChainDataOutputFileName);
  m_parser->getOption<std::string >(m_option_filteredChain_dataOutputFileType,           m_filteredChainDataOutputFileType);
  m_parser->getOption<bool        >(m_option_filteredChain_dataOutputAllowAll,           m_filteredChainDataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_filteredChain_dataOutputAllowedSet,         m_filteredChainDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool        >(m_option_filteredChain_computeStats,                 m_filteredChain_computeStats);
#endif
  m_parser->getOption<bool        >(m_option_displayCandidates,                          m_displayCandidates);
  m_parser->getOption<bool        >(m_option_putOutOfBoundsInChain,                      m_putOutOfBoundsInChain);
  m_parser->getOption<bool        >(m_option_tk_useLocalHessian,                         m_tkUseLocalHessian);
  m_parser->getOption<bool        >(m_option_tk_useNewtonComponent,                      m_tkUseNewtonComponent);
  m_parser->getOption<unsigned int>(m_option_dr_maxNumExtraStages,                       m_drMaxNumExtraStages);
  m_parser->getOption<std::vector<double> >(m_option_dr_listOfScalesForExtraStages,              m_drScalesForExtraStages);
  m_parser->getOption<bool        >(m_option_dr_duringAmNonAdaptiveInt,                  m_drDuringAmNonAdaptiveInt);
  m_parser->getOption<bool        >(m_option_am_keepInitialMatrix,                       m_amKeepInitialMatrix);
  m_parser->getOption<unsigned int>(m_option_am_initialNonAdaptInterval,                 m_amInitialNonAdaptInterval);
  m_parser->getOption<unsigned int>(m_option_am_adaptInterval,                           m_amAdaptInterval);
  m_parser->getOption<unsigned int>(m_option_am_adaptedMatrices_dataOutputPeriod,        m_amAdaptedMatricesDataOutputPeriod);
  m_parser->getOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileName,      m_amAdaptedMatricesDataOutputFileName);
  m_parser->getOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileType,      m_amAdaptedMatricesDataOutputFileType);
  m_parser->getOption<bool        >(m_option_am_adaptedMatrices_dataOutputAllowAll,      m_amAdaptedMatricesDataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_am_adaptedMatrices_dataOutputAllowedSet,    m_amAdaptedMatricesDataOutputAllowedSet);
  m_parser->getOption<double      >(m_option_am_eta,                                     m_amEta);
  m_parser->getOption<double      >(m_option_am_epsilon,                                 m_amEpsilon);
  m_parser->getOption<unsigned int>(m_option_enableBrooksGelmanConvMonitor,              m_enableBrooksGelmanConvMonitor);
  m_parser->getOption<unsigned int>(m_option_BrooksGelmanLag,                            m_BrooksGelmanLag);
  m_parser->getOption<bool        >(m_option_outputLogLikelihood,                        m_outputLogLikelihood);
  m_parser->getOption<bool        >(m_option_outputLogTarget,                            m_outputLogTarget);
  m_parser->getOption<bool        >(m_option_doLogitTransform,                           m_doLogitTransform);
  m_parser->getOption<std::string >(m_option_algorithm,                                  m_algorithm);
  m_parser->getOption<std::string >(m_option_tk,                                         m_tk);
  m_parser->getOption<unsigned int>(m_option_updateInterval,                             m_updateInterval);
#else
  m_help = m_env->input()(m_option_help, UQ_MH_SG_HELP);
  m_dataOutputFileName = m_env->input()(m_option_dataOutputFileName, UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV);
  m_dataOutputAllowAll = m_env->input()(m_option_dataOutputAllowAll, UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV);

  // UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = m_env->input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

  m_totallyMute = m_env->input()(m_option_totallyMute, UQ_MH_SG_TOTALLY_MUTE_ODV);
  m_initialPositionDataInputFileName = m_env->input()(m_option_initialPosition_dataInputFileName, UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV);
  m_initialPositionDataInputFileType = m_env->input()(m_option_initialPosition_dataInputFileType, UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV);
  m_initialProposalCovMatrixDataInputFileName = m_env->input()(m_option_initialProposalCovMatrix_dataInputFileName, UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV);
  m_initialProposalCovMatrixDataInputFileType = m_env->input()(m_option_initialProposalCovMatrix_dataInputFileType, UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV);

  // UQ_MH_SG_LIST_OF_DISABLED_PARAMETERS_ODV is the empty set (string) by
  // default
  size = m_env->input().vector_variable_size(m_option_listOfDisabledParameters);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int disabled = m_env->input()(m_option_listOfDisabledParameters, i, i);
    m_parameterDisabledSet.insert(disabled);
  }

  m_rawChainDataInputFileName = m_env->input()(m_option_rawChain_dataInputFileName, UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV);
  m_rawChainDataInputFileType = m_env->input()(m_option_rawChain_dataInputFileType, UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV);
  m_rawChainSize = m_env->input()(m_option_rawChain_size, UQ_MH_SG_RAW_CHAIN_SIZE_ODV);
  m_rawChainGenerateExtra = m_env->input()(m_option_rawChain_generateExtra, UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV);
  m_rawChainDisplayPeriod = m_env->input()(m_option_rawChain_displayPeriod, UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV);
  m_rawChainMeasureRunTimes = m_env->input()(m_option_rawChain_measureRunTimes, UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV);
  m_rawChainDataOutputPeriod = m_env->input()(m_option_rawChain_dataOutputPeriod, UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV);
  m_rawChainDataOutputFileName = m_env->input()(m_option_rawChain_dataOutputFileName, UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV);
  m_rawChainDataOutputFileType = m_env->input()(m_option_rawChain_dataOutputFileType, UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV);
  m_rawChainDataOutputAllowAll = m_env->input()(m_option_rawChain_dataOutputAllowAll, UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV);

  // UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_rawChain_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_rawChain_dataOutputAllowedSet, i, i);
    m_rawChainDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChain_computeStats = m_env->input()(m_option_rawChain_computeStats, UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV);
#endif
  m_filteredChainGenerate = m_env->input()(m_option_filteredChain_generate, UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV);
  m_filteredChainDiscardedPortion = m_env->input()(m_option_filteredChain_discardedPortion, UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV);
  m_filteredChainLag = m_env->input()(m_option_filteredChain_lag, UQ_MH_SG_FILTERED_CHAIN_LAG_ODV);
  m_filteredChainDataOutputFileName = m_env->input()(m_option_filteredChain_dataOutputFileName, UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV);
  m_filteredChainDataOutputFileType = m_env->input()(m_option_filteredChain_dataOutputFileType, UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV);
  m_filteredChainDataOutputAllowAll = m_env->input()(m_option_filteredChain_dataOutputAllowAll, UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV);

  // UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_filteredChain_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_filteredChain_dataOutputAllowedSet, i, i);
    m_filteredChainDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_filteredChain_computeStats = m_env->input()(m_option_filteredChain_computeStats, UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV);
#endif
  m_displayCandidates = m_env->input()(m_option_displayCandidates, UQ_MH_SG_DISPLAY_CANDIDATES_ODV);
  m_putOutOfBoundsInChain = m_env->input()(m_option_putOutOfBoundsInChain, UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV);
  m_tkUseLocalHessian = m_env->input()(m_option_tk_useLocalHessian, UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV);
  m_tkUseNewtonComponent = m_env->input()(m_option_tk_useNewtonComponent, UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV);
  m_drMaxNumExtraStages = m_env->input()(m_option_dr_maxNumExtraStages, UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV);

  // UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_dr_listOfScalesForExtraStages);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_dr_listOfScalesForExtraStages, i, i);
    m_drScalesForExtraStages.push_back(allowed);
  }

  m_drDuringAmNonAdaptiveInt = m_env->input()(m_option_dr_duringAmNonAdaptiveInt, UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV);
  m_amKeepInitialMatrix = m_env->input()(m_option_am_keepInitialMatrix, UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV);
  m_amInitialNonAdaptInterval = m_env->input()(m_option_am_initialNonAdaptInterval, UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV);
  m_amAdaptInterval = m_env->input()(m_option_am_adaptInterval, UQ_MH_SG_AM_ADAPT_INTERVAL_ODV);
  m_amAdaptedMatricesDataOutputPeriod = m_env->input()(m_option_am_adaptedMatrices_dataOutputPeriod, UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV);
  m_amAdaptedMatricesDataOutputFileName = m_env->input()(m_option_am_adaptedMatrices_dataOutputFileName, UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV);
  m_amAdaptedMatricesDataOutputFileType = m_env->input()(m_option_am_adaptedMatrices_dataOutputFileType, UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV);
  m_amAdaptedMatricesDataOutputAllowAll = m_env->input()(m_option_am_adaptedMatrices_dataOutputAllowAll, UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV);

  // UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_am_adaptedMatrices_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_am_adaptedMatrices_dataOutputAllowedSet, i, i);
    m_amAdaptedMatricesDataOutputAllowedSet.insert(allowed);
  }

  m_amEta = m_env->input()(m_option_am_eta, UQ_MH_SG_AM_ETA_ODV);
  m_amEpsilon = m_env->input()(m_option_am_epsilon, UQ_MH_SG_AM_EPSILON_ODV);
  m_enableBrooksGelmanConvMonitor = m_env->input()(m_option_enableBrooksGelmanConvMonitor, UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR);
  m_BrooksGelmanLag = m_env->input()(m_option_BrooksGelmanLag, UQ_MH_SG_BROOKS_GELMAN_LAG);
  m_outputLogLikelihood = m_env->input()(m_option_outputLogLikelihood, UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD);
  m_outputLogTarget = m_env->input()(m_option_outputLogTarget, UQ_MH_SG_OUTPUT_LOG_TARGET);
  m_doLogitTransform = m_env->input()(m_option_doLogitTransform, UQ_MH_SG_DO_LOGIT_TRANSFORM);
  m_algorithm = m_env->input()(m_option_algorithm, UQ_MH_SG_ALGORITHM);
  m_tk = m_env->input()(m_option_tk, UQ_MH_SG_TK);
  m_updateInterval = m_env->input()(m_option_updateInterval, UQ_MH_SG_UPDATE_INTERVAL);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions(env);
}
// Copy constructor----------------------------------
MhOptionsValues::MhOptionsValues(const MhOptionsValues& src)
{
  this->copy(src);
}
// Destructor ---------------------------------------
MhOptionsValues::~MhOptionsValues()
{
}
// Set methods --------------------------------------
MhOptionsValues&
MhOptionsValues::operator=(const MhOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}
// Private methods-----------------------------------
void
MhOptionsValues::checkOptions(const BaseEnvironment * env)
{
  if (m_dataOutputAllowAll) {
    // So, ignore the 'set' option
    m_dataOutputAllowedSet.clear();
    m_dataOutputAllowedSet.insert(env->subId());
  }

  if (m_rawChainDataOutputAllowAll) {
    // Again, ignore the set
    m_rawChainDataOutputAllowedSet.clear();
    m_rawChainDataOutputAllowedSet.insert(env->subId());
  }

  if (m_filteredChainGenerate == true) {
    queso_require_greater_equal_msg(m_filteredChainLag, 2, "option `" << m_option_filteredChain_lag << "` must be at least 2");
  }

  if (m_filteredChainDataOutputAllowAll) {
    m_filteredChainDataOutputAllowedSet.clear();
    m_filteredChainDataOutputAllowedSet.insert(env->subId());
  }

  // If max is bigger than the list provided, then pad with ones
  if (m_drMaxNumExtraStages > 0) {
    unsigned int size = m_drScalesForExtraStages.size();
    if (m_drMaxNumExtraStages > size) {
      for (unsigned int i = size; i < m_drMaxNumExtraStages; i++) {
        double scale = 1.0;
        m_drScalesForExtraStages.push_back(scale);
      }
    }
  }
  else {
    m_drScalesForExtraStages.clear();
  }

  if (m_amAdaptedMatricesDataOutputAllowAll) {
    m_amAdaptedMatricesDataOutputAllowedSet.clear();
    m_amAdaptedMatricesDataOutputAllowedSet.insert(env->subId());
  }

  if ((m_tk == "random_walk") && (m_algorithm == "logit_random_walk")) {
      queso_error_msg("random_walk transition kernel and logit_random_walk algorithm are incompatible options");
  }

  if ((m_tk == "logit_random_walk") && (m_algorithm == "random_walk")) {
    queso_error_msg("logit_random_walk transition kernel and random_walk algorithm are incompatible options");
  }

  if (m_tk == "random_walk") {
    queso_require_equal_to_msg(
        m_doLogitTransform,
        0,
        "logit transform must be off to use random_walk");
    queso_require_equal_to_msg(
        m_tkUseLocalHessian,
        0,
        "local Hessian must be off to use random_walk");
  }

  if (m_tk == "logit_random_walk") {
    queso_require_equal_to_msg(
        m_doLogitTransform,
        1,
        "logit transform must be on to use logit_random_walk");
    queso_require_equal_to_msg(
        m_tkUseLocalHessian,
        0,
        "local Hessian must be off to use logit_random_walk");
  }

  if (m_tk == "stochastic_newton") {
    queso_require_equal_to_msg(
        m_doLogitTransform,
        0,
        "logit transform must be off to use stochastic_newton");
    queso_require_equal_to_msg(
        m_tkUseLocalHessian,
        1,
        "local Hessian must be on to use stochastic_newton");
  }

}

void
MhOptionsValues::copy(const MhOptionsValues& src)
{
  m_dataOutputFileName                        = src.m_dataOutputFileName;
  m_dataOutputAllowAll                        = src.m_dataOutputAllowAll;
  m_dataOutputAllowedSet                      = src.m_dataOutputAllowedSet;
  m_totallyMute                               = src.m_totallyMute;
  m_initialPositionDataInputFileName          = src.m_initialPositionDataInputFileName;
  m_initialPositionDataInputFileType          = src.m_initialPositionDataInputFileType;
  m_initialProposalCovMatrixDataInputFileName = src.m_initialProposalCovMatrixDataInputFileName;
  m_initialProposalCovMatrixDataInputFileType = src.m_initialProposalCovMatrixDataInputFileType;
  m_parameterDisabledSet                      = src.m_parameterDisabledSet;
  m_rawChainDataInputFileName                 = src.m_rawChainDataInputFileName;
  m_rawChainDataInputFileType                 = src.m_rawChainDataInputFileType;
  m_rawChainSize                              = src.m_rawChainSize;
  m_rawChainGenerateExtra                     = src.m_rawChainGenerateExtra;
  m_rawChainDisplayPeriod                     = src.m_rawChainDisplayPeriod;
  m_rawChainMeasureRunTimes                   = src.m_rawChainMeasureRunTimes;
  m_rawChainDataOutputPeriod                  = src.m_rawChainDataOutputPeriod;
  m_rawChainDataOutputFileName                = src.m_rawChainDataOutputFileName;
  m_rawChainDataOutputFileType                = src.m_rawChainDataOutputFileType;
  m_rawChainDataOutputAllowAll                = src.m_rawChainDataOutputAllowAll;
  m_rawChainDataOutputAllowedSet              = src.m_rawChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChainComputeStats                      = src.m_rawChainComputeStats;
#endif
//m_rawChainStatisticalOptionsObj             = src.m_rawChainStatisticalOptionsObj; // dakota
//m_rawChainStatOptsInstantiated              = src.m_rawChainStatOptsInstantiated; // dakota
  m_filteredChainGenerate                     = src.m_filteredChainGenerate;
  m_filteredChainDiscardedPortion             = src.m_filteredChainDiscardedPortion;
  m_filteredChainLag                          = src.m_filteredChainLag;
  m_filteredChainDataOutputFileName           = src.m_filteredChainDataOutputFileName;
  m_filteredChainDataOutputFileType           = src.m_filteredChainDataOutputFileType;
  m_filteredChainDataOutputAllowAll           = src.m_filteredChainDataOutputAllowAll;
  m_filteredChainDataOutputAllowedSet         = src.m_filteredChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_filteredChainComputeStats                 = src.m_filteredChainComputeStats;
#endif
//m_filteredChainStatisticalOptionsObj        = src.m_filteredChainStatisticalOptionsObj; // dakota
//m_filteredChainStatOptsInstantiated         = src.m_filteredChainStatOptsInstantiated; // dakota
  m_displayCandidates                         = src.m_displayCandidates;
  m_putOutOfBoundsInChain                     = src.m_putOutOfBoundsInChain;
  m_tkUseLocalHessian                         = src.m_tkUseLocalHessian;
  m_tkUseNewtonComponent                      = src.m_tkUseNewtonComponent;
  m_drMaxNumExtraStages                       = src.m_drMaxNumExtraStages;
  m_drScalesForExtraStages                    = src.m_drScalesForExtraStages;
  m_drDuringAmNonAdaptiveInt                  = src.m_drDuringAmNonAdaptiveInt;
  m_amKeepInitialMatrix                       = src.m_amKeepInitialMatrix;
  m_amInitialNonAdaptInterval                 = src.m_amInitialNonAdaptInterval;
  m_amAdaptInterval                           = src.m_amAdaptInterval;
  m_amAdaptedMatricesDataOutputPeriod         = src.m_amAdaptedMatricesDataOutputPeriod;
  m_amAdaptedMatricesDataOutputFileName       = src.m_amAdaptedMatricesDataOutputFileName;
  m_amAdaptedMatricesDataOutputFileType       = src.m_amAdaptedMatricesDataOutputFileType;
  m_amAdaptedMatricesDataOutputAllowAll       = src.m_amAdaptedMatricesDataOutputAllowAll;
  m_amAdaptedMatricesDataOutputAllowedSet     = src.m_amAdaptedMatricesDataOutputAllowedSet;
  m_amEta                                     = src.m_amEta;
  m_amEpsilon                                 = src.m_amEpsilon;
  m_enableBrooksGelmanConvMonitor             = src.m_enableBrooksGelmanConvMonitor;
  m_BrooksGelmanLag                           = src.m_BrooksGelmanLag;
  m_outputLogLikelihood                       = src.m_outputLogLikelihood;
  m_outputLogTarget                           = src.m_outputLogTarget;
  m_doLogitTransform                          = src.m_doLogitTransform;
  m_algorithm                                 = src.m_algorithm;
  m_tk                                        = src.m_tk;
  m_updateInterval                            = src.m_updateInterval;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_alternativeRawSsOptionsValues             = src.m_alternativeRawSsOptionsValues;
  m_alternativeFilteredSsOptionsValues        = src.m_alternativeFilteredSsOptionsValues;
#endif
  return;
}

std::ostream & operator<<(std::ostream & os, const MhOptionsValues & obj)
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  os <<         obj.m_option_dataOutputFileName                         << " = " << obj.m_dataOutputFileName
     << "\n" << obj.m_option_dataOutputAllowAll                         << " = " << obj.m_dataOutputAllowAll
     << "\n" << obj.m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_dataOutputAllowedSet.begin(); setIt != obj.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << obj.m_option_totallyMute                                << " = " << obj.m_totallyMute
     << "\n" << obj.m_option_initialPosition_dataInputFileName          << " = " << obj.m_initialPositionDataInputFileName
     << "\n" << obj.m_option_initialPosition_dataInputFileType          << " = " << obj.m_initialPositionDataInputFileType
     << "\n" << obj.m_option_initialProposalCovMatrix_dataInputFileName << " = " << obj.m_initialProposalCovMatrixDataInputFileName
     << "\n" << obj.m_option_initialProposalCovMatrix_dataInputFileType << " = " << obj.m_initialProposalCovMatrixDataInputFileType
     << "\n" << obj.m_option_listOfDisabledParameters                   << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_parameterDisabledSet.begin(); setIt != obj.m_parameterDisabledSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << obj.m_option_rawChain_dataInputFileName                 << " = " << obj.m_rawChainDataInputFileName
     << "\n" << obj.m_option_rawChain_dataInputFileType                 << " = " << obj.m_rawChainDataInputFileType
     << "\n" << obj.m_option_rawChain_size                              << " = " << obj.m_rawChainSize
     << "\n" << obj.m_option_rawChain_generateExtra                     << " = " << obj.m_rawChainGenerateExtra
     << "\n" << obj.m_option_rawChain_displayPeriod                     << " = " << obj.m_rawChainDisplayPeriod
     << "\n" << obj.m_option_rawChain_measureRunTimes                   << " = " << obj.m_rawChainMeasureRunTimes
     << "\n" << obj.m_option_rawChain_dataOutputPeriod                  << " = " << obj.m_rawChainDataOutputPeriod
     << "\n" << obj.m_option_rawChain_dataOutputFileName                << " = " << obj.m_rawChainDataOutputFileName
     << "\n" << obj.m_option_rawChain_dataOutputFileType                << " = " << obj.m_rawChainDataOutputFileType
     << "\n" << obj.m_option_rawChain_dataOutputAllowAll                << " = " << obj.m_rawChainDataOutputAllowAll
     << "\n" << obj.m_option_rawChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_rawChainDataOutputAllowedSet.begin(); setIt != obj.m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << obj.m_option_rawChain_computeStats                      << " = " << obj.m_rawChainComputeStats
#endif
     << "\n" << obj.m_option_filteredChain_generate                     << " = " << obj.m_filteredChainGenerate
     << "\n" << obj.m_option_filteredChain_discardedPortion             << " = " << obj.m_filteredChainDiscardedPortion
     << "\n" << obj.m_option_filteredChain_lag                          << " = " << obj.m_filteredChainLag
     << "\n" << obj.m_option_filteredChain_dataOutputFileName           << " = " << obj.m_filteredChainDataOutputFileName
     << "\n" << obj.m_option_filteredChain_dataOutputFileType           << " = " << obj.m_filteredChainDataOutputFileType
     << "\n" << obj.m_option_filteredChain_dataOutputAllowAll           << " = " << obj.m_filteredChainDataOutputAllowAll
     << "\n" << obj.m_option_filteredChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_filteredChainDataOutputAllowedSet.begin(); setIt != obj.m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << obj.m_option_filteredChain_computeStats                 << " = " << obj.m_filteredChainComputeStats
#endif
     << "\n" << obj.m_option_displayCandidates                          << " = " << obj.m_displayCandidates
     << "\n" << obj.m_option_putOutOfBoundsInChain                      << " = " << obj.m_putOutOfBoundsInChain
     << "\n" << obj.m_option_tk_useLocalHessian                         << " = " << obj.m_tkUseLocalHessian
     << "\n" << obj.m_option_tk_useNewtonComponent                      << " = " << obj.m_tkUseNewtonComponent
     << "\n" << obj.m_option_dr_maxNumExtraStages                       << " = " << obj.m_drMaxNumExtraStages
     << "\n" << obj.m_option_dr_listOfScalesForExtraStages << " = ";
  for (unsigned int i = 0; i < obj.m_drScalesForExtraStages.size(); ++i) {
    os << obj.m_drScalesForExtraStages[i] << " ";
  }
  os << "\n" << obj.m_option_dr_duringAmNonAdaptiveInt                  << " = " << obj.m_drDuringAmNonAdaptiveInt
     << "\n" << obj.m_option_am_keepInitialMatrix                       << " = " << obj.m_amKeepInitialMatrix
     << "\n" << obj.m_option_am_initialNonAdaptInterval                 << " = " << obj.m_amInitialNonAdaptInterval
     << "\n" << obj.m_option_am_adaptInterval                           << " = " << obj.m_amAdaptInterval
     << "\n" << obj.m_option_am_adaptedMatrices_dataOutputPeriod        << " = " << obj.m_amAdaptedMatricesDataOutputPeriod
     << "\n" << obj.m_option_am_adaptedMatrices_dataOutputFileName      << " = " << obj.m_amAdaptedMatricesDataOutputFileName
     << "\n" << obj.m_option_am_adaptedMatrices_dataOutputFileType      << " = " << obj.m_amAdaptedMatricesDataOutputFileType
     << "\n" << obj.m_option_am_adaptedMatrices_dataOutputAllowAll      << " = " << obj.m_amAdaptedMatricesDataOutputAllowAll
     << "\n" << obj.m_option_am_adaptedMatrices_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_amAdaptedMatricesDataOutputAllowedSet.begin(); setIt != obj.m_amAdaptedMatricesDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << obj.m_option_am_eta                                     << " = " << obj.m_amEta
     << "\n" << obj.m_option_am_epsilon                                 << " = " << obj.m_amEpsilon
     << "\n" << obj.m_option_enableBrooksGelmanConvMonitor              << " = " << obj.m_enableBrooksGelmanConvMonitor
     << "\n" << obj.m_option_BrooksGelmanLag                            << " = " << obj.m_BrooksGelmanLag
     << "\n" << obj.m_option_outputLogLikelihood                        << " = " << obj.m_outputLogLikelihood
     << "\n" << obj.m_option_outputLogTarget                            << " = " << obj.m_outputLogTarget
     << "\n" << obj.m_option_doLogitTransform                           << " = " << obj.m_doLogitTransform
     << "\n" << obj.m_option_algorithm                                  << " = " << obj.m_algorithm
     << "\n" << obj.m_option_tk                                         << " = " << obj.m_tk
     << "\n" << obj.m_option_updateInterval                             << " = " << obj.m_updateInterval
     << std::endl;

  return os;
}

//---------------------------------------------------
// MetropolisHastingsSGOptions ---------------
//---------------------------------------------------

// Default constructor -----------------------------
MetropolisHastingsSGOptions::MetropolisHastingsSGOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_ov                                               (NULL,NULL), // dakota
  m_rawChainStatisticalOptionsObj                    (NULL),
  m_rawChainStatOptsInstantiated                     (false),
  m_filteredChainStatisticalOptionsObj               (NULL),
  m_filteredChainStatOptsInstantiated                (false),
#else
  m_ov                                               (),
#endif
  m_prefix                                           ((std::string)(prefix) + "mh_"),
  m_env                                              (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                                      (new boost::program_options::options_description("Bayesian Metropolis-Hastings options")),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                                      (m_prefix + "help"                                      ),
  m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
  m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
  m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
  m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
  m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
  m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
  m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
  m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
  m_option_listOfDisabledParameters                  (m_prefix + "listOfDisabledParameters"                  ),
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
  m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "am_adaptedMatrices_dataOutputPeriod"       ),
  m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "am_adaptedMatrices_dataOutputFileName"     ),
  m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "am_adaptedMatrices_dataOutputFileType"     ),
  m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "am_adaptedMatrices_dataOutputAllowAll"     ),
  m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "am_adaptedMatrices_dataOutputAllowedSet"   ),
  m_option_am_eta                                    (m_prefix + "am_eta"                                    ),
  m_option_am_epsilon                                (m_prefix + "am_epsilon"                                ),
  m_option_enableBrooksGelmanConvMonitor             (m_prefix + "enableBrooksGelmanConvMonitor"             ),
  m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           ),
  m_option_outputLogLikelihood                       (m_prefix + "outputLogLikelihood"                       ),
  m_option_outputLogTarget                           (m_prefix + "outputLogTarget"                           ),
  m_option_doLogitTransform                          (m_prefix + "doLogitTransform"                          ),
  m_option_algorithm                                 (m_prefix + "algorithm"                                 ),
  m_option_tk                                        (m_prefix + "tk"                                        ),
  m_option_updateInterval                            (m_prefix + "updateInterval"                            )
{
  queso_deprecated();

  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the absence of an options input file"));
}
// Constructor 2------------------------------------
MetropolisHastingsSGOptions::MetropolisHastingsSGOptions(
  const BaseEnvironment& env,
  const char*                   prefix,
  const MhOptionsValues& alternativeOptionsValues)
  :
  m_ov                                               (alternativeOptionsValues),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChainStatisticalOptionsObj                    (NULL),
  m_rawChainStatOptsInstantiated                     (false),
  m_filteredChainStatisticalOptionsObj               (NULL),
  m_filteredChainStatOptsInstantiated                (false),
#endif
  m_prefix                                           ((std::string)(prefix) + "mh_"),
  m_env                                              (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                                      (),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                                      (m_prefix + "help"                                      ),
  m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
  m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
  m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
  m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
  m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
  m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
  m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
  m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
  m_option_listOfDisabledParameters                  (m_prefix + "listOfDisabledParameters"                  ),
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
  m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "am_adaptedMatrices_dataOutputPeriod"       ),
  m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "am_adaptedMatrices_dataOutputFileName"     ),
  m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "am_adaptedMatrices_dataOutputFileType"     ),
  m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "am_adaptedMatrices_dataOutputAllowAll"     ),
  m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "am_adaptedMatrices_dataOutputAllowedSet"   ),
  m_option_am_eta                                    (m_prefix + "am_eta"                                    ),
  m_option_am_epsilon                                (m_prefix + "am_epsilon"                                ),
  m_option_enableBrooksGelmanConvMonitor             (m_prefix + "enableBrooksGelmanConvMonitor"             ),
  m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           ),
  m_option_outputLogLikelihood                       (m_prefix + "outputLogLikelihood"                       ),
  m_option_outputLogTarget                           (m_prefix + "outputLogTarget"                           ),
  m_option_doLogitTransform                          (m_prefix + "doLogitTransform"                          ),
  m_option_algorithm                                 (m_prefix + "algorithm"                                 ),
  m_option_tk                                        (m_prefix + "tk"                                        ),
  m_option_updateInterval                            (m_prefix + "updateInterval"                            )
{
  queso_deprecated();

  queso_require_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the existence of an options input file"));

  if ((m_env.subDisplayFile() != NULL ) &&
      (m_ov.m_totallyMute     == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSGOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_ov.m_rawChainComputeStats) {
    m_rawChainStatisticalOptionsObj = new SequenceStatisticalOptions(m_env,m_prefix + "rawChain_",m_ov.m_alternativeRawSsOptionsValues);
    m_rawChainStatOptsInstantiated  = true;
  }
  if (m_ov.m_filteredChainComputeStats) {
    m_filteredChainStatisticalOptionsObj = new SequenceStatisticalOptions(m_env,m_prefix + "filteredChain_",m_ov.m_alternativeFilteredSsOptionsValues);
    m_filteredChainStatOptsInstantiated  = true;
  }
#endif
}
// Copy constructor---------------------------------
MetropolisHastingsSGOptions::MetropolisHastingsSGOptions(
  const MLSamplingLevelOptions& mlOptions)
  :
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_ov                                               (NULL,NULL), // dakota
  m_rawChainStatisticalOptionsObj                    (NULL),
  m_rawChainStatOptsInstantiated                     (false),
  m_filteredChainStatisticalOptionsObj               (NULL),
  m_filteredChainStatOptsInstantiated                (false),
#else
  m_ov                                               (),
#endif
  m_prefix                                           (mlOptions.m_prefix),
  m_env                                              (mlOptions.env()),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                                      (),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                                      (m_prefix + "help"                                      ),
  m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
  m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
  m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
  m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
  m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
  m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
  m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
  m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
  m_option_listOfDisabledParameters                  (m_prefix + "listOfDisabledParameters"                  ),
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
  m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "am_adaptedMatrices_dataOutputPeriod"       ),
  m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "am_adaptedMatrices_dataOutputFileName"     ),
  m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "am_adaptedMatrices_dataOutputFileType"     ),
  m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "am_adaptedMatrices_dataOutputAllowAll"     ),
  m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "am_adaptedMatrices_dataOutputAllowedSet"   ),
  m_option_am_eta                                    (m_prefix + "am_eta"                                    ),
  m_option_am_epsilon                                (m_prefix + "am_epsilon"                                ),
  m_option_enableBrooksGelmanConvMonitor             (m_prefix + "enableBrooksGelmanConvMonitor"             ),
  m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           ),
  m_option_outputLogLikelihood                       (m_prefix + "outputLogLikelihood"                       ),
  m_option_outputLogTarget                           (m_prefix + "outputLogTarget"                           ),
  m_option_doLogitTransform                          (m_prefix + "doLogitTransform"                          ),
  m_option_algorithm                                 (m_prefix + "algorithm"                                 ),
  m_option_tk                                        (m_prefix + "tk"                                        ),
  m_option_updateInterval                            (m_prefix + "updateInterval"                            )
{
  queso_deprecated();

  m_ov.m_dataOutputFileName                        = mlOptions.m_dataOutputFileName;
  m_ov.m_dataOutputAllowAll                        = mlOptions.m_dataOutputAllowAll;
  m_ov.m_dataOutputAllowedSet                      = mlOptions.m_dataOutputAllowedSet;
  m_ov.m_totallyMute                               = mlOptions.m_totallyMute;
  m_ov.m_initialPositionDataInputFileName          = mlOptions.m_initialPositionDataInputFileName;
  m_ov.m_initialPositionDataInputFileType          = mlOptions.m_initialPositionDataInputFileType;
  m_ov.m_initialProposalCovMatrixDataInputFileName = mlOptions.m_initialProposalCovMatrixDataInputFileName;
  m_ov.m_initialProposalCovMatrixDataInputFileType = mlOptions.m_initialProposalCovMatrixDataInputFileType;
  m_ov.m_parameterDisabledSet                      = mlOptions.m_parameterDisabledSet;
  m_ov.m_rawChainDataInputFileName                 = mlOptions.m_rawChainDataInputFileName;
  m_ov.m_rawChainDataInputFileType                 = mlOptions.m_rawChainDataInputFileType;
  m_ov.m_rawChainSize                              = mlOptions.m_rawChainSize;
  m_ov.m_rawChainGenerateExtra                     = mlOptions.m_rawChainGenerateExtra;
  m_ov.m_rawChainDisplayPeriod                     = mlOptions.m_rawChainDisplayPeriod;
  m_ov.m_rawChainMeasureRunTimes                   = mlOptions.m_rawChainMeasureRunTimes;
  m_ov.m_rawChainDataOutputPeriod                  = mlOptions.m_rawChainDataOutputPeriod;
  m_ov.m_rawChainDataOutputFileName                = mlOptions.m_rawChainDataOutputFileName;
  m_ov.m_rawChainDataOutputFileType                = mlOptions.m_rawChainDataOutputFileType;
  m_ov.m_rawChainDataOutputAllowAll                = mlOptions.m_rawChainDataOutputAllowAll;
  m_ov.m_rawChainDataOutputAllowedSet              = mlOptions.m_rawChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_ov.m_rawChainComputeStats                      = mlOptions.m_rawChainComputeStats;
#endif
  m_ov.m_filteredChainGenerate                     = mlOptions.m_filteredChainGenerate;
  m_ov.m_filteredChainDiscardedPortion             = mlOptions.m_filteredChainDiscardedPortion;
  m_ov.m_filteredChainLag                          = mlOptions.m_filteredChainLag;
  m_ov.m_filteredChainDataOutputFileName           = mlOptions.m_filteredChainDataOutputFileName;
  m_ov.m_filteredChainDataOutputFileType           = mlOptions.m_filteredChainDataOutputFileType;
  m_ov.m_filteredChainDataOutputAllowAll           = mlOptions.m_filteredChainDataOutputAllowAll;
  m_ov.m_filteredChainDataOutputAllowedSet         = mlOptions.m_filteredChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_ov.m_filteredChainComputeStats                 = mlOptions.m_filteredChainComputeStats;
#endif
  m_ov.m_displayCandidates                         = mlOptions.m_displayCandidates;
  m_ov.m_putOutOfBoundsInChain                     = mlOptions.m_putOutOfBoundsInChain;
  m_ov.m_tkUseLocalHessian                         = mlOptions.m_tkUseLocalHessian;
  m_ov.m_tkUseNewtonComponent                      = mlOptions.m_tkUseNewtonComponent;
  m_ov.m_drMaxNumExtraStages                       = mlOptions.m_drMaxNumExtraStages;
  m_ov.m_drScalesForExtraStages                    = mlOptions.m_drScalesForExtraStages;
  m_ov.m_drDuringAmNonAdaptiveInt                  = mlOptions.m_drDuringAmNonAdaptiveInt;
  m_ov.m_amKeepInitialMatrix                       = mlOptions.m_amKeepInitialMatrix;
  m_ov.m_amInitialNonAdaptInterval                 = mlOptions.m_amInitialNonAdaptInterval;
  m_ov.m_amAdaptInterval                           = mlOptions.m_amAdaptInterval;
  m_ov.m_amAdaptedMatricesDataOutputPeriod         = mlOptions.m_amAdaptedMatricesDataOutputPeriod;
  m_ov.m_amAdaptedMatricesDataOutputFileName       = mlOptions.m_amAdaptedMatricesDataOutputFileName;
  m_ov.m_amAdaptedMatricesDataOutputFileType       = mlOptions.m_amAdaptedMatricesDataOutputFileType;
  m_ov.m_amAdaptedMatricesDataOutputAllowAll       = mlOptions.m_amAdaptedMatricesDataOutputAllowAll;
  m_ov.m_amAdaptedMatricesDataOutputAllowedSet     = mlOptions.m_amAdaptedMatricesDataOutputAllowedSet;
  m_ov.m_amEta                                     = mlOptions.m_amEta;
  m_ov.m_amEpsilon                                 = mlOptions.m_amEpsilon;
  m_ov.m_enableBrooksGelmanConvMonitor             = UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR;
  m_ov.m_BrooksGelmanLag                           = UQ_MH_SG_BROOKS_GELMAN_LAG;
  m_ov.m_outputLogLikelihood                       = UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD;
  m_ov.m_outputLogTarget                           = UQ_MH_SG_OUTPUT_LOG_TARGET;
  m_ov.m_doLogitTransform                          = mlOptions.m_doLogitTransform;
  m_ov.m_algorithm                                 = mlOptions.m_algorithm;
  m_ov.m_tk                                        = mlOptions.m_tk;
  m_ov.m_updateInterval                            = mlOptions.m_updateInterval;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
//m_ov.m_alternativeRawSsOptionsValues             = mlOptions.; // dakota
//m_ov.m_alternativeFilteredSsOptionsValues        = mlOptions.; // dakota
#endif

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChainStatisticalOptionsObj                  = mlOptions.m_rawChainStatisticalOptionsObj; // dakota
  m_rawChainStatOptsInstantiated                   = false;
  m_filteredChainStatisticalOptionsObj             = mlOptions.m_filteredChainStatisticalOptionsObj; // dakota
  m_filteredChainStatOptsInstantiated              = false;
#endif
  if ((m_env.subDisplayFile() != NULL ) &&
      (m_ov.m_totallyMute     == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSGOptions::constructor(3)"
                            << ": after copying values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}
// Destructor --------------------------------------
MetropolisHastingsSGOptions::~MetropolisHastingsSGOptions()
{
  queso_deprecated();

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_filteredChainStatOptsInstantiated) delete m_filteredChainStatisticalOptionsObj;
  if (m_rawChainStatOptsInstantiated     ) delete m_rawChainStatisticalOptionsObj;
#endif
}

// I/O methods -------------------------------------
void
MetropolisHastingsSGOptions::scanOptionsValues()
{
  queso_deprecated();

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  if ((m_env.subDisplayFile() != NULL) &&
      (m_ov.m_totallyMute == false   )) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSGOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_ov.m_rawChainComputeStats) {
    m_rawChainStatisticalOptionsObj = new SequenceStatisticalOptions(m_env,m_prefix + "rawChain_");
    m_rawChainStatOptsInstantiated  = true;
  }
  if (m_ov.m_filteredChainComputeStats) {
    m_filteredChainStatisticalOptionsObj = new SequenceStatisticalOptions(m_env,m_prefix + "filteredChain_");
    m_filteredChainStatOptsInstantiated  = true;
  }
#endif

  return;
}
// -------------------------------------------------
void
MetropolisHastingsSGOptions::print(std::ostream& os) const
{
  queso_deprecated();

  os <<         m_option_dataOutputFileName                         << " = " << m_ov.m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowAll                         << " = " << m_ov.m_dataOutputAllowAll
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_totallyMute                                << " = " << m_ov.m_totallyMute
     << "\n" << m_option_initialPosition_dataInputFileName          << " = " << m_ov.m_initialPositionDataInputFileName
     << "\n" << m_option_initialPosition_dataInputFileType          << " = " << m_ov.m_initialPositionDataInputFileType
     << "\n" << m_option_initialProposalCovMatrix_dataInputFileName << " = " << m_ov.m_initialProposalCovMatrixDataInputFileName
     << "\n" << m_option_initialProposalCovMatrix_dataInputFileType << " = " << m_ov.m_initialProposalCovMatrixDataInputFileType
     << "\n" << m_option_listOfDisabledParameters                   << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_parameterDisabledSet.begin(); setIt != m_ov.m_parameterDisabledSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_rawChain_dataInputFileName                 << " = " << m_ov.m_rawChainDataInputFileName
     << "\n" << m_option_rawChain_dataInputFileType                 << " = " << m_ov.m_rawChainDataInputFileType
     << "\n" << m_option_rawChain_size                              << " = " << m_ov.m_rawChainSize
     << "\n" << m_option_rawChain_generateExtra                     << " = " << m_ov.m_rawChainGenerateExtra
     << "\n" << m_option_rawChain_displayPeriod                     << " = " << m_ov.m_rawChainDisplayPeriod
     << "\n" << m_option_rawChain_measureRunTimes                   << " = " << m_ov.m_rawChainMeasureRunTimes
     << "\n" << m_option_rawChain_dataOutputPeriod                  << " = " << m_ov.m_rawChainDataOutputPeriod
     << "\n" << m_option_rawChain_dataOutputFileName                << " = " << m_ov.m_rawChainDataOutputFileName
     << "\n" << m_option_rawChain_dataOutputFileType                << " = " << m_ov.m_rawChainDataOutputFileType
     << "\n" << m_option_rawChain_dataOutputAllowAll                << " = " << m_ov.m_rawChainDataOutputAllowAll
     << "\n" << m_option_rawChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_rawChainDataOutputAllowedSet.begin(); setIt != m_ov.m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << m_option_rawChain_computeStats                      << " = " << m_ov.m_rawChainComputeStats
#endif
     << "\n" << m_option_filteredChain_generate                     << " = " << m_ov.m_filteredChainGenerate
     << "\n" << m_option_filteredChain_discardedPortion             << " = " << m_ov.m_filteredChainDiscardedPortion
     << "\n" << m_option_filteredChain_lag                          << " = " << m_ov.m_filteredChainLag
     << "\n" << m_option_filteredChain_dataOutputFileName           << " = " << m_ov.m_filteredChainDataOutputFileName
     << "\n" << m_option_filteredChain_dataOutputFileType           << " = " << m_ov.m_filteredChainDataOutputFileType
     << "\n" << m_option_filteredChain_dataOutputAllowAll           << " = " << m_ov.m_filteredChainDataOutputAllowAll
     << "\n" << m_option_filteredChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_filteredChainDataOutputAllowedSet.begin(); setIt != m_ov.m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << m_option_filteredChain_computeStats                 << " = " << m_ov.m_filteredChainComputeStats
#endif
     << "\n" << m_option_displayCandidates                          << " = " << m_ov.m_displayCandidates
     << "\n" << m_option_putOutOfBoundsInChain                      << " = " << m_ov.m_putOutOfBoundsInChain
     << "\n" << m_option_tk_useLocalHessian                         << " = " << m_ov.m_tkUseLocalHessian
     << "\n" << m_option_tk_useNewtonComponent                      << " = " << m_ov.m_tkUseNewtonComponent
     << "\n" << m_option_dr_maxNumExtraStages                       << " = " << m_ov.m_drMaxNumExtraStages
     << "\n" << m_option_dr_listOfScalesForExtraStages << " = ";
  for (unsigned int i = 0; i < m_ov.m_drScalesForExtraStages.size(); ++i) {
    os << m_ov.m_drScalesForExtraStages[i] << " ";
  }
  os << "\n" << m_option_dr_duringAmNonAdaptiveInt                  << " = " << m_ov.m_drDuringAmNonAdaptiveInt
     << "\n" << m_option_am_keepInitialMatrix                       << " = " << m_ov.m_amKeepInitialMatrix
     << "\n" << m_option_am_initialNonAdaptInterval                 << " = " << m_ov.m_amInitialNonAdaptInterval
     << "\n" << m_option_am_adaptInterval                           << " = " << m_ov.m_amAdaptInterval
     << "\n" << m_option_am_adaptedMatrices_dataOutputPeriod        << " = " << m_ov.m_amAdaptedMatricesDataOutputPeriod
     << "\n" << m_option_am_adaptedMatrices_dataOutputFileName      << " = " << m_ov.m_amAdaptedMatricesDataOutputFileName
     << "\n" << m_option_am_adaptedMatrices_dataOutputFileType      << " = " << m_ov.m_amAdaptedMatricesDataOutputFileType
     << "\n" << m_option_am_adaptedMatrices_dataOutputAllowAll      << " = " << m_ov.m_amAdaptedMatricesDataOutputAllowAll
     << "\n" << m_option_am_adaptedMatrices_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_amAdaptedMatricesDataOutputAllowedSet.begin(); setIt != m_ov.m_amAdaptedMatricesDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_am_eta                                     << " = " << m_ov.m_amEta
     << "\n" << m_option_am_epsilon                                 << " = " << m_ov.m_amEpsilon
     << "\n" << m_option_enableBrooksGelmanConvMonitor              << " = " << m_ov.m_enableBrooksGelmanConvMonitor
     << "\n" << m_option_BrooksGelmanLag                            << " = " << m_ov.m_BrooksGelmanLag
     << "\n" << m_option_outputLogLikelihood                        << " = " << m_ov.m_outputLogLikelihood
     << "\n" << m_option_outputLogTarget                            << " = " << m_ov.m_outputLogTarget
     << "\n" << m_option_doLogitTransform                           << " = " << m_ov.m_doLogitTransform
     << "\n" << m_option_algorithm                                  << " = " << m_ov.m_algorithm
     << "\n" << m_option_tk                                         << " = " << m_ov.m_tk
     << "\n" << m_option_updateInterval                             << " = " << m_ov.m_updateInterval
     << std::endl;

  return;
}

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
// Private methods----------------------------------
void
MetropolisHastingsSGOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                                                "produce help msg for Bayesian Metropolis-Hastings"          )
    (m_option_dataOutputFileName.c_str(),                         boost::program_options::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                           ), "name of generic output file"                                )
    (m_option_dataOutputAllowAll.c_str(),                         boost::program_options::value<bool        >()->default_value(UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV                           ), "allow all subEnvs write to a generic output file"           )
    (m_option_dataOutputAllowedSet.c_str(),                       boost::program_options::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV                         ), "subEnvs that will write to generic output file"             )
    (m_option_totallyMute.c_str(),                                boost::program_options::value<bool        >()->default_value(UQ_MH_SG_TOTALLY_MUTE_ODV                                    ), "totally mute (no printout msg)"                             )
    (m_option_initialPosition_dataInputFileName.c_str(),          boost::program_options::value<std::string >()->default_value(UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV           ), "name of input file for raw chain "                          )
    (m_option_initialPosition_dataInputFileType.c_str(),          boost::program_options::value<std::string >()->default_value(UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV           ), "type of input file for raw chain "                          )
    (m_option_initialProposalCovMatrix_dataInputFileName.c_str(), boost::program_options::value<std::string >()->default_value(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV), "name of input file for raw chain "                          )
    (m_option_initialProposalCovMatrix_dataInputFileType.c_str(), boost::program_options::value<std::string >()->default_value(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV), "type of input file for raw chain "                          )
    (m_option_listOfDisabledParameters.c_str(),                   boost::program_options::value<std::string >()->default_value(UQ_MH_SG_LIST_OF_DISABLED_PARAMETERS_ODV                     ), "list of disabled parameters"                                )
    (m_option_rawChain_dataInputFileName.c_str(),                 boost::program_options::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV                  ), "name of input file for raw chain "                          )
    (m_option_rawChain_dataInputFileType.c_str(),                 boost::program_options::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV                  ), "type of input file for raw chain "                          )
    (m_option_rawChain_size.c_str(),                              boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_SIZE_ODV                                  ), "size of raw chain"                                          )
    (m_option_rawChain_generateExtra.c_str(),                     boost::program_options::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV                        ), "generate extra information about raw chain"                 )
    (m_option_rawChain_displayPeriod.c_str(),                     boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV                        ), "period of msg display during raw chain generation"          )
    (m_option_rawChain_measureRunTimes.c_str(),                   boost::program_options::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV                     ), "measure run times"                                          )
    (m_option_rawChain_dataOutputPeriod.c_str(),                  boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV                    ), "period of msg display during raw chain generation"          )
    (m_option_rawChain_dataOutputFileName.c_str(),                boost::program_options::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV                 ), "name of output file for raw chain "                         )
    (m_option_rawChain_dataOutputFileType.c_str(),                boost::program_options::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV                 ), "type of output file for raw chain "                         )
    (m_option_rawChain_dataOutputAllowAll.c_str(),                boost::program_options::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV                 ), "allow all subEnvs to write raw chain to an output file"     )
    (m_option_rawChain_dataOutputAllowedSet.c_str(),              boost::program_options::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV               ), "subEnvs that will write raw chain to output file"           )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    (m_option_rawChain_computeStats.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV                         ), "compute statistics on raw chain"                            )
#endif
    (m_option_filteredChain_generate.c_str(),                     boost::program_options::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV                         ), "generate filtered chain"                                    )
    (m_option_filteredChain_discardedPortion.c_str(),             boost::program_options::value<double      >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV                ), "initial discarded portion for chain filtering"              )
    (m_option_filteredChain_lag.c_str(),                          boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                              ), "spacing for chain filtering"                                )
    (m_option_filteredChain_dataOutputFileName.c_str(),           boost::program_options::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV            ), "name of output file for filtered chain"                     )
    (m_option_filteredChain_dataOutputFileType.c_str(),           boost::program_options::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV            ), "type of output file for filtered chain"                     )
    (m_option_filteredChain_dataOutputAllowAll.c_str(),           boost::program_options::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV            ), "allow all subEnvs to write filt chain to an output file"    )
    (m_option_filteredChain_dataOutputAllowedSet.c_str(),         boost::program_options::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV          ), "subEnvs that will write filt chain to output file"          )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    (m_option_filteredChain_computeStats.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV                    ), "compute statistics on filtered chain"                       )
#endif
    (m_option_displayCandidates.c_str(),                          boost::program_options::value<bool        >()->default_value(UQ_MH_SG_DISPLAY_CANDIDATES_ODV                              ), "display candidates in the core MH algorithm"                )
    (m_option_putOutOfBoundsInChain.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV                      ), "put 'out of bound' candidates in chain as well"             )
    (m_option_tk_useLocalHessian.c_str(),                         boost::program_options::value<bool        >()->default_value(UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                            ), "'proposal' use local Hessian"                               )
    (m_option_tk_useNewtonComponent.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV                         ), "'proposal' use Newton component"                            )
    (m_option_dr_maxNumExtraStages.c_str(),                       boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV                         ), "'dr' maximum number of extra stages"                        )
    (m_option_dr_listOfScalesForExtraStages.c_str(),              boost::program_options::value<std::string >()->default_value(UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV              ), "'dr' scales for prop cov matrices from 2nd stage on"        )
    (m_option_dr_duringAmNonAdaptiveInt.c_str(),                  boost::program_options::value<bool        >()->default_value(UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV                   ), "'dr' used during 'am' non adaptive interval"                )
    (m_option_am_keepInitialMatrix.c_str(),                       boost::program_options::value<bool        >()->default_value(UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV                          ), "'am' keep initial (given) matrix"                           )
    (m_option_am_initialNonAdaptInterval.c_str(),                 boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                           ), "'am' initial non adaptation interval"                       )
    (m_option_am_adaptInterval.c_str(),                           boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                               ), "'am' adaptation interval"                                   )
    (m_option_am_adaptedMatrices_dataOutputPeriod.c_str(),        boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV          ), "period for outputting 'am' adapted matrices"                )
    (m_option_am_adaptedMatrices_dataOutputFileName.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV       ), "name of output file for 'am' adapted matrices"              )
    (m_option_am_adaptedMatrices_dataOutputFileType.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV       ), "type of output file for 'am' adapted matrices"              )
    (m_option_am_adaptedMatrices_dataOutputAllowAll.c_str(),      boost::program_options::value<bool        >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV       ), "type of output file for 'am' adapted matrices"              )
    (m_option_am_adaptedMatrices_dataOutputAllowedSet.c_str(),    boost::program_options::value<std::string >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV     ), "type of output file for 'am' adapted matrices"              )
    (m_option_am_eta.c_str(),                                     boost::program_options::value<double      >()->default_value(UQ_MH_SG_AM_ETA_ODV                                          ), "'am' eta"                                                   )
    (m_option_am_epsilon.c_str(),                                 boost::program_options::value<double      >()->default_value(UQ_MH_SG_AM_EPSILON_ODV                                      ), "'am' epsilon"                                               )
    (m_option_enableBrooksGelmanConvMonitor.c_str(),              boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR                   ), "assess convergence using Brooks-Gelman metric"              )
    (m_option_BrooksGelmanLag.c_str(),                            boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_BROOKS_GELMAN_LAG                                   ), "number of chain positions before starting to compute metric")
    (m_option_outputLogLikelihood.c_str(),                        boost::program_options::value<bool        >()->default_value(UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD                               ), "flag to toggle output of log likelihood values"             )
    (m_option_outputLogTarget.c_str(),                            boost::program_options::value<bool        >()->default_value(UQ_MH_SG_OUTPUT_LOG_TARGET                                   ), "flag to toggle output of log target values"                 )
    (m_option_doLogitTransform.c_str(),                           boost::program_options::value<bool        >()->default_value(UQ_MH_SG_DO_LOGIT_TRANSFORM                                  ), "flag to toggle logit transform for bounded domains"         )
    (m_option_algorithm.c_str(),                                  boost::program_options::value<std::string >()->default_value(UQ_MH_SG_ALGORITHM                                           ), "which mcmc algorithm to use"                                )
    (m_option_tk.c_str(),                                         boost::program_options::value<std::string >()->default_value(UQ_MH_SG_TK                                                  ), "which mcmc tk to use"                                       )
    (m_option_updateInterval.c_str(),                             boost::program_options::value<unsigned int>()->default_value(UQ_MH_SG_UPDATE_INTERVAL                                     ), "how often to call updateTK method"                          )
  ;

  return;
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
// -------------------------------------------------
void
MetropolisHastingsSGOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

  if (m_env.allOptionsMap().count(m_option_help)) {
    if ((m_env.subDisplayFile()) &&
        (m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowAll.c_str())) {
    m_ov.m_dataOutputAllowAll = m_env.allOptionsMap()[m_option_dataOutputAllowAll].as<bool>();
  }

  if (m_ov.m_dataOutputAllowAll) {
    m_ov.m_dataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_ov.m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_totallyMute)) {
    m_ov.m_totallyMute = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_totallyMute]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_initialPosition_dataInputFileName)) {
    m_ov.m_initialPositionDataInputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_initialPosition_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_initialPosition_dataInputFileType)) {
    m_ov.m_initialPositionDataInputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_initialPosition_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_initialProposalCovMatrix_dataInputFileName)) {
    m_ov.m_initialProposalCovMatrixDataInputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_initialProposalCovMatrix_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_initialProposalCovMatrix_dataInputFileType)) {
    m_ov.m_initialProposalCovMatrixDataInputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_initialProposalCovMatrix_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_listOfDisabledParameters)) {
    m_ov.m_parameterDisabledSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_listOfDisabledParameters].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);
    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_parameterDisabledSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileName)) {
    m_ov.m_rawChainDataInputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileType)) {
    m_ov.m_rawChainDataInputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_size)) {
    m_ov.m_rawChainSize = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_size]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_displayPeriod)) {
    m_ov.m_rawChainDisplayPeriod = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_displayPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_measureRunTimes)) {
    m_ov.m_rawChainMeasureRunTimes = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_measureRunTimes]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputPeriod)) {
    m_ov.m_rawChainDataOutputPeriod = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileName)) {
    m_ov.m_rawChainDataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileType)) {
    m_ov.m_rawChainDataOutputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputAllowAll.c_str())) {
    m_ov.m_rawChainDataOutputAllowAll = m_env.allOptionsMap()[m_option_rawChain_dataOutputAllowAll].as<bool>();
  }

  if (m_ov.m_rawChainDataOutputAllowAll) {
    m_ov.m_rawChainDataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputAllowedSet)) {
    m_ov.m_rawChainDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_rawChain_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_rawChainDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_env.allOptionsMap().count(m_option_rawChain_computeStats)) {
    m_ov.m_rawChainComputeStats = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_computeStats]).as<bool>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_rawChain_generateExtra)) {
    m_ov.m_rawChainGenerateExtra = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_rawChain_generateExtra]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_generate)) {
    m_ov.m_filteredChainGenerate = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_generate]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_discardedPortion)) {
    m_ov.m_filteredChainDiscardedPortion = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_discardedPortion]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_lag)) {
    m_ov.m_filteredChainLag = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_lag]).as<unsigned int>();
  }
  if ((m_ov.m_filteredChainGenerate == true) &&
      (m_ov.m_filteredChainLag      < 2    )) {
    std::cerr << "WARNING In MetropolisHastingsSG<P_V,P_M>::getMyOptionsValues()"
              << ", worldRank "             << m_env.worldRank()
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_filteredChain_lag
              << "' from "                  << m_ov.m_filteredChainLag
              << " to "                     << 2
              << std::endl;
    m_ov.m_filteredChainLag = 2;
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputFileName)) {
    m_ov.m_filteredChainDataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputFileType)) {
    m_ov.m_filteredChainDataOutputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputAllowAll.c_str())) {
    m_ov.m_filteredChainDataOutputAllowAll = m_env.allOptionsMap()[m_option_filteredChain_dataOutputAllowAll].as<bool>();
  }

  if (m_ov.m_filteredChainDataOutputAllowAll) {
    m_ov.m_filteredChainDataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputAllowedSet)) {
    m_ov.m_filteredChainDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_filteredChain_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_filteredChainDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_env.allOptionsMap().count(m_option_filteredChain_computeStats)) {
    m_ov.m_filteredChainComputeStats = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_computeStats]).as<bool>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_displayCandidates)) {
    m_ov.m_displayCandidates = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_displayCandidates]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_putOutOfBoundsInChain)) {
    m_ov.m_putOutOfBoundsInChain = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_putOutOfBoundsInChain]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useLocalHessian)) {
    m_ov.m_tkUseLocalHessian = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_tk_useLocalHessian]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useNewtonComponent)) {
    m_ov.m_tkUseNewtonComponent = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_tk_useNewtonComponent]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dr_maxNumExtraStages)) {
    m_ov.m_drMaxNumExtraStages = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dr_maxNumExtraStages]).as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count(m_option_dr_listOfScalesForExtraStages)) {
    std::string inputString = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dr_listOfScalesForExtraStages]).as<std::string>();
    MiscReadDoublesFromString(inputString,tmpScales);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::getMyOptionValues(): scales =";
    //  for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpScales[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}
  }

  if (m_ov.m_drMaxNumExtraStages > 0) {
    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();

    m_ov.m_drScalesForExtraStages.clear();
    m_ov.m_drScalesForExtraStages.resize(m_ov.m_drMaxNumExtraStages,1.);
    for (unsigned int i = 0; i < m_ov.m_drMaxNumExtraStages; ++i) {
      if (i < tmpSize) scale = tmpScales[i];
      m_ov.m_drScalesForExtraStages[i] = scale;
    }
    //updateTK();
  }

  if (m_env.allOptionsMap().count(m_option_dr_duringAmNonAdaptiveInt)) {
    m_ov.m_drDuringAmNonAdaptiveInt = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dr_duringAmNonAdaptiveInt]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_am_keepInitialMatrix)) {
    m_ov.m_amKeepInitialMatrix = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_keepInitialMatrix]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_am_initialNonAdaptInterval)) {
    m_ov.m_amInitialNonAdaptInterval = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_initialNonAdaptInterval]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptInterval)) {
    m_ov.m_amAdaptInterval = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_adaptInterval]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputPeriod)) {
    m_ov.m_amAdaptedMatricesDataOutputPeriod = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputFileName)) {
    m_ov.m_amAdaptedMatricesDataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputFileType)) {
    m_ov.m_amAdaptedMatricesDataOutputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputAllowAll.c_str())) {
    m_ov.m_amAdaptedMatricesDataOutputAllowAll = m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputAllowAll].as<bool>();
  }

  if (m_ov.m_amAdaptedMatricesDataOutputAllowAll) {
    m_ov.m_amAdaptedMatricesDataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputAllowedSet)) {
    m_ov.m_amAdaptedMatricesDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_amAdaptedMatricesDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_am_eta)) {
    m_ov.m_amEta = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_eta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_am_epsilon)) {
    m_ov.m_amEpsilon = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_am_epsilon]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_enableBrooksGelmanConvMonitor)) {
    m_ov.m_enableBrooksGelmanConvMonitor = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_enableBrooksGelmanConvMonitor]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_BrooksGelmanLag)) {
    m_ov.m_BrooksGelmanLag = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_BrooksGelmanLag]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_outputLogLikelihood)) {
    m_ov.m_outputLogLikelihood = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_outputLogLikelihood]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_outputLogTarget)) {
    m_ov.m_outputLogTarget = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_outputLogTarget]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_doLogitTransform)) {
    m_ov.m_doLogitTransform = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_doLogitTransform]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_algorithm)) {
    m_ov.m_algorithm = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_algorithm]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_tk)) {
    m_ov.m_tk = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_tk]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_updateInterval)) {
    m_ov.m_updateInterval = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_updateInterval]).as<unsigned int>();
  }
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

// --------------------------------------------------
// Operator declared outside class definition ------
// --------------------------------------------------

std::ostream& operator<<(std::ostream& os, const MetropolisHastingsSGOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
