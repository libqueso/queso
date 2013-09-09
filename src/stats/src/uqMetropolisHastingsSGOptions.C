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

#include <uqMetropolisHastingsSGOptions.h>
#include <uqMiscellaneous.h>

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
  m_dataOutputFileName                       (UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV),
  m_dataOutputAllowAll                       (UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV),
//m_dataOutputAllowedSet                     (),
  m_totallyMute                              (UQ_MH_SG_TOTALLY_MUTE_ODV),
  m_initialPositionDataInputFileName         (UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV),
  m_initialPositionDataInputFileType         (UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV),
  m_initialProposalCovMatrixDataInputFileName(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV),
  m_initialProposalCovMatrixDataInputFileType(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV),
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
  m_BrooksGelmanLag                          (UQ_MH_SG_BROOKS_GELMAN_LAG)
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  ,
  m_alternativeRawSsOptionsValues            (),
  m_alternativeFilteredSsOptionsValues       ()
#endif
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativeRawSsOptionsValues     ) m_alternativeRawSsOptionsValues      = *alternativeRawSsOptionsValues;
  if (alternativeFilteredSsOptionsValues) m_alternativeFilteredSsOptionsValues = *alternativeFilteredSsOptionsValues;
#endif
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

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_alternativeRawSsOptionsValues             = src.m_alternativeRawSsOptionsValues;
  m_alternativeFilteredSsOptionsValues        = src.m_alternativeFilteredSsOptionsValues;
#endif
  return;
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
  m_optionsDesc                                      (new po::options_description("Bayesian Metropolis-Hastings options")),
  m_option_help                                      (m_prefix + "help"                                       ),
  m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                         ),
  m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                         ),
  m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                       ),
  m_option_totallyMute                               (m_prefix + "totallyMute"                                ),
  m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"          ),
  m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"          ),
  m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName" ),
  m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType" ),
  m_option_rawChain_dataInputFileName                (m_prefix + "rawChain_dataInputFileName"                 ),
  m_option_rawChain_dataInputFileType                (m_prefix + "rawChain_dataInputFileType"                 ),
  m_option_rawChain_size                             (m_prefix + "rawChain_size"                              ),
  m_option_rawChain_generateExtra                    (m_prefix + "rawChain_generateExtra"                     ),
  m_option_rawChain_displayPeriod                    (m_prefix + "rawChain_displayPeriod"                     ),
  m_option_rawChain_measureRunTimes                  (m_prefix + "rawChain_measureRunTimes"                   ),
  m_option_rawChain_dataOutputPeriod                 (m_prefix + "rawChain_dataOutputPeriod"                  ),
  m_option_rawChain_dataOutputFileName               (m_prefix + "rawChain_dataOutputFileName"                ),
  m_option_rawChain_dataOutputFileType               (m_prefix + "rawChain_dataOutputFileType"                ),
  m_option_rawChain_dataOutputAllowAll               (m_prefix + "rawChain_dataOutputAllowAll"                ),
  m_option_rawChain_dataOutputAllowedSet             (m_prefix + "rawChain_dataOutputAllowedSet"              ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_rawChain_computeStats                     (m_prefix + "rawChain_computeStats"                      ),
#endif
  m_option_filteredChain_generate                    (m_prefix + "filteredChain_generate"                     ),
  m_option_filteredChain_discardedPortion            (m_prefix + "filteredChain_discardedPortion"             ),
  m_option_filteredChain_lag                         (m_prefix + "filteredChain_lag"                          ),
  m_option_filteredChain_dataOutputFileName          (m_prefix + "filteredChain_dataOutputFileName"           ),
  m_option_filteredChain_dataOutputFileType          (m_prefix + "filteredChain_dataOutputFileType"           ),
  m_option_filteredChain_dataOutputAllowAll          (m_prefix + "filteredChain_dataOutputAllowAll"           ),
  m_option_filteredChain_dataOutputAllowedSet        (m_prefix + "filteredChain_dataOutputAllowedSet"         ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_filteredChain_computeStats                (m_prefix + "filteredChain_computeStats"                 ),
#endif
  m_option_displayCandidates                         (m_prefix + "displayCandidates"                          ),
  m_option_putOutOfBoundsInChain                     (m_prefix + "putOutOfBoundsInChain"                      ),
  m_option_tk_useLocalHessian                        (m_prefix + "tk_useLocalHessian"                         ),
  m_option_tk_useNewtonComponent                     (m_prefix + "tk_useNewtonComponent"                      ),
  m_option_dr_maxNumExtraStages                      (m_prefix + "dr_maxNumExtraStages"                       ),
  m_option_dr_listOfScalesForExtraStages             (m_prefix + "dr_listOfScalesForExtraStages"              ),
  m_option_dr_duringAmNonAdaptiveInt                 (m_prefix + "dr_duringAmNonAdaptiveInt"                  ),
  m_option_am_keepInitialMatrix                      (m_prefix + "am_keepInitialMatrix"                       ),
  m_option_am_initialNonAdaptInterval                (m_prefix + "am_initialNonAdaptInterval"                 ),
  m_option_am_adaptInterval                          (m_prefix + "am_adaptInterval"                           ),
  m_option_am_adaptedMatrices_dataOutputPeriod       (m_prefix + "am_adaptedMatrices_dataOutputPeriod"        ),
  m_option_am_adaptedMatrices_dataOutputFileName     (m_prefix + "am_adaptedMatrices_dataOutputFileName"      ),
  m_option_am_adaptedMatrices_dataOutputFileType     (m_prefix + "am_adaptedMatrices_dataOutputFileType"      ),
  m_option_am_adaptedMatrices_dataOutputAllowAll     (m_prefix + "am_adaptedMatrices_dataOutputAllowAll"      ),
  m_option_am_adaptedMatrices_dataOutputAllowedSet   (m_prefix + "am_adaptedMatrices_dataOutputAllowedSet"    ),
  m_option_am_eta                                    (m_prefix + "am_eta"                                     ),
  m_option_am_epsilon                                (m_prefix + "am_epsilon"                                 ),
  m_option_enableBrooksGelmanConvMonitor             (m_prefix + "enableBrooksGelmanConvMonitor"              ),
  m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                            )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "MetropolisHastingsSGOptions::constructor(1)",
                      "this constructor is incompatible with the absence of an options input file");
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
  m_optionsDesc                                      (NULL),
  m_option_help                                      (m_prefix + "help"                                      ),
  m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
  m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
  m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
  m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
  m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
  m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
  m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
  m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
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
  m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "MetropolisHastingsSGOptions::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

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
  m_optionsDesc                                      (NULL),
  m_option_help                                      (m_prefix + "help"                                      ),
  m_option_dataOutputFileName                        (m_prefix + "dataOutputFileName"                        ),
  m_option_dataOutputAllowAll                        (m_prefix + "dataOutputAllowAll"                        ),
  m_option_dataOutputAllowedSet                      (m_prefix + "dataOutputAllowedSet"                      ),
  m_option_totallyMute                               (m_prefix + "totallyMute"                               ),
  m_option_initialPosition_dataInputFileName         (m_prefix + "initialPosition_dataInputFileName"         ),
  m_option_initialPosition_dataInputFileType         (m_prefix + "initialPosition_dataInputFileType"         ),
  m_option_initialProposalCovMatrix_dataInputFileName(m_prefix + "initialProposalCovMatrix_dataInputFileName"),
  m_option_initialProposalCovMatrix_dataInputFileType(m_prefix + "initialProposalCovMatrix_dataInputFileType"),
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
  m_option_BrooksGelmanLag                           (m_prefix + "BrooksGelmanLag"                           )
{
  m_ov.m_dataOutputFileName                        = mlOptions.m_dataOutputFileName;
  m_ov.m_dataOutputAllowAll                        = mlOptions.m_dataOutputAllowAll;
  m_ov.m_dataOutputAllowedSet                      = mlOptions.m_dataOutputAllowedSet;
  m_ov.m_totallyMute                               = mlOptions.m_totallyMute;
  m_ov.m_initialPositionDataInputFileName          = mlOptions.m_initialPositionDataInputFileName;
  m_ov.m_initialPositionDataInputFileType          = mlOptions.m_initialPositionDataInputFileType;
  m_ov.m_initialProposalCovMatrixDataInputFileName = mlOptions.m_initialProposalCovMatrixDataInputFileName;
  m_ov.m_initialProposalCovMatrixDataInputFileType = mlOptions.m_initialProposalCovMatrixDataInputFileType;
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
  if ((m_env.subDisplayFile()        != NULL ) &&
      (m_ov.m_totallyMute == false)) {
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
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_filteredChainStatOptsInstantiated) delete m_filteredChainStatisticalOptionsObj;
  if (m_rawChainStatOptsInstantiated     ) delete m_rawChainStatisticalOptionsObj;
#endif
  if (m_optionsDesc                      ) delete m_optionsDesc;
} 

// I/O methods -------------------------------------
void
MetropolisHastingsSGOptions::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "MetropolisHastingsSGOptions::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

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
     << "\n" << m_option_rawChain_dataInputFileName                 << " = " << m_ov.m_rawChainDataInputFileName
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
     << std::endl;

  return;
}

// Private methods----------------------------------
void
MetropolisHastingsSGOptions::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                                                                "produce help msg for Bayesian Metropolis-Hastings"          )
    (m_option_dataOutputFileName.c_str(),                         po::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                           ), "name of generic output file"                                )
    (m_option_dataOutputAllowAll.c_str(),                         po::value<bool        >()->default_value(UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV                           ), "allow all subEnvs write to a generic output file"           )
    (m_option_dataOutputAllowedSet.c_str(),                       po::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV                         ), "subEnvs that will write to generic output file"             )
    (m_option_totallyMute.c_str(),                                po::value<bool        >()->default_value(UQ_MH_SG_TOTALLY_MUTE_ODV                                    ), "totally mute (no printout msg)"                             )
    (m_option_initialPosition_dataInputFileName.c_str(),          po::value<std::string >()->default_value(UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV           ), "name of input file for raw chain "                          )
    (m_option_initialPosition_dataInputFileType.c_str(),          po::value<std::string >()->default_value(UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV           ), "type of input file for raw chain "                          )
    (m_option_initialProposalCovMatrix_dataInputFileName.c_str(), po::value<std::string >()->default_value(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV), "name of input file for raw chain "                          )
    (m_option_initialProposalCovMatrix_dataInputFileType.c_str(), po::value<std::string >()->default_value(UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV), "type of input file for raw chain "                          )
    (m_option_rawChain_dataInputFileName.c_str(),                 po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV                  ), "name of input file for raw chain "                          )
    (m_option_rawChain_dataInputFileType.c_str(),                 po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV                  ), "type of input file for raw chain "                          )
    (m_option_rawChain_size.c_str(),                              po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_SIZE_ODV                                  ), "size of raw chain"                                          )
    (m_option_rawChain_generateExtra.c_str(),                     po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV                        ), "generate extra information about raw chain"                 )
    (m_option_rawChain_displayPeriod.c_str(),                     po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV                        ), "period of msg display during raw chain generation"          )
    (m_option_rawChain_measureRunTimes.c_str(),                   po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV                     ), "measure run times"                                          )
    (m_option_rawChain_dataOutputPeriod.c_str(),                  po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV                    ), "period of msg display during raw chain generation"          )
    (m_option_rawChain_dataOutputFileName.c_str(),                po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV                 ), "name of output file for raw chain "                         )
    (m_option_rawChain_dataOutputFileType.c_str(),                po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV                 ), "type of output file for raw chain "                         )
    (m_option_rawChain_dataOutputAllowAll.c_str(),                po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV                 ), "allow all subEnvs to write raw chain to an output file"     )
    (m_option_rawChain_dataOutputAllowedSet.c_str(),              po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV               ), "subEnvs that will write raw chain to output file"           )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    (m_option_rawChain_computeStats.c_str(),                      po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV                         ), "compute statistics on raw chain"                            )
#endif
    (m_option_filteredChain_generate.c_str(),                     po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV                         ), "generate filtered chain"                                    )
    (m_option_filteredChain_discardedPortion.c_str(),             po::value<double      >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV                ), "initial discarded portion for chain filtering"              )
    (m_option_filteredChain_lag.c_str(),                          po::value<unsigned int>()->default_value(UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                              ), "spacing for chain filtering"                                )
    (m_option_filteredChain_dataOutputFileName.c_str(),           po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV            ), "name of output file for filtered chain"                     )
    (m_option_filteredChain_dataOutputFileType.c_str(),           po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV            ), "type of output file for filtered chain"                     )
    (m_option_filteredChain_dataOutputAllowAll.c_str(),           po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV            ), "allow all subEnvs to write filt chain to an output file"    )
    (m_option_filteredChain_dataOutputAllowedSet.c_str(),         po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV          ), "subEnvs that will write filt chain to output file"          )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    (m_option_filteredChain_computeStats.c_str(),                 po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV                    ), "compute statistics on filtered chain"                       )
#endif
    (m_option_displayCandidates.c_str(),                          po::value<bool        >()->default_value(UQ_MH_SG_DISPLAY_CANDIDATES_ODV                              ), "display candidates in the core MH algorithm"                )
    (m_option_putOutOfBoundsInChain.c_str(),                      po::value<bool        >()->default_value(UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV                      ), "put 'out of bound' candidates in chain as well"             )
    (m_option_tk_useLocalHessian.c_str(),                         po::value<bool        >()->default_value(UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                            ), "'proposal' use local Hessian"                               )
    (m_option_tk_useNewtonComponent.c_str(),                      po::value<bool        >()->default_value(UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV                         ), "'proposal' use Newton component"                            )
    (m_option_dr_maxNumExtraStages.c_str(),                       po::value<unsigned int>()->default_value(UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV                         ), "'dr' maximum number of extra stages"                        )
    (m_option_dr_listOfScalesForExtraStages.c_str(),              po::value<std::string >()->default_value(UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV              ), "'dr' scales for prop cov matrices from 2nd stage on"        )
    (m_option_dr_duringAmNonAdaptiveInt.c_str(),                  po::value<bool        >()->default_value(UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV                   ), "'dr' used during 'am' non adaptive interval"                )
    (m_option_am_keepInitialMatrix.c_str(),                       po::value<bool        >()->default_value(UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV                          ), "'am' keep initial (given) matrix"                           )
    (m_option_am_initialNonAdaptInterval.c_str(),                 po::value<unsigned int>()->default_value(UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                           ), "'am' initial non adaptation interval"                       )
    (m_option_am_adaptInterval.c_str(),                           po::value<unsigned int>()->default_value(UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                               ), "'am' adaptation interval"                                   )
    (m_option_am_adaptedMatrices_dataOutputPeriod.c_str(),        po::value<unsigned int>()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV          ), "period for outputting 'am' adapted matrices"                 )
    (m_option_am_adaptedMatrices_dataOutputFileName.c_str(),      po::value<std::string >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV       ), "name of output file for 'am' adapted matrices"              )
    (m_option_am_adaptedMatrices_dataOutputFileType.c_str(),      po::value<std::string >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV       ), "type of output file for 'am' adapted matrices"              )
    (m_option_am_adaptedMatrices_dataOutputAllowAll.c_str(),      po::value<bool        >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV       ), "type of output file for 'am' adapted matrices"              )
    (m_option_am_adaptedMatrices_dataOutputAllowedSet.c_str(),    po::value<std::string >()->default_value(UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV     ), "type of output file for 'am' adapted matrices"              )
    (m_option_am_eta.c_str(),                                     po::value<double      >()->default_value(UQ_MH_SG_AM_ETA_ODV                                          ), "'am' eta"                                                   )
    (m_option_am_epsilon.c_str(),                                 po::value<double      >()->default_value(UQ_MH_SG_AM_EPSILON_ODV                                      ), "'am' epsilon"                                               )
    (m_option_enableBrooksGelmanConvMonitor.c_str(),              po::value<unsigned int>()->default_value(UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR                   ), "assess convergence using Brooks-Gelman metric"              )
    (m_option_BrooksGelmanLag.c_str(),                            po::value<unsigned int>()->default_value(UQ_MH_SG_BROOKS_GELMAN_LAG                                   ), "number of chain positions before starting to compute metric")
  ;

  return;
}
// -------------------------------------------------
void
MetropolisHastingsSGOptions::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if ((m_env.subDisplayFile()) &&
        (m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
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
    m_ov.m_totallyMute = ((const po::variable_value&) m_env.allOptionsMap()[m_option_totallyMute]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_initialPosition_dataInputFileName)) {
    m_ov.m_initialPositionDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_initialPosition_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_initialPosition_dataInputFileType)) {
    m_ov.m_initialPositionDataInputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_initialPosition_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_initialProposalCovMatrix_dataInputFileName)) {
    m_ov.m_initialProposalCovMatrixDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_initialProposalCovMatrix_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_initialProposalCovMatrix_dataInputFileType)) {
    m_ov.m_initialProposalCovMatrixDataInputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_initialProposalCovMatrix_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileName)) {
    m_ov.m_rawChainDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileType)) {
    m_ov.m_rawChainDataInputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_size)) {
    m_ov.m_rawChainSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_size]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_displayPeriod)) {
    m_ov.m_rawChainDisplayPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_displayPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_measureRunTimes)) {
    m_ov.m_rawChainMeasureRunTimes = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_measureRunTimes]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputPeriod)) {
    m_ov.m_rawChainDataOutputPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileName)) {
    m_ov.m_rawChainDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileType)) {
    m_ov.m_rawChainDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileType]).as<std::string>();
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
    m_ov.m_rawChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_computeStats]).as<bool>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_rawChain_generateExtra)) {
    m_ov.m_rawChainGenerateExtra = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_generateExtra]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_generate)) {
    m_ov.m_filteredChainGenerate = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_generate]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_discardedPortion)) {
    m_ov.m_filteredChainDiscardedPortion = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_discardedPortion]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_lag)) {
    m_ov.m_filteredChainLag = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_lag]).as<unsigned int>();
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
    m_ov.m_filteredChainDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputFileType)) {
    m_ov.m_filteredChainDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileType]).as<std::string>();
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
    m_ov.m_filteredChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_computeStats]).as<bool>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_displayCandidates)) {
    m_ov.m_displayCandidates = ((const po::variable_value&) m_env.allOptionsMap()[m_option_displayCandidates]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_putOutOfBoundsInChain)) {
    m_ov.m_putOutOfBoundsInChain = ((const po::variable_value&) m_env.allOptionsMap()[m_option_putOutOfBoundsInChain]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useLocalHessian)) {
    m_ov.m_tkUseLocalHessian = ((const po::variable_value&) m_env.allOptionsMap()[m_option_tk_useLocalHessian]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useNewtonComponent)) {
    m_ov.m_tkUseNewtonComponent = ((const po::variable_value&) m_env.allOptionsMap()[m_option_tk_useNewtonComponent]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dr_maxNumExtraStages)) {
    m_ov.m_drMaxNumExtraStages = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_maxNumExtraStages]).as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count(m_option_dr_listOfScalesForExtraStages)) {
    std::string inputString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_listOfScalesForExtraStages]).as<std::string>();
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
    m_ov.m_drDuringAmNonAdaptiveInt = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_duringAmNonAdaptiveInt]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_am_keepInitialMatrix)) {
    m_ov.m_amKeepInitialMatrix = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_keepInitialMatrix]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_am_initialNonAdaptInterval)) {
    m_ov.m_amInitialNonAdaptInterval = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_initialNonAdaptInterval]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptInterval)) {
    m_ov.m_amAdaptInterval = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_adaptInterval]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputPeriod)) {
    m_ov.m_amAdaptedMatricesDataOutputPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputFileName)) {
    m_ov.m_amAdaptedMatricesDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptedMatrices_dataOutputFileType)) {
    m_ov.m_amAdaptedMatricesDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_adaptedMatrices_dataOutputFileType]).as<std::string>();
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
    m_ov.m_amEta = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_eta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_am_epsilon)) {
    m_ov.m_amEpsilon = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_epsilon]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_enableBrooksGelmanConvMonitor)) {
    m_ov.m_enableBrooksGelmanConvMonitor = ((const po::variable_value&) m_env.allOptionsMap()[m_option_enableBrooksGelmanConvMonitor]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_BrooksGelmanLag)) {
    m_ov.m_BrooksGelmanLag = ((const po::variable_value&) m_env.allOptionsMap()[m_option_BrooksGelmanLag]).as<unsigned int>();
  }

  return;
}

// --------------------------------------------------
// Operator declared outside class definition ------
// --------------------------------------------------

std::ostream& operator<<(std::ostream& os, const MetropolisHastingsSGOptions& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO
