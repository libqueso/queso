/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <uqMetropolisHastingsSGOptions.h>
#include <uqMiscellaneous.h>

uqMhOptionsValuesClass::uqMhOptionsValuesClass(
  const uqSsOptionsValuesClass* rawSsOptionsValues,
  const uqSsOptionsValuesClass* filteredSsOptionsValues)
  :
  m_dataOutputFileName                (UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV),
//m_dataOutputAllowedSet              (),
  m_totallyMute                       (UQ_MH_SG_TOTALLY_MUTE_ODV),
  m_rawChainDataInputFileName         (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV),
  m_rawChainDataInputFileType         (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV),
  m_rawChainSize                      (UQ_MH_SG_RAW_CHAIN_SIZE_ODV),
  m_rawChainGenerateExtra             (UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV),
  m_rawChainDisplayPeriod             (UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV),
  m_rawChainMeasureRunTimes           (UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV),
  m_rawChainDataOutputFileName        (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
  m_rawChainDataOutputFileType        (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
//m_rawChainDataOutputAllowedSet      (),
  m_rawChainComputeStats              (UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV),
  m_rawChainStatisticalOptionsValues  (),
  m_filteredChainGenerate             (UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV),
  m_filteredChainDiscardedPortion     (UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
  m_filteredChainLag                  (UQ_MH_SG_FILTERED_CHAIN_LAG_ODV),
  m_filteredChainDataOutputFileName   (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
  m_filteredChainDataOutputFileType   (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV),
//m_filteredChainDataOutputAllowedSet (),
  m_filteredChainComputeStats         (UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV),
  m_filteredChainStatisticalOptionsValues(),
  m_displayCandidates                 (UQ_MH_SG_DISPLAY_CANDIDATES_ODV),
  m_putOutOfBoundsInChain             (UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV),
  m_tkUseLocalHessian                 (UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV),
  m_tkUseNewtonComponent              (UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV),
  m_drMaxNumExtraStages               (UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_drScalesForExtraStages            (0),
  m_drDuringAmNonAdaptiveInt          (UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV),
  m_amKeepInitialMatrix               (UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV),
  m_amInitialNonAdaptInterval         (UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV),
  m_amAdaptInterval                   (UQ_MH_SG_AM_ADAPT_INTERVAL_ODV),
  m_amEta                             (UQ_MH_SG_AM_ETA_ODV),
  m_amEpsilon                         (UQ_MH_SG_AM_EPSILON_ODV),
  m_enableBrooksGelmanConvMonitor     (UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR),
  m_BrooksGelmanLag                   (UQ_MH_SG_BROOKS_GELMAN_LAG)
{
  if (rawSsOptionsValues     ) m_rawChainStatisticalOptionsValues      = *rawSsOptionsValues;
  if (filteredSsOptionsValues) m_filteredChainStatisticalOptionsValues = *filteredSsOptionsValues;
}

uqMhOptionsValuesClass::~uqMhOptionsValuesClass()
{
}

uqMhOptionsValuesClass::uqMhOptionsValuesClass(const uqMhOptionsValuesClass& src)
{
  this->copy(src);
}

uqMhOptionsValuesClass&
uqMhOptionsValuesClass::operator=(const uqMhOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}

void
uqMhOptionsValuesClass::copy(const uqMhOptionsValuesClass& src)
{
  m_dataOutputFileName                 = src.m_dataOutputFileName;
  m_dataOutputAllowedSet               = src.m_dataOutputAllowedSet;
  m_totallyMute                        = src.m_totallyMute;
  m_rawChainDataInputFileName          = src.m_rawChainDataInputFileName;
  m_rawChainDataInputFileType          = src.m_rawChainDataInputFileType;
  m_rawChainSize                       = src.m_rawChainSize;
  m_rawChainGenerateExtra              = src.m_rawChainGenerateExtra;
  m_rawChainDisplayPeriod              = src.m_rawChainDisplayPeriod;
  m_rawChainMeasureRunTimes            = src.m_rawChainMeasureRunTimes;
  m_rawChainDataOutputFileName         = src.m_rawChainDataOutputFileName;
  m_rawChainDataOutputFileType         = src.m_rawChainDataOutputFileType;
  m_rawChainDataOutputAllowedSet       = src.m_rawChainDataOutputAllowedSet;
  m_rawChainComputeStats               = src.m_rawChainComputeStats;
  m_rawChainStatisticalOptionsValues   = src.m_rawChainStatisticalOptionsValues;
  //m_rawChainStatisticalOptionsObj      = src.m_rawChainStatisticalOptionsObj; // dakota
  //m_rawChainStatOptsInstantiated       = src.m_rawChainStatOptsInstantiated; // dakota
  m_filteredChainGenerate              = src.m_filteredChainGenerate;
  m_filteredChainDiscardedPortion      = src.m_filteredChainDiscardedPortion;
  m_filteredChainLag                   = src.m_filteredChainLag;
  m_filteredChainDataOutputFileName    = src.m_filteredChainDataOutputFileName;
  m_filteredChainDataOutputFileType    = src.m_filteredChainDataOutputFileType;
  m_filteredChainDataOutputAllowedSet  = src.m_filteredChainDataOutputAllowedSet;
  m_filteredChainComputeStats          = src.m_filteredChainComputeStats;
  m_filteredChainStatisticalOptionsValues = src.m_rawChainStatisticalOptionsValues;
  //m_filteredChainStatisticalOptionsObj = src.m_filteredChainStatisticalOptionsObj; // dakota
  //m_filteredChainStatOptsInstantiated  = src.m_filteredChainStatOptsInstantiated; // dakota
  m_displayCandidates                  = src.m_displayCandidates;
  m_putOutOfBoundsInChain              = src.m_putOutOfBoundsInChain;
  m_tkUseLocalHessian                  = src.m_tkUseLocalHessian;
  m_tkUseNewtonComponent               = src.m_tkUseNewtonComponent;
  m_drMaxNumExtraStages                = src.m_drMaxNumExtraStages;
  m_drScalesForExtraStages             = src.m_drScalesForExtraStages;
  m_drDuringAmNonAdaptiveInt           = src.m_drDuringAmNonAdaptiveInt;
  m_amKeepInitialMatrix                = src.m_amKeepInitialMatrix;
  m_amInitialNonAdaptInterval          = src.m_amInitialNonAdaptInterval;
  m_amAdaptInterval                    = src.m_amAdaptInterval;
  m_amEta                              = src.m_amEta;
  m_amEpsilon                          = src.m_amEpsilon;
  m_enableBrooksGelmanConvMonitor      = src.m_enableBrooksGelmanConvMonitor;
  m_BrooksGelmanLag                    = src.m_BrooksGelmanLag;

  return;
}

uqMetropolisHastingsSGOptionsClass::uqMetropolisHastingsSGOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_optionsValues                            (NULL,NULL), // dakota
  m_rawChainStatisticalOptionsObj            (NULL),
  m_rawChainStatOptsInstantiated             (false),
  m_filteredChainStatisticalOptionsObj       (NULL),
  m_filteredChainStatOptsInstantiated        (false),
  m_prefix                                   ((std::string)(prefix) + "mh_"),
  m_env                                      (env),
  m_optionsDesc                              (new po::options_description("Bayesian Metropolis-Hastings options")),
  m_option_help                              (m_prefix + "help"                              ),
  m_option_dataOutputFileName                (m_prefix + "dataOutputFileName"                ),
  m_option_dataOutputAllowedSet              (m_prefix + "dataOutputAllowedSet"              ),
  m_option_totallyMute                       (m_prefix + "totallyMute"                       ),
  m_option_rawChain_dataInputFileName        (m_prefix + "rawChain_dataInputFileName"        ),
  m_option_rawChain_dataInputFileType        (m_prefix + "rawChain_dataInputFileType"        ),
  m_option_rawChain_size                     (m_prefix + "rawChain_size"                     ),
  m_option_rawChain_generateExtra            (m_prefix + "rawChain_generateExtra"            ),
  m_option_rawChain_displayPeriod            (m_prefix + "rawChain_displayPeriod"            ),
  m_option_rawChain_measureRunTimes          (m_prefix + "rawChain_measureRunTimes"          ),
  m_option_rawChain_dataOutputFileName       (m_prefix + "rawChain_dataOutputFileName"       ),
  m_option_rawChain_dataOutputFileType       (m_prefix + "rawChain_dataOutputFileType"       ),
  m_option_rawChain_dataOutputAllowedSet     (m_prefix + "rawChain_dataOutputAllowedSet"     ),
  m_option_rawChain_computeStats             (m_prefix + "rawChain_computeStats"             ),
  m_option_filteredChain_generate            (m_prefix + "filteredChain_generate"            ),
  m_option_filteredChain_discardedPortion    (m_prefix + "filteredChain_discardedPortion"    ),
  m_option_filteredChain_lag                 (m_prefix + "filteredChain_lag"                 ),
  m_option_filteredChain_dataOutputFileName  (m_prefix + "filteredChain_dataOutputFileName"  ),
  m_option_filteredChain_dataOutputFileType  (m_prefix + "filteredChain_dataOutputFileType"  ),
  m_option_filteredChain_dataOutputAllowedSet(m_prefix + "filteredChain_dataOutputAllowedSet"),
  m_option_filteredChain_computeStats        (m_prefix + "filteredChain_computeStats"        ),
  m_option_displayCandidates                 (m_prefix + "displayCandidates"                 ),
  m_option_putOutOfBoundsInChain             (m_prefix + "putOutOfBoundsInChain"             ),
  m_option_tk_useLocalHessian                (m_prefix + "tk_useLocalHessian"                ),
  m_option_tk_useNewtonComponent             (m_prefix + "tk_useNewtonComponent"             ),
  m_option_dr_maxNumExtraStages              (m_prefix + "dr_maxNumExtraStages"              ),
  m_option_dr_listOfScalesForExtraStages     (m_prefix + "dr_listOfScalesForExtraStages"     ),
  m_option_dr_duringAmNonAdaptiveInt         (m_prefix + "dr_duringAmNonAdaptiveInt"         ),
  m_option_am_keepInitialMatrix              (m_prefix + "am_keepInitialMatrix"              ),
  m_option_am_initialNonAdaptInterval        (m_prefix + "am_initialNonAdaptInterval"        ),
  m_option_am_adaptInterval                  (m_prefix + "am_adaptInterval"                  ),
  m_option_am_eta                            (m_prefix + "am_eta"                            ),
  m_option_am_epsilon                        (m_prefix + "am_epsilon"                        ),
  m_option_enableBrooksGelmanConvMonitor     (m_prefix + "enableBrooksGelmanConvMonitor"     ),
  m_option_BrooksGelmanLag                   (m_prefix + "BrooksGelmanLag"                   )
{
}

uqMetropolisHastingsSGOptionsClass::uqMetropolisHastingsSGOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
  const uqMhOptionsValuesClass& optionsValues)
  :
  m_optionsValues                            (optionsValues),
  m_rawChainStatisticalOptionsObj            (NULL),
  m_rawChainStatOptsInstantiated             (false),
  m_filteredChainStatisticalOptionsObj       (NULL),
  m_filteredChainStatOptsInstantiated        (false),
  m_prefix                                   ((std::string)(prefix) + "mh_"),
  m_env                                      (env),
  m_optionsDesc                              (NULL),
  m_option_help                              (m_prefix + "help"                              ),
  m_option_dataOutputFileName                (m_prefix + "dataOutputFileName"                ),
  m_option_dataOutputAllowedSet              (m_prefix + "dataOutputAllowedSet"              ),
  m_option_totallyMute                       (m_prefix + "totallyMute"                       ),
  m_option_rawChain_dataInputFileName        (m_prefix + "rawChain_dataInputFileName"        ),
  m_option_rawChain_dataInputFileType        (m_prefix + "rawChain_dataInputFileType"        ),
  m_option_rawChain_size                     (m_prefix + "rawChain_size"                     ),
  m_option_rawChain_generateExtra            (m_prefix + "rawChain_generateExtra"            ),
  m_option_rawChain_displayPeriod            (m_prefix + "rawChain_displayPeriod"            ),
  m_option_rawChain_measureRunTimes          (m_prefix + "rawChain_measureRunTimes"          ),
  m_option_rawChain_dataOutputFileName       (m_prefix + "rawChain_dataOutputFileName"       ),
  m_option_rawChain_dataOutputFileType       (m_prefix + "rawChain_dataOutputFileType"       ),
  m_option_rawChain_dataOutputAllowedSet     (m_prefix + "rawChain_dataOutputAllowedSet"     ),
  m_option_rawChain_computeStats             (m_prefix + "rawChain_computeStats"             ),
  m_option_filteredChain_generate            (m_prefix + "filteredChain_generate"            ),
  m_option_filteredChain_discardedPortion    (m_prefix + "filteredChain_discardedPortion"    ),
  m_option_filteredChain_lag                 (m_prefix + "filteredChain_lag"                 ),
  m_option_filteredChain_dataOutputFileName  (m_prefix + "filteredChain_dataOutputFileName"  ),
  m_option_filteredChain_dataOutputFileType  (m_prefix + "filteredChain_dataOutputFileType"  ),
  m_option_filteredChain_dataOutputAllowedSet(m_prefix + "filteredChain_dataOutputAllowedSet"),
  m_option_filteredChain_computeStats        (m_prefix + "filteredChain_computeStats"        ),
  m_option_displayCandidates                 (m_prefix + "displayCandidates"                 ),
  m_option_putOutOfBoundsInChain             (m_prefix + "putOutOfBoundsInChain"             ),
  m_option_tk_useLocalHessian                (m_prefix + "tk_useLocalHessian"                ),
  m_option_tk_useNewtonComponent             (m_prefix + "tk_useNewtonComponent"             ),
  m_option_dr_maxNumExtraStages              (m_prefix + "dr_maxNumExtraStages"              ),
  m_option_dr_listOfScalesForExtraStages     (m_prefix + "dr_listOfScalesForExtraStages"     ),
  m_option_dr_duringAmNonAdaptiveInt         (m_prefix + "dr_duringAmNonAdaptiveInt"         ),
  m_option_am_keepInitialMatrix              (m_prefix + "am_keepInitialMatrix"              ),
  m_option_am_initialNonAdaptInterval        (m_prefix + "am_initialNonAdaptInterval"        ),
  m_option_am_adaptInterval                  (m_prefix + "am_adaptInterval"                  ),
  m_option_am_eta                            (m_prefix + "am_eta"                            ),
  m_option_am_epsilon                        (m_prefix + "am_epsilon"                        ),
  m_option_enableBrooksGelmanConvMonitor     (m_prefix + "enableBrooksGelmanConvMonitor"     ),
  m_option_BrooksGelmanLag                   (m_prefix + "BrooksGelmanLag"                   )
{
  // dakota
}

uqMetropolisHastingsSGOptionsClass::uqMetropolisHastingsSGOptionsClass(
  const uqMLSamplingLevelOptionsClass& inputOptions)
  :
  m_optionsValues                            (NULL,NULL), // dakota
  m_rawChainStatisticalOptionsObj            (NULL),
  m_rawChainStatOptsInstantiated             (false),
  m_filteredChainStatisticalOptionsObj       (NULL),
  m_filteredChainStatOptsInstantiated        (false),
  m_prefix                                   (inputOptions.m_prefix),
  m_env                                      (inputOptions.env()),
  m_optionsDesc                              (NULL),
  m_option_help                              (m_prefix + "help"                              ),
  m_option_dataOutputFileName                (m_prefix + "dataOutputFileName"                ),
  m_option_dataOutputAllowedSet              (m_prefix + "dataOutputAllowedSet"              ),
  m_option_totallyMute                       (m_prefix + "totallyMute"                       ),
  m_option_rawChain_dataInputFileName        (m_prefix + "rawChain_dataInputFileName"        ),
  m_option_rawChain_dataInputFileType        (m_prefix + "rawChain_dataInputFileType"        ),
  m_option_rawChain_size                     (m_prefix + "rawChain_size"                     ),
  m_option_rawChain_generateExtra            (m_prefix + "rawChain_generateExtra"            ),
  m_option_rawChain_displayPeriod            (m_prefix + "rawChain_displayPeriod"            ),
  m_option_rawChain_measureRunTimes          (m_prefix + "rawChain_measureRunTimes"          ),
  m_option_rawChain_dataOutputFileName       (m_prefix + "rawChain_dataOutputFileName"       ),
  m_option_rawChain_dataOutputFileType       (m_prefix + "rawChain_dataOutputFileType"       ),
  m_option_rawChain_dataOutputAllowedSet     (m_prefix + "rawChain_dataOutputAllowedSet"     ),
  m_option_rawChain_computeStats             (m_prefix + "rawChain_computeStats"             ),
  m_option_filteredChain_generate            (m_prefix + "filteredChain_generate"            ),
  m_option_filteredChain_discardedPortion    (m_prefix + "filteredChain_discardedPortion"    ),
  m_option_filteredChain_lag                 (m_prefix + "filteredChain_lag"                 ),
  m_option_filteredChain_dataOutputFileName  (m_prefix + "filteredChain_dataOutputFileName"  ),
  m_option_filteredChain_dataOutputFileType  (m_prefix + "filteredChain_dataOutputFileType"  ),
  m_option_filteredChain_dataOutputAllowedSet(m_prefix + "filteredChain_dataOutputAllowedSet"),
  m_option_filteredChain_computeStats        (m_prefix + "filteredChain_computeStats"        ),
  m_option_displayCandidates                 (m_prefix + "displayCandidates"                 ),
  m_option_putOutOfBoundsInChain             (m_prefix + "putOutOfBoundsInChain"             ),
  m_option_tk_useLocalHessian                (m_prefix + "tk_useLocalHessian"                ),
  m_option_tk_useNewtonComponent             (m_prefix + "tk_useNewtonComponent"             ),
  m_option_dr_maxNumExtraStages              (m_prefix + "dr_maxNumExtraStages"              ),
  m_option_dr_listOfScalesForExtraStages     (m_prefix + "dr_listOfScalesForExtraStages"     ),
  m_option_dr_duringAmNonAdaptiveInt         (m_prefix + "dr_duringAmNonAdaptiveInt"         ),
  m_option_am_keepInitialMatrix              (m_prefix + "am_keepInitialMatrix"              ),
  m_option_am_initialNonAdaptInterval        (m_prefix + "am_initialNonAdaptInterval"        ),
  m_option_am_adaptInterval                  (m_prefix + "am_adaptInterval"                  ),
  m_option_am_eta                            (m_prefix + "am_eta"                            ),
  m_option_am_epsilon                        (m_prefix + "am_epsilon"                        ),
  m_option_enableBrooksGelmanConvMonitor     (m_prefix + "enableBrooksGelmanConvMonitor"     ),
  m_option_BrooksGelmanLag                   (m_prefix + "BrooksGelmanLag"                   )
{
  m_optionsValues.m_dataOutputFileName                 = inputOptions.m_dataOutputFileName;
  m_optionsValues.m_dataOutputAllowedSet               = inputOptions.m_dataOutputAllowedSet;
  m_optionsValues.m_totallyMute                        = inputOptions.m_totallyMute;
  m_optionsValues.m_rawChainDataInputFileName          = inputOptions.m_rawChainDataInputFileName;
  m_optionsValues.m_rawChainDataInputFileType          = inputOptions.m_rawChainDataInputFileType;
  m_optionsValues.m_rawChainSize                       = inputOptions.m_rawChainSize;
  m_optionsValues.m_rawChainGenerateExtra              = inputOptions.m_rawChainGenerateExtra;
  m_optionsValues.m_rawChainDisplayPeriod              = inputOptions.m_rawChainDisplayPeriod;
  m_optionsValues.m_rawChainMeasureRunTimes            = inputOptions.m_rawChainMeasureRunTimes;
  m_optionsValues.m_rawChainDataOutputFileName         = inputOptions.m_rawChainDataOutputFileName;
  m_optionsValues.m_rawChainDataOutputFileType         = inputOptions.m_rawChainDataOutputFileType;
  m_optionsValues.m_rawChainDataOutputAllowedSet       = inputOptions.m_rawChainDataOutputAllowedSet;
  m_optionsValues.m_rawChainComputeStats               = inputOptions.m_rawChainComputeStats;
  //m_optionsValues.m_rawChainStatisticalOptionsValues   = inputOptions.;  // dakota
  m_optionsValues.m_filteredChainGenerate              = inputOptions.m_filteredChainGenerate;
  m_optionsValues.m_filteredChainDiscardedPortion      = inputOptions.m_filteredChainDiscardedPortion;
  m_optionsValues.m_filteredChainLag                   = inputOptions.m_filteredChainLag;
  m_optionsValues.m_filteredChainDataOutputFileName    = inputOptions.m_filteredChainDataOutputFileName;
  m_optionsValues.m_filteredChainDataOutputFileType    = inputOptions.m_filteredChainDataOutputFileType;
  m_optionsValues.m_filteredChainDataOutputAllowedSet  = inputOptions.m_filteredChainDataOutputAllowedSet;
  m_optionsValues.m_filteredChainComputeStats          = inputOptions.m_filteredChainComputeStats;
  //m_optionsValues.m_filteredChainStatisticalOptionsValues = inputOptions.;  // dakota
  m_optionsValues.m_displayCandidates                  = inputOptions.m_displayCandidates;
  m_optionsValues.m_putOutOfBoundsInChain              = inputOptions.m_putOutOfBoundsInChain;
  m_optionsValues.m_tkUseLocalHessian                  = inputOptions.m_tkUseLocalHessian;
  m_optionsValues.m_tkUseNewtonComponent               = inputOptions.m_tkUseNewtonComponent;
  m_optionsValues.m_drMaxNumExtraStages                = inputOptions.m_drMaxNumExtraStages;
  m_optionsValues.m_drScalesForExtraStages             = inputOptions.m_drScalesForExtraStages;
  m_optionsValues.m_drDuringAmNonAdaptiveInt           = inputOptions.m_drDuringAmNonAdaptiveInt;
  m_optionsValues.m_amKeepInitialMatrix                = inputOptions.m_amKeepInitialMatrix;
  m_optionsValues.m_amInitialNonAdaptInterval          = inputOptions.m_amInitialNonAdaptInterval;
  m_optionsValues.m_amAdaptInterval                    = inputOptions.m_amAdaptInterval;
  m_optionsValues.m_amEta                              = inputOptions.m_amEta;
  m_optionsValues.m_amEpsilon                          = inputOptions.m_amEpsilon;
  m_optionsValues.m_enableBrooksGelmanConvMonitor      = UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR;
  m_optionsValues.m_BrooksGelmanLag                    = UQ_MH_SG_BROOKS_GELMAN_LAG;

  m_rawChainStatisticalOptionsObj      = inputOptions.m_rawChainStatisticalOptionsObj; // dakota
  m_rawChainStatOptsInstantiated       = false;
  m_filteredChainStatisticalOptionsObj = inputOptions.m_filteredChainStatisticalOptionsObj; // dakota
  m_filteredChainStatOptsInstantiated  = false;

  if ((m_env.subDisplayFile()        != NULL ) &&
      (m_optionsValues.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGOptionsClass::constructor(2)"
                            << ": after copying values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

uqMetropolisHastingsSGOptionsClass::~uqMetropolisHastingsSGOptionsClass()
{
  if (m_filteredChainStatOptsInstantiated) delete m_filteredChainStatisticalOptionsObj;
  if (m_rawChainStatOptsInstantiated     ) delete m_rawChainStatisticalOptionsObj;
  if (m_optionsDesc                      ) delete m_optionsDesc;
} 

void
uqMetropolisHastingsSGOptionsClass::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "uqMetropolisHastingsSGOptionsClass::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if ((m_env.subDisplayFile() != NULL) &&
      (m_optionsValues.m_totallyMute == false        )) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

  if (m_optionsValues.m_rawChainComputeStats) {
    m_rawChainStatisticalOptionsObj = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "rawChain_");
    m_rawChainStatOptsInstantiated  = true;
  }
  if (m_optionsValues.m_filteredChainComputeStats) {
    m_filteredChainStatisticalOptionsObj = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "filteredChain_");
    m_filteredChainStatOptsInstantiated  = true;
  }

  return;
}

void
uqMetropolisHastingsSGOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                                              "produce help msg for Bayesian Metropolis-Hastings"          )
    (m_option_dataOutputFileName.c_str(),                 po::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                 ), "name of generic output file"                                )
    (m_option_dataOutputAllowedSet.c_str(),               po::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV               ), "subEnvs that will write to generic output file"             )
    (m_option_totallyMute.c_str(),                        po::value<bool        >()->default_value(UQ_MH_SG_TOTALLY_MUTE_ODV                          ), "totally mute (no printout msg)"                             )
    (m_option_rawChain_dataInputFileName.c_str(),         po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV        ), "name of input file for raw chain "                          )
    (m_option_rawChain_dataInputFileType.c_str(),         po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV        ), "type of input file for raw chain "                          )
    (m_option_rawChain_size.c_str(),                      po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_SIZE_ODV                        ), "size of raw chain"                                          )
    (m_option_rawChain_generateExtra.c_str(),             po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV              ), "generate extra information about raw chain"                 )
    (m_option_rawChain_displayPeriod.c_str(),             po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV              ), "period of msg display during raw chain generation"          )
    (m_option_rawChain_measureRunTimes.c_str(),           po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV           ), "measure run times"                                          )
    (m_option_rawChain_dataOutputFileName.c_str(),        po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV       ), "name of output file for raw chain "                         )
    (m_option_rawChain_dataOutputFileType.c_str(),        po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV       ), "type of output file for raw chain "                         )
    (m_option_rawChain_dataOutputAllowedSet.c_str(),      po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV     ), "subEnvs that will write raw chain to output file"           )
    (m_option_rawChain_computeStats.c_str(),              po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV               ), "compute statistics on raw chain"                            )
    (m_option_filteredChain_generate.c_str(),             po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV               ), "generate filtered chain"                                    )
    (m_option_filteredChain_discardedPortion.c_str(),     po::value<double      >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV      ), "initial discarded portion for chain filtering"              )
    (m_option_filteredChain_lag.c_str(),                  po::value<unsigned int>()->default_value(UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                    ), "spacing for chain filtering"                                )
    (m_option_filteredChain_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV  ), "name of output file for filtered chain"                     )
    (m_option_filteredChain_dataOutputFileType.c_str(),   po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV  ), "type of output file for filtered chain"                     )
    (m_option_filteredChain_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write filt chain to output file"          )
    (m_option_filteredChain_computeStats.c_str(),         po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV          ), "compute statistics on filtered chain"                       )
    (m_option_displayCandidates.c_str(),                  po::value<bool        >()->default_value(UQ_MH_SG_DISPLAY_CANDIDATES_ODV                    ), "display candidates in the core MH algorithm"                )
    (m_option_putOutOfBoundsInChain.c_str(),              po::value<bool        >()->default_value(UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV            ), "put 'out of bound' candidates in chain as well"             )
    (m_option_tk_useLocalHessian.c_str(),                 po::value<bool        >()->default_value(UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                  ), "'proposal' use local Hessian"                               )
    (m_option_tk_useNewtonComponent.c_str(),              po::value<bool        >()->default_value(UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV               ), "'proposal' use Newton component"                            )
    (m_option_dr_maxNumExtraStages.c_str(),               po::value<unsigned int>()->default_value(UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV               ), "'dr' maximum number of extra stages"                        )
    (m_option_dr_listOfScalesForExtraStages.c_str(),      po::value<std::string >()->default_value(UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV    ), "'dr' scales for prop cov matrices from 2nd stage on"        )
    (m_option_dr_duringAmNonAdaptiveInt.c_str(),          po::value<bool        >()->default_value(UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV         ), "'dr' used during 'am' non adaptive interval"                )
    (m_option_am_keepInitialMatrix.c_str(),               po::value<bool        >()->default_value(UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV                ), "'am' keep initial (given) matrix"                           )
    (m_option_am_initialNonAdaptInterval.c_str(),         po::value<unsigned int>()->default_value(UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                 ), "'am' initial non adaptation interval"                       )
    (m_option_am_adaptInterval.c_str(),                   po::value<unsigned int>()->default_value(UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                     ), "'am' adaptation interval"                                   )
    (m_option_am_eta.c_str(),                             po::value<double      >()->default_value(UQ_MH_SG_AM_ETA_ODV                                ), "'am' eta"                                                   )
    (m_option_am_epsilon.c_str(),                         po::value<double      >()->default_value(UQ_MH_SG_AM_EPSILON_ODV                            ), "'am' epsilon"                                               )
    (m_option_enableBrooksGelmanConvMonitor.c_str(),      po::value<unsigned int>()->default_value(UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR         ), "assess convergence using Brooks-Gelman metric"              )
    (m_option_BrooksGelmanLag.c_str(),                    po::value<unsigned int>()->default_value(UQ_MH_SG_BROOKS_GELMAN_LAG                         ), "number of chain positions before starting to compute metric")
  ;

  return;
}

void
uqMetropolisHastingsSGOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if ((m_env.subDisplayFile()) &&
        (m_optionsValues.m_totallyMute == false)) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_optionsValues.m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_optionsValues.m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_optionsValues.m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_totallyMute)) {
    m_optionsValues.m_totallyMute = ((const po::variable_value&) m_env.allOptionsMap()[m_option_totallyMute]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileName)) {
    m_optionsValues.m_rawChainDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileType)) {
    m_optionsValues.m_rawChainDataInputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_size)) {
    m_optionsValues.m_rawChainSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_size]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_displayPeriod)) {
    m_optionsValues.m_rawChainDisplayPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_displayPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_measureRunTimes)) {
    m_optionsValues.m_rawChainMeasureRunTimes = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_measureRunTimes]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileName)) {
    m_optionsValues.m_rawChainDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileType)) {
    m_optionsValues.m_rawChainDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputAllowedSet)) {
    m_optionsValues.m_rawChainDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_rawChain_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_optionsValues.m_rawChainDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_computeStats)) {
    m_optionsValues.m_rawChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_computeStats]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_generateExtra)) {
    m_optionsValues.m_rawChainGenerateExtra = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_generateExtra]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_generate)) {
    m_optionsValues.m_filteredChainGenerate = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_generate]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_discardedPortion)) {
    m_optionsValues.m_filteredChainDiscardedPortion = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_discardedPortion]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_lag)) {
    m_optionsValues.m_filteredChainLag = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_lag]).as<unsigned int>();
  }
  if ((m_optionsValues.m_filteredChainGenerate == true) &&
      (m_optionsValues.m_filteredChainLag      < 2    )) {
    std::cerr << "WARNING In uqMetropolisHastingsSGClass<P_V,P_M>::getMyOptionsValues()"
              << ", worldRank "             << m_env.worldRank()
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_filteredChain_lag
              << "' from "                  << m_optionsValues.m_filteredChainLag
              << " to "                     << 2
              << std::endl;
    m_optionsValues.m_filteredChainLag = 2;
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputFileName)) {
    m_optionsValues.m_filteredChainDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputFileType)) {
    m_optionsValues.m_filteredChainDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputAllowedSet)) {
    m_optionsValues.m_filteredChainDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_filteredChain_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_optionsValues.m_filteredChainDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_computeStats)) {
    m_optionsValues.m_filteredChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_computeStats]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_displayCandidates)) {
    m_optionsValues.m_displayCandidates = ((const po::variable_value&) m_env.allOptionsMap()[m_option_displayCandidates]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_putOutOfBoundsInChain)) {
    m_optionsValues.m_putOutOfBoundsInChain = ((const po::variable_value&) m_env.allOptionsMap()[m_option_putOutOfBoundsInChain]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useLocalHessian)) {
    m_optionsValues.m_tkUseLocalHessian = ((const po::variable_value&) m_env.allOptionsMap()[m_option_tk_useLocalHessian]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useNewtonComponent)) {
    m_optionsValues.m_tkUseNewtonComponent = ((const po::variable_value&) m_env.allOptionsMap()[m_option_tk_useNewtonComponent]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dr_maxNumExtraStages)) {
    m_optionsValues.m_drMaxNumExtraStages = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_maxNumExtraStages]).as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count(m_option_dr_listOfScalesForExtraStages)) {
    std::string inputString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_listOfScalesForExtraStages]).as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpScales);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::getMyOptionValues(): scales =";
    //  for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpScales[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}
  }

  if (m_optionsValues.m_drMaxNumExtraStages > 0) {
    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();

    m_optionsValues.m_drScalesForExtraStages.clear();
    m_optionsValues.m_drScalesForExtraStages.resize(m_optionsValues.m_drMaxNumExtraStages,1.);
    for (unsigned int i = 0; i < m_optionsValues.m_drMaxNumExtraStages; ++i) {
      if (i < tmpSize) scale = tmpScales[i];
      m_optionsValues.m_drScalesForExtraStages[i] = scale;
    }
    //updateTK();
  }

  if (m_env.allOptionsMap().count(m_option_dr_duringAmNonAdaptiveInt)) {
    m_optionsValues.m_drDuringAmNonAdaptiveInt = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_duringAmNonAdaptiveInt]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_am_keepInitialMatrix)) {
    m_optionsValues.m_amKeepInitialMatrix = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_keepInitialMatrix]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_am_initialNonAdaptInterval)) {
    m_optionsValues.m_amInitialNonAdaptInterval = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_initialNonAdaptInterval]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptInterval)) {
    m_optionsValues.m_amAdaptInterval = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_adaptInterval]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_eta)) {
    m_optionsValues.m_amEta = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_eta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_am_epsilon)) {
    m_optionsValues.m_amEpsilon = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_epsilon]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_enableBrooksGelmanConvMonitor)) {
    m_optionsValues.m_enableBrooksGelmanConvMonitor = ((const po::variable_value&) m_env.allOptionsMap()[m_option_enableBrooksGelmanConvMonitor]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_BrooksGelmanLag)) {
    m_optionsValues.m_BrooksGelmanLag = ((const po::variable_value&) m_env.allOptionsMap()[m_option_BrooksGelmanLag]).as<unsigned int>();
  }
  return;
}

void
uqMetropolisHastingsSGOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_dataOutputFileName   << " = " << m_optionsValues.m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_optionsValues.m_dataOutputAllowedSet.begin(); setIt != m_optionsValues.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_totallyMute                   << " = " << m_optionsValues.m_totallyMute
     << "\n" << m_option_rawChain_dataInputFileName    << " = " << m_optionsValues.m_rawChainDataInputFileName
     << "\n" << m_option_rawChain_dataInputFileType    << " = " << m_optionsValues.m_rawChainDataInputFileType
     << "\n" << m_option_rawChain_size                 << " = " << m_optionsValues.m_rawChainSize
     << "\n" << m_option_rawChain_generateExtra        << " = " << m_optionsValues.m_rawChainGenerateExtra
     << "\n" << m_option_rawChain_displayPeriod        << " = " << m_optionsValues.m_rawChainDisplayPeriod
     << "\n" << m_option_rawChain_measureRunTimes      << " = " << m_optionsValues.m_rawChainMeasureRunTimes
     << "\n" << m_option_rawChain_dataOutputFileName   << " = " << m_optionsValues.m_rawChainDataOutputFileName
     << "\n" << m_option_rawChain_dataOutputFileType   << " = " << m_optionsValues.m_rawChainDataOutputFileType
     << "\n" << m_option_rawChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_optionsValues.m_rawChainDataOutputAllowedSet.begin(); setIt != m_optionsValues.m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_rawChain_computeStats              << " = " << m_optionsValues.m_rawChainComputeStats
     << "\n" << m_option_filteredChain_generate             << " = " << m_optionsValues.m_filteredChainGenerate
     << "\n" << m_option_filteredChain_discardedPortion     << " = " << m_optionsValues.m_filteredChainDiscardedPortion
     << "\n" << m_option_filteredChain_lag                  << " = " << m_optionsValues.m_filteredChainLag
     << "\n" << m_option_filteredChain_dataOutputFileName   << " = " << m_optionsValues.m_filteredChainDataOutputFileName
     << "\n" << m_option_filteredChain_dataOutputFileType   << " = " << m_optionsValues.m_filteredChainDataOutputFileType
     << "\n" << m_option_filteredChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_optionsValues.m_filteredChainDataOutputAllowedSet.begin(); setIt != m_optionsValues.m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_filteredChain_computeStats    << " = " << m_optionsValues.m_filteredChainComputeStats
     << "\n" << m_option_displayCandidates             << " = " << m_optionsValues.m_displayCandidates
     << "\n" << m_option_putOutOfBoundsInChain         << " = " << m_optionsValues.m_putOutOfBoundsInChain
     << "\n" << m_option_tk_useLocalHessian            << " = " << m_optionsValues.m_tkUseLocalHessian
     << "\n" << m_option_tk_useNewtonComponent         << " = " << m_optionsValues.m_tkUseNewtonComponent
     << "\n" << m_option_dr_maxNumExtraStages          << " = " << m_optionsValues.m_drMaxNumExtraStages
     << "\n" << m_option_dr_listOfScalesForExtraStages << " = ";
  for (unsigned int i = 0; i < m_optionsValues.m_drScalesForExtraStages.size(); ++i) {
    os << m_optionsValues.m_drScalesForExtraStages[i] << " ";
  }
  os << "\n" << m_option_dr_duringAmNonAdaptiveInt     << " = " << m_optionsValues.m_drDuringAmNonAdaptiveInt
     << "\n" << m_option_am_keepInitialMatrix          << " = " << m_optionsValues.m_amKeepInitialMatrix
     << "\n" << m_option_am_initialNonAdaptInterval    << " = " << m_optionsValues.m_amInitialNonAdaptInterval
     << "\n" << m_option_am_adaptInterval              << " = " << m_optionsValues.m_amAdaptInterval
     << "\n" << m_option_am_eta                        << " = " << m_optionsValues.m_amEta
     << "\n" << m_option_am_epsilon                    << " = " << m_optionsValues.m_amEpsilon
     << "\n" << m_option_enableBrooksGelmanConvMonitor << " = " << m_optionsValues.m_enableBrooksGelmanConvMonitor
     << "\n" << m_option_BrooksGelmanLag               << " = " << m_optionsValues.m_BrooksGelmanLag
     << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqMetropolisHastingsSGOptionsClass& obj)
{
  obj.print(os);

  return os;
}
