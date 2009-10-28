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

#include <uqMLSamplingLevelOptions.h>
#include <uqMiscellaneous.h>

uqMLSamplingLevelOptionsClass::uqMLSamplingLevelOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_prefix                                   ((std::string)(prefix) + ""),
  m_dataOutputFileName                       (UQ_ML_SAMPLING_L_DATA_OUTPUT_FILE_NAME_ODV),
//m_dataOutputAllowedSet                     (),
  m_str1                                     (""),
  m_minEffectiveSizeRatio                    (UQ_ML_SAMPLING_L_MIN_EFFECTIVE_SIZE_RATIO_ODV),
  m_maxEffectiveSizeRatio                    (UQ_ML_SAMPLING_L_MAX_EFFECTIVE_SIZE_RATIO_ODV),
  m_scaleCovMatrix                           (UQ_ML_SAMPLING_L_SCALE_COV_MATRIX_ODV),
  m_minRejectionRate                         (UQ_ML_SAMPLING_L_MIN_REJECTION_RATE_ODV),
  m_maxRejectionRate                         (UQ_ML_SAMPLING_L_MAX_REJECTION_RATE_ODV),
  m_covRejectionRate                         (UQ_ML_SAMPLING_L_COV_REJECTION_RATE_ODV),
  m_totallyMute                              (UQ_ML_SAMPLING_L_TOTALLY_MUTE_ODV),
  m_rawChainDataInputFileName                (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV),
  m_rawChainSize                             (UQ_ML_SAMPLING_L_RAW_CHAIN_SIZE_ODV),
  m_rawChainGenerateExtra                    (UQ_ML_SAMPLING_L_RAW_CHAIN_GENERATE_EXTRA_ODV),
  m_rawChainDisplayPeriod                    (UQ_ML_SAMPLING_L_RAW_CHAIN_DISPLAY_PERIOD_ODV),
  m_rawChainMeasureRunTimes                  (UQ_ML_SAMPLING_L_RAW_CHAIN_MEASURE_RUN_TIMES_ODV),
  m_rawChainDataOutputFileName               (UQ_ML_SAMPLING_L_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
//m_rawChainDataOutputAllowedSet             (),
  m_str2                                     (""),
  m_rawChainComputeStats                     (UQ_ML_SAMPLING_L_RAW_CHAIN_COMPUTE_STATS_ODV),
  m_rawChainStatisticalOptions               (NULL),
  m_rawChainStatOptsInstantiated             (false),
  m_filteredChainGenerate                    (UQ_ML_SAMPLING_L_FILTERED_CHAIN_GENERATE_ODV),
  m_filteredChainDiscardedPortion            (UQ_ML_SAMPLING_L_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
  m_filteredChainLag                         (UQ_ML_SAMPLING_L_FILTERED_CHAIN_LAG_ODV),
  m_filteredChainDataOutputFileName          (UQ_ML_SAMPLING_L_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
//m_filteredChainDataOutputAllowedSet        (),
  m_str3                                     (""),
  m_filteredChainComputeStats                (UQ_ML_SAMPLING_L_FILTERED_CHAIN_COMPUTE_STATS_ODV),
  m_filteredChainStatisticalOptions          (NULL),
  m_filteredChainStatOptsInstantiated        (false),
  m_displayCandidates                        (UQ_ML_SAMPLING_L_DISPLAY_CANDIDATES_ODV),
  m_putOutOfBoundsInChain                    (UQ_ML_SAMPLING_L_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV),
  m_tkUseLocalHessian                        (UQ_ML_SAMPLING_L_TK_USE_LOCAL_HESSIAN_ODV),
  m_tkUseNewtonComponent                     (UQ_ML_SAMPLING_L_TK_USE_NEWTON_COMPONENT_ODV),
  m_drMaxNumExtraStages                      (UQ_ML_SAMPLING_L_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_drScalesForExtraStages                   (0),
  m_str4                                     ("1. "),
  m_amInitialNonAdaptInterval                (UQ_ML_SAMPLING_L_AM_INIT_NON_ADAPT_INT_ODV),
  m_amAdaptInterval                          (UQ_ML_SAMPLING_L_AM_ADAPT_INTERVAL_ODV),
  m_amEta                                    (UQ_ML_SAMPLING_L_AM_ETA_ODV),
  m_amEpsilon                                (UQ_ML_SAMPLING_L_AM_EPSILON_ODV),
  m_env                                      (env),
  m_optionsDesc                              (new po::options_description("Multilevel sampling level options")),
  m_option_help                              (m_prefix + "help"                              ),
  m_option_dataOutputFileName                (m_prefix + "dataOutputFileName"                ),
  m_option_dataOutputAllowedSet              (m_prefix + "dataOutputAllowedSet"              ),
  m_option_minEffectiveSizeRatio             (m_prefix + "minEffectiveSizeRatio"             ),
  m_option_maxEffectiveSizeRatio             (m_prefix + "maxEffectiveSizeRatio"             ),
  m_option_scaleCovMatrix                    (m_prefix + "scaleCovMatrix"                    ),
  m_option_minRejectionRate                  (m_prefix + "minRejectionRate"                  ),
  m_option_maxRejectionRate                  (m_prefix + "maxRejectionRate"                  ),
  m_option_covRejectionRate                  (m_prefix + "covRejectionRate"                  ),
  m_option_totallyMute                       (m_prefix + "totallyMute"                       ),
  m_option_rawChain_dataInputFileName        (m_prefix + "rawChain_dataInputFileName"        ),
  m_option_rawChain_size                     (m_prefix + "rawChain_size"                     ),
  m_option_rawChain_generateExtra            (m_prefix + "rawChain_generateExtra"            ),
  m_option_rawChain_displayPeriod            (m_prefix + "rawChain_displayPeriod"            ),
  m_option_rawChain_measureRunTimes          (m_prefix + "rawChain_measureRunTimes"          ),
  m_option_rawChain_dataOutputFileName       (m_prefix + "rawChain_dataOutputFileName"       ),
  m_option_rawChain_dataOutputAllowedSet     (m_prefix + "rawChain_dataOutputAllowedSet"     ),
  m_option_rawChain_computeStats             (m_prefix + "rawChain_computeStats"             ),
  m_option_filteredChain_generate            (m_prefix + "filteredChain_generate"            ),
  m_option_filteredChain_discardedPortion    (m_prefix + "filteredChain_discardedPortion"    ),
  m_option_filteredChain_lag                 (m_prefix + "filteredChain_lag"                 ),
  m_option_filteredChain_dataOutputFileName  (m_prefix + "filteredChain_dataOutputFileName"  ),
  m_option_filteredChain_dataOutputAllowedSet(m_prefix + "filteredChain_dataOutputAllowedSet"),
  m_option_filteredChain_computeStats        (m_prefix + "filteredChain_computeStats"        ),
  m_option_displayCandidates                 (m_prefix + "displayCandidates"                 ),
  m_option_putOutOfBoundsInChain             (m_prefix + "putOutOfBoundsInChain"             ),
  m_option_tk_useLocalHessian                (m_prefix + "tk_useLocalHessian"                ),
  m_option_tk_useNewtonComponent             (m_prefix + "tk_useNewtonComponent"             ),
  m_option_dr_maxNumExtraStages              (m_prefix + "dr_maxNumExtraStages"              ),
  m_option_dr_listOfScalesForExtraStages     (m_prefix + "dr_listOfScalesForExtraStages"     ),
  m_option_am_initialNonAdaptInterval        (m_prefix + "am_initialNonAdaptInterval"        ),
  m_option_am_adaptInterval                  (m_prefix + "am_adaptInterval"                  ),
  m_option_am_eta                            (m_prefix + "am_eta"                            ),
  m_option_am_epsilon                        (m_prefix + "am_epsilon"                        )
{
}

void
uqMLSamplingLevelOptionsClass::copyOptionsValues(const uqMLSamplingLevelOptionsClass& srcOptions)
{
  m_dataOutputFileName                = srcOptions.m_dataOutputFileName;
  m_dataOutputAllowedSet              = srcOptions.m_dataOutputAllowedSet;
  m_str1                              = srcOptions.m_str1;
  m_minEffectiveSizeRatio             = srcOptions.m_minEffectiveSizeRatio;
  m_maxEffectiveSizeRatio             = srcOptions.m_maxEffectiveSizeRatio;
  m_scaleCovMatrix                    = srcOptions.m_scaleCovMatrix;
  m_minRejectionRate                  = srcOptions.m_minRejectionRate;
  m_maxRejectionRate                  = srcOptions.m_maxRejectionRate;
  m_covRejectionRate                  = srcOptions.m_covRejectionRate;
  m_totallyMute                       = srcOptions.m_totallyMute;
  m_rawChainDataInputFileName         = srcOptions.m_rawChainDataInputFileName;
  m_rawChainSize                      = srcOptions.m_rawChainSize;
//std::cout << "In copy(), rawChainSize = " << m_rawChainSize << std::endl;
  m_rawChainGenerateExtra             = srcOptions.m_rawChainGenerateExtra;
  m_rawChainDisplayPeriod             = srcOptions.m_rawChainDisplayPeriod;
  m_rawChainMeasureRunTimes           = srcOptions.m_rawChainMeasureRunTimes;
  m_rawChainDataOutputFileName        = srcOptions.m_rawChainDataOutputFileName;
  m_rawChainDataOutputAllowedSet      = srcOptions.m_rawChainDataOutputAllowedSet;
  m_str2                              = srcOptions.m_str2;
  m_rawChainComputeStats              = srcOptions.m_rawChainComputeStats;
  m_rawChainStatisticalOptions        = NULL; // Yes, 'NULL'
  m_rawChainStatOptsInstantiated      = false;
  m_filteredChainGenerate             = srcOptions.m_filteredChainGenerate;
  m_filteredChainDiscardedPortion     = srcOptions.m_filteredChainDiscardedPortion;
  m_filteredChainLag                  = srcOptions.m_filteredChainLag;
  m_filteredChainDataOutputFileName   = srcOptions.m_filteredChainDataOutputFileName;
  m_filteredChainDataOutputAllowedSet = srcOptions.m_filteredChainDataOutputAllowedSet;
  m_str3                              = srcOptions.m_str3;
  m_filteredChainComputeStats         = srcOptions.m_filteredChainComputeStats;
  m_filteredChainStatisticalOptions   = NULL; // Yes, 'NULL'
  m_filteredChainStatOptsInstantiated = false;
  m_displayCandidates                 = srcOptions.m_displayCandidates;
  m_putOutOfBoundsInChain             = srcOptions.m_putOutOfBoundsInChain;
  m_tkUseLocalHessian                 = srcOptions.m_tkUseLocalHessian;
  m_tkUseNewtonComponent              = srcOptions.m_tkUseNewtonComponent;
  m_drMaxNumExtraStages               = srcOptions.m_drMaxNumExtraStages;
  m_drScalesForExtraStages            = srcOptions.m_drScalesForExtraStages;
  m_str4                              = srcOptions.m_str4;
  m_amInitialNonAdaptInterval         = srcOptions.m_amInitialNonAdaptInterval;
  m_amAdaptInterval                   = srcOptions.m_amAdaptInterval;
  m_amEta                             = srcOptions.m_amEta;
  m_amEpsilon                         = srcOptions.m_amEpsilon;

  return;
}

uqMLSamplingLevelOptionsClass::~uqMLSamplingLevelOptionsClass()
{
  //std::cout << "In uqMLSamplingLevelOptionsClass::destructor()"
  //          << ": m_filteredChainStatOptsInstantiated = " << m_filteredChainStatOptsInstantiated
  //          << ", m_rawChainStatOptsInstantiated = "      << m_rawChainStatOptsInstantiated
  //          << std::endl;
  //sleep(1);
  if (m_filteredChainStatOptsInstantiated) delete m_filteredChainStatisticalOptions;
  if (m_rawChainStatOptsInstantiated     ) delete m_rawChainStatisticalOptions;
  if (m_optionsDesc                      ) delete m_optionsDesc;
} 

void
uqMLSamplingLevelOptionsClass::scanOptionsValues(const uqMLSamplingLevelOptionsClass* defaultOptions)
{
  if (m_optionsDesc == NULL) m_optionsDesc = new po::options_description("Multilevel sampling level options");
  if (defaultOptions) this->copyOptionsValues(*defaultOptions);

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if ((m_env.subDisplayFile() != NULL ) &&
      (1)) { //m_totallyMute          == false)) {
    *m_env.subDisplayFile() << "In uqMLSamplingLevelOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

  if (m_rawChainComputeStats) {
    m_rawChainStatisticalOptions   = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "rawChain_");
    m_rawChainStatOptsInstantiated = true;
  }
  if (m_filteredChainComputeStats) {
    m_filteredChainStatisticalOptions   = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "filteredChain_");
    m_filteredChainStatOptsInstantiated = true;
  }

  return;
};

void
uqMLSamplingLevelOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                                                      "produce help message for Bayesian Markov chain distr. calculator")
    (m_option_dataOutputFileName.c_str(),                 po::value<std::string >()->default_value(m_dataOutputFileName               ), "name of generic output file"                                     )
    (m_option_dataOutputAllowedSet.c_str(),               po::value<std::string >()->default_value(m_str1                             ), "subEnvs that will write to generic output file"                  )
    (m_option_minEffectiveSizeRatio.c_str(),              po::value<double      >()->default_value(m_minEffectiveSizeRatio            ), "minimum allowed effective size ratio wrt previous level"         )
    (m_option_maxEffectiveSizeRatio.c_str(),              po::value<double      >()->default_value(m_maxEffectiveSizeRatio            ), "maximum allowed effective size ratio wrt previous level"         )
    (m_option_scaleCovMatrix.c_str(),                     po::value<bool        >()->default_value(m_scaleCovMatrix                   ), "scale proposal covariance matrix"                                )
    (m_option_minRejectionRate.c_str(),                   po::value<double      >()->default_value(m_minRejectionRate                 ), "minimum allowed attempted rejection rate at current level"       )
    (m_option_maxRejectionRate.c_str(),                   po::value<double      >()->default_value(m_maxRejectionRate                 ), "maximum allowed attempted rejection rate at current level"       )
    (m_option_covRejectionRate.c_str(),                   po::value<double      >()->default_value(m_covRejectionRate                 ), "c.o.v. for judging attempted rejection rate at current level"    )
    (m_option_totallyMute.c_str(),                        po::value<bool        >()->default_value(m_totallyMute                      ), "totally mute (no printout message)"                              )
    (m_option_rawChain_dataInputFileName.c_str(),         po::value<std::string >()->default_value(m_rawChainDataInputFileName        ), "name of input file for raw chain "                               )
    (m_option_rawChain_size.c_str(),                      po::value<unsigned int>()->default_value(m_rawChainSize                     ), "size of raw chain"                                               )
    (m_option_rawChain_generateExtra.c_str(),             po::value<bool        >()->default_value(m_rawChainGenerateExtra            ), "generate extra information about raw chain"                      )
    (m_option_rawChain_displayPeriod.c_str(),             po::value<unsigned int>()->default_value(m_rawChainDisplayPeriod            ), "period of message display during raw chain generation"           )
    (m_option_rawChain_measureRunTimes.c_str(),           po::value<bool        >()->default_value(m_rawChainMeasureRunTimes          ), "measure run times"                                               )
    (m_option_rawChain_dataOutputFileName.c_str(),        po::value<std::string >()->default_value(m_rawChainDataOutputFileName       ), "name of output file for raw chain "                              )
    (m_option_rawChain_dataOutputAllowedSet.c_str(),      po::value<std::string >()->default_value(m_str2                             ), "subEnvs that will write to output file for raw chain"            )
    (m_option_rawChain_computeStats.c_str(),              po::value<bool        >()->default_value(m_rawChainComputeStats             ), "compute statistics on raw chain"                                 )
    (m_option_filteredChain_generate.c_str(),             po::value<bool        >()->default_value(m_filteredChainGenerate            ), "generate filtered chain"                                         )
    (m_option_filteredChain_discardedPortion.c_str(),     po::value<double      >()->default_value(m_filteredChainDiscardedPortion    ), "initial discarded portion for chain filtering"                   )
    (m_option_filteredChain_lag.c_str(),                  po::value<unsigned int>()->default_value(m_filteredChainLag                 ), "spacing for chain filtering"                                     )
    (m_option_filteredChain_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(m_filteredChainDataOutputFileName  ), "name of output file for filtered chain"                          )
    (m_option_filteredChain_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(m_str3                             ), "subEnvs that will write to output file for filtered chain"       )
    (m_option_filteredChain_computeStats.c_str(),         po::value<bool        >()->default_value(m_filteredChainComputeStats        ), "compute statistics on filtered chain"                            )
    (m_option_displayCandidates.c_str(),                  po::value<bool        >()->default_value(m_displayCandidates                ), "display candidates generated in the core MH algorithm"           )
    (m_option_putOutOfBoundsInChain.c_str(),              po::value<bool        >()->default_value(m_putOutOfBoundsInChain            ), "put 'out of bound' candidates in chain as well"                  )
    (m_option_tk_useLocalHessian.c_str(),                 po::value<bool        >()->default_value(m_tkUseLocalHessian                ), "'proposal' use local Hessian"                                    )
    (m_option_tk_useNewtonComponent.c_str(),              po::value<bool        >()->default_value(m_tkUseNewtonComponent             ), "'proposal' use Newton component"                                 )
    (m_option_dr_maxNumExtraStages.c_str(),               po::value<unsigned int>()->default_value(m_drMaxNumExtraStages              ), "'dr' maximum number of extra stages"                             )
    (m_option_dr_listOfScalesForExtraStages.c_str(),      po::value<std::string >()->default_value(m_str4                             ), "'dr' list of scales for proposal cov matrices from 2nd stage on" )
    (m_option_am_initialNonAdaptInterval.c_str(),         po::value<unsigned int>()->default_value(m_amInitialNonAdaptInterval        ), "'am' initial non adaptation interval"                            )
    (m_option_am_adaptInterval.c_str(),                   po::value<unsigned int>()->default_value(m_amAdaptInterval                  ), "'am' adaptation interval"                                        )
    (m_option_am_eta.c_str(),                             po::value<double      >()->default_value(m_amEta                            ), "'am' eta"                                                        )
    (m_option_am_epsilon.c_str(),                         po::value<double      >()->default_value(m_amEpsilon                        ), "'am' epsilon"                                                    )
  ;

  return;
}

void
uqMLSamplingLevelOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  char tmpStr[64];

  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName.c_str())) {
    m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet.c_str())) {
    m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }
  m_str1.clear();
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",*setIt);
    m_str1 += tmpStr;
    m_str1 += " ";
  }

  if (m_env.allOptionsMap().count(m_option_minEffectiveSizeRatio.c_str())) {
    m_minEffectiveSizeRatio = ((const po::variable_value&) m_env.allOptionsMap()[m_option_minEffectiveSizeRatio.c_str()]).as<double>();
  }
  if (m_minEffectiveSizeRatio >= 1.) {
    std::cerr << "WARNING In uqMLSamplingLevelOptionsClass::getMyOptionsValues()"
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_minEffectiveSizeRatio.c_str()
              << "' from "                  << m_minEffectiveSizeRatio
              << " to "                     << .5
              << std::endl;
    m_minEffectiveSizeRatio = .5;
  }

  if (m_env.allOptionsMap().count(m_option_maxEffectiveSizeRatio.c_str())) {
    m_maxEffectiveSizeRatio = ((const po::variable_value&) m_env.allOptionsMap()[m_option_maxEffectiveSizeRatio.c_str()]).as<double>();
  }
  if (m_maxEffectiveSizeRatio >= 1.) {
    std::cerr << "WARNING In uqMLSamplingLevelOptionsClass::getMyOptionsValues()"
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_maxEffectiveSizeRatio.c_str()
              << "' from "                  << m_maxEffectiveSizeRatio
              << " to "                     << .5
              << std::endl;
    m_maxEffectiveSizeRatio = .5;
  }

  if (m_env.allOptionsMap().count(m_option_scaleCovMatrix.c_str())) {
    m_scaleCovMatrix = ((const po::variable_value&) m_env.allOptionsMap()[m_option_scaleCovMatrix.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_minRejectionRate.c_str())) {
    m_minRejectionRate = ((const po::variable_value&) m_env.allOptionsMap()[m_option_minRejectionRate.c_str()]).as<double>();
  }
  if (m_minRejectionRate >= 1.) {
    std::cerr << "WARNING In uqMLSamplingLevelOptionsClass::getMyOptionsValues()"
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_minRejectionRate.c_str()
              << "' from "                  << m_minRejectionRate
              << " to "                     << .5
              << std::endl;
    m_minRejectionRate = .5;
  }

  if (m_env.allOptionsMap().count(m_option_maxRejectionRate.c_str())) {
    m_maxRejectionRate = ((const po::variable_value&) m_env.allOptionsMap()[m_option_maxRejectionRate.c_str()]).as<double>();
  }
  if (m_maxRejectionRate >= 1.) {
    std::cerr << "WARNING In uqMLSamplingLevelOptionsClass::getMyOptionsValues()"
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_maxRejectionRate.c_str()
              << "' from "                  << m_maxRejectionRate
              << " to "                     << .5
              << std::endl;
    m_maxRejectionRate = .5;
  }

  if (m_env.allOptionsMap().count(m_option_covRejectionRate.c_str())) {
    m_covRejectionRate = ((const po::variable_value&) m_env.allOptionsMap()[m_option_covRejectionRate.c_str()]).as<double>();
  }
  if (m_covRejectionRate >= 1.) {
    std::cerr << "WARNING In uqMLSamplingLevelOptionsClass::getMyOptionsValues()"
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_covRejectionRate.c_str()
              << "' from "                  << m_covRejectionRate
              << " to "                     << .5
              << std::endl;
    m_covRejectionRate = .5;
  }

  if (m_env.allOptionsMap().count(m_option_totallyMute.c_str())) {
    m_totallyMute = ((const po::variable_value&) m_env.allOptionsMap()[m_option_totallyMute.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileName.c_str())) {
    m_rawChainDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_size.c_str())) {
    //std::cout << "In count()=true, rawChainSize = " << m_rawChainSize << std::endl;
    m_rawChainSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_size.c_str()]).as<unsigned int>();
  }
//std::cout << "After count(), rawChainSize = " << m_rawChainSize << std::endl;

  if (m_env.allOptionsMap().count(m_option_rawChain_displayPeriod.c_str())) {
    m_rawChainDisplayPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_displayPeriod.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_measureRunTimes.c_str())) {
    m_rawChainMeasureRunTimes = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_measureRunTimes.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputFileName.c_str())) {
    m_rawChainDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataOutputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataOutputAllowedSet.c_str())) {
    m_rawChainDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_rawChain_dataOutputAllowedSet.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_rawChainDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }
  m_str2.clear();
  for (std::set<unsigned int>::iterator setIt = m_rawChainDataOutputAllowedSet.begin(); setIt != m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",*setIt);
    m_str2 += tmpStr;
    m_str2 += " ";
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_computeStats.c_str())) {
    m_rawChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_computeStats.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_generateExtra.c_str())) {
    m_rawChainGenerateExtra = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_generateExtra.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_generate.c_str())) {
    m_filteredChainGenerate = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_generate.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_discardedPortion.c_str())) {
    m_filteredChainDiscardedPortion = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_discardedPortion.c_str()]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_lag.c_str())) {
    m_filteredChainLag = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_lag.c_str()]).as<unsigned int>();
  }
  if ((m_filteredChainGenerate == true) &&
      (m_filteredChainLag      < 2    )) {
    std::cerr << "WARNING In uqMLSamplingLevelOptionsClass::getMyOptionsValues()"
              << ", fullRank "              << m_env.fullRank()
              << ", subEnvironment "        << m_env.subId()
              << ", subRank "               << m_env.subRank()
              << ", inter0Rank "            << m_env.inter0Rank()
              << ": forcing the value of '" << m_option_filteredChain_lag.c_str()
              << "' from "                  << m_filteredChainLag
              << " to "                     << 2
              << std::endl;
    m_filteredChainLag = 2;
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputFileName.c_str())) {
    m_filteredChainDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_dataOutputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_dataOutputAllowedSet.c_str())) {
    m_filteredChainDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_filteredChain_dataOutputAllowedSet.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_filteredChainDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }
  m_str3.clear();
  for (std::set<unsigned int>::iterator setIt = m_filteredChainDataOutputAllowedSet.begin(); setIt != m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    sprintf(tmpStr,"%d",*setIt);
    m_str3 += tmpStr;
    m_str3 += " ";
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_computeStats.c_str())) {
    m_filteredChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_computeStats.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_displayCandidates.c_str())) {
    m_displayCandidates = ((const po::variable_value&) m_env.allOptionsMap()[m_option_displayCandidates.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_putOutOfBoundsInChain.c_str())) {
    m_putOutOfBoundsInChain = ((const po::variable_value&) m_env.allOptionsMap()[m_option_putOutOfBoundsInChain.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useLocalHessian.c_str())) {
    m_tkUseLocalHessian = ((const po::variable_value&) m_env.allOptionsMap()[m_option_tk_useLocalHessian.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_tk_useNewtonComponent.c_str())) {
    m_tkUseNewtonComponent = ((const po::variable_value&) m_env.allOptionsMap()[m_option_tk_useNewtonComponent.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dr_maxNumExtraStages.c_str())) {
    m_drMaxNumExtraStages = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_maxNumExtraStages.c_str()]).as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count(m_option_dr_listOfScalesForExtraStages.c_str())) {
    std::string inputString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dr_listOfScalesForExtraStages.c_str()]).as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpScales);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqMLSamplingLevelOptionsClass::getMyOptionValues(): scales =";
    //  for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpScales[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}
  }

  if (m_drMaxNumExtraStages > 0) {
    m_drScalesForExtraStages.clear();

    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();

    m_drScalesForExtraStages.resize(m_drMaxNumExtraStages,1.);

    for (unsigned int i = 0; i < m_drMaxNumExtraStages; ++i) {
      if (i < tmpSize) scale = tmpScales[i];
      m_drScalesForExtraStages[i] = scale;
    }
    //updateTK();
  }

  m_str4.clear();
  for (unsigned int i = 0; i < m_drScalesForExtraStages.size(); ++i) {
    sprintf(tmpStr,"%e",m_drScalesForExtraStages[i]);
    m_str4 += tmpStr;
    m_str4 += " ";
  }
//std::cout << "m_str4 = " << m_str4 << std::endl;

  if (m_env.allOptionsMap().count(m_option_am_initialNonAdaptInterval.c_str())) {
    m_amInitialNonAdaptInterval = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_initialNonAdaptInterval.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptInterval.c_str())) {
    m_amAdaptInterval = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_adaptInterval.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_eta.c_str())) {
    m_amEta = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_eta.c_str()]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_am_epsilon.c_str())) {
    m_amEpsilon = ((const po::variable_value&) m_env.allOptionsMap()[m_option_am_epsilon.c_str()]).as<double>();
  }

  return;
}

void
uqMLSamplingLevelOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_dataOutputFileName   << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_minEffectiveSizeRatio          << " = " << m_minEffectiveSizeRatio
     << "\n" << m_option_maxEffectiveSizeRatio          << " = " << m_maxEffectiveSizeRatio
     << "\n" << m_option_scaleCovMatrix                 << " = " << m_scaleCovMatrix
     << "\n" << m_option_minRejectionRate               << " = " << m_minRejectionRate
     << "\n" << m_option_maxRejectionRate               << " = " << m_maxRejectionRate
     << "\n" << m_option_covRejectionRate               << " = " << m_covRejectionRate
     << "\n" << m_option_totallyMute                    << " = " << m_totallyMute
     << "\n" << m_option_rawChain_dataInputFileName     << " = " << m_rawChainDataInputFileName
     << "\n" << m_option_rawChain_size                  << " = " << m_rawChainSize
     << "\n" << m_option_rawChain_generateExtra         << " = " << m_rawChainGenerateExtra
     << "\n" << m_option_rawChain_displayPeriod         << " = " << m_rawChainDisplayPeriod
     << "\n" << m_option_rawChain_measureRunTimes       << " = " << m_rawChainMeasureRunTimes
     << "\n" << m_option_rawChain_dataOutputFileName    << " = " << m_rawChainDataOutputFileName
     << "\n" << m_option_rawChain_dataOutputAllowedSet  << " = ";
  for (std::set<unsigned int>::iterator setIt = m_rawChainDataOutputAllowedSet.begin(); setIt != m_rawChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_rawChain_computeStats              << " = " << m_rawChainComputeStats
     << "\n" << m_option_filteredChain_generate             << " = " << m_filteredChainGenerate
     << "\n" << m_option_filteredChain_discardedPortion     << " = " << m_filteredChainDiscardedPortion
     << "\n" << m_option_filteredChain_lag                  << " = " << m_filteredChainLag
     << "\n" << m_option_filteredChain_dataOutputFileName   << " = " << m_filteredChainDataOutputFileName
     << "\n" << m_option_filteredChain_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_filteredChainDataOutputAllowedSet.begin(); setIt != m_filteredChainDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_filteredChain_computeStats    << " = " << m_filteredChainComputeStats
     << "\n" << m_option_displayCandidates             << " = " << m_displayCandidates
     << "\n" << m_option_putOutOfBoundsInChain         << " = " << m_putOutOfBoundsInChain
     << "\n" << m_option_tk_useLocalHessian            << " = " << m_tkUseLocalHessian
     << "\n" << m_option_tk_useNewtonComponent         << " = " << m_tkUseNewtonComponent
     << "\n" << m_option_dr_maxNumExtraStages          << " = " << m_drMaxNumExtraStages
     << "\n" << m_option_dr_listOfScalesForExtraStages << " = ";
  for (unsigned int i = 0; i < m_drScalesForExtraStages.size(); ++i) {
    os << m_drScalesForExtraStages[i] << " ";
  }
  os << "\n" << m_option_am_initialNonAdaptInterval << " = " << m_amInitialNonAdaptInterval
     << "\n" << m_option_am_adaptInterval           << " = " << m_amAdaptInterval
     << "\n" << m_option_am_eta                     << " = " << m_amEta
     << "\n" << m_option_am_epsilon                 << " = " << m_amEpsilon
     << std::endl;

  return;
}

const uqBaseEnvironmentClass& 
uqMLSamplingLevelOptionsClass::env() const
{
  return m_env;
}

std::ostream& operator<<(std::ostream& os, const uqMLSamplingLevelOptionsClass& obj)
{
  obj.print(os);

  return os;
}
