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

uqMetropolisHastingsSGOptionsClass::uqMetropolisHastingsSGOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char* prefix)
  :
  m_prefix                                   ((std::string)(prefix) + "mc_"),
  m_dataOutputFileName                       (UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV),
//m_dataOutputAllowedSet                     (),
  m_totallyMute                              (UQ_MH_SG_TOTALLY_MUTE_ODV),
  m_rawChainType                             (UQ_MH_SG_RAW_CHAIN_TYPE_ODV),
  m_rawChainDataInputFileName                (UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV),
  m_rawChainSize                             (UQ_MH_SG_RAW_CHAIN_SIZE_ODV),
  m_rawChainGenerateExtra                    (UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV),
  m_rawChainDisplayPeriod                    (UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV),
  m_rawChainMeasureRunTimes                  (UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV),
  m_rawChainDataOutputFileName               (UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
//m_rawChainDataOutputAllowedSet             (),
  m_rawChainComputeStats                     (UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV),
  m_rawChainStatisticalOptions               (NULL),
  m_rawChainStatOptsInstantiated             (false),
  m_filteredChainGenerate                    (UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV),
  m_filteredChainDiscardedPortion            (UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
  m_filteredChainLag                         (UQ_MH_SG_FILTERED_CHAIN_LAG_ODV),
  m_filteredChainDataOutputFileName          (UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV),
//m_filteredChainDataOutputAllowedSet        (),
  m_filteredChainComputeStats                (UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV),
  m_filteredChainStatisticalOptions          (NULL),
  m_filteredChainStatOptsInstantiated        (false),
  m_mhDisplayCandidates                      (UQ_MH_SG_MH_DISPLAY_CANDIDATES_ODV),
  m_mhPutOutOfBoundsInChain                  (UQ_MH_SG_MH_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV),
  m_tkUseLocalHessian                        (UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV),
  m_tkUseNewtonComponent                     (UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV),
  m_drMaxNumExtraStages                      (UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV),
#ifdef UQ_READS_ONLY_EXTRA_STAGES
  m_drScalesForExtraStages                   (0),
#else
  m_drScalesForExtraStages                   (1,1.),
#endif
  m_amInitialNonAdaptInterval                (UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV),
  m_amAdaptInterval                          (UQ_MH_SG_AM_ADAPT_INTERVAL_ODV),
  m_amEta                                    (UQ_MH_SG_AM_ETA_ODV),
  m_amEpsilon                                (UQ_MH_SG_AM_EPSILON_ODV),
  m_env                                      (env),
  m_optionsDesc                              (new po::options_description("Bayesian Metropolis-Hastings options")),
  m_option_help                              (m_prefix + "help"                              ),
  m_option_dataOutputFileName                (m_prefix + "dataOutputFileName"                ),
  m_option_dataOutputAllowedSet              (m_prefix + "dataOutputAllowedSet"              ),
  m_option_totallyMute                       (m_prefix + "totallyMute"                       ),
  m_option_rawChain_type                     (m_prefix + "rawChain_type"                     ),
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
  m_option_mh_displayCandidates              (m_prefix + "mh_displayCandidates"              ),
  m_option_mh_putOutOfBoundsInChain          (m_prefix + "mh_putOutOfBoundsInChain"          ),
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

uqMetropolisHastingsSGOptionsClass::~uqMetropolisHastingsSGOptionsClass()
{
  if (m_filteredChainStatOptsInstantiated) delete m_filteredChainStatisticalOptions;
  if (m_rawChainStatOptsInstantiated     ) delete m_rawChainStatisticalOptions;
  if (m_optionsDesc                      ) delete m_optionsDesc;
} 

void
uqMetropolisHastingsSGOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if ((m_env.subDisplayFile() != NULL) &&
      (m_totallyMute == false        )) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGOptionsClass::scanOptionsValues()"
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
}

void
uqMetropolisHastingsSGOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                                               "produce help message for Bayesian Metropolis-Hastings")
    (m_option_dataOutputFileName.c_str(),                 po::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                 ), "name of generic output file"                                     )
    (m_option_dataOutputAllowedSet.c_str(),               po::value<std::string >()->default_value(UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV               ), "subEnvs that will write to generic output file"                  )
    (m_option_totallyMute.c_str(),                        po::value<bool        >()->default_value(UQ_MH_SG_TOTALLY_MUTE_ODV                          ), "totally mute (no printout message)"                              )
    (m_option_rawChain_type.c_str(),                      po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_TYPE_ODV                        ), "type of raw chain (1=Markov, 2=White noise)"                     )
    (m_option_rawChain_dataInputFileName.c_str(),         po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV        ), "name of input file for raw chain "                               )
    (m_option_rawChain_size.c_str(),                      po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_SIZE_ODV                        ), "size of raw chain"                                               )
    (m_option_rawChain_generateExtra.c_str(),             po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV              ), "generate extra information about raw chain"                      )
    (m_option_rawChain_displayPeriod.c_str(),             po::value<unsigned int>()->default_value(UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV              ), "period of message display during raw chain generation"           )
    (m_option_rawChain_measureRunTimes.c_str(),           po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV           ), "measure run times"                                               )
    (m_option_rawChain_dataOutputFileName.c_str(),        po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV       ), "name of output file for raw chain "                              )
    (m_option_rawChain_dataOutputAllowedSet.c_str(),      po::value<std::string >()->default_value(UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV     ), "subEnvs that will write to output file for raw chain"            )
    (m_option_rawChain_computeStats.c_str(),              po::value<bool        >()->default_value(UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV               ), "compute statistics on raw chain"                                 )
    (m_option_filteredChain_generate.c_str(),             po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV               ), "generate filtered chain"                                         )
    (m_option_filteredChain_discardedPortion.c_str(),     po::value<double      >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV      ), "initial discarded portion for chain filtering"                   )
    (m_option_filteredChain_lag.c_str(),                  po::value<unsigned int>()->default_value(UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                    ), "spacing for chain filtering"                                     )
    (m_option_filteredChain_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV  ), "name of output file for filtered chain"                          )
    (m_option_filteredChain_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to output file for filtered chain"       )
    (m_option_filteredChain_computeStats.c_str(),         po::value<bool        >()->default_value(UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV          ), "compute statistics on filtered chain"                            )
    (m_option_mh_displayCandidates.c_str(),               po::value<bool        >()->default_value(UQ_MH_SG_MH_DISPLAY_CANDIDATES_ODV                 ), "display candidates generated in the core MH algorithm"           )
    (m_option_mh_putOutOfBoundsInChain.c_str(),           po::value<bool        >()->default_value(UQ_MH_SG_MH_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV         ), "put 'out of bound' candidates in chain as well"                  )
    (m_option_tk_useLocalHessian.c_str(),                 po::value<bool        >()->default_value(UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                  ), "'proposal' use local Hessian"                                    )
    (m_option_tk_useNewtonComponent.c_str(),              po::value<bool        >()->default_value(UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV               ), "'proposal' use Newton component"                                 )
    (m_option_dr_maxNumExtraStages.c_str(),               po::value<unsigned int>()->default_value(UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV               ), "'dr' maximum number of extra stages"                             )
    (m_option_dr_listOfScalesForExtraStages.c_str(),      po::value<std::string >()->default_value(UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV    ), "'dr' list of scales for proposal cov matrices from 2nd stage on" )
    (m_option_am_initialNonAdaptInterval.c_str(),         po::value<unsigned int>()->default_value(UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                 ), "'am' initial non adaptation interval"                            )
    (m_option_am_adaptInterval.c_str(),                   po::value<unsigned int>()->default_value(UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                     ), "'am' adaptation interval"                                        )
    (m_option_am_eta.c_str(),                             po::value<double      >()->default_value(UQ_MH_SG_AM_ETA_ODV                                ), "'am' eta"                                                        )
    (m_option_am_epsilon.c_str(),                         po::value<double      >()->default_value(UQ_MH_SG_AM_EPSILON_ODV                            ), "'am' epsilon"                                                    )
  ;

  return;
}

void
uqMetropolisHastingsSGOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if ((m_env.subDisplayFile()) &&
        (m_totallyMute == false)) {
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

  if (m_env.allOptionsMap().count(m_option_totallyMute.c_str())) {
    m_totallyMute = ((const po::variable_value&) m_env.allOptionsMap()[m_option_totallyMute.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_type.c_str())) {
    m_rawChainType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_type.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_dataInputFileName.c_str())) {
    m_rawChainDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_dataInputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_rawChain_size.c_str())) {
    m_rawChainSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rawChain_size.c_str()]).as<unsigned int>();
  }

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
    std::cerr << "WARNING In uqMetropolisHastingsSGClass<P_V,P_M>::getMyOptionsValues()"
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

  if (m_env.allOptionsMap().count(m_option_filteredChain_computeStats.c_str())) {
    m_filteredChainComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_filteredChain_computeStats.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_displayCandidates.c_str())) {
    m_mhDisplayCandidates = ((const po::variable_value&) m_env.allOptionsMap()[m_option_mh_displayCandidates.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_putOutOfBoundsInChain.c_str())) {
    m_mhPutOutOfBoundsInChain = ((const po::variable_value&) m_env.allOptionsMap()[m_option_mh_putOutOfBoundsInChain.c_str()]).as<bool>();
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
    //  *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::getMyOptionValues(): scales =";
    //  for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpScales[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}
  }

  if (m_drMaxNumExtraStages > 0) {
    m_drScalesForExtraStages.clear();
#ifdef UQ_USES_TK_CLASS
#else
    m_lowerCholProposalCovMatrices.clear();
    m_proposalCovMatrices.clear();
#endif

    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();

#ifdef UQ_READS_ONLY_EXTRA_STAGES
    m_drScalesForExtraStages.resize(m_drMaxNumExtraStages,1.);
#else
    m_drScalesForExtraStages.resize(m_drMaxNumExtraStages+1,1.);
#endif
#ifdef UQ_USES_TK_CLASS
#else
    m_lowerCholProposalCovMatrices.resize(m_drMaxNumExtraStages+1,NULL);
    m_proposalCovMatrices.resize         (m_drMaxNumExtraStages+1,NULL);
#endif

#ifdef UQ_READS_ONLY_EXTRA_STAGES
    for (unsigned int i = 0; i < m_drMaxNumExtraStages; ++i) {
      if (i < tmpSize) scale = tmpScales[i];
#else
    for (unsigned int i = 1; i < (m_drMaxNumExtraStages+1); ++i) {
      if (i <= tmpSize) scale = tmpScales[i-1];
#endif
      m_drScalesForExtraStages[i] = scale;
    }
    //updateTK();
  }

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
uqMetropolisHastingsSGOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_dataOutputFileName   << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_totallyMute                   << " = " << m_totallyMute
     << "\n" << m_option_rawChain_type                 << " = " << m_rawChainType
     << "\n" << m_option_rawChain_dataInputFileName    << " = " << m_rawChainDataInputFileName
     << "\n" << m_option_rawChain_size                 << " = " << m_rawChainSize
     << "\n" << m_option_rawChain_generateExtra        << " = " << m_rawChainGenerateExtra
     << "\n" << m_option_rawChain_displayPeriod        << " = " << m_rawChainDisplayPeriod
     << "\n" << m_option_rawChain_measureRunTimes      << " = " << m_rawChainMeasureRunTimes
     << "\n" << m_option_rawChain_dataOutputFileName   << " = " << m_rawChainDataOutputFileName
     << "\n" << m_option_rawChain_dataOutputAllowedSet << " = ";
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
     << "\n" << m_option_mh_displayCandidates          << " = " << m_mhDisplayCandidates
     << "\n" << m_option_mh_putOutOfBoundsInChain      << " = " << m_mhPutOutOfBoundsInChain
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

std::ostream& operator<<(std::ostream& os, const uqMetropolisHastingsSGOptionsClass& obj)
{
  obj.print(os);

  return os;
}
