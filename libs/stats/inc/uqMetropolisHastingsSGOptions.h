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

#ifndef __UQ_MH_SG_OPTIONS_H__
#define __UQ_MH_SG_OPTIONS_H__

#include <uqEnvironment.h>
#include <uqMLSamplingLevelOptions.h>
#include <uqSequenceStatisticalOptions.h>

#undef  UQ_MH_SG_REQUIRES_INVERTED_COV_MATRICES
#define UQ_NOTHING_JUST_FOR_TEST_OF_SVN_ID 1

#define UQ_MH_SG_FILENAME_FOR_NO_FILE   "."

// _ODV = option default value
#define UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                  UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV                ""

#define UQ_MH_SG_TOTALLY_MUTE_ODV                           0
#define UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV         UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_RAW_CHAIN_SIZE_ODV                         100
#define UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV               0
#define UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV               500
#define UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV            0
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV        UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV      ""
#define UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV                0
#define UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV                0
#define UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV       0.
#define UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                     1
#define UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV   UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV ""
#define UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV           0
#define UQ_MH_SG_DISPLAY_CANDIDATES_ODV                     0
#define UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV             1
#define UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                   0
#define UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV                1
#define UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV                0
#define UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV     ""
#define UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                  0
#define UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                      0
#define UQ_MH_SG_AM_ETA_ODV                                 1.
#define UQ_MH_SG_AM_EPSILON_ODV                             1.e-5

class uqMetropolisHastingsSGOptionsClass
{
public:
  uqMetropolisHastingsSGOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
  uqMetropolisHastingsSGOptionsClass(const uqMLSamplingLevelOptionsClass& inputOptions);
 ~uqMetropolisHastingsSGOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  std::string                        m_prefix;

  std::string                        m_dataOutputFileName;
  std::set<unsigned int>             m_dataOutputAllowedSet;

  bool                               m_totallyMute;
  std::string                        m_rawChainDataInputFileName;
  unsigned int                       m_rawChainSize;
  bool                               m_rawChainGenerateExtra;
  unsigned int                       m_rawChainDisplayPeriod;
  bool                               m_rawChainMeasureRunTimes;
  std::string                        m_rawChainDataOutputFileName;
  std::set<unsigned int>             m_rawChainDataOutputAllowedSet;
  bool                               m_rawChainComputeStats;
  uqSequenceStatisticalOptionsClass* m_rawChainStatisticalOptions;
  bool                               m_rawChainStatOptsInstantiated;

  bool                               m_filteredChainGenerate;
  double                             m_filteredChainDiscardedPortion; // input or set during run time
  unsigned int                       m_filteredChainLag;              // input or set during run time
  std::string                        m_filteredChainDataOutputFileName;
  std::set<unsigned int>             m_filteredChainDataOutputAllowedSet;
  bool                               m_filteredChainComputeStats;
  uqSequenceStatisticalOptionsClass* m_filteredChainStatisticalOptions;
  bool                               m_filteredChainStatOptsInstantiated;

  bool                               m_displayCandidates;
  bool                               m_putOutOfBoundsInChain;
  bool                               m_tkUseLocalHessian;
  bool                               m_tkUseNewtonComponent;
  unsigned int                       m_drMaxNumExtraStages;
  std::vector<double>                m_drScalesForExtraStages;
  unsigned int                       m_amInitialNonAdaptInterval;
  unsigned int                       m_amAdaptInterval;
  double                             m_amEta;
  double                             m_amEpsilon;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;

  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;

  std::string                   m_option_totallyMute;
  std::string                   m_option_rawChain_dataInputFileName;
  std::string                   m_option_rawChain_size;
  std::string                   m_option_rawChain_generateExtra;
  std::string                   m_option_rawChain_displayPeriod;
  std::string                   m_option_rawChain_measureRunTimes;
  std::string                   m_option_rawChain_dataOutputFileName;
  std::string                   m_option_rawChain_dataOutputAllowedSet;
  std::string                   m_option_rawChain_computeStats;
  std::string                   m_option_filteredChain_generate;
  std::string                   m_option_filteredChain_discardedPortion;
  std::string                   m_option_filteredChain_lag;
  std::string                   m_option_filteredChain_dataOutputFileName;
  std::string                   m_option_filteredChain_dataOutputAllowedSet;
  std::string                   m_option_filteredChain_computeStats;
  std::string                   m_option_displayCandidates;
  std::string                   m_option_putOutOfBoundsInChain;
  std::string                   m_option_tk_useLocalHessian;
  std::string                   m_option_tk_useNewtonComponent;
  std::string                   m_option_dr_maxNumExtraStages;
  std::string                   m_option_dr_listOfScalesForExtraStages;
  std::string                   m_option_am_initialNonAdaptInterval;
  std::string                   m_option_am_adaptInterval;
  std::string                   m_option_am_eta;
  std::string                   m_option_am_epsilon;
};

std::ostream& operator<<(std::ostream& os, const uqMetropolisHastingsSGOptionsClass& obj);
#endif // __UQ_MH_SG_OPTIONS_H__
