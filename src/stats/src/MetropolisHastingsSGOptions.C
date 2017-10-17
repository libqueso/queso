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

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#define GETPOT_NAMESPACE QUESO // So we don't clash with other getpots
#include <queso/getpot.h>
#undef GETPOT_NAMESPACE
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

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
  m_env(NULL)
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  , m_parser(new BoostInputOptionsParser())
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativeRawSsOptionsValues     ) m_alternativeRawSsOptionsValues      = *alternativeRawSsOptionsValues;
  if (alternativeFilteredSsOptionsValues) m_alternativeFilteredSsOptionsValues = *alternativeFilteredSsOptionsValues;
#endif

  this->set_defaults();
  this->set_prefix("");
}

MhOptionsValues::MhOptionsValues(
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  const SsOptionsValues* alternativeRawSsOptionsValues,
  const SsOptionsValues* alternativeFilteredSsOptionsValues,
#endif
  const BaseEnvironment * env,
  const char * prefix
  )
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativeRawSsOptionsValues     ) m_alternativeRawSsOptionsValues      = *alternativeRawSsOptionsValues;
  if (alternativeFilteredSsOptionsValues) m_alternativeFilteredSsOptionsValues = *alternativeFilteredSsOptionsValues;
#endif

  this->set_defaults();
  this->parse(*env, prefix);
}

// Copy constructor----------------------------------
MhOptionsValues::MhOptionsValues(const MhOptionsValues& src)
{
  this->copy(src);
}

MhOptionsValues::MhOptionsValues(
  const MLSamplingLevelOptions& mlOptions)
  :
  m_prefix                                           (mlOptions.m_prefix),
  m_env                                              (&mlOptions.env()),
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

  m_dataOutputFileName                        = mlOptions.m_dataOutputFileName;
  m_dataOutputAllowAll                        = mlOptions.m_dataOutputAllowAll;
  m_dataOutputAllowedSet                      = mlOptions.m_dataOutputAllowedSet;
  m_totallyMute                               = mlOptions.m_totallyMute;
  m_initialPositionDataInputFileName          = mlOptions.m_initialPositionDataInputFileName;
  m_initialPositionDataInputFileType          = mlOptions.m_initialPositionDataInputFileType;
  m_initialProposalCovMatrixDataInputFileName = mlOptions.m_initialProposalCovMatrixDataInputFileName;
  m_initialProposalCovMatrixDataInputFileType = mlOptions.m_initialProposalCovMatrixDataInputFileType;
  m_parameterDisabledSet                      = mlOptions.m_parameterDisabledSet;
  m_rawChainDataInputFileName                 = mlOptions.m_rawChainDataInputFileName;
  m_rawChainDataInputFileType                 = mlOptions.m_rawChainDataInputFileType;
  m_rawChainSize                              = mlOptions.m_rawChainSize;
  m_rawChainGenerateExtra                     = mlOptions.m_rawChainGenerateExtra;
  m_rawChainDisplayPeriod                     = mlOptions.m_rawChainDisplayPeriod;
  m_rawChainMeasureRunTimes                   = mlOptions.m_rawChainMeasureRunTimes;
  m_rawChainDataOutputPeriod                  = mlOptions.m_rawChainDataOutputPeriod;
  m_rawChainDataOutputFileName                = mlOptions.m_rawChainDataOutputFileName;
  m_rawChainDataOutputFileType                = mlOptions.m_rawChainDataOutputFileType;
  m_rawChainDataOutputAllowAll                = mlOptions.m_rawChainDataOutputAllowAll;
  m_rawChainDataOutputAllowedSet              = mlOptions.m_rawChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChainComputeStats                      = mlOptions.m_rawChainComputeStats;
#endif
  m_filteredChainGenerate                     = mlOptions.m_filteredChainGenerate;
  m_filteredChainDiscardedPortion             = mlOptions.m_filteredChainDiscardedPortion;
  m_filteredChainLag                          = mlOptions.m_filteredChainLag;
  m_filteredChainDataOutputFileName           = mlOptions.m_filteredChainDataOutputFileName;
  m_filteredChainDataOutputFileType           = mlOptions.m_filteredChainDataOutputFileType;
  m_filteredChainDataOutputAllowAll           = mlOptions.m_filteredChainDataOutputAllowAll;
  m_filteredChainDataOutputAllowedSet         = mlOptions.m_filteredChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_filteredChainComputeStats                 = mlOptions.m_filteredChainComputeStats;
#endif
  m_displayCandidates                         = mlOptions.m_displayCandidates;
  m_putOutOfBoundsInChain                     = mlOptions.m_putOutOfBoundsInChain;
  m_tkUseLocalHessian                         = mlOptions.m_tkUseLocalHessian;
  m_tkUseNewtonComponent                      = mlOptions.m_tkUseNewtonComponent;
  m_drMaxNumExtraStages                       = mlOptions.m_drMaxNumExtraStages;
  m_drScalesForExtraStages                    = mlOptions.m_drScalesForExtraStages;
  m_drDuringAmNonAdaptiveInt                  = mlOptions.m_drDuringAmNonAdaptiveInt;
  m_amKeepInitialMatrix                       = mlOptions.m_amKeepInitialMatrix;
  m_amInitialNonAdaptInterval                 = mlOptions.m_amInitialNonAdaptInterval;
  m_amAdaptInterval                           = mlOptions.m_amAdaptInterval;
  m_amAdaptedMatricesDataOutputPeriod         = mlOptions.m_amAdaptedMatricesDataOutputPeriod;
  m_amAdaptedMatricesDataOutputFileName       = mlOptions.m_amAdaptedMatricesDataOutputFileName;
  m_amAdaptedMatricesDataOutputFileType       = mlOptions.m_amAdaptedMatricesDataOutputFileType;
  m_amAdaptedMatricesDataOutputAllowAll       = mlOptions.m_amAdaptedMatricesDataOutputAllowAll;
  m_amAdaptedMatricesDataOutputAllowedSet     = mlOptions.m_amAdaptedMatricesDataOutputAllowedSet;
  m_amEta                                     = mlOptions.m_amEta;
  m_amEpsilon                                 = mlOptions.m_amEpsilon;
  m_enableBrooksGelmanConvMonitor             = UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR;
  m_BrooksGelmanLag                           = UQ_MH_SG_BROOKS_GELMAN_LAG;
  m_outputLogLikelihood                       = UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD;
  m_outputLogTarget                           = UQ_MH_SG_OUTPUT_LOG_TARGET;
  m_doLogitTransform                          = mlOptions.m_doLogitTransform;
  m_algorithm                                 = mlOptions.m_algorithm;
  m_tk                                        = mlOptions.m_tk;
  m_updateInterval                            = mlOptions.m_updateInterval;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
//m_alternativeRawSsOptionsValues             = mlOptions.; // dakota
//m_alternativeFilteredSsOptionsValues        = mlOptions.; // dakota
#endif

  if ((m_env->subDisplayFile() != NULL ) &&
      (m_totallyMute     == false)) {
    *(m_env->subDisplayFile()) << "In MhOptionsValues::constructor(3)"
			       << ": after copying values of options with prefix '" << m_prefix
			       << "', state of object is:"
			       << "\n" << *this
			       << std::endl;
  }
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
MhOptionsValues::checkOptions()
{
  if (m_dataOutputAllowAll) {
    // So, ignore the 'set' option
    m_dataOutputAllowedSet.clear();
    m_dataOutputAllowedSet.insert(m_env->subId());
  }

  if (m_rawChainDataOutputAllowAll) {
    // Again, ignore the set
    m_rawChainDataOutputAllowedSet.clear();
    m_rawChainDataOutputAllowedSet.insert(m_env->subId());
  }

  if (m_filteredChainGenerate == true) {
    queso_require_greater_equal_msg(m_filteredChainLag, 2, "option `" << m_option_filteredChain_lag << "` must be at least 2");
  }

  if (m_filteredChainDataOutputAllowAll) {
    m_filteredChainDataOutputAllowedSet.clear();
    m_filteredChainDataOutputAllowedSet.insert(m_env->subId());
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
    m_amAdaptedMatricesDataOutputAllowedSet.insert(m_env->subId());
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
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

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


void MhOptionsValues::set_prefix(const std::string& prefix)
{
  m_prefix = prefix + "mh_";

  m_option_help = m_prefix + "help";
  m_option_dataOutputFileName = m_prefix + "dataOutputFileName";
  m_option_dataOutputAllowAll = m_prefix + "dataOutputAllowAll";
  m_option_dataOutputAllowedSet = m_prefix + "dataOutputAllowedSet";
  m_option_totallyMute = m_prefix + "totallyMute";
  m_option_initialPosition_dataInputFileName = m_prefix + "initialPosition_dataInputFileName";
  m_option_initialPosition_dataInputFileType = m_prefix + "initialPosition_dataInputFileType";
  m_option_initialProposalCovMatrix_dataInputFileName = m_prefix + "initialProposalCovMatrix_dataInputFileName";
  m_option_initialProposalCovMatrix_dataInputFileType = m_prefix + "initialProposalCovMatrix_dataInputFileType";
  m_option_listOfDisabledParameters = m_prefix + "listOfDisabledParameters";
  m_option_rawChain_dataInputFileName = m_prefix + "rawChain_dataInputFileName";
  m_option_rawChain_dataInputFileType = m_prefix + "rawChain_dataInputFileType";
  m_option_rawChain_size = m_prefix + "rawChain_size";
  m_option_rawChain_generateExtra = m_prefix + "rawChain_generateExtra";
  m_option_rawChain_displayPeriod = m_prefix + "rawChain_displayPeriod";
  m_option_rawChain_measureRunTimes = m_prefix + "rawChain_measureRunTimes";
  m_option_rawChain_dataOutputPeriod = m_prefix + "rawChain_dataOutputPeriod";
  m_option_rawChain_dataOutputFileName = m_prefix + "rawChain_dataOutputFileName";
  m_option_rawChain_dataOutputFileType = m_prefix + "rawChain_dataOutputFileType";
  m_option_rawChain_dataOutputAllowAll = m_prefix + "rawChain_dataOutputAllowAll";
  m_option_rawChain_dataOutputAllowedSet = m_prefix + "rawChain_dataOutputAllowedSet";
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_rawChain_computeStats = m_prefix + "rawChain_computeStats";
#endif
  m_option_filteredChain_generate = m_prefix + "filteredChain_generate";
  m_option_filteredChain_discardedPortion = m_prefix + "filteredChain_discardedPortion";
  m_option_filteredChain_lag = m_prefix + "filteredChain_lag";
  m_option_filteredChain_dataOutputFileName = m_prefix + "filteredChain_dataOutputFileName";
  m_option_filteredChain_dataOutputFileType = m_prefix + "filteredChain_dataOutputFileType";
  m_option_filteredChain_dataOutputAllowAll = m_prefix + "filteredChain_dataOutputAllowAll";
  m_option_filteredChain_dataOutputAllowedSet = m_prefix + "filteredChain_dataOutputAllowedSet";
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_filteredChain_computeStats = m_prefix + "filteredChain_computeStats";
#endif
  m_option_displayCandidates = m_prefix + "displayCandidates";
  m_option_putOutOfBoundsInChain = m_prefix + "putOutOfBoundsInChain";
  m_option_tk_useLocalHessian = m_prefix + "tk_useLocalHessian";
  m_option_tk_useNewtonComponent = m_prefix + "tk_useNewtonComponent";
  m_option_dr_maxNumExtraStages = m_prefix + "dr_maxNumExtraStages";
  m_option_dr_listOfScalesForExtraStages = m_prefix + "dr_listOfScalesForExtraStages";
  m_option_dr_duringAmNonAdaptiveInt = m_prefix + "dr_duringAmNonAdaptiveInt";
  m_option_am_keepInitialMatrix = m_prefix + "am_keepInitialMatrix";
  m_option_am_initialNonAdaptInterval = m_prefix + "am_initialNonAdaptInterval";
  m_option_am_adaptInterval = m_prefix + "am_adaptInterval";
  m_option_am_adaptedMatrices_dataOutputPeriod = m_prefix + "am_adaptedMatrices_dataOutputPeriod";
  m_option_am_adaptedMatrices_dataOutputFileName = m_prefix + "am_adaptedMatrices_dataOutputFileName";
  m_option_am_adaptedMatrices_dataOutputFileType = m_prefix + "am_adaptedMatrices_dataOutputFileType";
  m_option_am_adaptedMatrices_dataOutputAllowAll = m_prefix + "am_adaptedMatrices_dataOutputAllowAll";
  m_option_am_adaptedMatrices_dataOutputAllowedSet = m_prefix + "am_adaptedMatrices_dataOutputAllowedSet";
  m_option_am_eta = m_prefix + "am_eta";
  m_option_am_epsilon = m_prefix + "am_epsilon";
  m_option_enableBrooksGelmanConvMonitor = m_prefix + "enableBrooksGelmanConvMonitor";
  m_option_BrooksGelmanLag = m_prefix + "BrooksGelmanLag";
  m_option_outputLogLikelihood = m_prefix + "outputLogLikelihood";
  m_option_outputLogTarget = m_prefix + "outputLogTarget";
  m_option_doLogitTransform = m_prefix + "doLogitTransform";
  m_option_algorithm = m_prefix + "algorithm";
  m_option_tk = m_prefix + "tk";
  m_option_updateInterval = m_prefix + "updateInterval";
}


void MhOptionsValues::set_defaults()
{
    m_help = UQ_MH_SG_HELP;
    m_dataOutputFileName = UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV;
    m_dataOutputAllowAll = UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV;
  //m_dataOutputAllowedSet                     (),
    m_totallyMute = UQ_MH_SG_TOTALLY_MUTE_ODV;
    m_initialPositionDataInputFileName = UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV;
    m_initialPositionDataInputFileType = UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV;
    m_initialProposalCovMatrixDataInputFileName = UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV;
    m_initialProposalCovMatrixDataInputFileType = UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV;
  //m_parameterDisabledSet                     (),
    m_rawChainDataInputFileName = UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV;
    m_rawChainDataInputFileType = UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV;
    m_rawChainSize = UQ_MH_SG_RAW_CHAIN_SIZE_ODV;
    m_rawChainGenerateExtra = UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV;
    m_rawChainDisplayPeriod = UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV;
    m_rawChainMeasureRunTimes = UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV;
    m_rawChainDataOutputPeriod = UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV;
    m_rawChainDataOutputFileName = UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV;
    m_rawChainDataOutputFileType = UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV;
    m_rawChainDataOutputAllowAll = UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV;
  //m_rawChainDataOutputAllowedSet             (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_rawChainComputeStats = UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV;
#endif
    m_filteredChainGenerate = UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV;
    m_filteredChainDiscardedPortion = UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV;
    m_filteredChainLag = UQ_MH_SG_FILTERED_CHAIN_LAG_ODV;
    m_filteredChainDataOutputFileName = UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV;
    m_filteredChainDataOutputFileType = UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV;
    m_filteredChainDataOutputAllowAll = UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV;
  //m_filteredChainDataOutputAllowedSet        (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_filteredChainComputeStats = UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV;
#endif
    m_displayCandidates = UQ_MH_SG_DISPLAY_CANDIDATES_ODV;
    m_putOutOfBoundsInChain = UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV;
    m_tkUseLocalHessian = UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV;
    m_tkUseNewtonComponent = UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV;
    m_drMaxNumExtraStages = UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV;
    m_drScalesForExtraStages.resize(0);
    m_drDuringAmNonAdaptiveInt = UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV;
    m_amKeepInitialMatrix = UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV;
    m_amInitialNonAdaptInterval = UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV;
    m_amAdaptInterval = UQ_MH_SG_AM_ADAPT_INTERVAL_ODV;
    m_amAdaptedMatricesDataOutputPeriod = UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV;
    m_amAdaptedMatricesDataOutputFileName = UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV;
    m_amAdaptedMatricesDataOutputFileType = UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV;
    m_amAdaptedMatricesDataOutputAllowAll = UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV;
  //m_amAdaptedMatricesDataOutputAllowedSet    (),
    m_amEta = UQ_MH_SG_AM_ETA_ODV;
    m_amEpsilon = UQ_MH_SG_AM_EPSILON_ODV;
    m_enableBrooksGelmanConvMonitor = UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR;
    m_BrooksGelmanLag = UQ_MH_SG_BROOKS_GELMAN_LAG;
    m_outputLogLikelihood = UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD;
    m_outputLogTarget = UQ_MH_SG_OUTPUT_LOG_TARGET;
    m_doLogitTransform = UQ_MH_SG_DO_LOGIT_TRANSFORM;
    m_algorithm = UQ_MH_SG_ALGORITHM;
    m_tk = UQ_MH_SG_TK;
    m_updateInterval = UQ_MH_SG_UPDATE_INTERVAL;
}

void
MhOptionsValues::parse(const BaseEnvironment& env, const std::string& prefix)
{
  m_env = &env;

  this->set_prefix(prefix);

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser.reset(new BoostInputOptionsParser(env.optionsInputFileName()));

  m_parser->registerOption<std::string >(m_option_help,                                       m_help,                                       "produce help msg for Bayesian Metropolis-Hastings"          );
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,                         m_dataOutputFileName,                         "name of generic output file"                                );
  m_parser->registerOption<bool        >(m_option_dataOutputAllowAll,                         m_dataOutputAllowAll,                         "allow all subEnvs write to a generic output file"           );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,                       container_to_string(m_dataOutputAllowedSet),  "subEnvs that will write to generic output file"             );
  m_parser->registerOption<bool        >(m_option_totallyMute,                                m_totallyMute,                                "totally mute (no printout msg)"                             );
  m_parser->registerOption<std::string >(m_option_initialPosition_dataInputFileName,          m_initialPositionDataInputFileName,          "name of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_initialPosition_dataInputFileType,          m_initialPositionDataInputFileType,          "type of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName, "name of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType, "type of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_listOfDisabledParameters,                   container_to_string(m_parameterDisabledSet), "list of disabled parameters"                                );
  m_parser->registerOption<std::string >(m_option_rawChain_dataInputFileName,                 m_rawChainDataInputFileName,                 "name of input file for raw chain "                          );
  m_parser->registerOption<std::string >(m_option_rawChain_dataInputFileType,                 m_rawChainDataInputFileType,                 "type of input file for raw chain "                          );
  m_parser->registerOption<unsigned int>(m_option_rawChain_size,                              m_rawChainSize,                              "size of raw chain"                                          );
  m_parser->registerOption<bool        >(m_option_rawChain_generateExtra,                     m_rawChainGenerateExtra,                     "generate extra information about raw chain"                 );
  m_parser->registerOption<unsigned int>(m_option_rawChain_displayPeriod,                     m_rawChainDisplayPeriod,                     "period of msg display during raw chain generation"          );
  m_parser->registerOption<bool        >(m_option_rawChain_measureRunTimes,                   m_rawChainMeasureRunTimes,                   "measure run times"                                          );
  m_parser->registerOption<unsigned int>(m_option_rawChain_dataOutputPeriod,                  m_rawChainDataOutputPeriod,                  "period of msg display during raw chain generation"          );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputFileName,                m_rawChainDataOutputFileName,                "name of output file for raw chain "                         );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputFileType,                m_rawChainDataOutputFileType,                "type of output file for raw chain "                         );
  m_parser->registerOption<bool        >(m_option_rawChain_dataOutputAllowAll,                m_rawChainDataOutputAllowAll,                "allow all subEnvs to write raw chain to an output file"     );
  m_parser->registerOption<std::string >(m_option_rawChain_dataOutputAllowedSet,              container_to_string(m_rawChainDataOutputAllowedSet), "subEnvs that will write raw chain to output file"           );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_rawChain_computeStats,                      m_rawChain_computeStats,                      "compute statistics on raw chain"                            );
#endif
  m_parser->registerOption<bool        >(m_option_filteredChain_generate,                     m_filteredChainGenerate,                     "generate filtered chain"                                    );
  m_parser->registerOption<double      >(m_option_filteredChain_discardedPortion,             m_filteredChainDiscardedPortion,             "initial discarded portion for chain filtering"              );
  m_parser->registerOption<unsigned int>(m_option_filteredChain_lag,                          m_filteredChainLag,                          "spacing for chain filtering"                                );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputFileName,           m_filteredChainDataOutputFileName,           "name of output file for filtered chain"                     );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputFileType,           m_filteredChainDataOutputFileType,           "type of output file for filtered chain"                     );
  m_parser->registerOption<bool        >(m_option_filteredChain_dataOutputAllowAll,           m_filteredChainDataOutputAllowAll,           "allow all subEnvs to write filt chain to an output file"    );
  m_parser->registerOption<std::string >(m_option_filteredChain_dataOutputAllowedSet,         container_to_string(m_filteredChainDataOutputAllowedSet), "subEnvs that will write filt chain to output file"          );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_filteredChain_computeStats,                 m_filteredChainComputeStats,                 "compute statistics on filtered chain"                       );
#endif
  m_parser->registerOption<bool        >(m_option_displayCandidates,                          m_displayCandidates,                          "display candidates in the core MH algorithm"                );
  m_parser->registerOption<bool        >(m_option_putOutOfBoundsInChain,                      m_putOutOfBoundsInChain,                      "put 'out of bound' candidates in chain as well"             );
  m_parser->registerOption<bool        >(m_option_tk_useLocalHessian,                         m_tkUseLocalHessian,                         "'proposal' use local Hessian"                               );
  m_parser->registerOption<bool        >(m_option_tk_useNewtonComponent,                      m_tkUseNewtonComponent,                      "'proposal' use Newton component"                            );
  m_parser->registerOption<unsigned int>(m_option_dr_maxNumExtraStages,                       m_drMaxNumExtraStages,                       "'dr' maximum number of extra stages"                        );
  m_parser->registerOption<std::string >(m_option_dr_listOfScalesForExtraStages,              container_to_string(m_drScalesForExtraStages), "'dr' scales for prop cov matrices from 2nd stage on"        );
  m_parser->registerOption<bool        >(m_option_dr_duringAmNonAdaptiveInt,                  m_drDuringAmNonAdaptiveInt,                  "'dr' used during 'am' non adaptive interval"                );
  m_parser->registerOption<bool        >(m_option_am_keepInitialMatrix,                       m_amKeepInitialMatrix,                       "'am' keep initial (given) matrix"                           );
  m_parser->registerOption<unsigned int>(m_option_am_initialNonAdaptInterval,                 m_amInitialNonAdaptInterval,                 "'am' initial non adaptation interval"                       );
  m_parser->registerOption<unsigned int>(m_option_am_adaptInterval,                           m_amAdaptInterval,                           "'am' adaptation interval"                                   );
  m_parser->registerOption<unsigned int>(m_option_am_adaptedMatrices_dataOutputPeriod,        m_amAdaptedMatricesDataOutputPeriod,        "period for outputting 'am' adapted matrices"                );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileName,      m_amAdaptedMatricesDataOutputFileName,      "name of output file for 'am' adapted matrices"              );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputFileType,      m_amAdaptedMatricesDataOutputFileType,      "type of output file for 'am' adapted matrices"              );
  m_parser->registerOption<bool        >(m_option_am_adaptedMatrices_dataOutputAllowAll,      m_amAdaptedMatricesDataOutputAllowAll,      "type of output file for 'am' adapted matrices"              );
  m_parser->registerOption<std::string >(m_option_am_adaptedMatrices_dataOutputAllowedSet,    container_to_string(m_amAdaptedMatricesDataOutputAllowedSet), "type of output file for 'am' adapted matrices"              );
  m_parser->registerOption<double      >(m_option_am_eta,                                     m_amEta,                                     "'am' eta"                                                   );
  m_parser->registerOption<double      >(m_option_am_epsilon,                                 m_amEpsilon,                                 "'am' epsilon"                                               );
  m_parser->registerOption<unsigned int>(m_option_enableBrooksGelmanConvMonitor,              m_enableBrooksGelmanConvMonitor,              "assess convergence using Brooks-Gelman metric"              );
  m_parser->registerOption<unsigned int>(m_option_BrooksGelmanLag,                            m_BrooksGelmanLag,                            "number of chain positions before starting to compute metric");
  m_parser->registerOption<bool        >(m_option_outputLogLikelihood,                        m_outputLogLikelihood,                        "flag to toggle output of log likelihood values"             );
  m_parser->registerOption<bool        >(m_option_outputLogTarget,                            m_outputLogTarget,                            "flag to toggle output of log target values"                 );
  m_parser->registerOption<bool        >(m_option_doLogitTransform,                           m_doLogitTransform,                           "flag to toggle logit transform for bounded domains"         );
  m_parser->registerOption<std::string >(m_option_algorithm,                                  m_algorithm,                                  "which MCMC algorithm to use"                                );
  m_parser->registerOption<std::string >(m_option_tk,                                         m_tk,                                         "which MCMC transition kernel to use"                        );
  m_parser->registerOption<unsigned int>(m_option_updateInterval,                             m_updateInterval,                             "how often to call updateTK method"                          );

  m_parser->scanInputFile();

  m_parser->getOption<std::string >(m_option_help,                                       m_help);
  m_parser->getOption<std::string >(m_option_dataOutputFileName,                         m_dataOutputFileName);
  m_parser->getOption<bool        >(m_option_dataOutputAllowAll,                         m_dataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,            m_dataOutputAllowedSet);
  m_parser->getOption<bool        >(m_option_totallyMute,                                m_totallyMute);
  m_parser->getOption<std::string >(m_option_initialPosition_dataInputFileName,          m_initialPositionDataInputFileName);
  m_parser->getOption<std::string >(m_option_initialPosition_dataInputFileType,          m_initialPositionDataInputFileType);
  m_parser->getOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName);
  m_parser->getOption<std::string >(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType);
  m_parser->getOption<std::set<unsigned int> >(m_option_listOfDisabledParameters,        m_parameterDisabledSet);
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
  m_parser->getOption<std::set<unsigned int> >(m_option_rawChain_dataOutputAllowedSet,   m_rawChainDataOutputAllowedSet);
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
  m_parser->getOption<std::vector<double> >(m_option_dr_listOfScalesForExtraStages,      m_drScalesForExtraStages);
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
  m_help = m_env->input()(m_option_help, m_help);
  m_dataOutputFileName = m_env->input()(m_option_dataOutputFileName, m_dataOutputFileName);
  m_dataOutputAllowAll = m_env->input()(m_option_dataOutputAllowAll, m_dataOutputAllowAll);

  // UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = m_env->input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

  m_totallyMute = m_env->input()(m_option_totallyMute, m_totallyMute);
  m_initialPositionDataInputFileName = m_env->input()(m_option_initialPosition_dataInputFileName, m_initialPositionDataInputFileName);
  m_initialPositionDataInputFileType = m_env->input()(m_option_initialPosition_dataInputFileType, m_initialPositionDataInputFileType);
  m_initialProposalCovMatrixDataInputFileName = m_env->input()(m_option_initialProposalCovMatrix_dataInputFileName, m_initialProposalCovMatrixDataInputFileName);
  m_initialProposalCovMatrixDataInputFileType = m_env->input()(m_option_initialProposalCovMatrix_dataInputFileType, m_initialProposalCovMatrixDataInputFileType);

  // UQ_MH_SG_LIST_OF_DISABLED_PARAMETERS_ODV is the empty set (string) by
  // default
  size = m_env->input().vector_variable_size(m_option_listOfDisabledParameters);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int disabled = m_env->input()(m_option_listOfDisabledParameters, i, i);
    m_parameterDisabledSet.insert(disabled);
  }

  m_rawChainDataInputFileName = m_env->input()(m_option_rawChain_dataInputFileName, m_rawChainDataInputFileName);
  m_rawChainDataInputFileType = m_env->input()(m_option_rawChain_dataInputFileType, m_rawChainDataInputFileType);
  m_rawChainSize = m_env->input()(m_option_rawChain_size, m_rawChainSize);
  m_rawChainGenerateExtra = m_env->input()(m_option_rawChain_generateExtra, m_rawChainGenerateExtra);
  m_rawChainDisplayPeriod = m_env->input()(m_option_rawChain_displayPeriod, m_rawChainDisplayPeriod);
  m_rawChainMeasureRunTimes = m_env->input()(m_option_rawChain_measureRunTimes, m_rawChainMeasureRunTimes);
  m_rawChainDataOutputPeriod = m_env->input()(m_option_rawChain_dataOutputPeriod, m_rawChainDataOutputPeriod);
  m_rawChainDataOutputFileName = m_env->input()(m_option_rawChain_dataOutputFileName, m_rawChainDataOutputFileName);
  m_rawChainDataOutputFileType = m_env->input()(m_option_rawChain_dataOutputFileType, m_rawChainDataOutputFileType);
  m_rawChainDataOutputAllowAll = m_env->input()(m_option_rawChain_dataOutputAllowAll, m_rawChainDataOutputAllowAll);

  // UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_rawChain_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_rawChain_dataOutputAllowedSet, i, i);
    m_rawChainDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_rawChain_computeStats = m_env->input()(m_option_rawChain_computeStats, m_rawChainComputeStats);
#endif
  m_filteredChainGenerate = m_env->input()(m_option_filteredChain_generate, m_filteredChainGenerate);
  m_filteredChainDiscardedPortion = m_env->input()(m_option_filteredChain_discardedPortion, m_filteredChainDiscardedPortion);
  m_filteredChainLag = m_env->input()(m_option_filteredChain_lag, m_filteredChainLag);
  m_filteredChainDataOutputFileName = m_env->input()(m_option_filteredChain_dataOutputFileName, m_filteredChainDataOutputFileName);
  m_filteredChainDataOutputFileType = m_env->input()(m_option_filteredChain_dataOutputFileType, m_filteredChainDataOutputFileType);
  m_filteredChainDataOutputAllowAll = m_env->input()(m_option_filteredChain_dataOutputAllowAll, m_filteredChainDataOutputAllowAll);

  // UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_filteredChain_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_filteredChain_dataOutputAllowedSet, i, i);
    m_filteredChainDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_filteredChain_computeStats = m_env->input()(m_option_filteredChain_computeStats, m_filteredChainComputeStats);
#endif
  m_displayCandidates = m_env->input()(m_option_displayCandidates, m_displayCandidates);
  m_putOutOfBoundsInChain = m_env->input()(m_option_putOutOfBoundsInChain, m_putOutOfBoundsInChain);
  m_tkUseLocalHessian = m_env->input()(m_option_tk_useLocalHessian, m_tkUseLocalHessian);
  m_tkUseNewtonComponent = m_env->input()(m_option_tk_useNewtonComponent, m_tkUseNewtonComponent);
  m_drMaxNumExtraStages = m_env->input()(m_option_dr_maxNumExtraStages, m_drMaxNumExtraStages);

  // UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_dr_listOfScalesForExtraStages);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_dr_listOfScalesForExtraStages, i, i);
    m_drScalesForExtraStages.push_back(allowed);
  }

  m_drDuringAmNonAdaptiveInt = m_env->input()(m_option_dr_duringAmNonAdaptiveInt, m_drDuringAmNonAdaptiveInt);
  m_amKeepInitialMatrix = m_env->input()(m_option_am_keepInitialMatrix, m_amKeepInitialMatrix);
  m_amInitialNonAdaptInterval = m_env->input()(m_option_am_initialNonAdaptInterval, m_amInitialNonAdaptInterval);
  m_amAdaptInterval = m_env->input()(m_option_am_adaptInterval, m_amAdaptInterval);
  m_amAdaptedMatricesDataOutputPeriod = m_env->input()(m_option_am_adaptedMatrices_dataOutputPeriod, m_amAdaptedMatricesDataOutputPeriod);
  m_amAdaptedMatricesDataOutputFileName = m_env->input()(m_option_am_adaptedMatrices_dataOutputFileName, m_amAdaptedMatricesDataOutputFileName);
  m_amAdaptedMatricesDataOutputFileType = m_env->input()(m_option_am_adaptedMatrices_dataOutputFileType, m_amAdaptedMatricesDataOutputFileType);
  m_amAdaptedMatricesDataOutputAllowAll = m_env->input()(m_option_am_adaptedMatrices_dataOutputAllowAll, m_amAdaptedMatricesDataOutputAllowAll);

  // UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = m_env->input().vector_variable_size(m_option_am_adaptedMatrices_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never used
    // here
    unsigned int allowed = m_env->input()(m_option_am_adaptedMatrices_dataOutputAllowedSet, i, i);
    m_amAdaptedMatricesDataOutputAllowedSet.insert(allowed);
  }

  m_amEta = m_env->input()(m_option_am_eta, m_amEta);
  m_amEpsilon = m_env->input()(m_option_am_epsilon, m_amEpsilon);
  m_enableBrooksGelmanConvMonitor = m_env->input()(m_option_enableBrooksGelmanConvMonitor, m_enableBrooksGelmanConvMonitor);
  m_BrooksGelmanLag = m_env->input()(m_option_BrooksGelmanLag, m_BrooksGelmanLag);
  m_outputLogLikelihood = m_env->input()(m_option_outputLogLikelihood, m_outputLogLikelihood);
  m_outputLogTarget = m_env->input()(m_option_outputLogTarget, m_outputLogTarget);
  m_doLogitTransform = m_env->input()(m_option_doLogitTransform, m_doLogitTransform);
  m_algorithm = m_env->input()(m_option_algorithm, m_algorithm);
  m_tk = m_env->input()(m_option_tk, m_tk);
  m_updateInterval = m_env->input()(m_option_updateInterval, m_updateInterval);
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();

}

}  // End namespace QUESO
