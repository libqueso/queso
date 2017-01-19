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

#include <queso/Defines.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/MonteCarloSGOptions.h>
#include <queso/Miscellaneous.h>

// -------------------------------------------------
// McOptionsValues --------------------------
// -------------------------------------------------

namespace QUESO {

// Default constructor -----------------------------
McOptionsValues::McOptionsValues(
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  const SsOptionsValues* alternativePSsOptionsValues,
  const SsOptionsValues* alternativeQSsOptionsValues
#endif
  )
  :
    m_prefix                          ("mc_"),
    m_help                       (UQ_MOC_SG_HELP),
    m_dataOutputFileName         (UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV     ),
  //m_dataOutputAllowedSet       (),
    m_pseqDataOutputPeriod       (UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV   ),
    m_pseqDataOutputFileName     (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV),
    m_pseqDataOutputFileType     (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV),
  //m_pseqDataOutputAllowedSet   (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_pseqComputeStats           (UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV        ),
#endif
    m_qseqDataInputFileName      (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV ),
    m_qseqDataInputFileType      (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV ),
    m_qseqSize                   (UQ_MOC_SG_QSEQ_SIZE_ODV                 ),
    m_qseqDisplayPeriod          (UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV       ),
    m_qseqMeasureRunTimes        (UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV    ),
    m_qseqDataOutputPeriod       (UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV   ),
    m_qseqDataOutputFileName     (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV),
    m_qseqDataOutputFileType     (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV),
  //m_qseqDataOutputAllowedSet   (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_qseqComputeStats           (UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV        ),
    m_alternativePSsOptionsValues(),
    m_alternativeQSsOptionsValues(),
#endif
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(NULL),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help                     (m_prefix + "help"                       ),
    m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"         ),
    m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"       ),
    m_option_pseq_dataOutputPeriod    (m_prefix + "pseq_dataOutputPeriod"      ),
    m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"    ),
    m_option_pseq_dataOutputFileType  (m_prefix + "pseq_dataOutputFileType"    ),
    m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"  ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"          ),
#endif
    m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"     ),
    m_option_qseq_dataInputFileType   (m_prefix + "qseq_dataInputFileType"     ),
    m_option_qseq_size                (m_prefix + "qseq_size"                  ),
    m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"         ),
    m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"       ),
    m_option_qseq_dataOutputPeriod    (m_prefix + "qseq_dataOutputPeriod"      ),
    m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"    ),
    m_option_qseq_dataOutputFileType  (m_prefix + "qseq_dataOutputFileType"    ),
    m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet"  )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    ,
    m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"          )
#endif
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativePSsOptionsValues) m_alternativePSsOptionsValues = *alternativePSsOptionsValues;
  if (alternativeQSsOptionsValues) m_alternativeQSsOptionsValues = *alternativeQSsOptionsValues;
#endif
}

McOptionsValues::McOptionsValues(
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  const SsOptionsValues* alternativePSsOptionsValues,
  const SsOptionsValues* alternativeQSsOptionsValues,
#endif
  const BaseEnvironment * env, const char * prefix
  )
  :
    m_prefix                          ((std::string)(prefix) + "mc_"),
    m_help                       (UQ_MOC_SG_HELP),
    m_dataOutputFileName         (UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV     ),
  //m_dataOutputAllowedSet       (),
    m_pseqDataOutputPeriod       (UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV   ),
    m_pseqDataOutputFileName     (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV),
    m_pseqDataOutputFileType     (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV),
  //m_pseqDataOutputAllowedSet   (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_pseqComputeStats           (UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV        ),
#endif
    m_qseqDataInputFileName      (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV ),
    m_qseqDataInputFileType      (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV ),
    m_qseqSize                   (UQ_MOC_SG_QSEQ_SIZE_ODV                 ),
    m_qseqDisplayPeriod          (UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV       ),
    m_qseqMeasureRunTimes        (UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV    ),
    m_qseqDataOutputPeriod       (UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV   ),
    m_qseqDataOutputFileName     (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV),
    m_qseqDataOutputFileType     (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV),
  //m_qseqDataOutputAllowedSet   (),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_qseqComputeStats           (UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV        ),
    m_alternativePSsOptionsValues(),
    m_alternativeQSsOptionsValues(),
#endif
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
    m_option_help                     (m_prefix + "help"                       ),
    m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"         ),
    m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"       ),
    m_option_pseq_dataOutputPeriod    (m_prefix + "pseq_dataOutputPeriod"      ),
    m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"    ),
    m_option_pseq_dataOutputFileType  (m_prefix + "pseq_dataOutputFileType"    ),
    m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"  ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"          ),
#endif
    m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"     ),
    m_option_qseq_dataInputFileType   (m_prefix + "qseq_dataInputFileType"     ),
    m_option_qseq_size                (m_prefix + "qseq_size"                  ),
    m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"         ),
    m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"       ),
    m_option_qseq_dataOutputPeriod    (m_prefix + "qseq_dataOutputPeriod"      ),
    m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"    ),
    m_option_qseq_dataOutputFileType  (m_prefix + "qseq_dataOutputFileType"    ),
    m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet"  )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    ,
    m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"          )
#endif
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativePSsOptionsValues) m_alternativePSsOptionsValues = *alternativePSsOptionsValues;
  if (alternativeQSsOptionsValues) m_alternativeQSsOptionsValues = *alternativeQSsOptionsValues;
#endif

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string >(m_option_help,                      UQ_MOC_SG_HELP                            , "produce help message for Monte Carlo distribution calculator");
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,        UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV       , "name of generic data output file"                            );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,      UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV     , "subEnvs that will write to generic data output file"         );
  m_parser->registerOption<unsigned int>(m_option_pseq_dataOutputPeriod,     UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV     , "period of message display during param sequence generation"  );
  m_parser->registerOption<std::string >(m_option_pseq_dataOutputFileName,   UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV  , "name of data output file for parameters"                     );
  m_parser->registerOption<std::string >(m_option_pseq_dataOutputFileType,   UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV  , "type of data output file for parameters"                     );
  m_parser->registerOption<std::string >(m_option_pseq_dataOutputAllowedSet, UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV, "subEnvs that will write to data output file for parameters"  );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_pseq_computeStats,         UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV          , "compute statistics on sequence of parameter"                 );
#endif
  m_parser->registerOption<std::string >(m_option_qseq_dataInputFileName,    UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV   , "name of data input file for qois"                            );
  m_parser->registerOption<std::string >(m_option_qseq_dataInputFileType,    UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV   , "type of data input file for qois"                            );
  m_parser->registerOption<unsigned int>(m_option_qseq_size,                 UQ_MOC_SG_QSEQ_SIZE_ODV                   , "size of qoi sequence"                                        );
  m_parser->registerOption<unsigned int>(m_option_qseq_displayPeriod,        UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV         , "period of message display during qoi sequence generation"    );
  m_parser->registerOption<bool        >(m_option_qseq_measureRunTimes,      UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV      , "measure run times"                                           );
  m_parser->registerOption<unsigned int>(m_option_qseq_dataOutputPeriod,     UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV     , "period of message display during qoi sequence generation"    );
  m_parser->registerOption<std::string >(m_option_qseq_dataOutputFileName,   UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV  , "name of data output file for qois"                           );
  m_parser->registerOption<std::string >(m_option_qseq_dataOutputFileType,   UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV  , "type of data output file for qois"                           );
  m_parser->registerOption<std::string >(m_option_qseq_dataOutputAllowedSet, UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV, "subEnvs that will write to data output file for qois"        );
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool        >(m_option_qseq_computeStats,         UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV          , "compute statistics on sequence of qoi"                       );
#endif

  m_parser->scanInputFile();

  m_parser->getOption<std::string >(m_option_help,        m_help);
  m_parser->getOption<std::string >(m_option_dataOutputFileName,        m_dataOutputFileName);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,      m_dataOutputAllowedSet);
  m_parser->getOption<unsigned int>(m_option_pseq_dataOutputPeriod,     m_pseqDataOutputPeriod);
  m_parser->getOption<std::string >(m_option_pseq_dataOutputFileName,   m_pseqDataOutputFileName);
  m_parser->getOption<std::string >(m_option_pseq_dataOutputFileType,   m_pseqDataOutputFileType);
  m_parser->getOption<std::set<unsigned int> >(m_option_pseq_dataOutputAllowedSet, m_pseqDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool        >(m_option_pseq_computeStats,         m_pseq_computeStats);
#endif
  m_parser->getOption<std::string >(m_option_qseq_dataInputFileName,    m_qseqDataInputFileName);
  m_parser->getOption<std::string >(m_option_qseq_dataInputFileType,    m_qseqDataInputFileType);
  m_parser->getOption<unsigned int>(m_option_qseq_size,                 m_qseqSize);
  m_parser->getOption<unsigned int>(m_option_qseq_displayPeriod,        m_qseqDisplayPeriod);
  m_parser->getOption<bool        >(m_option_qseq_measureRunTimes,      m_qseqMeasureRunTimes);
  m_parser->getOption<unsigned int>(m_option_qseq_dataOutputPeriod,     m_qseqDataOutputPeriod);
  m_parser->getOption<std::string >(m_option_qseq_dataOutputFileName,   m_qseqDataOutputFileName);
  m_parser->getOption<std::string >(m_option_qseq_dataOutputFileType,   m_qseqDataOutputFileType);
  m_parser->getOption<std::set<unsigned int> >(m_option_qseq_dataOutputAllowedSet, m_qseqDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool        >(m_option_qseq_computeStats,         m_qseq_computeStats);
#endif
#else
  m_help = env->input()(m_option_help, UQ_MOC_SG_HELP);
  m_dataOutputFileName = env->input()(m_option_dataOutputFileName, UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV);

  // UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = env->input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env->input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

  m_pseqDataOutputPeriod = env->input()(m_option_pseq_dataOutputPeriod, UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV);
  m_pseqDataOutputFileName = env->input()(m_option_pseq_dataOutputFileName, UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV);
  m_pseqDataOutputFileType = env->input()(m_option_pseq_dataOutputFileType, UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV);

  // UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = env->input().vector_variable_size(m_option_pseq_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env->input()(m_option_pseq_dataOutputAllowedSet, i, i);
    m_pseqDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_pseq_computeStats = env->input()(m_option_pseq_computeStats, UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV);
#endif
  m_qseqDataInputFileName = env->input()(m_option_qseq_dataInputFileName, UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV);
  m_qseqDataInputFileType = env->input()(m_option_qseq_dataInputFileType, UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV);
  m_qseqSize = env->input()(m_option_qseq_size, UQ_MOC_SG_QSEQ_SIZE_ODV);
  m_qseqDisplayPeriod = env->input()(m_option_qseq_displayPeriod, UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV);
  m_qseqMeasureRunTimes = env->input()(m_option_qseq_measureRunTimes, UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV);
  m_qseqDataOutputPeriod = env->input()(m_option_qseq_dataOutputPeriod, UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV);
  m_qseqDataOutputFileName = env->input()(m_option_qseq_dataOutputFileName, UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV);
  m_qseqDataOutputFileType = env->input()(m_option_qseq_dataOutputFileType, UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV);

  // UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = env->input().vector_variable_size(m_option_qseq_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env->input()(m_option_qseq_dataOutputAllowedSet, i, i);
    m_qseqDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_qseq_computeStats = env->input()(m_option_qseq_computeStats, UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV);
#endif
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
}

// Copy constructor --------------------------------
McOptionsValues::McOptionsValues(const McOptionsValues& src)
{
  this->copy(src);
}
// Destructor ---------------------------------------
McOptionsValues::~McOptionsValues()
{
}
// Set methods --------------------------------------
McOptionsValues&
McOptionsValues::operator=(const McOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}
// Private methods-----------------------------------
void
McOptionsValues::checkOptions()
{
  // Do nothing
}

void
McOptionsValues::copy(const McOptionsValues& src)
{
  m_dataOutputFileName          = src.m_dataOutputFileName;
  m_dataOutputAllowedSet        = src.m_dataOutputAllowedSet;
  m_pseqDataOutputPeriod        = src.m_pseqDataOutputPeriod;
  m_pseqDataOutputFileName      = src.m_pseqDataOutputFileName;
  m_pseqDataOutputFileType      = src.m_pseqDataOutputFileType;
  m_pseqDataOutputAllowedSet    = src.m_pseqDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_pseqComputeStats            = src.m_pseqComputeStats;
#endif
  m_qseqDataInputFileName       = src.m_qseqDataInputFileName;
  m_qseqDataInputFileType       = src.m_qseqDataInputFileType;
  m_qseqSize                    = src.m_qseqSize;
  m_qseqDisplayPeriod           = src.m_qseqDisplayPeriod;
  m_qseqMeasureRunTimes         = src.m_qseqMeasureRunTimes;
  m_qseqDataOutputPeriod        = src.m_qseqDataOutputPeriod;
  m_qseqDataOutputFileName      = src.m_qseqDataOutputFileName;
  m_qseqDataOutputFileType      = src.m_qseqDataOutputFileType;
  m_qseqDataOutputAllowedSet    = src.m_qseqDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_qseqComputeStats            = src.m_qseqComputeStats;
#endif

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_alternativePSsOptionsValues = src.m_alternativePSsOptionsValues;
  m_alternativeQSsOptionsValues = src.m_alternativeQSsOptionsValues;
#endif

  return;
}

std::ostream & operator<<(std::ostream & os, const McOptionsValues & obj)
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  os <<         obj.m_option_dataOutputFileName   << " = " << obj.m_dataOutputFileName
     << "\n" << obj.m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_dataOutputAllowedSet.begin(); setIt != obj.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << obj.m_option_pseq_dataOutputPeriod     << " = " << obj.m_pseqDataOutputPeriod
     << "\n" << obj.m_option_pseq_dataOutputFileName   << " = " << obj.m_pseqDataOutputFileName
     << "\n" << obj.m_option_pseq_dataOutputFileType   << " = " << obj.m_pseqDataOutputFileType
     << "\n" << obj.m_option_pseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_pseqDataOutputAllowedSet.begin(); setIt != obj.m_pseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << obj.m_option_pseq_computeStats         << " = " << obj.m_pseqComputeStats
#endif
     << "\n" << obj.m_option_qseq_dataInputFileName    << " = " << obj.m_qseqDataInputFileName
     << "\n" << obj.m_option_qseq_dataInputFileType    << " = " << obj.m_qseqDataInputFileType
     << "\n" << obj.m_option_qseq_size                 << " = " << obj.m_qseqSize
     << "\n" << obj.m_option_qseq_displayPeriod        << " = " << obj.m_qseqDisplayPeriod
     << "\n" << obj.m_option_qseq_measureRunTimes      << " = " << obj.m_qseqMeasureRunTimes
     << "\n" << obj.m_option_qseq_dataOutputPeriod     << " = " << obj.m_qseqDataOutputPeriod
     << "\n" << obj.m_option_qseq_dataOutputFileName   << " = " << obj.m_qseqDataOutputFileName
     << "\n" << obj.m_option_qseq_dataOutputFileType   << " = " << obj.m_qseqDataOutputFileType
     << "\n" << obj.m_option_qseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_qseqDataOutputAllowedSet.begin(); setIt != obj.m_qseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  os << "\n" << obj.m_option_qseq_computeStats << " = " << obj.m_qseqComputeStats;
#endif

  return os;
}

// --------------------------------------------------
//MonteCarloSGOptions ------------------------
// --------------------------------------------------

// Default constructor -----------------------------
MonteCarloSGOptions::MonteCarloSGOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_ov                              (NULL,NULL),
  m_pseqStatisticalOptionsObj       (NULL),
  m_qseqStatisticalOptionsObj       (NULL),
#else
  m_ov                              (),
#endif
  m_prefix                          ((std::string)(prefix) + "mc_"),
  m_env                             (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                     (new boost::program_options::options_description("Monte Carlo options")),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                     (m_prefix + "help"                       ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"         ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"       ),
  m_option_pseq_dataOutputPeriod    (m_prefix + "pseq_dataOutputPeriod"      ),
  m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"    ),
  m_option_pseq_dataOutputFileType  (m_prefix + "pseq_dataOutputFileType"    ),
  m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"  ),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"          ),
#endif
  m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"     ),
  m_option_qseq_dataInputFileType   (m_prefix + "qseq_dataInputFileType"     ),
  m_option_qseq_size                (m_prefix + "qseq_size"                  ),
  m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"         ),
  m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"       ),
  m_option_qseq_dataOutputPeriod    (m_prefix + "qseq_dataOutputPeriod"      ),
  m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"    ),
  m_option_qseq_dataOutputFileType  (m_prefix + "qseq_dataOutputFileType"    ),
  m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet"  )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  ,
  m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"          )
#endif
{
  queso_deprecated();
  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the absence of an options input file"));
}
// Constructor 2 -----------------------------------
MonteCarloSGOptions::MonteCarloSGOptions(
  const BaseEnvironment& env,
  const char*                   prefix,
  const McOptionsValues& alternativeOptionsValues)
  :
  m_ov                              (alternativeOptionsValues),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_pseqStatisticalOptionsObj       (NULL),
  m_qseqStatisticalOptionsObj       (NULL),
#endif
  m_prefix                          ((std::string)(prefix) + "mc_"),
  m_env                             (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                     (NULL),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                     (m_prefix + "help"                     ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"       ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"     ),
  m_option_pseq_dataOutputPeriod    (m_prefix + "pseq_dataOutputPeriod"    ),
  m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"  ),
  m_option_pseq_dataOutputFileType  (m_prefix + "pseq_dataOutputFileType"  ),
  m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"        ),
#endif
  m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"   ),
  m_option_qseq_dataInputFileType   (m_prefix + "qseq_dataInputFileType"   ),
  m_option_qseq_size                (m_prefix + "qseq_size"                ),
  m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"       ),
  m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"     ),
  m_option_qseq_dataOutputPeriod    (m_prefix + "qseq_dataOutputPeriod"    ),
  m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"  ),
  m_option_qseq_dataOutputFileType  (m_prefix + "qseq_dataOutputFileType"  ),
  m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet")
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  ,
  m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"        )
#endif
{
  queso_deprecated();
  queso_require_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the existence of an options input file"));

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In MonteCarloSGOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

  // dakota
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_ov.m_pseqComputeStats) m_pseqStatisticalOptionsObj =
    new SequenceStatisticalOptions(m_env,m_prefix + "pseq_",m_ov.m_alternativePSsOptionsValues);
  if (m_ov.m_qseqComputeStats) m_qseqStatisticalOptionsObj =
    new SequenceStatisticalOptions(m_env,m_prefix + "qseq_",m_ov.m_alternativeQSsOptionsValues);
#endif
}
// Destructor --------------------------------------
MonteCarloSGOptions::~MonteCarloSGOptions()
{
  queso_deprecated();

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_pseqStatisticalOptionsObj) delete m_pseqStatisticalOptionsObj; // dakota
  if (m_qseqStatisticalOptionsObj) delete m_qseqStatisticalOptionsObj; // dakota
#endif
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  if (m_optionsDesc              ) delete m_optionsDesc;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
}
// I/O methods -------------------------------------
void
MonteCarloSGOptions::scanOptionsValues()
{
  queso_deprecated();

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In MonteCarloSGOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

  // dakota
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_ov.m_pseqComputeStats) m_pseqStatisticalOptionsObj =
    new SequenceStatisticalOptions(m_env,m_prefix + "pseq_");
  if (m_ov.m_qseqComputeStats) m_qseqStatisticalOptionsObj =
    new SequenceStatisticalOptions(m_env,m_prefix + "qseq_");
#endif
  return;
}
// Private methods ---------------------------------
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
void
MonteCarloSGOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                            "produce help message for Monte Carlo distribution calculator")
    (m_option_dataOutputFileName.c_str(),        boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV       ), "name of generic data output file"                            )
    (m_option_dataOutputAllowedSet.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV     ), "subEnvs that will write to generic data output file"         )
    (m_option_pseq_dataOutputPeriod.c_str(),     boost::program_options::value<unsigned int>()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV     ), "period of message display during param sequence generation"  )
    (m_option_pseq_dataOutputFileName.c_str(),   boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file for parameters"                     )
    (m_option_pseq_dataOutputFileType.c_str(),   boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV  ), "type of data output file for parameters"                     )
    (m_option_pseq_dataOutputAllowedSet.c_str(), boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file for parameters"  )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    (m_option_pseq_computeStats.c_str(),         boost::program_options::value<bool        >()->default_value(UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV          ), "compute statistics on sequence of parameter"                 )
#endif
    (m_option_qseq_dataInputFileName.c_str(),    boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV   ), "name of data input file for qois"                            )
    (m_option_qseq_dataInputFileType.c_str(),    boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV   ), "type of data input file for qois"                            )
    (m_option_qseq_size.c_str(),                 boost::program_options::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_SIZE_ODV                   ), "size of qoi sequence"                                        )
    (m_option_qseq_displayPeriod.c_str(),        boost::program_options::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV         ), "period of message display during qoi sequence generation"    )
    (m_option_qseq_measureRunTimes.c_str(),      boost::program_options::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV      ), "measure run times"                                           )
    (m_option_qseq_dataOutputPeriod.c_str(),     boost::program_options::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV     ), "period of message display during qoi sequence generation"    )
    (m_option_qseq_dataOutputFileName.c_str(),   boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file for qois"                           )
    (m_option_qseq_dataOutputFileType.c_str(),   boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV  ), "type of data output file for qois"                           )
    (m_option_qseq_dataOutputAllowedSet.c_str(), boost::program_options::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file for qois"        )
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    (m_option_qseq_computeStats.c_str(),         boost::program_options::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV          ), "compute statistics on sequence of qoi"                       )
#endif
   ;

  return;
}

void
MonteCarloSGOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
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

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputPeriod)) {
    m_ov.m_pseqDataOutputPeriod = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_pseq_dataOutputPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputFileName)) {
    m_ov.m_pseqDataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_pseq_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputFileType)) {
    m_ov.m_pseqDataOutputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_pseq_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputAllowedSet)) {
    m_ov.m_pseqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_pseq_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_pseqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_env.allOptionsMap().count(m_option_pseq_computeStats)) {
    m_ov.m_pseqComputeStats = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_pseq_computeStats]).as<bool>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_qseq_dataInputFileName)) {
    m_ov.m_qseqDataInputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataInputFileType)) {
    m_ov.m_qseqDataInputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_size)) {
    m_ov.m_qseqSize = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_size]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_displayPeriod)) {
    m_ov.m_qseqDisplayPeriod = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_displayPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_measureRunTimes)) {
    m_ov.m_qseqMeasureRunTimes = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_measureRunTimes]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputPeriod)) {
    m_ov.m_qseqDataOutputPeriod = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataOutputPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputFileName)) {
    m_ov.m_qseqDataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputFileType)) {
    m_ov.m_qseqDataOutputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputAllowedSet)) {
    m_ov.m_qseqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_qseq_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_qseqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_env.allOptionsMap().count(m_option_qseq_computeStats)) {
    m_ov.m_qseqComputeStats = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_qseq_computeStats]).as<bool>();
  }
#endif
  return;
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

void
MonteCarloSGOptions::print(std::ostream& os) const
{
  queso_deprecated();

  os <<         m_option_dataOutputFileName   << " = " << m_ov.m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_dataOutputPeriod     << " = " << m_ov.m_pseqDataOutputPeriod
     << "\n" << m_option_pseq_dataOutputFileName   << " = " << m_ov.m_pseqDataOutputFileName
     << "\n" << m_option_pseq_dataOutputFileType   << " = " << m_ov.m_pseqDataOutputFileType
     << "\n" << m_option_pseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_pseqDataOutputAllowedSet.begin(); setIt != m_ov.m_pseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
     << "\n" << m_option_pseq_computeStats         << " = " << m_ov.m_pseqComputeStats
#endif
     << "\n" << m_option_qseq_dataInputFileName    << " = " << m_ov.m_qseqDataInputFileName
     << "\n" << m_option_qseq_dataInputFileType    << " = " << m_ov.m_qseqDataInputFileType
     << "\n" << m_option_qseq_size                 << " = " << m_ov.m_qseqSize
     << "\n" << m_option_qseq_displayPeriod        << " = " << m_ov.m_qseqDisplayPeriod
     << "\n" << m_option_qseq_measureRunTimes      << " = " << m_ov.m_qseqMeasureRunTimes
     << "\n" << m_option_qseq_dataOutputPeriod     << " = " << m_ov.m_qseqDataOutputPeriod
     << "\n" << m_option_qseq_dataOutputFileName   << " = " << m_ov.m_qseqDataOutputFileName
     << "\n" << m_option_qseq_dataOutputFileType   << " = " << m_ov.m_qseqDataOutputFileType
     << "\n" << m_option_qseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_qseqDataOutputAllowedSet.begin(); setIt != m_ov.m_qseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  os << "\n" << m_option_qseq_computeStats << " = " << m_ov.m_qseqComputeStats;
#endif

  return;
}

std::ostream& operator<<(std::ostream& os, const MonteCarloSGOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
