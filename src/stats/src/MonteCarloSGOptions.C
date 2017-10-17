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

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#define GETPOT_NAMESPACE QUESO // So we don't clash with other getpots
#include <queso/getpot.h>
#undef GETPOT_NAMESPACE
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

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
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  : m_parser(new BoostInputOptionsParser())
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativePSsOptionsValues) m_alternativePSsOptionsValues = *alternativePSsOptionsValues;
  if (alternativeQSsOptionsValues) m_alternativeQSsOptionsValues = *alternativeQSsOptionsValues;
#endif
  this->set_defaults();
  this->set_prefix("");
}

McOptionsValues::McOptionsValues(
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  const SsOptionsValues* alternativePSsOptionsValues,
  const SsOptionsValues* alternativeQSsOptionsValues,
#endif
  const BaseEnvironment * env, const char * prefix
  )
{
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (alternativePSsOptionsValues) m_alternativePSsOptionsValues = *alternativePSsOptionsValues;
  if (alternativeQSsOptionsValues) m_alternativeQSsOptionsValues = *alternativeQSsOptionsValues;
#endif

  this->set_defaults();
  this->parse(*env, prefix);
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
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

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


void McOptionsValues::set_defaults()
{
  m_help = UQ_MOC_SG_HELP;
  m_dataOutputFileName = UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV;
  //m_dataOutputAllowedSet
  m_pseqDataOutputPeriod = UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV;
  m_pseqDataOutputFileName = UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV;
  m_pseqDataOutputFileType = UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV;
  //m_pseqDataOutputAllowedSet
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_pseqComputeStats = UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV;
#endif
  m_qseqDataInputFileName = UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV;
  m_qseqDataInputFileType = UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV;
  m_qseqSize = UQ_MOC_SG_QSEQ_SIZE_ODV;
  m_qseqDisplayPeriod = UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV;
  m_qseqMeasureRunTimes = UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV;
  m_qseqDataOutputPeriod = UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV;
  m_qseqDataOutputFileName = UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV;
  m_qseqDataOutputFileType = UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV;
  //m_qseqDataOutputAllowedSet
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_qseqComputeStats            = UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV;
#endif
}

void McOptionsValues::set_prefix(const std::string& prefix)
{
  m_prefix = prefix + "mc_";

  m_option_help = m_prefix + "help";
  m_option_dataOutputFileName = m_prefix + "dataOutputFileName";
  m_option_dataOutputAllowedSet = m_prefix + "dataOutputAllowedSet";
  m_option_pseq_dataOutputPeriod = m_prefix + "pseq_dataOutputPeriod";
  m_option_pseq_dataOutputFileName = m_prefix + "pseq_dataOutputFileName";
  m_option_pseq_dataOutputFileType = m_prefix + "pseq_dataOutputFileType";
  m_option_pseq_dataOutputAllowedSet = m_prefix + "pseq_dataOutputAllowedSet";
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_pseq_computeStats = m_prefix + "pseq_computeStats";
#endif
  m_option_qseq_dataInputFileName = m_prefix + "qseq_dataInputFileName";
  m_option_qseq_dataInputFileType = m_prefix + "qseq_dataInputFileType";
  m_option_qseq_size = m_prefix + "qseq_size";
  m_option_qseq_displayPeriod = m_prefix + "qseq_displayPeriod";
  m_option_qseq_measureRunTimes = m_prefix + "qseq_measureRunTimes";
  m_option_qseq_dataOutputPeriod = m_prefix + "qseq_dataOutputPeriod";
  m_option_qseq_dataOutputFileName = m_prefix + "qseq_dataOutputFileName";
  m_option_qseq_dataOutputFileType = m_prefix + "qseq_dataOutputFileType";
  m_option_qseq_dataOutputAllowedSet = m_prefix + "qseq_dataOutputAllowedSet";
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_option_qseq_computeStats = m_prefix + "qseq_computeStats";
#endif
}

void McOptionsValues::parse(const BaseEnvironment& env, const std::string& prefix)
{
  this->set_prefix(prefix);

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  m_parser.reset(new BoostInputOptionsParser(env.optionsInputFileName()));

  m_parser->registerOption<std::string>
    (m_option_help, m_help,
     "produce help message for Monte Carlo distribution calculator");
  m_parser->registerOption<std::string>
    (m_option_dataOutputFileName, m_dataOutputFileName,
     "name of generic data output file");
  m_parser->registerOption<std::string>
    (m_option_dataOutputAllowedSet, container_to_string(m_dataOutputAllowedSet),
     "subEnvs that will write to generic data output file");
  m_parser->registerOption<unsigned int>
    (m_option_pseq_dataOutputPeriod, m_pseqDataOutputPeriod,
     "period of message display during param sequence generation");
  m_parser->registerOption<std::string>
    (m_option_pseq_dataOutputFileName, m_pseqDataOutputFileName,
     "name of data output file for parameters");
  m_parser->registerOption<std::string>
    (m_option_pseq_dataOutputFileType, m_pseqDataOutputFileType,
     "type of data output file for parameters");
  m_parser->registerOption<std::string>
    (m_option_pseq_dataOutputAllowedSet,
     container_to_string(m_pseqDataOutputAllowedSet),
     "subEnvs that will write to data output file for parameters");
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool>
    (m_option_pseq_computeStats, m_pseqComputeStats,
    "compute statistics on sequence of parameter");
#endif
  m_parser->registerOption<std::string>
    (m_option_qseq_dataInputFileName, m_qseqDataInputFileName,
    "name of data input file for qois");
  m_parser->registerOption<std::string>
    (m_option_qseq_dataInputFileType, m_qseqDataInputFileType,
    "type of data input file for qois");
  m_parser->registerOption<unsigned int>
    (m_option_qseq_size, m_qseqSize, "size of qoi sequence");
  m_parser->registerOption<unsigned int>
    (m_option_qseq_displayPeriod, m_qseqDisplayPeriod,
    "period of message display during qoi sequence generation");
  m_parser->registerOption<bool>
    (m_option_qseq_measureRunTimes, m_qseqMeasureRunTimes, "measure run times");
  m_parser->registerOption<unsigned int>
    (m_option_qseq_dataOutputPeriod, m_qseqDataOutputPeriod,
    "period of message display during qoi sequence generation");
  m_parser->registerOption<std::string>
    (m_option_qseq_dataOutputFileName, m_qseqDataOutputFileName,
    "name of data output file for qois");
  m_parser->registerOption<std::string>
    (m_option_qseq_dataOutputFileType, m_qseqDataOutputFileType,
    "type of data output file for qois");
  m_parser->registerOption<std::string>
    (m_option_qseq_dataOutputAllowedSet,
     container_to_string(m_qseqDataOutputAllowedSet),
    "subEnvs that will write to data output file for qois");
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->registerOption<bool>
    (m_option_qseq_computeStats, m_qseqComputeStats,
     "compute statistics on sequence of qoi");
#endif

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help, m_help);
  m_parser->getOption<std::string>(m_option_dataOutputFileName, m_dataOutputFileName);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet, m_dataOutputAllowedSet);
  m_parser->getOption<unsigned int>(m_option_pseq_dataOutputPeriod, m_pseqDataOutputPeriod);
  m_parser->getOption<std::string>(m_option_pseq_dataOutputFileName, m_pseqDataOutputFileName);
  m_parser->getOption<std::string>(m_option_pseq_dataOutputFileType, m_pseqDataOutputFileType);
  m_parser->getOption<std::set<unsigned int> >(m_option_pseq_dataOutputAllowedSet, m_pseqDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool>(m_option_pseq_computeStats, m_pseq_computeStats);
#endif
  m_parser->getOption<std::string>(m_option_qseq_dataInputFileName, m_qseqDataInputFileName);
  m_parser->getOption<std::string>(m_option_qseq_dataInputFileType, m_qseqDataInputFileType);
  m_parser->getOption<unsigned int>(m_option_qseq_size, m_qseqSize);
  m_parser->getOption<unsigned int>(m_option_qseq_displayPeriod, m_qseqDisplayPeriod);
  m_parser->getOption<bool>(m_option_qseq_measureRunTimes, m_qseqMeasureRunTimes);
  m_parser->getOption<unsigned int>(m_option_qseq_dataOutputPeriod, m_qseqDataOutputPeriod);
  m_parser->getOption<std::string>(m_option_qseq_dataOutputFileName, m_qseqDataOutputFileName);
  m_parser->getOption<std::string>(m_option_qseq_dataOutputFileType, m_qseqDataOutputFileType);
  m_parser->getOption<std::set<unsigned int> >(m_option_qseq_dataOutputAllowedSet, m_qseqDataOutputAllowedSet);
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_parser->getOption<bool>(m_option_qseq_computeStats, m_qseq_computeStats);
#endif
#else
  m_help = env.input()(m_option_help, m_help);
  m_dataOutputFileName = env.input()(m_option_dataOutputFileName, m_dataOutputFileName);

  // UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = env.input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env.input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

  m_pseqDataOutputPeriod = env.input()(m_option_pseq_dataOutputPeriod, m_pseqDataOutputPeriod);
  m_pseqDataOutputFileName = env.input()(m_option_pseq_dataOutputFileName, m_pseqDataOutputFileName);
  m_pseqDataOutputFileType = env.input()(m_option_pseq_dataOutputFileType, m_pseqDataOutputFileType);

  // UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = env.input().vector_variable_size(m_option_pseq_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env.input()(m_option_pseq_dataOutputAllowedSet, i, i);
    m_pseqDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_pseq_computeStats = env.input()(m_option_pseq_computeStats, m_pseq_computeStats);
#endif
  m_qseqDataInputFileName = env.input()(m_option_qseq_dataInputFileName, m_qseqDataInputFileName);
  m_qseqDataInputFileType = env.input()(m_option_qseq_dataInputFileType, m_qseqDataInputFileType);
  m_qseqSize = env.input()(m_option_qseq_size, m_qseqSize);
  m_qseqDisplayPeriod = env.input()(m_option_qseq_displayPeriod, m_qseqDisplayPeriod);
  m_qseqMeasureRunTimes = env.input()(m_option_qseq_measureRunTimes, m_qseqMeasureRunTimes);
  m_qseqDataOutputPeriod = env.input()(m_option_qseq_dataOutputPeriod, m_qseqDataOutputPeriod);
  m_qseqDataOutputFileName = env.input()(m_option_qseq_dataOutputFileName, m_qseqDataOutputFileName);
  m_qseqDataOutputFileType = env.input()(m_option_qseq_dataOutputFileType, m_qseqDataOutputFileType);

  // UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  size = env.input().vector_variable_size(m_option_qseq_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env.input()(m_option_qseq_dataOutputAllowedSet, i, i);
    m_qseqDataOutputAllowedSet.insert(allowed);
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_qseq_computeStats = env.input()(m_option_qseq_computeStats, m_qseq_computeStats);
#endif
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

}  // End namespace QUESO
