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

#include <uqMonteCarloSGOptions.h>
#include <uqMiscellaneous.h>

uqMcOptionsValuesClass::uqMcOptionsValuesClass(
  const uqSsOptionsValuesClass* alternativePSsOptionsValues,
  const uqSsOptionsValuesClass* alternativeQSsOptionsValues)
  :
  m_dataOutputFileName         (UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV     ),
//m_dataOutputAllowedSet       (),
  m_pseqDataOutputFileName     (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV),
  m_pseqDataOutputFileType     (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV),
//m_pseqDataOutputAllowedSet   (),
  m_pseqComputeStats           (UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV        ),
  m_qseqDataInputFileName      (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV ),
  m_qseqDataInputFileType      (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV ),
  m_qseqSize                   (UQ_MOC_SG_QSEQ_SIZE_ODV                 ),
  m_qseqDisplayPeriod          (UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV       ),
  m_qseqMeasureRunTimes        (UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV    ),
  m_qseqDataOutputFileName     (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV),
  m_qseqDataOutputFileType     (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV),
//m_qseqDataOutputAllowedSet   (),
  m_qseqComputeStats           (UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV        ),
  m_alternativePSsOptionsValues(),
  m_alternativeQSsOptionsValues()
{
  if (alternativePSsOptionsValues) m_alternativePSsOptionsValues = *alternativePSsOptionsValues;
  if (alternativeQSsOptionsValues) m_alternativeQSsOptionsValues = *alternativeQSsOptionsValues;
}

uqMcOptionsValuesClass::~uqMcOptionsValuesClass()
{
}

uqMcOptionsValuesClass::uqMcOptionsValuesClass(const uqMcOptionsValuesClass& src)
{
  this->copy(src);
}

uqMcOptionsValuesClass&
uqMcOptionsValuesClass::operator=(const uqMcOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}

void
uqMcOptionsValuesClass::copy(const uqMcOptionsValuesClass& src)
{
  m_dataOutputFileName          = src.m_dataOutputFileName;
  m_dataOutputAllowedSet        = src.m_dataOutputAllowedSet;
  m_pseqDataOutputFileName      = src.m_pseqDataOutputFileName;
  m_pseqDataOutputFileType      = src.m_pseqDataOutputFileType;
  m_pseqDataOutputAllowedSet    = src.m_pseqDataOutputAllowedSet; 
  m_pseqComputeStats            = src.m_pseqComputeStats;
  m_qseqDataInputFileName       = src.m_qseqDataInputFileName;
  m_qseqDataInputFileType       = src.m_qseqDataInputFileType;
  m_qseqSize                    = src.m_qseqSize;
  m_qseqDisplayPeriod           = src.m_qseqDisplayPeriod;
  m_qseqMeasureRunTimes         = src.m_qseqMeasureRunTimes;
  m_qseqDataOutputFileName      = src.m_qseqDataOutputFileName;
  m_qseqDataOutputFileType      = src.m_qseqDataOutputFileType;
  m_qseqDataOutputAllowedSet    = src.m_qseqDataOutputAllowedSet; 
  m_qseqComputeStats            = src.m_qseqComputeStats;

  m_alternativePSsOptionsValues = src.m_alternativePSsOptionsValues;
  m_alternativeQSsOptionsValues = src.m_alternativeQSsOptionsValues;

  return;
}

uqMonteCarloSGOptionsClass::uqMonteCarloSGOptionsClass(
  const uqBaseEnvironmentClass& env, 
  const char*                   prefix)
  :
  m_ov                              (NULL,NULL),
  m_pseqStatisticalOptionsObj       (NULL),
  m_qseqStatisticalOptionsObj       (NULL),
  m_prefix                          ((std::string)(prefix) + "mc_"),
  m_env                             (env),
  m_optionsDesc                     (new po::options_description("Monte Carlo options")),
  m_option_help                     (m_prefix + "help"                       ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"         ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"       ),
  m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"    ),
  m_option_pseq_dataOutputFileType  (m_prefix + "pseq_dataOutputFileType"    ),
  m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"  ),
  m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"          ),
  m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"     ),
  m_option_qseq_dataInputFileType   (m_prefix + "qseq_dataInputFileType"     ),
  m_option_qseq_size                (m_prefix + "qseq_size"                  ),
  m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"         ),
  m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"       ),
  m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"    ),
  m_option_qseq_dataOutputFileType  (m_prefix + "qseq_dataOutputFileType"    ),
  m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet"  ),
  m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"          )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "uqMonteCarloSGOptionsClass::constructor(1)",
                      "this constructor is incompatible with the abscense of an options input file");
}

uqMonteCarloSGOptionsClass::uqMonteCarloSGOptionsClass(
  const uqBaseEnvironmentClass& env, 
  const char*                   prefix,
  const uqMcOptionsValuesClass& alternativeOptionsValues)
  :
  m_ov                              (alternativeOptionsValues),
  m_pseqStatisticalOptionsObj       (NULL),
  m_qseqStatisticalOptionsObj       (NULL),
  m_prefix                          ((std::string)(prefix) + "mc_"),
  m_env                             (env),
  m_optionsDesc                     (NULL),
  m_option_help                     (m_prefix + "help"                     ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"       ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"     ),
  m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"  ),
  m_option_pseq_dataOutputFileType  (m_prefix + "pseq_dataOutputFileType"  ),
  m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"),
  m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"        ),
  m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"   ),
  m_option_qseq_dataInputFileType   (m_prefix + "qseq_dataInputFileType"   ),
  m_option_qseq_size                (m_prefix + "qseq_size"                ),
  m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"       ),
  m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"     ),
  m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"  ),
  m_option_qseq_dataOutputFileType  (m_prefix + "qseq_dataOutputFileType"  ),
  m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet"),
  m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"        )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "uqMonteCarloSGOptionsClass::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGOptionsClass::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

uqMonteCarloSGOptionsClass::~uqMonteCarloSGOptionsClass()
{
//if (m_filteredChainStatisticalOptionsObj) delete m_filteredChainStatisticalOptionsObj; // dakota
//if (m_rawChainStatisticalOptionsObj     ) delete m_rawChainStatisticalOptionsObj;      // dakota
  if (m_optionsDesc                       ) delete m_optionsDesc;
} 

void
uqMonteCarloSGOptionsClass::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "uqMonteCarloSGOptionsClass::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGOptionsClass::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

//if (m_rawChainComputeStats     ) m_rawChainStatisticalOptionsObj      = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "rawChain_"     ); // dakota
//if (m_filteredChainComputeStats) m_filteredChainStatisticalOptionsObj = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "filteredChain_"); // dakota

  return;
}

void
uqMonteCarloSGOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                            "produce help message for Monte Carlo distribution calculator")
    (m_option_dataOutputFileName.c_str(),        po::value<std::string >()->default_value(UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV       ), "name of generic data output file"                            )
    (m_option_dataOutputAllowedSet.c_str(),      po::value<std::string >()->default_value(UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV     ), "subEnvs that will write to generic data output file"         )
    (m_option_pseq_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file for parameters"                     )
    (m_option_pseq_dataOutputFileType.c_str(),   po::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV  ), "type of data output file for parameters"                     )
    (m_option_pseq_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file for parameters"  )
    (m_option_pseq_computeStats.c_str(),         po::value<bool        >()->default_value(UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV          ), "compute statistics on sequence of parameter"                 )
    (m_option_qseq_dataInputFileName.c_str(),    po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV   ), "name of data input file for qois"                            )
    (m_option_qseq_dataInputFileType.c_str(),    po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV   ), "type of data input file for qois"                            )
    (m_option_qseq_size.c_str(),                 po::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_SIZE_ODV                   ), "size of qoi sequence"                                        )
    (m_option_qseq_displayPeriod.c_str(),        po::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV         ), "period of message display during qoi sequence generation"    )
    (m_option_qseq_measureRunTimes.c_str(),      po::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV      ), "measure run times"                                           )
    (m_option_qseq_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file for qois"                           )
    (m_option_qseq_dataOutputFileType.c_str(),   po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV  ), "type of data output file for qois"                           )
    (m_option_qseq_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file for qois"        )
    (m_option_qseq_computeStats.c_str(),         po::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV          ), "compute statistics on sequence of qoi"                       )
   ;

  return;
}

void
uqMonteCarloSGOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_ov.m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputFileName)) {
    m_ov.m_pseqDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputFileType)) {
    m_ov.m_pseqDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputAllowedSet)) {
    m_ov.m_pseqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_pseq_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_pseqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_pseq_computeStats)) {
    m_ov.m_pseqComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_computeStats]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataInputFileName)) {
    m_ov.m_qseqDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataInputFileType)) {
    m_ov.m_qseqDataInputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataInputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_size)) {
    m_ov.m_qseqSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_size]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_displayPeriod)) {
    m_ov.m_qseqDisplayPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_displayPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_measureRunTimes)) {
    m_ov.m_qseqMeasureRunTimes = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_measureRunTimes]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputFileName)) {
    m_ov.m_qseqDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputFileType)) {
    m_ov.m_qseqDataOutputFileType = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputAllowedSet)) {
    m_ov.m_qseqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_qseq_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_qseqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_qseq_computeStats)) {
    m_ov.m_qseqComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_computeStats]).as<bool>();
  }

  return;
}

void
uqMonteCarloSGOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_dataOutputFileName   << " = " << m_ov.m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_dataOutputFileName   << " = " << m_ov.m_pseqDataOutputFileName
     << "\n" << m_option_pseq_dataOutputFileType   << " = " << m_ov.m_pseqDataOutputFileType
     << "\n" << m_option_pseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_pseqDataOutputAllowedSet.begin(); setIt != m_ov.m_pseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_computeStats         << " = " << m_ov.m_pseqComputeStats
     << "\n" << m_option_qseq_dataInputFileName    << " = " << m_ov.m_qseqDataInputFileName
     << "\n" << m_option_qseq_dataInputFileType    << " = " << m_ov.m_qseqDataInputFileType
     << "\n" << m_option_qseq_size                 << " = " << m_ov.m_qseqSize
     << "\n" << m_option_qseq_displayPeriod        << " = " << m_ov.m_qseqDisplayPeriod
     << "\n" << m_option_qseq_measureRunTimes      << " = " << m_ov.m_qseqMeasureRunTimes
     << "\n" << m_option_qseq_dataOutputFileName   << " = " << m_ov.m_qseqDataOutputFileName
     << "\n" << m_option_qseq_dataOutputFileType   << " = " << m_ov.m_qseqDataOutputFileType
     << "\n" << m_option_qseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_qseqDataOutputAllowedSet.begin(); setIt != m_ov.m_qseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_qseq_computeStats << " = " << m_ov.m_qseqComputeStats;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGOptionsClass& obj)
{
  obj.print(os);

  return os;
}
