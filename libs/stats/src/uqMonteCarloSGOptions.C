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

uqMonteCarloSGOptionsClass::uqMonteCarloSGOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_prefix                          ((std::string)(prefix) + "mc_"),
  m_dataOutputFileName              (UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV     ),
//m_dataOutputAllowedSet            (),
  m_pseqDataOutputFileName          (UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV),
//m_pseqDataOutputAllowedSet        (),
  m_pseqComputeStats                (UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV        ),
  m_pseqStatisticalOptions          (NULL),
  m_qseqDataInputFileName           (UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV ),
  m_qseqSize                        (UQ_MOC_SG_QSEQ_SIZE_ODV                 ),
  m_qseqDisplayPeriod               (UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV       ),
  m_qseqMeasureRunTimes             (UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV    ),
  m_qseqDataOutputFileName          (UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV),
//m_qseqDataOutputAllowedSet        (),
  m_qseqComputeStats                (UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV        ),
  m_qseqStatisticalOptions          (NULL),
  m_env                             (env),
  m_optionsDesc                     (new po::options_description("Monte Carlo options")),
  m_option_help                     (m_prefix + "help"                       ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"         ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"       ),
  m_option_pseq_dataOutputFileName  (m_prefix + "pseq_dataOutputFileName"    ),
  m_option_pseq_dataOutputAllowedSet(m_prefix + "pseq_dataOutputAllowedSet"  ),
  m_option_pseq_computeStats        (m_prefix + "pseq_computeStats"          ),
  m_option_qseq_dataInputFileName   (m_prefix + "qseq_dataInputFileName"     ),
  m_option_qseq_size                (m_prefix + "qseq_size"                  ),
  m_option_qseq_displayPeriod       (m_prefix + "qseq_displayPeriod"         ),
  m_option_qseq_measureRunTimes     (m_prefix + "qseq_measureRunTimes"       ),
  m_option_qseq_dataOutputFileName  (m_prefix + "qseq_dataOutputFileName"    ),
  m_option_qseq_dataOutputAllowedSet(m_prefix + "qseq_dataOutputAllowedSet"  ),
  m_option_qseq_computeStats        (m_prefix + "qseq_computeStats"          )
{
}

uqMonteCarloSGOptionsClass::~uqMonteCarloSGOptionsClass()
{
  //if (m_filteredChainStatisticalOptions) delete m_filteredChainStatisticalOptions;
  //if (m_rawChainStatisticalOptions     ) delete m_rawChainStatisticalOptions;
  if (m_optionsDesc                    ) delete m_optionsDesc;
} 

void
uqMonteCarloSGOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  //if (m_rawChainComputeStats     ) m_rawChainStatisticalOptions      = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "rawChain_"     );
  //if (m_filteredChainComputeStats) m_filteredChainStatisticalOptions = new uqSequenceStatisticalOptionsClass(m_env,m_prefix + "filteredChain_");

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
    (m_option_pseq_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file for parameters"  )
    (m_option_pseq_computeStats.c_str(),         po::value<bool        >()->default_value(UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV          ), "compute statistics on sequence of parameter"                 )
    (m_option_qseq_dataInputFileName.c_str(),    po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV   ), "name of data input file for qois"                            )
    (m_option_qseq_size.c_str(),                 po::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_SIZE_ODV                   ), "size of qoi sequence"                                        )
    (m_option_qseq_displayPeriod.c_str(),        po::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV         ), "period of message display during qoi sequence generation"    )
    (m_option_qseq_measureRunTimes.c_str(),      po::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV      ), "measure run times"                                           )
    (m_option_qseq_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file for qois"                           )
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
    m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputFileName)) {
    m_pseqDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_dataOutputAllowedSet)) {
    m_pseqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_pseq_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_pseqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_pseq_computeStats)) {
    m_pseqComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_computeStats]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataInputFileName)) {
    m_qseqDataInputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataInputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_size)) {
    m_qseqSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_size]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_displayPeriod)) {
    m_qseqDisplayPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_displayPeriod]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_measureRunTimes)) {
    m_qseqMeasureRunTimes = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_measureRunTimes]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputFileName)) {
    m_qseqDataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_dataOutputAllowedSet)) {
    m_qseqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_qseq_dataOutputAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_qseqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_qseq_computeStats)) {
    m_qseqComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_computeStats]).as<bool>();
  }

  return;
}

void
uqMonteCarloSGOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_dataOutputFileName   << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_dataOutputFileName   << " = " << m_pseqDataOutputFileName
     << "\n" << m_option_pseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_pseqDataOutputAllowedSet.begin(); setIt != m_pseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_computeStats         << " = " << m_pseqComputeStats
     << "\n" << m_option_qseq_dataInputFileName    << " = " << m_qseqDataInputFileName
     << "\n" << m_option_qseq_size                 << " = " << m_qseqSize
     << "\n" << m_option_qseq_displayPeriod        << " = " << m_qseqDisplayPeriod
     << "\n" << m_option_qseq_measureRunTimes      << " = " << m_qseqMeasureRunTimes
     << "\n" << m_option_qseq_dataOutputFileName   << " = " << m_qseqDataOutputFileName
     << "\n" << m_option_qseq_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_qseqDataOutputAllowedSet.begin(); setIt != m_qseqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_qseq_computeStats << " = " << m_qseqComputeStats;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGOptionsClass& obj)
{
  obj.print(os);

  return os;
}
