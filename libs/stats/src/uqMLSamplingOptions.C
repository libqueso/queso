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

#include <uqMLSamplingOptions.h>
#include <uqMiscellaneous.h>

uqMLSamplingOptionsClass::uqMLSamplingOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_prefix                         ((std::string)(prefix) + "ml_"                 ),
  m_checkpointLevel                (UQ_ML_SAMPLING_CHECKPOINT_LEVEL_ODV           ),
  m_checkpointOutputFileName       (UQ_ML_SAMPLING_CHECKPOINT_OUTPUT_FILE_NAME_ODV),
  m_restartLevel                   (UQ_ML_SAMPLING_RESTART_LEVEL_ODV              ),
  m_restartInputFileName           (UQ_ML_SAMPLING_RESTART_INPUT_FILE_NAME_ODV    ),
  m_dataOutputFileName             (UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV      ),
//m_dataOutputAllowedSet           (),
  m_env                            (env),
  m_optionsDesc                    (new po::options_description("Multilevel sampling options")),
  m_option_help                    (m_prefix + "help"                    ),
  m_option_checkpointLevel         (m_prefix + "checkpointLevel"         ),
  m_option_checkpointOutputFileName(m_prefix + "checkpointOutputFileName"),
  m_option_restartLevel            (m_prefix + "restartLevel"            ),
  m_option_restartInputFileName    (m_prefix + "restartInputFileName"    ),
  m_option_dataOutputFileName      (m_prefix + "dataOutputFileName"      ),
  m_option_dataOutputAllowedSet    (m_prefix + "dataOutputAllowedSet"    )
{
}

uqMLSamplingOptionsClass::~uqMLSamplingOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
} 

void
uqMLSamplingOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqMLSamplingOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
uqMLSamplingOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  // ernesto
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                       "produce help message for multilevel sampling options")
    (m_option_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  ), "name of generic output file"                         )
    (m_option_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to generic output file"      )
  ;

  return;
}

void
uqMLSamplingOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  // ernesto

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

  return;
}

void
uqMLSamplingOptionsClass::print(std::ostream& os) const
{
  // ernesto
  os <<         m_option_dataOutputFileName   << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n";

  return;
}

std::ostream& operator<<(std::ostream& os, const uqMLSamplingOptionsClass& obj)
{
  obj.print(os);

  return os;
}
