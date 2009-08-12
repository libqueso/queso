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
  m_prefix                     ((std::string)(prefix) + "ml_"),
  m_dataOutputFileName         (UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV),
//m_dataOutputAllowedSet       (),
//m_initialExponent            (UQ_ML_SAMPLING_INITIAL_EXPONENT_ODV),
  m_maxNumberOfLevels          (UQ_ML_SAMPLING_MAX_NUMBER_OF_LEVELS_ODV),
  m_levelOptions               (0),
  m_env                        (env),
  m_optionsDesc                (new po::options_description("Multilevel sampling options")),
  m_option_help                (m_prefix + "help"                ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
//m_option_initialExponent     (m_prefix + "initialExponent"     ),
  m_option_maxNumberOfLevels   (m_prefix + "maxNumberOfLevels"   )
{
}

uqMLSamplingOptionsClass::~uqMLSamplingOptionsClass()
{
  for (unsigned int i = 0; i < m_levelOptions.size(); ++i) {
    if (m_levelOptions[i] != NULL) {
      delete m_levelOptions[i];
      m_levelOptions[i] = NULL;
    }
  }
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
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  char tmpSufix[256];
  m_levelOptions.resize(m_maxNumberOfLevels,NULL);
  for (unsigned int i = 0; i < m_levelOptions.size(); ++i) {
    sprintf(tmpSufix,"%d_",i+1); // Yes, '+1'
    m_levelOptions[i] = new uqMLSamplingLevelOptionsClass(m_env,(m_prefix + tmpSufix).c_str());
    m_levelOptions[i]->scanOptionsValues();
  }

  return;
};

void
uqMLSamplingOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                       "produce help message for multilevel sampling options")
    (m_option_dataOutputFileName.c_str(),   po::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_FILE_NAME_ODV  ), "name of generic output file"                         )
    (m_option_dataOutputAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_ML_SAMPLING_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to generic output file"      )
  //(m_option_initialExponent.c_str(),      po::value<double      >()->default_value(UQ_ML_SAMPLING_INITIAL_EXPONENT_ODV       ), "initial exponent"                                    )
    (m_option_maxNumberOfLevels.c_str(),    po::value<unsigned int>()->default_value(UQ_ML_SAMPLING_MAX_NUMBER_OF_LEVELS_ODV   ), "maximum number of levels"                            )
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

//if (m_env.allOptionsMap().count(m_option_initialExponent.c_str())) {
//  m_initialExponent = ((const po::variable_value&) m_env.allOptionsMap()[m_option_initialExponent.c_str()]).as<double>();
//}

  if (m_env.allOptionsMap().count(m_option_maxNumberOfLevels.c_str())) {
    m_maxNumberOfLevels = ((const po::variable_value&) m_env.allOptionsMap()[m_option_maxNumberOfLevels.c_str()]).as<unsigned int>();
  }

  return;
}

void
uqMLSamplingOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_dataOutputFileName   << " = " << m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
//os << "\n" << m_option_initialExponent   << " = " << m_initialExponent
  os << "\n" << m_option_maxNumberOfLevels << " = " << m_maxNumberOfLevels
     << "\n";

  for (unsigned int i = 0; i < m_levelOptions.size(); ++i) {
    os << "\n" << *(m_levelOptions[i]);
  }

  return;
}

std::ostream& operator<<(std::ostream& os, const uqMLSamplingOptionsClass& obj)
{
  obj.print(os);

  return os;
}
