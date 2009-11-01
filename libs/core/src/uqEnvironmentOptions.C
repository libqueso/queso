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

#include <uqEnvironmentOptions.h>
#include <uqMiscellaneous.h>

uqEnvironmentOptionsClass::uqEnvironmentOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_numSubEnvironments         (UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
  m_subDisplayFileName         (UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),
  m_subDisplayAllowAll         (UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),
//m_subDisplayAllowedSet       (),
  m_displayVerbosity           (UQ_ENV_DISPLAY_VERBOSITY_ODV),
  m_syncVerbosity              (UQ_ENV_SYNC_VERBOSITY_ODV),
  m_seed                       (UQ_ENV_SEED_ODV),
  m_numDebugParams             (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams                (m_numDebugParams,0.),
  m_env                        (env),
  m_prefix                     ((std::string)(prefix) + "env_"),
  m_optionsDesc                (new po::options_description("Environment options")),
  m_option_help                (m_prefix + "help"                ),
  m_option_numSubEnvironments  (m_prefix + "numSubEnvironments"  ),
  m_option_subDisplayFileName  (m_prefix + "subDisplayFileName"  ),
  m_option_subDisplayAllowAll  (m_prefix + "subDisplayAllowAll"  ),
  m_option_subDisplayAllowedSet(m_prefix + "subDisplayAllowedSet"),
  m_option_displayVerbosity    (m_prefix + "displayVerbosity"    ),
  m_option_syncVerbosity       (m_prefix + "syncVerbosity"       ),
  m_option_seed                (m_prefix + "seed"                )
{
}

uqEnvironmentOptionsClass::~uqEnvironmentOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
}

void
uqEnvironmentOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  // 'm_subDisplayOutputFile' is still not available at this moment. Use 'std::cout'
  //if (m_env.subScreenFile() != NULL) {
  //  *m_env.subScreenFile()
  if ((m_env.fullRank() == 0) && (m_env.displayVerbosity() >= 3)) {
    std::cout << "In uqEnvironmentOptionsClass::scanOptionsValues()"
              << ": after getting values of options with prefix '" << m_prefix
              << "', state of object is:"
              << "\n" << *this
              << std::endl;
  }

  return;
}

void
uqEnvironmentOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                               "produce help message for uq environment"  )
    (m_option_numSubEnvironments.c_str(),   po::value<unsigned int>()->default_value(UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),    "number of subEnvironments"                )
    (m_option_subDisplayFileName.c_str(),   po::value<std::string >()->default_value(UQ_ENV_SUB_DISPLAY_FILE_NAME_ODV),   "output filename for subscreen writing"    )
    (m_option_subDisplayAllowAll.c_str(),   po::value<bool        >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOW_ALL_ODV),   "Allow all subEnvs to write to output file")
    (m_option_subDisplayAllowedSet.c_str(), po::value<std::string >()->default_value(UQ_ENV_SUB_DISPLAY_ALLOWED_SET_ODV), "subEnvs that will write to output file"   )
    (m_option_displayVerbosity.c_str(),     po::value<unsigned int>()->default_value(UQ_ENV_DISPLAY_VERBOSITY_ODV),       "set verbosity"                            )
    (m_option_syncVerbosity.c_str(),        po::value<unsigned int>()->default_value(UQ_ENV_SYNC_VERBOSITY_ODV),          "set sync verbosity"                       )
    (m_option_seed.c_str(),                 po::value<int         >()->default_value(UQ_ENV_SEED_ODV),                    "set seed"                                 )
  //(m_option_numDebugParams.c_str(),       po::value<unsigned int>()->default_value(UQ_ENV_NUM_DEBUG_PARAMS_ODV),        "set number of debug parameters"           )
  ;

  return;
}

void
uqEnvironmentOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    // 'm_subDisplayOutputFile' is still not available at this moment. Use 'std::cout'
    if (m_env.fullRank() == 0) std::cout << optionsDesc
                                         << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_numSubEnvironments.c_str())) {
    m_numSubEnvironments = m_env.allOptionsMap()[m_option_numSubEnvironments].as<unsigned int>();
  }
  if ((m_env.fullComm().NumProc()%m_numSubEnvironments) != 0) {
    std::cerr << "In uqBaseEnvironmentClass::getMyOptionValues()"
              << ": m_env.fullComm().NumProc() = " << m_env.fullComm().NumProc()
              << ", m_numSubEnvironments = "       << m_numSubEnvironments
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_env.fullComm().NumProc()%m_numSubEnvironments) != 0,
                      m_env.fullRank(),
                      "uqBaseEnvironmentClass::getMyOptionValues()",
                      "total number of processors in environment must be multiple of the specified number of subEnvironments");

  if (m_env.allOptionsMap().count(m_option_subDisplayFileName.c_str())) {
    m_subDisplayFileName = m_env.allOptionsMap()[m_option_subDisplayFileName].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_subDisplayAllowAll.c_str())) {
    m_subDisplayAllowAll = m_env.allOptionsMap()[m_option_subDisplayAllowAll].as<bool>();
  }

  if (m_subDisplayAllowAll) {
    m_subDisplayAllowedSet.clear();
    // The line below is commented because 'm_subId' is not set at this point yet
    //m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }
  else if (m_env.allOptionsMap().count(m_option_subDisplayAllowedSet.c_str())) {
    m_subDisplayAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_subDisplayAllowedSet].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);
    //if (m_subDisplayOutputFile) {
    //  *m_subDisplayOutputFile << "In uqEnvironmentOptionsClass::getMyOptionValues(): allow = ";
    //  for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
    //    *m_subDisplayOutputFile << " " << tmpAllow[i];
    //  }
    //  *m_subDisplayOutputFile << std::endl;
    //}

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_subDisplayAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_displayVerbosity.c_str())) {
    m_displayVerbosity = m_env.allOptionsMap()[m_option_displayVerbosity].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_syncVerbosity.c_str())) {
    m_syncVerbosity = m_env.allOptionsMap()[m_option_syncVerbosity].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_seed.c_str())) {
    m_seed = m_env.allOptionsMap()[m_option_seed].as<int>();
  }

  //if (m_env.allOptionsMap().count(m_option_numDebugParams.c_str())) {
  //  m_seed = m_env.allOptionsMap()[m_option_numDebugParams].as<unsigned int>();
  //}

  return;
}

void
uqEnvironmentOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_numSubEnvironments   << " = " << m_numSubEnvironments
     << "\n" << m_option_subDisplayFileName   << " = " << m_subDisplayFileName
     << "\n" << m_option_subDisplayAllowAll   << " = " << m_subDisplayAllowAll
     << "\n" << m_option_subDisplayAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_subDisplayAllowedSet.begin(); setIt != m_subDisplayAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_displayVerbosity  << " = " << m_displayVerbosity
     << "\n" << m_option_syncVerbosity     << " = " << m_syncVerbosity
     << "\n" << m_option_seed              << " = " << m_seed
   //<< "\n" << m_option_numDebugParams    << " = " << m_numDebugParams
     << std::endl;
  return;
}

std::ostream& operator<<(std::ostream& os, const uqEnvironmentOptionsClass& obj)
{
  obj.print(os);

  return os;
}
