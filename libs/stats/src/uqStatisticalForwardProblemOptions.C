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

#include <uqStatisticalForwardProblemOptions.h>
#include <uqMiscellaneous.h>

uqStatisticalForwardProblemOptionsClass::uqStatisticalForwardProblemOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_prefix                     ((std::string)(prefix) + "fp_"   ),
  m_computeSolution            (UQ_SFP_COMPUTE_SOLUTION_ODV     ),
  m_computeCovariances         (UQ_SFP_COMPUTE_COVARIANCES_ODV  ),
  m_computeCorrelations        (UQ_SFP_COMPUTE_CORRELATIONS_ODV ),
  m_dataOutputFileName         (UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV),
//m_dataOutputAllowedSet       (),
#ifdef UQ_SFP_READS_SOLVER_OPTION
  m_solverString               (UQ_SFP_SOLVER_ODV               ),
#endif
  m_env                        (env),
  m_optionsDesc                (new po::options_description("Statistical Forward Problem options")),
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_computeCovariances  (m_prefix + "computeCovariances"  ),
  m_option_computeCorrelations (m_prefix + "computeCorrelations" ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet")
#ifdef UQ_SFP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              )
#endif
{
}

uqStatisticalForwardProblemOptionsClass::~uqStatisticalForwardProblemOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
} 

void
uqStatisticalForwardProblemOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqStatisticalForwardProblemOptionsClass::scanOptionsValues()"
                            << ": after getting values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
uqStatisticalForwardProblemOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                              "produce help message for statistical forward problem")
    (m_option_computeSolution.c_str(),      po::value<bool       >()->default_value(UQ_SFP_COMPUTE_SOLUTION_ODV       ), "compute solution process"                            )
    (m_option_computeCovariances.c_str(),   po::value<bool       >()->default_value(UQ_SFP_COMPUTE_COVARIANCES_ODV    ), "compute pq covariances"                              )
    (m_option_computeCorrelations.c_str(),  po::value<bool       >()->default_value(UQ_SFP_COMPUTE_CORRELATIONS_ODV   ), "compute pq correlations"                             )
    (m_option_dataOutputFileName.c_str(),   po::value<std::string>()->default_value(UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file"                            )
    (m_option_dataOutputAllowedSet.c_str(), po::value<std::string>()->default_value(UQ_SFP_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file"         )
#ifdef UQ_SFP_READS_SOLVER_OPTION
    (m_option_solver.c_str(),               po::value<std::string>()->default_value(UQ_SFP_SOLVER_ODV                 ), "algorithm for propagation"                           )
#endif
  ;

  return;
}

void
uqStatisticalForwardProblemOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_computeSolution)) {
    m_computeSolution = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeSolution]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeCovariances)) {
    m_computeCovariances = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeCovariances]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeCorrelations)) {
    m_computeCorrelations = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeCorrelations]).as<bool>();
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

#ifdef UQ_SFP_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver)) {
    m_solverString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_solver]).as<std::string>();
  }
#endif

  return;
}

void
uqStatisticalForwardProblemOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_computeSolution      << " = " << m_computeSolution
     << "\n" << m_option_computeCovariances   << " = " << m_computeCovariances
     << "\n" << m_option_computeCorrelations  << " = " << m_computeCorrelations
     << "\n" << m_option_dataOutputFileName   << " = " << m_dataOutputFileName;
  os << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef UQ_SFP_READS_SOLVER_OPTION
       << "\n" << m_option_solver << " = " << m_solverString
#endif
  os << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqStatisticalForwardProblemOptionsClass& obj)
{
  obj.print(os);

  return os;
}
