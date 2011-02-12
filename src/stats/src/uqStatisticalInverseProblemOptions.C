//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqStatisticalInverseProblemOptions.h>
#include <uqMiscellaneous.h>

uqSipOptionsValuesClass::uqSipOptionsValuesClass()
  :
  m_computeSolution     (UQ_SIP_COMPUTE_SOLUTION_ODV     ),
  m_dataOutputFileName  (UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV)
//m_dataOutputAllowedSet(),
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_solverString        (UQ_SIP_SOLVER_ODV),
#endif
{
}

uqSipOptionsValuesClass::~uqSipOptionsValuesClass()
{
}

uqSipOptionsValuesClass::uqSipOptionsValuesClass(const uqSipOptionsValuesClass& src)
{
  this->copy(src);
}

uqSipOptionsValuesClass&
uqSipOptionsValuesClass::operator=(const uqSipOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}

void
uqSipOptionsValuesClass::copy(const uqSipOptionsValuesClass& src)
{
  m_computeSolution      = src.m_computeSolution;
  m_dataOutputFileName   = src.m_dataOutputFileName;
  m_dataOutputAllowedSet = src.m_dataOutputAllowedSet;
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_solverString         = src.m_solverString;
#endif
//m_mhOptionsValues      = src.m_mhOptionsValues;

  return;
}

uqStatisticalInverseProblemOptionsClass::uqStatisticalInverseProblemOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_ov                         (),
  m_prefix                     ((std::string)(prefix) + "ip_"),
  m_env                        (env),
  m_optionsDesc                (new po::options_description("Statistical Inverse Problem options")),
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet")
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              )
#endif
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemOptionsClass::constructor(1)",
                      "this constructor is incompatible with the abscense of an options input file");
}

uqStatisticalInverseProblemOptionsClass::uqStatisticalInverseProblemOptionsClass(
  const uqBaseEnvironmentClass&  env,
  const char*                    prefix,
  const uqSipOptionsValuesClass& alternativeOptionsValues)
  :
  m_ov                         (alternativeOptionsValues),
  m_prefix                     ((std::string)(prefix) + "ip_"),
  m_env                        (env),
  m_optionsDesc                (NULL),
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet")
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              )
#endif
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemOptionsClass::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqStatisticalInverseProblemOptionsClass::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

uqStatisticalInverseProblemOptionsClass::~uqStatisticalInverseProblemOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
} 

void
uqStatisticalInverseProblemOptionsClass::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemOptionsClass::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  //std::cout << "scan 000\n"
  //          << std::endl;
  getMyOptionValues              (*m_optionsDesc);
  //std::cout << "scan 001\n"
  //          << std::endl;

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqStatisticalInverseProblemOptionsClass::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
uqStatisticalInverseProblemOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                              "produce help message for statistical inverse problem")
    (m_option_computeSolution.c_str(),      po::value<bool       >()->default_value(UQ_SIP_COMPUTE_SOLUTION_ODV       ), "compute solution process"                            )
    (m_option_dataOutputFileName.c_str(),   po::value<std::string>()->default_value(UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file"                            )
    (m_option_dataOutputAllowedSet.c_str(), po::value<std::string>()->default_value(UQ_SIP_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file"         )
#ifdef UQ_SIP_READS_SOLVER_OPTION
    (m_option_solver.c_str(),               po::value<std::string>()->default_value(UQ_SIP_SOLVER_ODV                 ), "algorithm for calibration"                           )
#endif
  ;

  return;
}

void
uqStatisticalInverseProblemOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_computeSolution)) {
    m_ov.m_computeSolution = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeSolution]).as<bool>();
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

#ifdef UQ_SIP_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver)) {
    m_ov.m_solverString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_solver]).as<std::string>();
  }
#endif

  return;
}

void
uqStatisticalInverseProblemOptionsClass::print(std::ostream& os) const
{
  os << "\n" << m_option_computeSolution      << " = " << m_ov.m_computeSolution
     << "\n" << m_option_dataOutputFileName   << " = " << m_ov.m_dataOutputFileName;
  os << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef UQ_SIP_READS_SOLVER_OPTION
     << "\n" << m_option_solver << " = " << m_ov.m_solverString
#endif
  os << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqStatisticalInverseProblemOptionsClass& obj)
{
  obj.print(os);

  return os;
}
