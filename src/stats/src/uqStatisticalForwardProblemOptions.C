//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <uqStatisticalForwardProblemOptions.h>
#include <uqMiscellaneous.h>

namespace QUESO {

// --------------------------------------------------
// SfpOptionsValuesClass---------------------------
// --------------------------------------------------

// Default constructor -----------------------------
SfpOptionsValuesClass::SfpOptionsValuesClass()
  :
  m_computeSolution     (UQ_SFP_COMPUTE_SOLUTION_ODV     ),
  m_computeCovariances  (UQ_SFP_COMPUTE_COVARIANCES_ODV  ),
  m_computeCorrelations (UQ_SFP_COMPUTE_CORRELATIONS_ODV ),
  m_dataOutputFileName  (UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV)
//m_dataOutputAllowedSet(),
#ifdef UQ_SFP_READS_SOLVER_OPTION
  m_solverString        (UQ_SFP_SOLVER_ODV               ),
#endif
{
}
// Copy constructor----------------------------------
SfpOptionsValuesClass::SfpOptionsValuesClass(const SfpOptionsValuesClass& src)
{
  this->copy(src);
}
// Destructor ---------------------------------------
SfpOptionsValuesClass::~SfpOptionsValuesClass()
{
}
// Set methods --------------------------------------
SfpOptionsValuesClass&
SfpOptionsValuesClass::operator=(const SfpOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}
// Private methods-----------------------------------
void
SfpOptionsValuesClass::copy(const SfpOptionsValuesClass& src)
{
  m_computeSolution      = src.m_computeSolution;
  m_computeCovariances   = src.m_computeCovariances;
  m_computeCorrelations  = src.m_computeCorrelations;
  m_dataOutputFileName   = src.m_dataOutputFileName;
  m_dataOutputAllowedSet = src.m_dataOutputAllowedSet;
#ifdef UQ_SFP_READS_SOLVER_OPTION
  m_solverString         = src.m_solverString;
#endif

  //m_mcOptionsValues      = src.m_mcOptionsValues;

  return;
}


// --------------------------------------------------
// StatisticalForwardProblemOptionsClass-----------
// --------------------------------------------------

// Default constructor -----------------------------
StatisticalForwardProblemOptionsClass::StatisticalForwardProblemOptionsClass(
  const BaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_ov                         (),
  m_prefix                     ((std::string)(prefix) + "fp_"   ),
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
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "StatisticalForwardProblemOptionsClass::constructor(1)",
                      "this constructor is incompatible with the absence of an options input file");
}

StatisticalForwardProblemOptionsClass::StatisticalForwardProblemOptionsClass(
  const BaseEnvironmentClass& env, 
  const char*                   prefix,
  const SfpOptionsValuesClass& alternativeOptionsValues)
  :
  m_ov                         (alternativeOptionsValues         ),
  m_prefix                     ((std::string)(prefix) + "fp_"    ),
  m_env                        (env),
  m_optionsDesc                (NULL),
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
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "StatisticalForwardProblemOptionsClass::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In StatisticalForwardProblemOptionsClass::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}
// Destructor --------------------------------------
StatisticalForwardProblemOptionsClass::~StatisticalForwardProblemOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
} 

// I/O methods -------------------------------------
void
StatisticalForwardProblemOptionsClass::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "StatisticalForwardProblemOptionsClass::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In StatisticalForwardProblemOptionsClass::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}
//--------------------------------------------------
void
StatisticalForwardProblemOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_computeSolution      << " = " << m_ov.m_computeSolution
     << "\n" << m_option_computeCovariances   << " = " << m_ov.m_computeCovariances
     << "\n" << m_option_computeCorrelations  << " = " << m_ov.m_computeCorrelations
     << "\n" << m_option_dataOutputFileName   << " = " << m_ov.m_dataOutputFileName;
  os << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef UQ_SFP_READS_SOLVER_OPTION
       << "\n" << m_option_solver << " = " << m_ov.m_solverString
#endif
  os << std::endl;

  return;
}

// Private methods ---------------------------------
void
StatisticalForwardProblemOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
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
//--------------------------------------------------
void
StatisticalForwardProblemOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
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

  if (m_env.allOptionsMap().count(m_option_computeCovariances)) {
    m_ov.m_computeCovariances = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeCovariances]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeCorrelations)) {
    m_ov.m_computeCorrelations = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeCorrelations]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
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

#ifdef UQ_SFP_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver)) {
    m_ov.m_solverString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_solver]).as<std::string>();
  }
#endif

  return;
}

// --------------------------------------------------
// Operator declared outside class definition ------
// --------------------------------------------------

std::ostream& operator<<(std::ostream& os, const StatisticalForwardProblemOptionsClass& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO
