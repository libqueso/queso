//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
#include <queso/StatisticalForwardProblemOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

// --------------------------------------------------
// SfpOptionsValues---------------------------
// --------------------------------------------------

// Default constructor -----------------------------
SfpOptionsValues::SfpOptionsValues()
  :
    BaseInputOptions(),
    m_prefix                     ("fp_"),
    m_computeSolution     (UQ_SFP_COMPUTE_SOLUTION_ODV     ),
    m_computeCovariances  (UQ_SFP_COMPUTE_COVARIANCES_ODV  ),
    m_computeCorrelations (UQ_SFP_COMPUTE_CORRELATIONS_ODV ),
    m_dataOutputFileName  (UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV),
    //m_dataOutputAllowedSet(),
    m_option_help                (m_prefix + "help"                ),
    m_option_computeSolution     (m_prefix + "computeSolution"     ),
    m_option_computeCovariances  (m_prefix + "computeCovariances"  ),
    m_option_computeCorrelations (m_prefix + "computeCorrelations" ),
    m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
    m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet")
#ifdef UQ_SFP_READS_SOLVER_OPTION
    m_option_solver              (m_prefix + "solver"              ),
    m_solverString        (UQ_SFP_SOLVER_ODV               )
#endif
{
}

SfpOptionsValues::SfpOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
    BaseInputOptions(env),
    m_prefix                     ((std::string)(prefix) + "fp_"),
    m_computeSolution     (UQ_SFP_COMPUTE_SOLUTION_ODV     ),
    m_computeCovariances  (UQ_SFP_COMPUTE_COVARIANCES_ODV  ),
    m_computeCorrelations (UQ_SFP_COMPUTE_CORRELATIONS_ODV ),
    m_dataOutputFileName  (UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV),
    //m_dataOutputAllowedSet(),
    m_option_help                (m_prefix + "help"                ),
    m_option_computeSolution     (m_prefix + "computeSolution"     ),
    m_option_computeCovariances  (m_prefix + "computeCovariances"  ),
    m_option_computeCorrelations (m_prefix + "computeCorrelations" ),
    m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
    m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet")
#ifdef UQ_SFP_READS_SOLVER_OPTION
    m_option_solver              (m_prefix + "solver"              ),
    m_solverString        (UQ_SFP_SOLVER_ODV               )
#endif
{
}

// Copy constructor----------------------------------
SfpOptionsValues::SfpOptionsValues(const SfpOptionsValues& src)
{
  this->copy(src);
}
// Destructor ---------------------------------------
SfpOptionsValues::~SfpOptionsValues()
{
}
// Set methods --------------------------------------
SfpOptionsValues&
SfpOptionsValues::operator=(const SfpOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}

void
SfpOptionsValues::defineOptions()
{
  (*m_optionsDescription).add_options()
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
}

void
SfpOptionsValues::getOptionValues()
{
  if (m_env->allOptionsMap().count(m_option_help)) {
    if (m_env->subDisplayFile()) {
      *m_env->subDisplayFile() << (*m_optionsDescription)
                              << std::endl;
    }
  }

  if (m_env->allOptionsMap().count(m_option_computeSolution)) {
    m_computeSolution = ((const po::variable_value&) m_env->allOptionsMap()[m_option_computeSolution]).as<bool>();
  }

  if (m_env->allOptionsMap().count(m_option_computeCovariances)) {
    m_computeCovariances = ((const po::variable_value&) m_env->allOptionsMap()[m_option_computeCovariances]).as<bool>();
  }

  if (m_env->allOptionsMap().count(m_option_computeCorrelations)) {
    m_computeCorrelations = ((const po::variable_value&) m_env->allOptionsMap()[m_option_computeCorrelations]).as<bool>();
  }

  if (m_env->allOptionsMap().count(m_option_dataOutputFileName)) {
    m_dataOutputFileName = ((const po::variable_value&) m_env->allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env->allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env->allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

#ifdef UQ_SFP_READS_SOLVER_OPTION
  if (m_env->allOptionsMap().count(m_option_solver)) {
    m_solverString = ((const po::variable_value&) m_env->allOptionsMap()[m_option_solver]).as<std::string>();
  }
#endif
}

// Private methods-----------------------------------
void
SfpOptionsValues::copy(const SfpOptionsValues& src)
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
// StatisticalForwardProblemOptions-----------
// --------------------------------------------------

// Default constructor -----------------------------
StatisticalForwardProblemOptions::StatisticalForwardProblemOptions(
  const BaseEnvironment& env,
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
  queso_deprecated();

  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "StatisticalForwardProblemOptions::constructor(1)",
                      "this constructor is incompatible with the absence of an options input file");
}

StatisticalForwardProblemOptions::StatisticalForwardProblemOptions(
  const BaseEnvironment& env,
  const char*                   prefix,
  const SfpOptionsValues& alternativeOptionsValues)
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
  queso_deprecated();

  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "StatisticalForwardProblemOptions::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In StatisticalForwardProblemOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}
// Destructor --------------------------------------
StatisticalForwardProblemOptions::~StatisticalForwardProblemOptions()
{
  queso_deprecated();

  if (m_optionsDesc) delete m_optionsDesc;
}

// I/O methods -------------------------------------
void
StatisticalForwardProblemOptions::scanOptionsValues()
{
  queso_deprecated();

  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "StatisticalForwardProblemOptions::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In StatisticalForwardProblemOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}
//--------------------------------------------------
void
StatisticalForwardProblemOptions::print(std::ostream& os) const
{
  queso_deprecated();

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
StatisticalForwardProblemOptions::defineMyOptions(po::options_description& optionsDesc) const
{
  queso_deprecated();

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
StatisticalForwardProblemOptions::getMyOptionValues(po::options_description& optionsDesc)
{
  queso_deprecated();

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

std::ostream& operator<<(std::ostream& os, const StatisticalForwardProblemOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
