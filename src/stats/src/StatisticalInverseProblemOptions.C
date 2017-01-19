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

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/StatisticalInverseProblemOptions.h>
#include <queso/Miscellaneous.h>

// -------------------------------------------------
// SipOptionsValues--------------------------
// -------------------------------------------------

namespace QUESO {

// Default constructor -----------------------------
SipOptionsValues::SipOptionsValues()
  :
  m_prefix("ip_"),
  m_help(UQ_SIP_HELP),
  m_computeSolution     (UQ_SIP_COMPUTE_SOLUTION_ODV     ),
  m_dataOutputFileName  (UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV),
  m_seedWithMAPEstimator(UQ_SIP_SEEDWITHMAPESTIMATOR),
  m_useOptimizerMonitor(UQ_SIP_USEOPTIMIZERMONITOR),
//m_dataOutputAllowedSet(),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser(),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              ),
  m_solverString        (UQ_SIP_SOLVER_ODV)
#endif
  m_option_seedWithMAPEstimator(m_prefix + "seedWithMAPEstimator"),
  m_option_useOptimizerMonitor(m_prefix + "useOptimizerMonitor")
{
}

SipOptionsValues::SipOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
  m_prefix((std::string)(prefix) + "ip_"),
  m_help(UQ_SIP_HELP),
  m_computeSolution     (UQ_SIP_COMPUTE_SOLUTION_ODV     ),
  m_dataOutputFileName  (UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV),
  m_seedWithMAPEstimator(UQ_SIP_SEEDWITHMAPESTIMATOR),
  m_useOptimizerMonitor(UQ_SIP_USEOPTIMIZERMONITOR),
//m_dataOutputAllowedSet(),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              ),
  m_solverString        (UQ_SIP_SOLVER_ODV)
#endif
  m_option_seedWithMAPEstimator(m_prefix + "seedWithMAPEstimator"),
  m_option_useOptimizerMonitor(m_prefix + "useOptimizerMonitor")
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser->registerOption<std::string>(m_option_help, UQ_SIP_HELP, "produce help message for statistical inverse problem");
  m_parser->registerOption<bool       >(m_option_computeSolution,      UQ_SIP_COMPUTE_SOLUTION_ODV       , "compute solution process"                            );
  m_parser->registerOption<std::string>(m_option_dataOutputFileName,   UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV  , "name of data output file"                            );
  m_parser->registerOption<std::string>(m_option_dataOutputAllowedSet, UQ_SIP_DATA_OUTPUT_ALLOWED_SET_ODV, "subEnvs that will write to data output file"         );
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_parser->registerOption<std::string>(m_option_solver,               UQ_SIP_SOLVER_ODV                 , "algorithm for calibration"                           );
#endif
  m_parser->registerOption<bool       >(m_option_seedWithMAPEstimator, UQ_SIP_SEEDWITHMAPESTIMATOR       , "toggle for seeding chain at MAP"                     );
  m_parser->registerOption<bool       >(m_option_useOptimizerMonitor,  UQ_SIP_USEOPTIMIZERMONITOR        , "toggle for using optimizer monitor (prints diagnostics");

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help, m_help);
  m_parser->getOption<bool       >(m_option_computeSolution,      m_computeSolution);
  m_parser->getOption<std::string>(m_option_dataOutputFileName,   m_dataOutputFileName);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet, m_dataOutputAllowedSet);
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_parser->getOption<std::string>(m_option_solver,               m_solver);
#endif
  m_parser->getOption<bool       >(m_option_seedWithMAPEstimator, m_seedWithMAPEstimator);
  m_parser->getOption<bool       >(m_option_useOptimizerMonitor,  m_useOptimizerMonitor);
#else
  m_help = env->input()(m_option_help, UQ_SIP_HELP);
  m_computeSolution = env->input()(m_option_computeSolution, UQ_SIP_COMPUTE_SOLUTION_ODV);
  m_dataOutputFileName = env->input()(m_option_dataOutputFileName, UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV);

  // UQ_SIP_DATA_OUTPUT_ALLOWED_SET_ODV is the empty set (string) by default
  unsigned int size = env->input().vector_variable_size(m_option_dataOutputAllowedSet);
  for (unsigned int i = 0; i < size; i++) {
    // We default to empty set, so the default values are actually never
    // used here
    unsigned int allowed = env->input()(m_option_dataOutputAllowedSet, i, i);
    m_dataOutputAllowedSet.insert(allowed);
  }

#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_solver = env->input()(m_option_solver, UQ_SIP_SOLVER_ODV);
#endif
  m_seedWithMAPEstimator = env->input()(m_option_seedWithMAPEstimator, UQ_SIP_SEEDWITHMAPESTIMATOR);
  m_useOptimizerMonitor = env->input()(m_option_useOptimizerMonitor,  UQ_SIP_USEOPTIMIZERMONITOR);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}
// Copy constructor - -----------------------------
SipOptionsValues::SipOptionsValues(const SipOptionsValues& src)
{
  this->copy(src);
}
// Destructor ---------------------------------------
SipOptionsValues::~SipOptionsValues()
{
}

// Set methods --------------------------------------
SipOptionsValues&
SipOptionsValues::operator=(const SipOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}
// Private methods-----------------------------------
void
SipOptionsValues::checkOptions()
{
  // Nothing
}

void
SipOptionsValues::copy(const SipOptionsValues& src)
{
  m_computeSolution      = src.m_computeSolution;
  m_dataOutputFileName   = src.m_dataOutputFileName;
  m_dataOutputAllowedSet = src.m_dataOutputAllowedSet;
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_solverString         = src.m_solverString;
#endif
  m_seedWithMAPEstimator = src.m_seedWithMAPEstimator;
  m_useOptimizerMonitor = src.m_useOptimizerMonitor;
}

std::ostream &
operator<<(std::ostream& os, const SipOptionsValues & obj)
{
  os << "\n" << obj.m_option_computeSolution      << " = " << obj.m_computeSolution
     << "\n" << obj.m_option_dataOutputFileName   << " = " << obj.m_dataOutputFileName;
  os << "\n" << obj.m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = obj.m_dataOutputAllowedSet.begin(); setIt != obj.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef UQ_SIP_READS_SOLVER_OPTION
     << "\n" << obj.m_option_solver << " = " << obj.m_solverString
#endif
  os << std::endl;
  return os;
}

// -------------------------------------------------
// StatisticalInverseProblemOptions----------
// -------------------------------------------------

// Default constructor -----------------------------
StatisticalInverseProblemOptions::StatisticalInverseProblemOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
:
  m_ov                         (),
  m_prefix                     ((std::string)(prefix) + "ip_"),
  m_env                        (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                (new boost::program_options::options_description("Statistical Inverse Problem options")),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              )
#endif
  m_option_seedWithMAPEstimator(m_prefix + "seedWithMAPEstimator"),
  m_option_useOptimizerMonitor(m_prefix + "useOptimizerMonitor")
{
  queso_deprecated();
  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the absence of an options input file"));
}

// Constructor 2------------------------------------
StatisticalInverseProblemOptions::StatisticalInverseProblemOptions(
  const BaseEnvironment&  env,
  const char*                    prefix,
  const SipOptionsValues& alternativeOptionsValues)
:
  m_ov                         (alternativeOptionsValues),
  m_prefix                     ((std::string)(prefix) + "ip_"),
  m_env                        (env),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_optionsDesc                (),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
#ifdef UQ_SIP_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              )
#endif
  m_option_seedWithMAPEstimator(m_prefix + "seedWithMAPEstimator"),
  m_option_useOptimizerMonitor(m_prefix + "useOptimizerMonitor")
{
  queso_deprecated();

  queso_require_equal_to_msg(m_env.optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the existence of an options input file"));

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In StatisticalInverseProblemOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}
// Destructor --------------------------------------
StatisticalInverseProblemOptions::~StatisticalInverseProblemOptions()
{
  queso_deprecated();
}

// I/O methods --------------------------------------
void
StatisticalInverseProblemOptions::scanOptionsValues()
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In StatisticalInverseProblemOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

// --------------------------------------------------
void
StatisticalInverseProblemOptions::print(std::ostream& os) const
{
  queso_deprecated();

  os << "\n" << m_option_computeSolution      << " = " << m_ov.m_computeSolution
     << "\n" << m_option_dataOutputFileName   << " = " << m_ov.m_dataOutputFileName;
  os << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef UQ_SIP_READS_SOLVER_OPTION
     << "\n" << m_option_solver << " = " << m_ov.m_solverString
#endif
  os << "\n" << m_option_seedWithMAPEstimator << " = " << m_ov.m_seedWithMAPEstimator
     << "\n" << m_option_useOptimizerMonitor << " = " << m_ov.m_useOptimizerMonitor;
  os << std::endl;

  return;
}

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
void
StatisticalInverseProblemOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                              "produce help message for statistical inverse problem")
    (m_option_computeSolution.c_str(),      boost::program_options::value<bool       >()->default_value(UQ_SIP_COMPUTE_SOLUTION_ODV       ), "compute solution process"                            )
    (m_option_dataOutputFileName.c_str(),   boost::program_options::value<std::string>()->default_value(UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV  ), "name of data output file"                            )
    (m_option_dataOutputAllowedSet.c_str(), boost::program_options::value<std::string>()->default_value(UQ_SIP_DATA_OUTPUT_ALLOWED_SET_ODV), "subEnvs that will write to data output file"         )
#ifdef UQ_SIP_READS_SOLVER_OPTION
    (m_option_solver.c_str(),               boost::program_options::value<std::string>()->default_value(UQ_SIP_SOLVER_ODV                 ), "algorithm for calibration"                           )
#endif
    (m_option_seedWithMAPEstimator.c_str(), boost::program_options::value<bool>()->default_value(UQ_SIP_SEEDWITHMAPESTIMATOR), "toggle optimisation of posterior pdf before sampling")
    (m_option_useOptimizerMonitor.c_str(), boost::program_options::value<bool>()->default_value(UQ_SIP_USEOPTIMIZERMONITOR), "toggle printing of optimisation progress when seed with MAP estimator")
  ;
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
//--------------------------------------------------
void
StatisticalInverseProblemOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_computeSolution)) {
    m_ov.m_computeSolution = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_computeSolution]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
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

#ifdef UQ_SIP_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver)) {
    m_ov.m_solverString = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_solver]).as<std::string>();
  }
#endif

  if (m_env.allOptionsMap().count(m_option_seedWithMAPEstimator)) {
    m_ov.m_seedWithMAPEstimator = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_seedWithMAPEstimator]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_useOptimizerMonitor)) {
    m_ov.m_useOptimizerMonitor = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_useOptimizerMonitor]).as<bool>();
  }
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

// --------------------------------------------------
// Operator declared outside class definition ------
// --------------------------------------------------

std::ostream& operator<<(std::ostream& os, const StatisticalInverseProblemOptions& obj)
{
  queso_deprecated();
  obj.print(os);

  return os;
}

}  // End namespace QUESO
