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

#include <queso/config_queso.h>
#include <queso/InfiniteDimensionalMCMCSamplerOptions.h>

// ODV = option default value
#define UQ_INF_DATA_OUTPUT_DIR_NAME_ODV "chain"
#define UQ_INF_DATA_OUTPUT_FILE_NAME_ODV "out.h5"
#define UQ_INF_NUM_ITERS_ODV 1000
#define UQ_INF_SAVE_FREQ_ODV 1
#define UQ_INF_RWMH_STEP_ODV 1e-2

namespace QUESO {

InfiniteDimensionalMCMCSamplerOptions::InfiniteDimensionalMCMCSamplerOptions()
  :
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser(new BoostInputOptionsParser()),
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  m_env(NULL)
{
  this->set_defaults();
  this->set_prefix("");
}

InfiniteDimensionalMCMCSamplerOptions::InfiniteDimensionalMCMCSamplerOptions(
    const BaseEnvironment& env,
    const char * prefix)
{
  queso_require_not_equal_to_msg(m_env->optionsInputFileName(), std::string(""), std::string("this constructor is incompatible with the abscense of an options input file"));

  this->set_defaults();
  this->parse(env, prefix);
}

void
InfiniteDimensionalMCMCSamplerOptions::checkOptions()
{
  queso_require_equal_to_msg(m_num_iters % m_save_freq, 0, "save frequency must divide number of iterations");
  queso_require_greater_msg(m_rwmh_step, 0, "random-walk Metropolis step size must be positive");
}

InfiniteDimensionalMCMCSamplerOptions::~InfiniteDimensionalMCMCSamplerOptions()
{
}

void InfiniteDimensionalMCMCSamplerOptions::print(std::ostream & os) const
{
  os << "\n" << this->m_option_dataOutputDirName << " = " << this->m_dataOutputDirName;
  os << "\n" << this->m_option_dataOutputFileName << " = " << this->m_dataOutputFileName;
  os << "\n" << this->m_option_num_iters << " = " << this->m_num_iters;
  os << "\n" << this->m_option_save_freq << " = " << this->m_save_freq;
  os << "\n" << this->m_option_rwmh_step << " = " << this->m_rwmh_step;
  os << std::endl;
  return;
}

std::ostream & operator<<(std::ostream & os,
    const InfiniteDimensionalMCMCSamplerOptions & obj)
{
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  obj.print(os);
  return os;
}

const BaseEnvironment&
InfiniteDimensionalMCMCSamplerOptions::env() const
{
  return *(this->m_env);
}


void InfiniteDimensionalMCMCSamplerOptions::set_defaults()
{
  m_dataOutputDirName = UQ_INF_DATA_OUTPUT_DIR_NAME_ODV;
  m_dataOutputFileName = UQ_INF_DATA_OUTPUT_FILE_NAME_ODV;
  m_num_iters = UQ_INF_NUM_ITERS_ODV;
  m_save_freq = UQ_INF_SAVE_FREQ_ODV;
  m_rwmh_step = UQ_INF_RWMH_STEP_ODV;
}


void InfiniteDimensionalMCMCSamplerOptions::set_prefix(const std::string& prefix)
{
  m_prefix = prefix+ "infmcmc_";

  m_option_help = m_prefix + "help";
  m_option_dataOutputDirName = m_prefix + "dataOutputDirName";
  m_option_dataOutputFileName = m_prefix + "dataOutputFileName";
  m_option_num_iters = m_prefix + "num_iters";
  m_option_save_freq = m_prefix + "save_freq";
  m_option_rwmh_step = m_prefix + "rwmh_step";
}


void InfiniteDimensionalMCMCSamplerOptions::parse(const BaseEnvironment& env, const std::string& prefix)
{
  m_env = &env;

  this->set_prefix(prefix);

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  m_parser.reset(new BoostInputOptionsParser(env.optionsInputFileName()));

  m_parser->registerOption(m_option_help, "produce help message for infinite dimensional sampler");
  m_parser->registerOption<std::string>(m_option_dataOutputDirName, m_dataOutputDirName, "name of data output dir");
  m_parser->registerOption<std::string>(m_option_dataOutputFileName, m_dataOutputFileName, "name of data output file (HDF5)");
  m_parser->registerOption<unsigned int>(m_option_num_iters, m_num_iters, "number of mcmc iterations to do");
  m_parser->registerOption<unsigned int>(m_option_save_freq, m_save_freq, "the frequency at which to save the chain state");
  m_parser->registerOption<double      >(m_option_rwmh_step, m_rwmh_step, "the step-size in the random-walk Metropolis proposal");;

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_dataOutputDirName,  m_dataOutputDirName);
  m_parser->getOption<std::string>(m_option_dataOutputFileName, m_dataOutputFileName);
  m_parser->getOption<unsigned int>(m_option_num_iters,         m_num_iters);
  m_parser->getOption<unsigned int>(m_option_save_freq,         m_save_freq);
  m_parser->getOption<double     >(m_option_rwmh_step,          m_rwmh_step);
#else
  m_dataOutputDirName = m_env->input()(m_option_dataOutputDirName, m_dataOutputDirName);
  m_dataOutputFileName = m_env->input()(m_option_dataOutputFileName, m_dataOutputFileName);
  m_num_iters = m_env->input()(m_option_num_iters, m_num_iters);
  m_save_freq = m_env->input()(m_option_save_freq, m_save_freq);
  m_rwmh_step = m_env->input()(m_option_rwmh_step, m_rwmh_step);
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

}  // End namespace QUESO
