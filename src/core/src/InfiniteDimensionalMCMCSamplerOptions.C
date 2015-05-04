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

#include <queso/InfiniteDimensionalMCMCSamplerOptions.h>

// ODV = option default value
#define UQ_INF_DATA_OUTPUT_DIR_NAME_ODV "chain"
#define UQ_INF_DATA_OUTPUT_FILE_NAME_ODV "out.h5"
#define UQ_INF_NUM_ITERS_ODV 1000
#define UQ_INF_SAVE_FREQ_ODV 1
#define UQ_INF_RWMH_STEP_ODV 1e-2

namespace QUESO {

InfiniteDimensionalMCMCSamplerOptions::InfiniteDimensionalMCMCSamplerOptions(
    const BaseEnvironment& env,
    const char * prefix)
  : BaseInputOptions(&env),
    m_prefix((std::string)(prefix) + "infmcmc_"),
    m_dataOutputDirName(UQ_INF_DATA_OUTPUT_DIR_NAME_ODV),
    m_dataOutputFileName(UQ_INF_DATA_OUTPUT_FILE_NAME_ODV),
    m_num_iters(UQ_INF_NUM_ITERS_ODV),
    m_save_freq(UQ_INF_SAVE_FREQ_ODV),
    m_rwmh_step(UQ_INF_RWMH_STEP_ODV),
    m_env(env),
    m_option_help(m_prefix + "help"),
    m_option_dataOutputDirName(m_prefix + "dataOutputDirName"),
    m_option_dataOutputFileName(m_prefix + "dataOutputFileName"),
    m_option_num_iters(m_prefix + "num_iters"),
    m_option_save_freq(m_prefix + "save_freq"),
    m_option_rwmh_step(m_prefix + "rwmh_step")
{
  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the abscense of an options input file");
}

InfiniteDimensionalMCMCSamplerOptions::~InfiniteDimensionalMCMCSamplerOptions()
{
}

void InfiniteDimensionalMCMCSamplerOptions::defineOptions()
{
  (*m_optionsDescription).add_options()
    (m_option_help.c_str(), "produce help message for infinite dimensional sampler")
    (m_option_dataOutputDirName.c_str(), po::value<std::string>()->default_value(UQ_INF_DATA_OUTPUT_DIR_NAME_ODV), "name of data output dir")
    (m_option_dataOutputFileName.c_str(), po::value<std::string>()->default_value(UQ_INF_DATA_OUTPUT_FILE_NAME_ODV), "name of data output file (HDF5)")
    (m_option_num_iters.c_str(), po::value<int>()->default_value(UQ_INF_NUM_ITERS_ODV), "number of mcmc iterations to do")
    (m_option_save_freq.c_str(), po::value<int>()->default_value(UQ_INF_SAVE_FREQ_ODV), "the frequency at which to save the chain state")
    (m_option_rwmh_step.c_str(), po::value<double>()->default_value(UQ_INF_RWMH_STEP_ODV), "the step-size in the random-walk Metropolis proposal");
}

void InfiniteDimensionalMCMCSamplerOptions::getOptionValues()
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*m_optionsDescription)
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputDirName)) {
    this->m_dataOutputDirName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputDirName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    this->m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_num_iters)) {
    this->m_num_iters = ((const po::variable_value&) m_env.allOptionsMap()[m_option_num_iters]).as<int>();
  }

  if (m_env.allOptionsMap().count(m_option_save_freq)) {
    this->m_save_freq = ((const po::variable_value&) m_env.allOptionsMap()[m_option_save_freq]).as<int>();
  }

  if (m_env.allOptionsMap().count(m_option_rwmh_step)) {
    this->m_rwmh_step = ((const po::variable_value&) m_env.allOptionsMap()[m_option_rwmh_step]).as<double>();
  }

  return;
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
  obj.print(os);
  return os;
}

const BaseEnvironment&
InfiniteDimensionalMCMCSamplerOptions::env() const
{
  return this->m_env;
}

}  // End namespace QUESO
