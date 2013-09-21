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
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __UQ_INFMCMC_OPTIONS_H__
#define __UQ_INFMCMC_OPTIONS_H__

#include <queso/Environment.h>

namespace QUESO {

class uqInfiniteDimensionalMCMCSamplerOptions
{
public:
  /*!
   * Given prefix, read the input file for parameters named prefix_*
   */
  uqInfiniteDimensionalMCMCSamplerOptions(const uqBaseEnvironmentClass& env, const char* prefix);

  /*
   * Destructor
   */
 ~uqInfiniteDimensionalMCMCSamplerOptions();

  /*
   * Scans the input file for options prefixed with \c prefix
   */
  void scanOptionsValues();

  /*
   * Prints \c this to \c os
   */
  void print(std::ostream& os) const;

  /*
   * The prefix to look for in the input file
   */
  std::string m_prefix;

  // Name of the output file
  std::string m_dataOutputFileName;

  // The total number of iterations to do
  unsigned int m_num_iters;

  // The frequency at which to save the state of the chain
  unsigned int m_save_freq;

  // The proposal step size
  double m_rwmh_step;

  /*!
   * Returns the QUESO environment
   */
  const uqBaseEnvironmentClass & env() const;

private:
  void defineMyOptions(po::options_description& optionsDesc) const;
  void getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass & m_env;

  po::options_description* m_optionsDesc;
  std::string m_option_help;
  std::string m_option_dataOutputFileName;
  std::string m_option_num_iters;
  std::string m_option_save_freq;
  std::string m_option_rwmh_step;
};

std::ostream& operator<<(std::ostream& os, const uqInfiniteDimensionalMCMCSamplerOptions & opts);

}  // End namespace QUESO

#endif // __UQ_INFMCMC_OPTIONS_H__
