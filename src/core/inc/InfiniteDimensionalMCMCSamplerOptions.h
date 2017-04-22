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

#ifndef UQ_INFMCMC_OPTIONS_H
#define UQ_INFMCMC_OPTIONS_H

#include <queso/Environment.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

namespace QUESO {

/*!
 * \file InfiniteDimensionalMCMCSamplerOptions.h
 * \brief This class defines the options that specify the behaviour of the MCMC sampler
 *
 * \class InfiniteDimensionalMCMCSamplerOptions
 * \brief This class defines the options that specify the behaviour of the MCMC sampler
 */

class InfiniteDimensionalMCMCSamplerOptions
{
public:
  //! Given prefix, read the input file for parameters named prefix_*
  InfiniteDimensionalMCMCSamplerOptions(const BaseEnvironment& env, const char* prefix);

  //! Destructor
  virtual ~InfiniteDimensionalMCMCSamplerOptions();

  //! Prints \c this to \c os
  void print(std::ostream& os) const;

  //! The prefix to look for in the input file
  std::string m_prefix;

  //! Name of the output dir to save infinite dimensional output files to.
  /*!
   * For example, if \c m_dataOutputDirName is set to 'outputData/chain' then
   * a series of folders will be created, one for each QUESO subenvironment,
   * called 'outputData/chain0', 'outputData/chain1', etc.  Inside these
   * folders, output data from the infinite dimensional chain will be saved as
   * an HDF5 file.
   */
  std::string m_dataOutputDirName;

  //! Name of the HDF5 output file to store chain statistics.
  std::string m_dataOutputFileName;

  //! The total number of iterations to do
  unsigned int m_num_iters;

  //! The frequency at which to save the state of the chain
  unsigned int m_save_freq;

  //! The proposal step size
  double m_rwmh_step;

  //! Returns the QUESO environment
  const BaseEnvironment& env() const;

private:
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  BoostInputOptionsParser * m_parser;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  const BaseEnvironment& m_env;

  std::string m_option_help;
  std::string m_option_dataOutputDirName;
  std::string m_option_dataOutputFileName;
  std::string m_option_num_iters;
  std::string m_option_save_freq;
  std::string m_option_rwmh_step;

  void checkOptions();

  friend std::ostream & operator<<(std::ostream & os,
      const InfiniteDimensionalMCMCSamplerOptions & opts);
};

}  // End namespace QUESO

#endif // UQ_INFMCMC_OPTIONS_H
