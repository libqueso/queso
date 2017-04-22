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

#ifdef QUESO_HAVE_HDF5

#ifndef QUESO_INFINITE_SAMPLER_H
#define QUESO_INFINITE_SAMPLER_H

// Queso includes
#include <queso/SharedPtr.h>
#include <queso/Environment.h>
#include <queso/InfiniteDimensionalMeasureBase.h>
#include <queso/InfiniteDimensionalLikelihoodBase.h>
#include <queso/InfiniteDimensionalMCMCSamplerOptions.h>

// HDF5 includes
#include <hdf5.h>

namespace QUESO {

/*!
 * \file InfiniteDimensionalMCMCSampler.h
 * \brief Class representing the infinite dimensional Markov chain Monte Carlo sampler
 * \class InfiniteDimensionalMCMCSampler
 * \brief Class representing the infinite dimensional Markov chain Monte Carlo sampler
 */

class InfiniteDimensionalMCMCSampler
{
public:
  //! Construct an infinite dimensional MCMC chain with given \c prior, likelihood \c llhd, and options \c ov.
  InfiniteDimensionalMCMCSampler(
      const BaseEnvironment& env,
      InfiniteDimensionalMeasureBase & prior,
      InfiniteDimensionalLikelihoodBase & llhd,
      InfiniteDimensionalMCMCSamplerOptions * ov);

  //! Destructor
  ~InfiniteDimensionalMCMCSampler();

  //! Do one iteration of the Markov chain
  void step();

  //! Get the current value of the llhd
  double llhd_val() const;

  //! Returns current acceptance probability
  double acc_prob();

  //! Returns current average acceptance probability
  double avg_acc_prob();

  //! Returns the current iteration number
  unsigned int iteration() const;

  //! Returns a pointer to new sampler, with all the moments reset.
  SharedPtr<InfiniteDimensionalMCMCSampler>::Type clone_and_reset() const;

private:
  // Current iteration
  unsigned int _iteration;

  // The current value of the negative log-likelihood functional
  double _llhd_val;

  // The current acceptance probability
  double _acc_prob;

  // The current running mean of the acceptance probability
  double _avg_acc_prob;

  // The prior measure from which to draw
  InfiniteDimensionalMeasureBase & prior;

  // The negative log-likelihood functional.  Operates on functions.
  // Should be const?
  InfiniteDimensionalLikelihoodBase & llhd;

  // Aggregate options object
  InfiniteDimensionalMCMCSamplerOptions * m_ov;

  // The QUESO environment
  const BaseEnvironment& m_env;

  // Pointer to the current physical state
  SharedPtr<FunctionBase>::Type current_physical_state;

  // Pointer to the current proposed state
  SharedPtr<FunctionBase>::Type proposed_physical_state;

  // Pointer to the current physical mean
  SharedPtr<FunctionBase>::Type current_physical_mean;

  // Pointer to the current physical variance
  SharedPtr<FunctionBase>::Type current_physical_var;

  // Stores the differences from the mean
  SharedPtr<FunctionBase>::Type _delta;

  // Stores a running sum-of-squares (kinda)
  SharedPtr<FunctionBase>::Type _M2;

  // HDF5 file identifier
  hid_t _outfile;

  // HDF5 datasets
  hid_t _acc_dset;
  hid_t _avg_acc_dset;
  hid_t _neg_log_llhd_dset;
  hid_t _L2_norm_samples_dset;
  hid_t _L2_norm_mean_dset;
  hid_t _L2_norm_var_dset;

  /*!
   * Make a proposal from the prior using a standard random walk
   */
  void _propose();

  // Do a metropolis random walk step
  void _metropolis_hastings();

  // Increments the current iteration and updates the running mean and variance.
  void _update_moments();

  // Write the current state of the chain to disk
  void _write_state();

  // Creates a scalar dataset called \c name in \c _outfile file and returns it
  hid_t _create_scalar_dataset(const std::string & name);

  // Appends to a scalar dataset in the hdf5 file
  void _append_scalar_dataset(hid_t dset, double data);
};

}  // End namespace QUESO

#endif // QUESO_INFINITE_SAMPLER_H

#endif  // QUESO_HAVE_HDF5
