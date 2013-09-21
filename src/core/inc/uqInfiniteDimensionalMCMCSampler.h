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
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_INFINITE_SAMPLER__
#define __QUESO_INFINITE_SAMPLER__

// Boost includes
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

// Queso includes
#include <queso/uqLibMeshFunction.h>

// GSL includes
#include <gsl/gsl_rng.h>

// HDF5 includes
#include <H5Cpp.h>

class BaseEnvironmentClass;
class uqInfiniteDimensionalMeasureBase;
class uqInfiniteDimensionalLikelihoodBase;
class uqInfiniteDimensionalMCMCSamplerOptions;
namespace QUESO {

class uqInfiniteDimensionalMCMCSampler
{
public:
  /*!
   * Constructor
   */
  uqInfiniteDimensionalMCMCSampler(
      const BaseEnvironmentClass & env,
      const uqInfiniteDimensionalMeasureBase & prior,
      uqInfiniteDimensionalLikelihoodBase & llhd,
      uqInfiniteDimensionalMCMCSamplerOptions * ov);

  /*!
   * Destructor
   */
  ~uqInfiniteDimensionalMCMCSampler();

  /*!
   * Do one iteration of the Markov chain
   */
  void step();

  /*! 
   * Get the current value of the llhd
   */
  double llhd_val() const;

  /*!
   * Returns current acceptance probability
   */
  double acc_prob();

  /*!
   * Returns current average acceptance probability
   */
  double avg_acc_prob();

  /*!
   * Returns the current iteration number
   */
  unsigned int iteration() const;

  /*!
   * Returns a pointer to new sampler, with all the moments reset.
   */
  boost::shared_ptr<uqInfiniteDimensionalMCMCSampler> clone_and_reset() const;

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
  const uqInfiniteDimensionalMeasureBase & prior;

  // The negative log-likelihood functional.  Operates on functions.
  // Should be const?
  uqInfiniteDimensionalLikelihoodBase & llhd;

  // Aggregate options object
  uqInfiniteDimensionalMCMCSamplerOptions * m_ov;

  // The QUESO environment
  const BaseEnvironmentClass & m_env;

  // Pointer to the current physical state
  boost::shared_ptr<uqFunctionBase> current_physical_state;

  // Pointer to the current proposed state
  boost::shared_ptr<uqFunctionBase> proposed_physical_state;

  // Pointer to the current physical mean
  boost::shared_ptr<uqFunctionBase> current_physical_mean;

  // Pointer to the current physical variance
  boost::shared_ptr<uqFunctionBase> current_physical_var;

  // Stores the differences from the mean
  boost::shared_ptr<uqFunctionBase> _delta;

  // Stores a running sum-of-squares (kinda)
  boost::shared_ptr<uqFunctionBase> _M2;

  // A pointer to the random number generator to use.
  // Should probably use the one in the queso environment.
  gsl_rng *r;

  // Is the output file open?
  bool _outfile_open;
  boost::scoped_ptr<H5::H5File> _outfile;

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

  // Creates a scalar dataset to the hdf5 file
  void _create_scalar_dataset(const std::string & name);

  // Appends to a scalar dataset in the hdf5 file
  void _append_scalar_dataset(const std::string & name, double data);
};

}  // End namespace QUESO

#endif // __QUESO_INFINITE_SAMPLER__
