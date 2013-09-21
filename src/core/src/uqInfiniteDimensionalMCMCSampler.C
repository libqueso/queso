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

#include <cmath>
#include <algorithm>
#include <boost/math/constants/constants.hpp>

// QUESO includes
// #include <uqLibMeshFunction.h>
#include <queso/uqInfiniteDimensionalMeasureBase.h>

#include <queso/uqInfiniteDimensionalLikelihoodBase.h>
// #include <forward_solver/forward_solver.h>
#include <queso/uqInfiniteDimensionalMCMCSampler.h>
#include <queso/uqInfiniteDimensionalMCMCSamplerOptions.h>

// libmesh includes
// #include <libmesh/equation_systems.h>
// #include <libmesh/linear_implicit_system.h>

uqInfiniteDimensionalMCMCSampler::uqInfiniteDimensionalMCMCSampler(
    const uqBaseEnvironmentClass & env,
    const uqInfiniteDimensionalMeasureBase & prior,
    uqInfiniteDimensionalLikelihoodBase & llhd,
    uqInfiniteDimensionalMCMCSamplerOptions * ov)
  : prior(prior),
    llhd(llhd),
    m_env(env),
    current_physical_state(prior.draw()),
    proposed_physical_state(prior.draw()),
    current_physical_mean(prior.draw()),
    current_physical_var(prior.draw()),
    _delta(prior.draw()),
    _M2(prior.draw())
{
  if (ov != NULL) {
    this->m_ov = ov;
  }

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering uqInfiniteDimensionalMCMCSampler class" << std::endl;
#endif
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqInfiniteDimensionalMCMCSampler::constructor()"
                            << ": prefix = " << this->m_ov->m_prefix
                            << ", m_env.optionsInputFileName() = " << this->m_env.optionsInputFileName()
                            << std::endl;
  }

  if (m_env.optionsInputFileName() != "") {
    this->m_ov->scanOptionsValues();
  }

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In uqInfiniteDimensionalMCMCSamplerOptions,"
            << " finished scanning options" << std::endl;
#endif

  this->_iteration = 0;
  this->_acc_prob = 0.0;
  this->_avg_acc_prob = 0.0;
  r = gsl_rng_alloc(gsl_rng_taus2);

  // Outfile is not yet open
  this->_outfile_open = false;

  // Zero out these guys.  There's probably a better way of doing this than
  // calling prior.draw() at the start, but creation of a Sampler object
  // should only be done a O(1) times anyway.
  this->_delta->zero();
  this->_M2->zero();
  this->current_physical_mean->zero();
  this->current_physical_var->zero();

  // uqLibMeshFunction & p = libmesh_cast_ref<uqLibMeshFunction &>(*(this->current_physical_state)); 

  // Seed the sampler at the truth
  // for (unsigned int ii = 0; ii < 513; ii++) {
  //   p.equation_systems->get_system<ExplicitSystem>("Function").solution->set(ii, std::cos(2.0 * boost::math::constants::pi<double>() * ii / 512.0));
  // }
  // p.equation_systems->get_system<ExplicitSystem>("Function").solution->close();
  // std::string seedname("seedfn");
  // p.save_function(seedname);

  // Initialise cost function
  this->_llhd_val = this->llhd.evaluate(*(this->current_physical_state));

  std::cout << "got to sampler ctor" << std::endl;
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(15);
}

uqInfiniteDimensionalMCMCSampler::~uqInfiniteDimensionalMCMCSampler()
{
  std::cout << "got to sampler dtor" << std::endl;
}

void uqInfiniteDimensionalMCMCSampler::_propose()
{
  const double rwmh_step_sq = (this->m_ov->m_rwmh_step * this->m_ov->m_rwmh_step);
  const double coeff = std::sqrt(1.0 - rwmh_step_sq);

  boost::shared_ptr<uqFunctionBase> p(prior.draw());

  this->proposed_physical_state->zero();
  this->proposed_physical_state->add(coeff, *(this->current_physical_state));
  this->proposed_physical_state->add(this->m_ov->m_rwmh_step, *p);
}

void uqInfiniteDimensionalMCMCSampler::_metropolis_hastings()
{
  // Downcast since we know we're dealing with a libmesh function
  // uqLibMeshFunction & p = static_cast<uqLibMeshFunction &>(*(this->proposed_physical_state));
  // uqLibMeshFunction & q = static_cast<uqLibMeshFunction &>(*(this->current_physical_state));
  // Evaluate the likelihood at the proposed state
  double proposed_llhd = this->llhd.evaluate(*(this->proposed_physical_state));
  // double current_llhd = this->llhd.evaluate(q);
  // std::cout << "llhd of  current is: " << current_llhd << std::endl;
  // std::cout << "llhd of proposal is: " << proposed_llhd << std::endl;

  double diff = this->_llhd_val - proposed_llhd;
  double alpha = std::min(1.0, std::exp(diff));
  double rand = gsl_rng_uniform(this->r);
  if (rand < alpha) {
    // Accept
    // std::cout << "accepted" << std::endl;
    this->current_physical_state.swap(this->proposed_physical_state);
    this->proposed_physical_state->zero();
    this->_acc_prob = alpha;
    this->_llhd_val = proposed_llhd;
  }
  else {
    // std::cout << "rejected" << std::endl;
    this->_acc_prob = 0.0;
  }
}

void uqInfiniteDimensionalMCMCSampler::_update_moments()
{
  // Increment the current iteration number and update the running mean and
  // variance.
  this->_iteration += 1;

  // Make _delta contain the difference from the mean
  this->_delta->zero();
  this->_delta->add(1.0, *(this->current_physical_state));
  this->_delta->add(-1.0, *(this->current_physical_mean));

  // Scale _delta to update the mean field
  this->current_physical_mean->add(1.0 / this->iteration(), *(this->_delta));

  // Update running sum-of-squares
  boost::shared_ptr<uqFunctionBase> temp_ptr(this->_delta->zero_clone());
  // uqLibMeshFunction & temp = static_cast<uqLibMeshFunction &>(*temp_ptr); 

  temp_ptr->pointwise_mult(*(this->_delta), *(this->current_physical_state));
  this->_M2->add(1.0, *temp_ptr);
  temp_ptr->pointwise_mult(*(this->_delta), *(this->current_physical_mean));
  this->_M2->add(-1.0, *temp_ptr);

  if (this->iteration() > 1) {
    this->current_physical_var->zero();
    this->current_physical_var->add(1.0 / (this->iteration() - 1), *(this->_M2));
  }

  // Update acceptance rate
  double delta_acc = this->acc_prob() - this->avg_acc_prob();
  this->_avg_acc_prob += delta_acc / this->iteration();
}

void uqInfiniteDimensionalMCMCSampler::step()
{
  this->_propose();
  this->_metropolis_hastings();
  this->_update_moments();
  
  if (this->_iteration % this->m_ov->m_save_freq == 0) {
    this->_write_state();
  }
}

void uqInfiniteDimensionalMCMCSampler::_create_scalar_dataset(const std::string & name)
{
  hsize_t      dims[1]  = {1};  // dataset dimensions at creation
  hsize_t      maxdims[1];
  maxdims[0] = this->m_ov->m_num_iters / this->m_ov->m_save_freq;
  const int rank = 1;
  H5::DataSpace mspace1(rank, dims, maxdims);

  H5::DSetCreatPropList cparms;

  hsize_t      chunk_dims[1] ={1};
  cparms.setChunk( rank, chunk_dims );

  double fill_val = 0.0;
  cparms.setFillValue(H5::PredType::NATIVE_DOUBLE, &fill_val);

  const H5std_string DATASET_NAME(name);
  H5::DataSet dataset2 = this->_outfile->createDataSet(DATASET_NAME, H5::PredType::NATIVE_DOUBLE, mspace1, cparms);
}

void uqInfiniteDimensionalMCMCSampler::_append_scalar_dataset(const std::string & name, double data)
{
  const H5std_string DATASET_NAME(name);
  H5::DataSet dataset = this->_outfile->openDataSet(DATASET_NAME);

  hsize_t      dims[1]  = {1};  // dataset dimensions at creation
  hsize_t      maxdims[1];
  maxdims[0] = this->m_ov->m_num_iters / this->m_ov->m_save_freq;
  const int rank = 1;
  H5::DataSpace mspace1(rank, dims, maxdims);

  hsize_t      size[1];
  size[0]   = this->_iteration / this->m_ov->m_save_freq;
  dataset.extend(size);

  H5::DataSpace fspace1 = dataset.getSpace();
  hsize_t     offset[1];
  offset[0] = (this->_iteration / this->m_ov->m_save_freq) - 1;
  hsize_t      dims1[1] = {1};            /* data1 dimensions */
  fspace1.selectHyperslab(H5S_SELECT_SET, dims1, offset);

  dataset.write(&data, H5::PredType::NATIVE_DOUBLE, mspace1, fspace1);
}

void uqInfiniteDimensionalMCMCSampler::_write_state()
{
  if (!this->_outfile_open) {
    // std::cout << "opening new file" << std::endl;
    const H5std_string file_name(this->m_ov->m_dataOutputFileName.c_str());
    this->_outfile.reset(new H5::H5File(file_name, H5F_ACC_TRUNC));
    // std::cout << "opened new file" << std::endl;
    this->_outfile_open = true;
    
    this->_create_scalar_dataset("acc");
    this->_create_scalar_dataset("avg_acc");
    this->_create_scalar_dataset("neg_log_llhd");
    this->_create_scalar_dataset("L2_norm_samples");
    this->_create_scalar_dataset("L2_norm_mean");
    this->_create_scalar_dataset("L2_norm_var");
  }
  // else {
  //   std::cout << "file already open" << std::endl;
  // }

  this->_append_scalar_dataset("acc", this->acc_prob());
  this->_append_scalar_dataset("avg_acc", this->avg_acc_prob());
  this->_append_scalar_dataset("neg_log_llhd", this->_llhd_val);
  this->_append_scalar_dataset("L2_norm_samples", this->current_physical_state->L2_norm());
  this->_append_scalar_dataset("L2_norm_mean", this->current_physical_mean->L2_norm());
  this->_append_scalar_dataset("L2_norm_var", this->current_physical_var->L2_norm());

  // Now to write the fields.  It appears to be a pain in the arse to write a
  // method to spit this out to HDF5 format.  Also, this won't scale to
  // non-uniform finite element meshes.  Therefore, I'm going to spit out an
  // ExodusII file for each sample, average and variance.  Got disk space?
  std::ostringstream curr_iter;
  curr_iter << this->iteration();

  std::string sample_name("sample_");
  sample_name += curr_iter.str();
  sample_name += ".vtk";
  this->current_physical_state->save_function(sample_name);

  std::string mean_name("mean_");
  mean_name += curr_iter.str();
  mean_name += ".vtk";
  this->current_physical_mean->save_function(mean_name);

  std::string var_name("var_");
  var_name += curr_iter.str();
  var_name += ".vtk";
  this->current_physical_var->save_function(var_name);

  // std::string qoi_name("qoi_");
  // qoi_name += curr_iter.str();
  // qoi_name += ".vtk";
  // this->llhd.get_qoi().save_function(qoi_name);
}

boost::shared_ptr<uqInfiniteDimensionalMCMCSampler> uqInfiniteDimensionalMCMCSampler::clone_and_reset() const
{
  // Set up a clone
  boost::shared_ptr<uqInfiniteDimensionalMCMCSampler> clone(new uqInfiniteDimensionalMCMCSampler(this->m_env, this->prior, this->llhd, this->m_ov));

  // Copy the state.
  clone->current_physical_state = this->current_physical_state;
  clone->proposed_physical_state = this->proposed_physical_state;
  clone->_llhd_val = this->_llhd_val;
  clone->_acc_prob = this->_acc_prob;
  clone->_avg_acc_prob = 0.0;

  return clone;
}

double uqInfiniteDimensionalMCMCSampler::acc_prob()
{
  return this->_acc_prob;
}

double uqInfiniteDimensionalMCMCSampler::avg_acc_prob()
{
  return this->_avg_acc_prob;
}

double uqInfiniteDimensionalMCMCSampler::llhd_val() const
{
  return this->_llhd_val;
}

unsigned int uqInfiniteDimensionalMCMCSampler::iteration() const
{
  return this->_iteration;
}
