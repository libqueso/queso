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

#include <queso/Defines.h>

#ifdef QUESO_HAVE_HDF5

#include <cmath>
#include <algorithm>
#include <sstream>

// QUESO includes
#include <queso/Miscellaneous.h>
#include <queso/InfiniteDimensionalMeasureBase.h>

#include <queso/InfiniteDimensionalLikelihoodBase.h>
#include <queso/InfiniteDimensionalMCMCSampler.h>
#include <queso/InfiniteDimensionalMCMCSamplerOptions.h>

namespace QUESO {

InfiniteDimensionalMCMCSampler::InfiniteDimensionalMCMCSampler(
    const BaseEnvironment& env,
    InfiniteDimensionalMeasureBase & prior,
    InfiniteDimensionalLikelihoodBase & llhd,
    InfiniteDimensionalMCMCSamplerOptions * ov)
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
  std::cout << "Entering InfiniteDimensionalMCMCSampler class" << std::endl;
#endif
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering InfiniteDimensionalMCMCSampler::constructor()"
                            << ": prefix = " << this->m_ov->m_prefix
                            << ", m_env.optionsInputFileName() = " << this->m_env.optionsInputFileName()
                            << std::endl;
  }

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In InfiniteDimensionalMCMCSamplerOptions,"
            << " finished scanning options" << std::endl;
#endif

  // Verify parent directory exists (for cases when a user specifies a
  // relative path for the desired output file).
  // Only subprocess of subrank 0 creates the output file
  if ((this->m_env).subRank() == 0) {
    int irtrn = CheckFilePath((this->m_ov->m_dataOutputDirName +
                               (this->m_env).subIdString() +
                               "/test.txt").c_str());
    queso_require_greater_equal_msg(irtrn, 0, "unable to verify output path");
  }

  // Only subprocess of subrank 0 creates the output file
  if ((this->m_env).subRank() == 0) {
    this->_outfile = H5Fcreate((this->m_ov->m_dataOutputDirName +
          (this->m_env).subIdString() + "/" +
          this->m_ov->m_dataOutputFileName).c_str(),
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  this->_acc_dset = this->_create_scalar_dataset("acc");
  this->_avg_acc_dset = this->_create_scalar_dataset("avg_acc");
  this->_neg_log_llhd_dset = this->_create_scalar_dataset("neg_log_llhd");
  this->_L2_norm_samples_dset = this->_create_scalar_dataset("L2_norm_samples");
  this->_L2_norm_mean_dset = this->_create_scalar_dataset("L2_norm_mean");
  this->_L2_norm_var_dset = this->_create_scalar_dataset("L2_norm_var");

  // Ensure that we created the path and the output files
  ((this->m_env).fullComm()).Barrier();

  this->_iteration = 0;
  this->_acc_prob = 0.0;
  this->_avg_acc_prob = 0.0;

  // Zero out these guys.  There's probably a better way of doing this than
  // calling prior.draw() at the start, but creation of a Sampler object
  // should only be done a O(1) times anyway.
  this->_delta->zero();
  this->_M2->zero();
  this->current_physical_mean->zero();
  this->current_physical_var->zero();

  // LibMeshFunction & p = libmesh_cast_ref<LibMeshFunction &>(*(this->current_physical_state));

  // Seed the sampler at the truth
  // for (unsigned int ii = 0; ii < 513; ii++) {
  //   p.equation_systems->get_system<ExplicitSystem>("Function").solution->set(ii, std::cos(2.0 * boost::math::constants::pi<double>() * ii / 512.0));
  // }
  // p.equation_systems->get_system<ExplicitSystem>("Function").solution->close();
  // std::string seedname("seedfn");
  // p.save_function(seedname);

  // Initialise cost function
  this->_llhd_val = this->llhd.evaluate(*(this->current_physical_state));

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(15);
}

InfiniteDimensionalMCMCSampler::~InfiniteDimensionalMCMCSampler()
{
  if ((this->m_env).subRank() == 0) {
    H5Dclose(this->_acc_dset);
    H5Dclose(this->_avg_acc_dset);
    H5Dclose(this->_neg_log_llhd_dset);
    H5Dclose(this->_L2_norm_samples_dset);
    H5Dclose(this->_L2_norm_mean_dset);
    H5Dclose(this->_L2_norm_var_dset);
    H5Fclose(this->_outfile);
  }
}

void InfiniteDimensionalMCMCSampler::_propose()
{
  const double rwmh_step_sq = (this->m_ov->m_rwmh_step * this->m_ov->m_rwmh_step);
  const double coeff = std::sqrt(1.0 - rwmh_step_sq);

  SharedPtr<FunctionBase>::Type p(prior.draw());

  this->proposed_physical_state->zero();
  this->proposed_physical_state->add(coeff, *(this->current_physical_state));
  this->proposed_physical_state->add(this->m_ov->m_rwmh_step, *p);
}

void InfiniteDimensionalMCMCSampler::_metropolis_hastings()
{
  // Downcast since we know we're dealing with a libmesh function
  // LibMeshFunction & p = static_cast<LibMeshFunction &>(*(this->proposed_physical_state));
  // LibMeshFunction & q = static_cast<LibMeshFunction &>(*(this->current_physical_state));
  // Evaluate the likelihood at the proposed state
  double proposed_llhd = this->llhd.evaluate(*(this->proposed_physical_state));
  // double current_llhd = this->llhd.evaluate(q);
  // std::cout << "llhd of  current is: " << current_llhd << std::endl;
  // std::cout << "llhd of proposal is: " << proposed_llhd << std::endl;

  double diff = this->_llhd_val - proposed_llhd;
  double alpha = std::min(1.0, std::exp(diff));
  double rand = this->m_env.rngObject()->uniformSample();
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

void InfiniteDimensionalMCMCSampler::_update_moments()
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
  SharedPtr<FunctionBase>::Type temp_ptr(this->_delta->zero_clone());
  // LibMeshFunction & temp = static_cast<LibMeshFunction &>(*temp_ptr);

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

void InfiniteDimensionalMCMCSampler::step()
{
  this->_propose();
  this->_metropolis_hastings();
  this->_update_moments();

  // We never save the 0th iteration
  if (this->_iteration % this->m_ov->m_save_freq == 0) {
    this->_write_state();
  }
}

hid_t InfiniteDimensionalMCMCSampler::_create_scalar_dataset(const std::string & name)
{
  // Only subprocess with rank 0 manipulates the output file
  if ((this->m_env).subRank() == 0) {
    // Create a 1D dataspace.  Unlimited size.  Initially set to 0.  We will
    // extend it later
    const int ndims = 1;
    hsize_t dims[ndims] = {0};  // dataset dimensions at creation
    hsize_t maxdims[ndims] = {H5S_UNLIMITED};

    hid_t file_space = H5Screate_simple(ndims, dims, maxdims);

    // Create dataset creation property list.  Unlimited datasets must be
    // chunked.  Choosing the chunk size is an issue, here we set it to 1
    hsize_t chunk_dims[ndims] = {1};
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_layout(plist, H5D_CHUNKED);
    H5Pset_chunk(plist, ndims, chunk_dims);

    // Create the dataset
    hid_t dset = H5Dcreate(this->_outfile, name.c_str(), H5T_NATIVE_DOUBLE,
        file_space, H5P_DEFAULT, plist, H5P_DEFAULT);

    // We don't need the property list anymore.  We also don't need the file
    // dataspace anymore because we'll extend it later, making this one
    // invalild anyway.
    H5Pclose(plist);
    H5Sclose(file_space);

    return dset;
  }

  hid_t dummy = -1;
  return dummy;
}

void InfiniteDimensionalMCMCSampler::_append_scalar_dataset(hid_t dset, double data)
{
  // Only subprocess with rank 0 manipulates the output file
  if ((this->m_env).subRank() == 0) {
    int err;
    // Create a memory dataspace for data to append
    const int ndims = 1;
    hsize_t dims[ndims] = { 1 };  // Only writing one double
    hid_t mem_space = H5Screate_simple(ndims, dims, NULL);

    // Extend the dataset
    // Set dims to be the *new* dimension of the extended dataset
    dims[0] = this->_iteration / this->m_ov->m_save_freq;
    err = H5Dset_extent(dset, dims);
    queso_require_greater_equal_msg(err, 0, "H5DSet_extent(dset, dims) failed");

    // Select hyperslab on file dataset
    hid_t file_space = H5Dget_space(dset);
    hsize_t start[1] = {(this->_iteration / this->m_ov->m_save_freq) - 1};
    hsize_t count[1] = {1};

    err = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    queso_require_greater_equal_msg(err, 0, "H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL) failed");

    // hsize_t      size[1];
    // size[0]   = this->_iteration / this->m_ov->m_save_freq;

    // Write the data
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &data);

    // Close a bunch of stuff
    H5Sclose(file_space);
    H5Sclose(mem_space);
  }
}

void InfiniteDimensionalMCMCSampler::_write_state()
{
  this->_append_scalar_dataset(this->_acc_dset, this->acc_prob());
  this->_append_scalar_dataset(this->_avg_acc_dset, this->avg_acc_prob());
  this->_append_scalar_dataset(this->_neg_log_llhd_dset, this->_llhd_val);
  this->_append_scalar_dataset(this->_L2_norm_samples_dset, this->current_physical_state->L2_norm());
  this->_append_scalar_dataset(this->_L2_norm_mean_dset, this->current_physical_mean->L2_norm());
  this->_append_scalar_dataset(this->_L2_norm_var_dset, this->current_physical_var->L2_norm());

  // Now to write the fields.  It appears to be a pain in the arse to write a
  // method to spit this out to HDF5 format.  Also, this won't scale to
  // non-uniform finite element meshes.  Therefore, I'm going to spit out an
  // ExodusII file for each sample, average and variance.  Got disk space?
  std::stringstream basename;
  basename << this->m_ov->m_dataOutputDirName;
  basename << (this->m_env).subIdString();
  basename << "/";  // Sigh

  std::ostringstream curr_iter;
  curr_iter << this->iteration();

  std::string sample_name(basename.str() + "sample.e-s.");
  sample_name += curr_iter.str();
  this->current_physical_state->save_function(sample_name, this->iteration());

  std::string mean_name(basename.str() + "mean.e-s.");
  mean_name += curr_iter.str();
  this->current_physical_mean->save_function(mean_name, this->iteration());

  std::string var_name(basename.str() + "var.e-s.");
  var_name += curr_iter.str();
  this->current_physical_var->save_function(var_name, this->iteration());
}

SharedPtr<InfiniteDimensionalMCMCSampler>::Type InfiniteDimensionalMCMCSampler::clone_and_reset() const
{
  // Set up a clone
  SharedPtr<InfiniteDimensionalMCMCSampler>::Type clone(new InfiniteDimensionalMCMCSampler(this->m_env, this->prior, this->llhd, this->m_ov));

  // Copy the state.
  clone->current_physical_state = this->current_physical_state;
  clone->proposed_physical_state = this->proposed_physical_state;
  clone->_llhd_val = this->_llhd_val;
  clone->_acc_prob = this->_acc_prob;
  clone->_avg_acc_prob = 0.0;

  return clone;
}

double InfiniteDimensionalMCMCSampler::acc_prob()
{
  return this->_acc_prob;
}

double InfiniteDimensionalMCMCSampler::avg_acc_prob()
{
  return this->_avg_acc_prob;
}

double InfiniteDimensionalMCMCSampler::llhd_val() const
{
  return this->_llhd_val;
}

unsigned int InfiniteDimensionalMCMCSampler::iteration() const
{
  return this->_iteration;
}

}  // End namespace QUESO

#endif  // QUESO_HAVE_HDF5
