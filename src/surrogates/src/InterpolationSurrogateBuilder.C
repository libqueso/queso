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

// This class
#include <queso/InterpolationSurrogateBuilder.h>

// QUESO
#include <queso/MpiComm.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/MultiDimensionalIndexing.h>
#include <queso/StreamUtilities.h>

// C++
#include <numeric>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateBuilder<V,M>::InterpolationSurrogateBuilder( InterpolationSurrogateDataSet<V,M>& data )
    : SurrogateBuilderBase<V>(),
    m_data(data),
    m_njobs(this->get_default_data().get_paramDomain().env().numSubEnvironments(), 0)
  {
    this->partition_work();
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::partition_work()
  {
    // Convenience
    unsigned int n_values = this->get_default_data().n_values();
    unsigned int n_workers = this->get_default_data().get_paramDomain().env().numSubEnvironments();

    unsigned int n_jobs = n_values/n_workers;
    unsigned int n_leftover = this->get_default_data().n_values() % n_workers;

    /* If the number of values is evenly divisible over all workers,
       then everyone gets the same amount work */
    if( n_leftover  == 0 )
      {
        for(unsigned int n = 0; n < n_workers; n++)
          this->m_njobs[n] = n_jobs;
      }
    /* Otherwise, some workers get more work than others*/
    else
      {
        for(unsigned int n = 0; n < n_workers; n++)
          {
            if( n < n_leftover )
              this->m_njobs[n] = n_jobs+1;
            else
              this->m_njobs[n] = n_jobs;
          }
      }

    // Sanity check
    queso_assert_equal_to( (int)n_values, std::accumulate( m_njobs.begin(), m_njobs.end(), 0 ) );
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::build_values()
  {
    unsigned int n_begin, n_end;
    this->set_work_bounds( n_begin, n_end );

    // Cache each processors work, then we only need to do 1 Allgather
    std::vector<unsigned int> local_n(n_end-n_begin);

    // We need to cache (n_end-n_begin) values for each dataset,
    std::vector<std::vector<double> > local_values(this->m_data.size());
    for( std::vector<std::vector<double> >::iterator it = local_values.begin();
         it != local_values.end(); ++it )
      it->resize(n_end-n_begin);

    unsigned int count = 0;

    // vector to store current domain value
    V domain_vector(this->get_default_data().get_paramDomain().vectorSpace().zeroVector());

    // vector to store values evaluated at the current domain_vector
    std::vector<double> values(this->m_data.size());

    for( unsigned int n = n_begin; n < n_end; n++ )
      {
        this->set_domain_vector( n, domain_vector );

        this->evaluate_model( domain_vector, values );

        local_n[count] = n;

        for( unsigned int s = 0; s < this->m_data.size(); s++ )
          local_values[s][count] = values[s];

        count += 1;
      }

    /* Sync all the locally computed values between the subenvironments
       so all processes have all the computed values. We need to sync
       values for every data set. */
    for( unsigned int s = 0; s < this->m_data.size(); s++ )
      this->sync_data( local_n, local_values[s], this->m_data.get_dataset(s) );
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::set_work_bounds( unsigned int& n_begin, unsigned int& n_end ) const
  {
    unsigned int my_subid = this->get_default_data().get_paramDomain().env().subId();

    /* Starting index will be the sum of the all the previous num jobs */
    n_begin = 0;
    for( unsigned int n = 0; n < my_subid; n++ )
      n_begin += m_njobs[n];

    n_end = n_begin + m_njobs[my_subid];
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::sync_data( std::vector<unsigned int>& local_n,
                                                      std::vector<double>& local_values,
                                                      InterpolationSurrogateData<V,M>& data )
  {
    // Only members of the inter0comm will do the communication of the local values
    unsigned int my_subrank = data.get_paramDomain().env().subRank();

    if( my_subrank == 0 )
      {
        std::vector<double> all_values(data.n_values());

        std::vector<unsigned int> all_indices(data.n_values());

        std::vector<int> strides;
        this->compute_strides( strides );

        const MpiComm& inter0comm = data.get_paramDomain().env().inter0Comm();

        /*! \todo Would be more efficient to pack local_n and local_values
            togethers and do Gatherv only once. */
        inter0comm.template Gatherv<unsigned int>(&local_n[0], local_n.size(),
            &all_indices[0], &m_njobs[0], &strides[0],
            0 /*root*/, "InterpolationSurrogateBuilder::sync_data()",
            "MpiComm::gatherv() failed!");

        inter0comm.template Gatherv<double>(&local_values[0],
            local_values.size(), &all_values[0], &m_njobs[0], &strides[0],
            0 /*root*/, "InterpolationSurrogateBuilder::sync_data()",
            "MpiComm::gatherv() failed!");

        // Now set the values.
        /* PB: Although we are guaranteed per-rank ordering of the data we gathered,
           I'm not sure we can assume the same continuity of the inter0 ranks, i.e.
           I'm not sure how QUESO ordered the inter0 ranks. So, we go ahead and
           manually set the values. */
        if( data.get_paramDomain().env().subRank() == 0 )
          {
            for( unsigned int n = 0; n < data.n_values(); n++ )
              data.set_value( all_indices[n], all_values[n] );
          }
      }

    // Now broadcast the values data to all other processes
    data.sync_values( 0 /*root*/);
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::set_domain_vector( unsigned int n, V& domain_vector ) const
  {
    // Convert global index n to local coordinates in each dimension
    std::vector<unsigned int> indices(this->get_default_data().dim());
    MultiDimensionalIndexing::globalToCoord( n, this->get_default_data().get_n_points(), indices );

    // Use indices to get x coordinates and populate domain_vector
    for( unsigned int d = 0; d < this->get_default_data().dim(); d++ )
      {
        domain_vector[d] = this->get_default_data().get_x( d, indices[d] );
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::compute_strides( std::vector<int>& strides ) const
  {
    unsigned int n_subenvs = this->get_default_data().get_paramDomain().env().numSubEnvironments();

    strides.resize(n_subenvs);

    // Don't stride the first entry
    strides[0] = 0;

    int stride = 0;
    for( unsigned int n = 1; n < n_subenvs; n++ )
      {
        // The stride is measured agaisnt the beginning of the buffer
        // We want things packed tightly together so just stride
        // by the number of entries from the previous group.
        stride += this->m_njobs[n-1];
        strides[n] = stride;
      }
  }

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateBuilder<QUESO::GslVector,QUESO::GslMatrix>;
