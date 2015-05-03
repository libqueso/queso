//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/InterpolationSurrogateHelper.h>
#include <queso/StreamUtilities.h>

// C++
#include <numeric>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateBuilder<V,M>::InterpolationSurrogateBuilder( InterpolationSurrogateData<V,M>& data )
    : SurrogateBuilderBase<V>(),
    m_data(data),
    m_njobs(data.get_paramDomain().env().numSubEnvironments(), 0)
  {
    this->check_process_config();

    this->partition_work();
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::check_process_config()
  {
    /* If fullComm() > 1 and n_subenvironments == 1 this is strange and means
       the user is doing redundant work. So, let's make that an error.
       This could happen if, for example the user forgets to change the number
       of subenvironments to > 1 and runs with multiple MPI processes. */
    unsigned int full_comm_size = this->m_data.get_paramDomain().env().fullComm().NumProc();
    unsigned int n_subenvs = this->m_data.get_paramDomain().env().numSubEnvironments();

    if( (full_comm_size > 1) &&  (n_subenvs == 1) )
      {
        std::string error = "ERROR: fullComm() size is greater than 1 and the number\n";
        error += "       of subenvrionments = 1. InterpolationSurrogateBuilder\n";
        error += "       is not compatible with this configuration.\n";
        error += "       Did you forget to change the number of subenvironments?\n";

        queso_error_msg(error);
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::partition_work()
  {
    // Convenience
    unsigned int n_values = this->m_data.n_values();
    unsigned int n_workers = this->m_data.get_paramDomain().env().numSubEnvironments();

    unsigned int n_jobs = n_values/n_workers;
    unsigned int n_leftover = this->m_data.n_values() % n_workers;

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
    queso_assert_equal_to( n_values, std::accumulate( m_njobs.begin(), m_njobs.end(), 0 ) );
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::build_values()
  {
    unsigned int n_begin, n_end;
    this->set_work_bounds( n_begin, n_end );

    // Cache each processors work, then we only need to do 1 Allgather
    std::vector<unsigned int> local_n(n_end-n_begin);
    std::vector<double> local_values(n_end-n_begin);
    unsigned int count = 0;

    // vector to store current domain value
    V domain_vector(this->m_data.get_paramDomain().vectorSpace().zeroVector());

    for( unsigned int n = n_begin; n < n_end; n++ )
      {
        this->set_domain_vector( n, domain_vector );

        double value = this->evaluate_model( domain_vector );

        local_n[count] = n;
        local_values[count] = value;
        count += 1;
      }

    /* Sync all the locally computed values between the subenvironments
       so all processes have all the computed values. */
    this->sync_data( local_n, local_values );
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::set_work_bounds( unsigned int& n_begin, unsigned int& n_end ) const
  {
    unsigned int my_subid = this->m_data.get_paramDomain().env().subId();

    /* Starting index will be the sum of the all the previous num jobs */
    n_begin = 0;
    for( unsigned int n = 0; n < my_subid; n++ )
      n_begin += m_njobs[n];

    n_end = n_begin + m_njobs[my_subid];
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::sync_data( std::vector<unsigned int>& local_n,
                                                      std::vector<double>& local_values  )
  {
    // Only members of the inter0comm will do the communication of the local values
    unsigned int my_subrank = this->m_data.get_paramDomain().env().subRank();

    if( my_subrank == 0 )
      {
        std::vector<double> all_values(this->m_data.n_values());

        std::vector<unsigned int> all_indices(this->m_data.n_values());

        std::vector<int> strides;
        this->compute_strides( strides );

        const MpiComm& inter0comm = this->m_data.get_paramDomain().env().inter0Comm();

        /*! \todo Would be more efficient to pack local_n and local_values
            togethers and go Gatherv only once. */
        inter0comm.Gatherv( &local_n[0], local_n.size(), MPI_UNSIGNED,
                            &all_indices[0], &m_njobs[0], &strides[0], MPI_UNSIGNED,
                            0 /*root*/,
                            "InterpolationSurrogateBuilder::sync_data()",
                            "MpiComm::gatherv() failed!" );

        inter0comm.Gatherv( &local_values[0], local_values.size(), MPI_DOUBLE,
                            &all_values[0], &m_njobs[0], &strides[0], MPI_DOUBLE,
                            0 /*root*/,
                            "InterpolationSurrogateBuilder::sync_data()",
                            "MpiComm::gatherv() failed!" );

        // Now set the values.
        /* PB: Although we are guaranteed per-rank ordering of the data we gathered,
           I'm not sure we can assume the same continuity of the inter0 ranks, i.e.
           I'm not sure how QUESO ordered the inter0 ranks. So, we go ahead and
           manually set the values. */
        if( this->m_data.get_paramDomain().env().subRank() == 0 )
          {
            for( unsigned int n = 0; n < this->m_data.n_values(); n++ )
              this->m_data.set_value( all_indices[n], all_values[n] );
          }
      }

    // Now broadcast the values data to all other processes
    const MpiComm& fullcomm = this->m_data.get_paramDomain().env().fullComm();

    std::vector<double>& values = this->m_data.get_values();
    fullcomm.Bcast( &values[0], values.size(), MPI_DOUBLE,
                    0 /*root*/,
                    "InterpolationSurrogateBuilder::sync_data()",
                    "MpiComm::Bcast() failed!" );
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::set_domain_vector( unsigned int n, V& domain_vector ) const
  {
    // Convert global index n to local coordinates in each dimension
    std::vector<unsigned int> indices(this->m_data.dim());
    InterpolationSurrogateHelper::globalToCoord( n, this->m_data.get_n_points(), indices );

    // Use indices to get x coordinates and populate domain_vector
    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
      {
        domain_vector[d] = this->m_data.get_x( d, indices[d] );
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::compute_strides( std::vector<int>& strides ) const
  {
    unsigned int n_subenvs = this->m_data.get_paramDomain().env().numSubEnvironments();

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

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::read( std::istream& input )
  {
    // skip the header
    StreamUtilities::skip_comment_lines(input, '#');

    // Read in dimension and confirm it matches with m_data
    unsigned int dim;
    input >> dim;

    this->check_parsed_dim( dim, this->m_data.get_paramDomain().vectorSpace().dimGlobal() );


    // Read in n_points in each dimension and confirm it matches with m_data
    {
      std::vector<unsigned int> n_points(dim);
      for( unsigned int d = 0; d < dim; d++ )
        {
          // Skip any comments
          StreamUtilities::skip_comment_lines(input, '#');

          // Make sure file is still good
          if( !input.good() )
            queso_error_msg("ERROR: Found unexpected end-of-file");

          input >> n_points[d];
        }

      this->check_parsed_points(n_points, this->m_data.get_n_points() );
    }

    // Read in bounds for each dimension and confirm it matches with m_data
    {
      // Cache
      std::vector<double> x_min(dim);
      std::vector<double> x_max(dim);
      std::vector<double> param_xmin(dim);
      std::vector<double> param_xmax(dim);

      for( unsigned int d = 0; d < dim; d++ )
        {
          // Instead of adding new API for the vectors, we just cache them here
          param_xmin[d] = this->m_data.x_min(d);
          param_xmax[d] = this->m_data.x_max(d);

          // Skip any comments
          StreamUtilities::skip_comment_lines(input, '#');

          // Make sure file is still good
          if( !input.good() )
            queso_error_msg("ERROR: Found unexpected end-of-file");

          input >> x_min[d] >> x_max[d];
        }

      this->check_parsed_bounds( x_min, x_max, param_xmin, param_xmax );
    }

    // Now read in values and set them in the data
    for( unsigned int n = 0; n < this->m_data.n_values(); n++ )
      {
        // Skip any comments
          StreamUtilities::skip_comment_lines(input, '#');

          // Make sure file is still good
          if( !input.good() )
            queso_error_msg("ERROR: Found unexpected end-of-file");

          double value;
          input >> value;

          this->m_data.set_value( n, value );
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::check_parsed_dim( unsigned int parsed_dim,
                                                             unsigned int param_dim )
  {
    if( parsed_dim != param_dim )
      {
        std::stringstream dim_str;
        dim_str << parsed_dim;

        std::stringstream param_dim_str;
        param_dim_str << param_dim;

        std::string error = "ERROR: Mismatch in parameter dimension and parsed dimension.\n";
        error += " parsed dimension   = "+dim_str.str()+"\n";
        error += " paramter dimension = "+param_dim_str.str()+"\n";
        queso_error_msg(error);
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::check_parsed_points( const std::vector<unsigned int>& parsed_points,
                                                                const std::vector<unsigned int>& param_points )
  {
    queso_assert_equal_to( parsed_points.size(), param_points.size() );

    bool mismatch = false;

    for( unsigned int d = 0; d < parsed_points.size(); d++ )
      {
        if( parsed_points[d] != param_points[d] )
          mismatch = true;
      }

    if( mismatch )
      {
        std::string error = "ERROR: Mismatch in parameter n_points and parsed n_points.\n";
        error += "parsed n_points    param n_points\n";
        for( unsigned int d = 0; d < parsed_points.size(); d++ )
          {
            std::stringstream parsed_str;
            parsed_str << parsed_points[d];

            std::stringstream param_str;
            param_str << param_points[d];

            error += parsed_str.str()+"    "+param_str.str()+"\n";
          }

        queso_error_msg(error);
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::check_parsed_bounds( const std::vector<double>& parsed_xmin,
                                                                const std::vector<double>& parsed_xmax,
                                                                const std::vector<double>& param_xmin,
                                                                const std::vector<double>& param_xmax )
  {
    queso_assert_equal_to( parsed_xmin.size(), parsed_xmax.size() );
    queso_assert_equal_to( parsed_xmin.size(), param_xmin.size() );
    queso_assert_equal_to( parsed_xmax.size(), param_xmax.size() );

    bool mismatch = false;

    for( unsigned int d = 0; d < parsed_xmin.size(); d++ )
      {
        double tol = std::numeric_limits<double>::epsilon();

        if( ( std::abs(parsed_xmin[d]-param_xmin[d]) > tol ) ||
            ( std::abs(parsed_xmin[d]-param_xmin[d]) > tol ) )
          mismatch = true;
      }

    if( mismatch )
      {
        std::string error = "ERROR: Mismatch in parameter bounds and parsed bounds.\n";
        error += "params bounds\n";
        for( unsigned int d = 0; d < parsed_xmin.size(); d++ )
          {
            std::stringstream xmin_str;
            xmin_str << param_xmin[d];

            std::stringstream xmax_str;
            xmax_str << param_xmax[d];

            error += xmin_str.str()+" "+xmax_str.str()+"\n";
          }
        error += "\n";

        error += "parsed bounds\n";
        for( unsigned int d = 0; d < parsed_xmin.size(); d++ )
          {
            std::stringstream xmin_str;
            xmin_str << parsed_xmin[d];

            std::stringstream xmax_str;
            xmax_str << parsed_xmax[d];

            error += xmin_str.str()+" "+xmax_str.str()+"\n";
          }

        queso_error_msg(error);
      }

  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::write( std::ostream& output ) const
  {
    // Write simpler header comments
    std::string header = "# Data for interpolation surrogate\n";
    header += "# Format is as follows:\n";
    header += "# dimension (unsigned int)\n";
    header += "# n_points in each dimension\n";
    header += "# x_min, x_max pairs for each dimension\n";
    header += "# values for each point in parameter space\n";
    header += "# values musted ordered in structured format.\n";
    output << header;

    // Write dimension
    unsigned int dim = this->m_data.get_paramDomain().vectorSpace().dimGlobal();
    output << dim << std::endl;

    // Write n_points
    output << "# n_points" << std::endl;
    for( unsigned int d = 0; d < dim; d++ )
      {
        output << this->m_data.get_n_points()[d] << std::endl;
      }

    // Set precision for double output
    output << std::scientific << std::setprecision(16);

    // Write domain bounds
    output << "# domain bounds" << std::endl;
    for( unsigned int d = 0; d < dim; d++ )
      {
        output << this->m_data.x_min(d) << " "
               << this->m_data.x_max(d) << std::endl;
      }

    // Write values
    output << "# values" << std::endl;
    for( unsigned int n = 0; n < this->m_data.n_values(); n++ )
      {
        output << this->m_data.get_value(n) << std::endl;
      }
  }

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateBuilder<QUESO::GslVector,QUESO::GslMatrix>;
