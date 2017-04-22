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
#include <queso/InterpolationSurrogateIOASCII.h>

// QUESO
#include <queso/MpiComm.h>
#include <queso/StreamUtilities.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

// C++
#include <fstream>

namespace QUESO
{

  template<class V, class M>
  InterpolationSurrogateIOASCII<V,M>::InterpolationSurrogateIOASCII()
    : InterpolationSurrogateIOBase<V,M>()
  {}

  template<class V, class M>
  void InterpolationSurrogateIOASCII<V,M>::read( const std::string& filename,
                                                 const FullEnvironment& env,
                                                 const std::string& vector_space_prefix,
                                                 int reading_rank )
  {
    // Root processor
    int root = reading_rank;

    MpiComm full_comm = env.fullComm();

    std::ifstream input;

    unsigned int dim;

    // Only processor 0 does the reading.
    // We'll broadcast the data as needed
    if( env.fullRank() == root )
      {
        input.open( filename.c_str() );

        // skip the header
        StreamUtilities::skip_comment_lines(input, '#');

        // Read in dimension
        input >> dim;
      }

    // Broadcast the parsed dimension
    full_comm.Bcast( &dim, 1, RawValue_MPI_UNSIGNED, root,
                     "InterpolationSurrogateIOASCII::read()",
                     "MpiComm::Bcast() failed!" );

    // Construct vector space
    this->m_vector_space.reset( new VectorSpace<V,M>(env,
                                                     vector_space_prefix.c_str(),
                                                     dim,
                                                     NULL) );

    // Read in n_points in each dimension
    this->m_n_points.resize(dim);

    if( env.fullRank() == root )
      {
        for( unsigned int d = 0; d < dim; d++ )
            {
              // Skip any comments
              StreamUtilities::skip_comment_lines(input, '#');

              // Make sure file is still good
              if( !input.good() )
                queso_error_msg("ERROR: Found unexpected end-of-file");

              input >> this->m_n_points[d];
            }
      }

    // Broadcast m_n_points
    full_comm.Bcast( &this->m_n_points[0], dim, RawValue_MPI_UNSIGNED, root,
                     "InterpolationSurrogateIOASCII::read()",
                     "MpiComm::Bcast() failed!" );

    // Read parameter bounds
    std::vector<double> param_mins(dim);
    std::vector<double> param_maxs(dim);

    if( env.fullRank() == root )
      {
        for( unsigned int d = 0; d < dim; d++ )
          {
            // Skip any comments
            StreamUtilities::skip_comment_lines(input, '#');

            // Make sure file is still good
            if( !input.good() )
              queso_error_msg("ERROR: Found unexpected end-of-file");

            input >> param_mins[d] >> param_maxs[d];
          }
      }

    // Broadcast the bounds
    full_comm.Bcast( &param_mins[0], dim, RawValue_MPI_DOUBLE, root,
                     "InterpolationSurrogateIOASCII::read()",
                     "MpiComm::Bcast() failed!" );

    full_comm.Bcast( &param_maxs[0], dim, RawValue_MPI_DOUBLE, root,
                     "InterpolationSurrogateIOASCII::read()",
                     "MpiComm::Bcast() failed!" );

    // Construct parameter domain
    /* BoxSubset copies the incoming paramMins/paramMaxs so we don't
       need to cache these copies, they can die. */
    QUESO::GslVector paramMins(this->m_vector_space->zeroVector());
    QUESO::GslVector paramMaxs(this->m_vector_space->zeroVector());

    for( unsigned int d = 0; d < dim; d++ )
      {
        paramMins[d] = param_mins[d];
        paramMaxs[d] = param_maxs[d];
      }

    this->m_domain.reset( new BoxSubset<V,M>(vector_space_prefix.c_str(),
                                             *(this->m_vector_space.get()),
                                             paramMins,
                                             paramMaxs) );

    // Construct data object
    this->m_data.reset( new InterpolationSurrogateData<V,M>(*(this->m_domain.get()),
                                                            this->m_n_points) );

    // Now read in the values
    if( env.fullRank() == root )
      {
        for( unsigned int n = 0; n < this->m_data->n_values(); n++ )
          {
            // Skip any comments
            StreamUtilities::skip_comment_lines(input, '#');

            // Make sure file is still good
            if( !input.good() )
              queso_error_msg("ERROR: Found unexpected end-of-file");

            double value;
            input >> value;

            this->m_data->set_value( n, value );
          }

        // We are done parsing now, so close the file
        input.close();
      }

    // Broadcast the values
    this->m_data->sync_values(root);

    // Fin
  }

  template<class V, class M>
  void InterpolationSurrogateIOASCII<V,M>::write( const std::string& filename,
                                                  const InterpolationSurrogateData<V,M>& data,
                                                  int writing_rank ) const
  {
    // Make sure there are values in the data. If not the user didn't populate the data
    if( !(data.n_values() > 0) )
      {
        std::string error = "ERROR: No values found in InterpolationSurrogateData.\n";
        error += "Cannot write data without values.\n";
        error += "Use InterpolationSurrogateBuilder or the read method to populate\n";
        error += "data values.\n";

        queso_error_msg(error);
      }

    std::ofstream output;

    // Only processor 0 does the writing
    if( data.get_paramDomain().env().fullRank() == writing_rank )
      {
        output.open( filename.c_str() );

        // Write simpler header comments
        std::string header = "# Data for interpolation surrogate\n";
        header += "# Format is as follows:\n";
        header += "# dimension (unsigned int)\n";
        header += "# n_points in each dimension\n";
        header += "# x_min, x_max pairs for each dimension\n";
        header += "# values for each point in parameter space\n";
        header += "# values must be ordered in structured format.\n";
        output << header;

        // Write dimension
        unsigned int dim = data.get_paramDomain().vectorSpace().dimGlobal();
        output << dim << std::endl;

        // Write n_points
        output << "# n_points" << std::endl;
        for( unsigned int d = 0; d < dim; d++ )
          {
            output << data.get_n_points()[d] << std::endl;
          }

        // Set precision for double output
        output << std::scientific << std::setprecision(16);

        // Write domain bounds
        output << "# domain bounds" << std::endl;
        for( unsigned int d = 0; d < dim; d++ )
          {
            output << data.x_min(d) << " "
                   << data.x_max(d) << std::endl;
          }

        // Write values
        output << "# values" << std::endl;
        for( unsigned int n = 0; n < data.n_values(); n++ )
          {
            output << data.get_value(n) << std::endl;
          }

        // All done
        output.close();

      } // data.get_paramDomain().env().fullRank() == writing_rank
  }

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateIOASCII<QUESO::GslVector,QUESO::GslMatrix>;
