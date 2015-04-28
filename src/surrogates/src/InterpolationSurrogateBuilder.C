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
    this->partition_work();
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
    // vector to store current domain value
    V* domain_vector = this->m_data.get_paramDomain().vectorSpace().newVector();

    unsigned int n_begin, n_end;
    this->set_work_bounds( n_begin, n_end );

    for( unsigned int n = n_begin; n < n_end; n++ )
      {
        this->set_domain_vector( n, *domain_vector );

        double value = this->evaluate_model( *domain_vector );

        this->set_value( n, value );
      }
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::set_work_bounds( unsigned int& n_begin, unsigned int& n_end ) const
  {
    n_begin = 0;
    n_end = this->m_data.n_values();
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::set_value( unsigned int n, double value )
  {
    this->m_data.set_value(n, value);
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

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateBuilder<QUESO::GslVector,QUESO::GslMatrix>;
