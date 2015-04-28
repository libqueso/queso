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

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateBuilder<V,M>::InterpolationSurrogateBuilder( InterpolationSurrogateData<V,M>& data )
    : SurrogateBuilderBase<V>(),
    m_data(data)
  {}

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
