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
#include <queso/InterpolationSurrogateBase.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

// C++
#include <sstream>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateBase<V,M>::InterpolationSurrogateBase(const BoxSubset<V,M> & domain,
                                                              const std::vector<unsigned int>& n_points,
                                                              const std::vector<double>& values)
    : SurrogateBase<V>(),
    m_domain(domain),
    m_n_points(n_points),
    m_values(values)
    {
      // This checks that the dimension of n_points and the domain are consistent
      this->check_dim_consistency();
    }

  template<class V, class M>
  unsigned int InterpolationSurrogateBase<V,M>::coordToGlobal( const std::vector<unsigned int>& coord_indices,
                                                               const std::vector<unsigned int>& n_points ) const
  {
    queso_assert_equal_to( coord_indices.size(), this->dim() );
    queso_assert_equal_to( n_points.size(), this->dim() );

    /* Mapping is: i + j*n_i + k*n_i*n_j + l*n_i*n_j*n_k + ...
       Initialize global_index to "i".
       Then loop and build up each term.*/
    queso_assert_less( coord_indices[0], n_points[0] );

    unsigned int global_index = coord_indices[0];

    for( unsigned int d = 1; d < this->dim(); d++ )
      {
        queso_assert_less( coord_indices[d], n_points[d] );

        // Accumulate the current term
        unsigned int idx = coord_indices[d];

        for( int local_d = d-1; local_d >=0; local_d -= 1)
          {
            idx *= n_points[local_d];
          }

        global_index += idx;
      }

    return global_index;
  }

  template<class V, class M>
  void InterpolationSurrogateBase<V,M>::check_dim_consistency() const
  {
    if( this->dim() != this->m_n_points.size() )
      {
        std::stringstream vspace_dim;
        vspace_dim << this->m_domain.vectorSpace().dimGlobal();

        std::stringstream n_points_dim;
        n_points_dim << this->m_n_points.size();

        std::string error = "ERROR: Mismatch between dimension of parameter space and number of points\n.";
        error += "        domain dimension = " + vspace_dim.str() + "\n";
        error += "        points dimension = " + n_points_dim.str() + "\n";

        queso_error_msg(error);
      }
  }

  template<class V, class M>
  double InterpolationSurrogateBase<V,M>::spacing( unsigned int dim ) const
  {
    unsigned int n_intervals = this->m_n_points[dim]-1;
    double x_min = this->x_min(dim);
    double x_max = this->x_max(dim);

    return (x_max-x_min)/n_intervals;
  }

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateBase<QUESO::GslVector,QUESO::GslMatrix>;
