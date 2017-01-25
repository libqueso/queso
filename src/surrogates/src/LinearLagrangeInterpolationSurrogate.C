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
#include <queso/LinearLagrangeInterpolationSurrogate.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/MultiDimensionalIndexing.h>
#include <queso/InterpolationSurrogateData.h>

namespace QUESO
{
  template<class V, class M>
  LinearLagrangeInterpolationSurrogate<V,M>::LinearLagrangeInterpolationSurrogate(const InterpolationSurrogateData<V,M>& data)
    : InterpolationSurrogateBase<V,M>(data)
  {}

  template<class V, class M>
  double LinearLagrangeInterpolationSurrogate<V,M>::evaluate(const V & domainVector) const
  {
    /* Populate indices. These are the lower bound global indices for the
       "element" containing the domainVector */
    std::vector<unsigned int> indices(this->m_data.dim());
    this->compute_interval_indices(domainVector, indices);

    // lower bound coordinates
    std::vector<double> x_min(this->m_data.dim());

    // upper bound coordinates
    std::vector<double> x_max(this->m_data.dim());

    // Function value at each of the "nodes" for the "element"
    std::vector<double> values(this->n_coeffs());

    // Use the indices to populate the x_min, etc. structures
    this->compute_interval_values( indices, x_min, x_max, values );

    // Now evaluate the interpolant
    return this->eval_interpolant( x_min, x_max, values, domainVector );
  }

  template<class V, class M>
  void LinearLagrangeInterpolationSurrogate<V,M>::compute_interval_indices(const V & domainVector,
                                                                           std::vector<unsigned int>& indices) const
  {
    queso_assert_equal_to( domainVector.sizeGlobal(), this->m_data.dim() );
    queso_assert_equal_to( indices.size(), this->m_data.dim() );

    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
      {
        double spacing = this->m_data.spacing(d);
        indices[d] = std::floor( (domainVector[d] - this->m_data.x_min(d))/spacing );

        // Index should be less than the number of point along this dimension
        queso_assert_less( indices[d], this->m_data.get_n_points()[d] );
      }
  }

  template<class V, class M>
  void LinearLagrangeInterpolationSurrogate<V,M>::compute_interval_values( const std::vector<unsigned int>& indices,
                                                                           std::vector<double>& x_min,
                                                                           std::vector<double>& x_max,
                                                                           std::vector<double>& values ) const
  {
    queso_assert_equal_to( x_min.size(), this->m_data.dim() );
    queso_assert_equal_to( x_max.size(), this->m_data.dim() );
    queso_assert_equal_to( values.size(), this->n_coeffs() );
    queso_assert_equal_to( indices.size(), this->m_data.dim() );

    // First, use the lower bound (global) indices to populate x_min, x_max
    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
      {
        x_min[d] = this->m_data.get_x( d, indices[d] );
        x_max[d] = x_min[d] + this->m_data.spacing( d );
      }

    // Now populate the values.
    std::vector<unsigned int> local_indices(this->m_data.dim());
    std::vector<unsigned int> global_indices(this->m_data.dim());
    for( unsigned int n = 0 ; n < this->n_coeffs(); n++ )
      {
        // Figure out local indices (each coordinate is 0 or 1)
        this->singleToCoords(n,local_indices);

        /* For each dimension, if the local index is 0, use the lower
           bound of the global index. If the local index is 1, use the
           "upper" global index.

           In 2-D:
           The lower left corner of the element will have local_indices = [0,0];
           Then, global_indices = [ indices[0], indices[1] ];
           The upper right corner will have local_indices = [1,1];
           Then, global_indices = [ indices[0]+1, indices[1]+1 ]; */
        for( unsigned int d = 0; d < this->m_data.dim(); d++ )
          {
            if( local_indices[d] == 0 )
              global_indices[d] = indices[d];

            else if( local_indices[d] == 1 )
              global_indices[d] = indices[d]+1;

            // This shouldn't happen
            else
              queso_error();
          }

        /* Now that we have the global indices for each coordinate,
           we get the "global" index. This is the index into the global
           values array */
        unsigned int global = MultiDimensionalIndexing::coordToGlobal( global_indices, this->m_data.get_n_points() );
        values[n] = this->m_data.get_value(global);
      }
  }

  template<class V, class M>
  double LinearLagrangeInterpolationSurrogate<V,M>::eval_interpolant( const std::vector<double>& x_min,
                                                                      const std::vector<double>& x_max,
                                                                      const std::vector<double>& values,
                                                                      const V & domainVector ) const
  {
    queso_assert_equal_to( x_min.size(), this->m_data.dim() );
    queso_assert_equal_to( x_max.size(), this->m_data.dim() );
    queso_assert_equal_to( values.size(), this->n_coeffs() );
    queso_assert_equal_to( domainVector.sizeGlobal(), this->m_data.dim() );

    double interp_value = 0.0;

    std::vector<unsigned int> indices( this->m_data.dim() );

    for( unsigned int n = 0; n < this->n_coeffs(); n++ )
      {
        this->singleToCoords( n, indices );

        double shape_fn = this->tensor_product_lagrange( x_min, x_max, indices, domainVector );

        interp_value += values[n]*shape_fn;
      }

    return interp_value;
  }

  template<class V, class M>
  unsigned int LinearLagrangeInterpolationSurrogate<V,M>::coordsToSingle( const std::vector<unsigned int>& indices ) const
  {
    unsigned int local_dim = 2;
    std::vector<unsigned int> n_points( this->m_data.dim(), local_dim );

    /* We're abusing this function as it does what we need to with local_dim = 2
       for all entries of n_points. */
    return MultiDimensionalIndexing::coordToGlobal( indices, n_points );
  }

  template<class V, class M>
  void LinearLagrangeInterpolationSurrogate<V,M>::singleToCoords( unsigned int global, std::vector<unsigned int>& indices ) const
  {
    unsigned int local_dim = 2;
    std::vector<unsigned int> n_points( this->m_data.dim(), local_dim );

    /* We're abusing this function as it does what we need to with local_dim = 2
       for all entries of n_points. */
    return MultiDimensionalIndexing::globalToCoord( global, n_points, indices );
  }

  template<class V, class M>
  double LinearLagrangeInterpolationSurrogate<V,M>::tensor_product_lagrange( const std::vector<double>& x_min,
                                                                             const std::vector<double>& x_max,
                                                                             const std::vector<unsigned int>& indices,
                                                                             const V & domainVector ) const
  {
    queso_assert_equal_to( x_min.size(), this->m_data.dim() );
    queso_assert_equal_to( x_max.size(), this->m_data.dim() );
    queso_assert_equal_to( indices.size(), this->m_data.dim() );
    queso_assert_equal_to( domainVector.sizeGlobal(), this->m_data.dim() );

    double value = 1.0;

    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
      {
        value *= this->lagrange_poly( x_min[d], x_max[d], domainVector[d], indices[d] );
      }

    return value;
  }

  template<class V, class M>
  double LinearLagrangeInterpolationSurrogate<V,M>::lagrange_poly( double x0, double x1, double x, unsigned int index ) const
  {
    // Make sure we're trying to interpolate between the given points
    queso_assert_less_equal( x, x1 );
    queso_assert_greater_equal( x, x0 );

    // Make sure index is either 0 or 1
    queso_assert( (index == 0) || (index == 1) );

    double value = 0.0;

    if( index == 0 )
      value = (x-x1)/(x0-x1);
    else
      value = (x-x0)/(x1-x0);

    return value;
  }

} // end namespace QUESO

// Instantiate
template class QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>;
