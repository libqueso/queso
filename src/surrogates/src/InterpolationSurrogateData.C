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
#include <queso/InterpolationSurrogateData.h>

// QUESO
#include <queso/MpiComm.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/math_macros.h>

// C++
#include <sstream>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateData<V,M>::InterpolationSurrogateData(const BoxSubset<V,M> & domain,
                                                              const std::vector<unsigned int>& n_points )
    : m_domain(domain),
      m_n_points(n_points)
  {
    // This checks that the dimension of n_points and the domain are consistent
    this->check_dim_consistency();

    // Size m_values
    this->init_values(this->m_n_points);

    // Check the domain is bounded
    this->check_domain_bounds();
  }

  template<class V, class M>
  void InterpolationSurrogateData<V,M>::check_dim_consistency() const
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
  void InterpolationSurrogateData<V,M>::init_values( const std::vector<unsigned int>& n_points )
  {
    unsigned int n_total_points = 1.0;
    for( std::vector<unsigned int>::const_iterator it = n_points.begin();
         it != n_points.end(); ++it )
      {
        n_total_points *= *it;
      }

    this->m_values.resize(n_total_points);
  }

  template<class V, class M>
  void InterpolationSurrogateData<V,M>::check_domain_bounds() const
  {
    for (unsigned int i = 0; i < m_domain.vectorSpace().dimLocal(); i++) {
      queso_require_msg(queso_isfinite(m_domain.minValues()[i]),
          "Interpolation with an unbounded domain is unsupported");
      queso_require_msg(queso_isfinite(m_domain.maxValues()[i]),
          "Interpolation with an unbounded domain is unsupported");
    }
  }

  template<class V, class M>
  void InterpolationSurrogateData<V,M>::set_values( std::vector<double>& values )
  {
    queso_assert_equal_to( values.size(), m_values.size() );

    this->m_values = values;
  }

  template<class V, class M>
  void InterpolationSurrogateData<V,M>::set_value( unsigned int n, double value )
  {
    queso_assert_less( n, m_values.size() );

    this->m_values[n] = value;
  }

  template<class V, class M>
  double InterpolationSurrogateData<V,M>::spacing( unsigned int dim ) const
  {
    unsigned int n_intervals = this->m_n_points[dim]-1;
    double x_min = this->x_min(dim);
    double x_max = this->x_max(dim);

    return (x_max-x_min)/n_intervals;
  }

  template<class V, class M>
  double InterpolationSurrogateData<V,M>::get_x( unsigned int dim, unsigned int index ) const
  {
    double x_min = this->x_min(dim);
    double spacing = this->spacing(dim);

    return x_min + spacing*index;
  }

  template<class V, class M>
  void InterpolationSurrogateData<V,M>::sync_values( unsigned int root )
  {
    MpiComm full_comm = this->m_domain.env().fullComm();

    full_comm.Bcast( &this->m_values[0], this->n_values(),
                     RawValue_MPI_DOUBLE, root,
                     "InterpolationSurrogateData::sync_values()",
                     "MpiComm::Bcast() failed!" );
  }


} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateData<QUESO::GslVector,QUESO::GslMatrix>;
