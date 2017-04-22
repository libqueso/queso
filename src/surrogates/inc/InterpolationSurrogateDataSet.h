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

#ifndef UQ_INTERPOLATION_SURROGATE_DATA_SET_H
#define UQ_INTERPOLATION_SURROGATE_DATA_SET_H

#include <queso/InterpolationSurrogateData.h>

namespace QUESO
{
  class GslVector;
  class GslMatrix;

  //! Container class for multiple, consistent InterpolationSurrogateData objects
  /*! This is mostly a convenience object to build and hold multiple
      InterpolationSurrogateData objects that will have the domian and
      n_points in each dimension. This is useful for building the surrogates
      using the InterpolationSurrogateBuilder. */
  template<class V = GslVector, class M = GslMatrix>
  class InterpolationSurrogateDataSet
  {
  public:

    InterpolationSurrogateDataSet( const BoxSubset<V,M>& domain,
                                   const std::vector<unsigned int>& n_points,
                                   unsigned int n_datasets );

    ~InterpolationSurrogateDataSet();

    const InterpolationSurrogateData<V,M>& get_dataset( unsigned int s ) const;

    InterpolationSurrogateData<V,M>& get_dataset( unsigned int s );

    unsigned int size() const
    { return m_datasets.size(); }

  private:

    InterpolationSurrogateDataSet();

    //! Data structure to hold all data sets
    std::vector<InterpolationSurrogateData<V,M>*> m_datasets;

  };

  template<class V, class M>
  inline
  const InterpolationSurrogateData<V,M>& InterpolationSurrogateDataSet<V,M>::get_dataset( unsigned int s ) const
  {
    queso_require_less( s, m_datasets.size() );
    queso_assert( m_datasets[s] );

    return *(m_datasets[s]);
  }

  template<class V, class M>
  inline
  InterpolationSurrogateData<V,M>& InterpolationSurrogateDataSet<V,M>::get_dataset( unsigned int s )
  {
    queso_require_less( s, m_datasets.size() );
    queso_assert( m_datasets[s] );

    return *(m_datasets[s]);
  }

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_DATA_SET_H
