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

#ifndef UQ_INTERPOLATION_SURROGATE_DATA_H
#define UQ_INTERPOLATION_SURROGATE_DATA_H

#include <queso/BoxSubset.h>

namespace QUESO
{
  class GslVector;
  class GslMatrix;

  template<class V = GslVector, class M = GslMatrix>
  class InterpolationSurrogateData
  {
  public:

    InterpolationSurrogateData( const BoxSubset<V,M>& domain,
                                const std::vector<unsigned int>& n_points );

    ~InterpolationSurrogateData(){};

    const BoxSubset<V,M>& get_paramDomain() const
    { return this->m_domain; };

    const std::vector<unsigned int>& get_n_points() const
    { return this->m_n_points; };

    const std::vector<double>& get_values() const
    { return this->m_values; };

    std::vector<double>& get_values()
    { return this->m_values; };

    double get_value( unsigned int n ) const
    { queso_assert_less(n,this->m_values.size());
      return this->m_values[n]; };

    unsigned int n_values() const
    { return this->m_values.size(); };

    //! Set all values. Dimension must be consistent with internal m_values.
    /*! This does a full copy of the values vector. This is mainly for testing,
        users are encouraged to use the InterpolationSurrogateBuilder. */
    void set_values( std::vector<double>& values );

    void set_value( unsigned int n, double value );

    //! Dimension of parameter space
    unsigned int dim() const
    { return this->m_domain.vectorSpace().dimGlobal(); };

    //! Lower bound of domain along dimension dim
    double x_min( unsigned int dim ) const
    { return this->m_domain.minValues()[dim]; };

    //! Upper bound of domain along dimension dim
    double x_max( unsigned int dim ) const
    { return this->m_domain.maxValues()[dim]; };

    //! Spacing between points along dimension dim
    double spacing( unsigned int dim ) const;

    //! Get spatial coordinate value at the node index along coordinate direction dim
    double get_x( unsigned int dim, unsigned int index ) const;

    //! Sync values across all processors from root processor
    /*! m_values may only be set on one processor, so this method
        is used for MPI_Bcast'ing those values to all processors.
        This calls MPI_Bcast, which is collective, so this *MUST*
        be called from all processors in env.fullComm(). */
    void sync_values( unsigned int root );

  private:

    InterpolationSurrogateData();

    //! Helper function for constructor
    void check_dim_consistency() const;

    //! Helper function for sizing m_values
    void init_values( const std::vector<unsigned int>& n_points );

    //! Helper function for checking domain extent
    void check_domain_bounds() const;

    //! Parameter domain over which we use surrogate
    const BoxSubset<V,M>& m_domain;

    //! vector to store number of points in each coordinate direction
    /*! We assume that the spacing in each coordinate direction is constant
        so we only need to know the number of points. */
    const std::vector<unsigned int>& m_n_points;

    //! vector to store values to be interpolated
    /*! These will be stored in a particular ordering. Helper functions will provide
        mapping between grid coordinates and a global coordinate; the global coordinate
        will be used to index into the values.
        \todo We currently store all values reside on all processes. Generalization would
              be to partition values across processes allocated for the subenvironment. */
    std::vector<double> m_values;
  };

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_DATA_H
