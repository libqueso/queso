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

#ifndef UQ_INTERPOLATION_SURROGATE_BUILDER_H
#define UQ_INTERPOLATION_SURROGATE_BUILDER_H

#include <queso/SurrogateBuilderBase.h>
#include <queso/InterpolationSurrogateData.h>

namespace QUESO
{
  //! Build interpolation-based surrogate
  /*! Interpolation surrogates assume a structured grid. So, given the domain
      and the number of equally space points desired in each dimension, this
      class will handle calling the user's model to populate the values needed
      by the surrogate objects. User should subclass this object and implement
      the evaluate_model method. */
  template<class V, class M>
  class InterpolationSurrogateBuilder : public SurrogateBuilderBase<V>
  {
  public:

    //! Constructor
    /*! We do not take a const& to the data because we want to compute and
        set the values directly. */
    InterpolationSurrogateBuilder( InterpolationSurrogateData<V,M>& data );

    virtual ~InterpolationSurrogateBuilder(){};

    //! Execute the user's model and populate m_values for the given n_points
    void build_values();

  protected:

    InterpolationSurrogateData<V,M>& m_data;

    void set_work_bounds( unsigned int& n_begin, unsigned int& n_end ) const;

    void set_value( unsigned int n, double value );

    void set_domain_vector( unsigned int n, V& domain_vector ) const;

  private:

    InterpolationSurrogateBuilder();

  };

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_BUILDER_H
