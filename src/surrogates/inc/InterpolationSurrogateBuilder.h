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
#include <queso/BoxSubset.h>

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

    InterpolationSurrogateBuilder( const BoxSubset<V,M> & domain,
                                   const std::vector<unsigned int>& n_points );

    virtual ~InterpolationSurrogateBuilder(){};

    const BoxSubset<V,M>& paramDomain() const
    { return this->m_domain; };

    const std::vector<unsigned int>& n_points() const
    { return this->m_n_points; };

    const std::vector<double>& values() const
    { return this->m_values; };

  protected:

    //! Parameter domain over which we use surrogate
    const BoxSubset<V,M>& m_domain;

    //! vector to store number of points in each coordinate direction
    /*! We assume that the spacing in each coordinate direction is constant
        so we only need to know the number of points. Then we can use the coordToGlobal
        function to map coordinate indices to global index for access m_values. */
    const std::vector<unsigned int>& m_n_points;

    //! vector to store values to be interpolated
    std::vector<double> m_values;

    void init_values( const std::vector<unsigned int>& n_points );

  private:

    InterpolationSurrogateBuilder();

  };

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_BUILDER_H
