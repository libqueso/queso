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

#ifndef UQ_INTERPOLATION_SURROGATE_BASE_H
#define UQ_INTERPOLATION_SURROGATE_BASE_H

#include <queso/SurrogateBase.h>

namespace QUESO
{
  //! Base class for interpolation-based surrogates
  /*! This class is used for surrogoate approximations of a model using interpolation.
      Subclasses will define behavior of interpolant, but common assumptions
      are: 1. Bounded domain; future work may extend behavior to unbounded domains.
           2. Structured grids on the parameter domain. Subclasses will determine behavior
              of spacing (uniform vs. directionally-uniform, etc.). */
  template<class V, class M>
  class InterpolationSurrogateBase : public SurrogateBase<V,M>
  {
  public:

    InterpolationSurrogateBase(const VectorSet<V,M> & domain )
      : SurrogateBase<V,M>(domain)
    {};

    virtual ~InterpolationSurrogateBase(){};

    //! Map coordinate indices to a singal global index.
    /*! e.g. in 3-D, you pass in i,j,k and the n_points in each of the 3 directions.
        This function will return the corresponding global index to which the value
        at i,j,k can be indexed. Ordering must be consistent between coord_indices
        and n_points. Must be ordered in increasing dimension. e.g. x,y,z,...
        The user shouldn't need to call this method, this is public mainly to facilitate
        testing. */
    unsigned int coordToGlobal( const std::vector<unsigned int>& coord_indices,
                                const std::vector<unsigned int>& n_points ) const;

  protected:

    //! vector to store values to be interpolated
    /*! These will be stored in a particular ordering. Subclasses
        will provide a helper function to give the index into this
        vector.
        \todo We currently store all values reside on all processes. Generalization would
              be to partition values across processes allocated for the subenvironment. */
    std::vector<double> m_values;

  private:

    InterpolationSurrogateBase();

  };

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_BASE_H
