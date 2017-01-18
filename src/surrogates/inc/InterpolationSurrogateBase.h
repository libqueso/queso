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

#ifndef UQ_INTERPOLATION_SURROGATE_BASE_H
#define UQ_INTERPOLATION_SURROGATE_BASE_H

#include <queso/SurrogateBase.h>

namespace QUESO
{
  // Forward declarations
  template<typename V, typename M>
  class InterpolationSurrogateData;
  class GslVector;
  class GslMatrix;

  //! Base class for interpolation-based surrogates
  /*! This class is used for surrogoate approximations of a model using interpolation.
      These surrogates map an \f$ n\f$ dimensional parameter space to the reals.
      That is \f$ f: \mathbb{R}^n \rightarrow \mathbb{R} \f$.
      Subclasses will define behavior of interpolant, but common assumptions
      are:
           -# Bounded domain; future work may extend behavior to unbounded domains.
           -# Structured grids on the parameter domain. Subclasses will determine behavior
              of spacing (uniform vs. directionally-uniform, etc.).

      For the structured grid, we think of referencing each "node" in the box by its
      index coordinates (i,j,k,...), where i runs from (0, n_points[0]-1),
      j run from (0,n_points[1]-1), etc. We use this indexing to build maps. */
  template<class V = GslVector, class M = GslMatrix>
  class InterpolationSurrogateBase : public SurrogateBase<V>
  {
  public:

    //! Constructor
    /*! The data object should be already fully populated when constructing
        this object. The user can directly set values or use the
        InterpolationSurrogateBuilder object to populate the data. */
    InterpolationSurrogateBase(const InterpolationSurrogateData<V,M>& data);

    virtual ~InterpolationSurrogateBase(){};

  protected:

    const InterpolationSurrogateData<V,M>& m_data;

  private:

    InterpolationSurrogateBase();

  };

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_BASE_H
