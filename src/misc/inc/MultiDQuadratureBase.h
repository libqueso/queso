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

#ifndef UQ_MULTI_D_QUADRATURE_BASE_H
#define UQ_MULTI_D_QUADRATURE_BASE_H

#include <queso/BaseQuadrature.h>
#include <queso/VectorSubset.h>
#include <queso/SharedPtr.h>

namespace QUESO
{
  // Forward declarations
  class GslVector;
  class GslMatrix;

  //! Base class for multi-dimensional quadrature rules
  /*! Base class for numerical integration of functions of multiple variables
   *  using quadrature rules. Subclasses will implement the actual quadrature rule
   *  and should populate the positions and weights. The integration will be performed
   *  over the user-supplied domain.
   *
   *  The intended use of these rules is that the constructor builds the rule and populates
   *  the positions and weights and then the user can perform the sum \f$ \sum_{i=1}^N f(x_i) w_i \f$.
   */
  template <class V = GslVector, class M = GslMatrix>
  class MultiDQuadratureBase : public BaseQuadrature
  {
  public:

    MultiDQuadratureBase( const VectorSubset<V,M> & domain )
      : BaseQuadrature(),
        m_domain(domain)
    {}

    //! Pure virtual destructor, forcing this to be an abstract object.
    virtual ~MultiDQuadratureBase() =0;

    const std::vector<typename QUESO::SharedPtr<V>::Type> & positions() const
    { queso_assert(!m_positions.empty());
      return m_positions; }

    const VectorSubset<V,M> & getDomain() const
    { return m_domain; }

  protected:

    //! Locations of quadrature points
    /*! \todo We should make this a vector of std::unique_ptr when we go full C++11. */
    std::vector<typename QUESO::SharedPtr<V>::Type> m_positions;

    //! Domain over which the quadrature will be performed
    const VectorSubset<V,M> & m_domain;
  };

  template <class V, class M>
  inline
  MultiDQuadratureBase<V,M>::~MultiDQuadratureBase(){}

} // end namespace QUESO

#endif // UQ_MULTI_D_QUADRATURE_BASE_H
