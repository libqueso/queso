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

#ifndef UQ_TENSOR_PRODUCT_QUADRATURE_H
#define UQ_TENSOR_PRODUCT_QUADRATURE_H

#include <queso/MultiDQuadratureBase.h>
#include <queso/SharedPtr.h>

namespace QUESO
{
  // Forward declarations
  class GslVector;
  class GslMatrix;
  class Base1DQuadrature;

  //! Numerical quadrature using a tensor product of Base1DQuadrature rules.
  /*!
   *  The user constructs a vector of Base1DQuadrature rules from which this class
   *  will construct a multi-dimensional quadrature rule using a tensor product of
   *  the one-dimensional rules. Thus, the dimension of this tensor product rule
   *  will be the size of the vector of one-dimensional rules provided by the user and,
   *  therefore, should match the dimension of the domain also supplied by the user.
   */
  template <class V = GslVector, class M = GslMatrix>
  class TensorProductQuadrature : public MultiDQuadratureBase<V,M>
  {
  public:

    TensorProductQuadrature( const VectorSubset<V,M> & domain,
                             const std::vector<QUESO::SharedPtr<Base1DQuadrature>::Type> & q_rules );

    virtual ~TensorProductQuadrature(){}

  };

} // end namespace QUESO

#endif // UQ_TENSOR_PRODUCT_QUADRATURE_H
