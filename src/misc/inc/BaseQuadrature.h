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

#ifndef UQ_BASE_QUADRATURE_H
#define UQ_BASE_QUADRATURE_H

#include <vector>

#include <queso/asserts.h>

namespace QUESO
{
  //! Base class for quadrature rules
  /*! This serves as the base class for approximating integrals using numerical quadrature. Namely,
   *   \f[ \int_{\Omega} f(x) \; dx \approx \sum_{i=1}^N f(x_i) w_i \f]
   * This class contains the weights. Subclasses will have the positions since the position can be one-dimensional
   * (Base1DQuadrature) or multi-dimensional (MultiDQuadratureBase).
   */
  class BaseQuadrature
  {
  public:

    BaseQuadrature(){}

    //! Pure virtual destructor, forcing this to be an abstract object.
    virtual ~BaseQuadrature() =0;

    //! Array of the weights used in the numerical integration.
    const std::vector<double> & weights() const
    { queso_assert(!m_weights.empty());
      return m_weights; }

  protected:

    std::vector<double> m_weights;
  };

  inline
  BaseQuadrature::~BaseQuadrature(){}

} // end namespace QUESO

#endif // UQ_BASE_QUADRATURE_H
