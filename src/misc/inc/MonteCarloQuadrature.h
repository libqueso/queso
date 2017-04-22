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

#ifndef UQ_MONTE_CARLO_QUADRATURE_H
#define UQ_MONTE_CARLO_QUADRATURE_H

#include <queso/MultiDQuadratureBase.h>

namespace QUESO
{
  // Forward declarations
  class GslVector;
  class GslMatrix;

  //! Numerical integration using Monte Carlo
  /*!
   * Monte Carlo integration approximates an integral by randomly sampling the function. Namely, if
   * \f$ <f> = 1/|\Omega| \int_{\Omega} f(x) \; dx \f$ is the average, then
   * \f[ \int_{\Omega} f(x) \; dx = |\Omega| <f> \approx |\Omega| \frac{1}{N} \sum_{i=1}^N f(x_i) \f]
   * where \f$N\f$ is the number of samples. This class, then, uniformly samples over the user-specified
   * domain for the given number of samples. Then, to adhere to the quadrature interface, the positions
   * are accessible in the positions vector and the factor \f$|\Omega|/N \f$ is cached in the weights
   * vector.
   */
  template <class V = GslVector, class M = GslMatrix>
  class MonteCarloQuadrature : public MultiDQuadratureBase<V,M>
  {
  public:

    MonteCarloQuadrature( const VectorSubset<V,M> & domain,
                          unsigned int n_samples );

    virtual ~MonteCarloQuadrature(){}

  };

} // end namespace QUESO

#endif // UQ_MONTE_CARLO_QUADRATURE_H
