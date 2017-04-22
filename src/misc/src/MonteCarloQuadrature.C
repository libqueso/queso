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

#include <queso/MonteCarloQuadrature.h>
#include <queso/UniformVectorRV.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO
{
  template <class V, class M>
  MonteCarloQuadrature<V,M>::MonteCarloQuadrature( const VectorSubset<V,M> & domain,
                                                   unsigned int n_samples )
    : MultiDQuadratureBase<V,M>(domain)
  {
    // We'll just uniformly sample over the domain to retrieve the quadrature
    // points.
    UniformVectorRV<V,M> uniform_rv("MonteCarloQuadrature_", //prefix
                                    this->m_domain);

    const BaseVectorRealizer< V, M > & realizer = uniform_rv.realizer();

    this->m_positions.resize(n_samples,SharedPtr<GslVector>::Type());

    for( unsigned int i = 0; i < n_samples; i++ )
      {
        // Create new vector from our vector space
        // We are taking ownership of this pointer.
        typename QUESO::SharedPtr<V>::Type domain_vec(domain.vectorSpace().newVector());

        // Draw the new point
        realizer.realization(*domain_vec);

        // This is now the quadrature point position
        this->m_positions[i] = domain_vec;
      }

    // The weights are just the scaling factor for the MC sum
    this->m_weights.resize(n_samples,domain.volume()/(double)(n_samples));
  }

  // Instantiate
  template class MonteCarloQuadrature<GslVector,GslMatrix>;

} // end namespace QUESO
