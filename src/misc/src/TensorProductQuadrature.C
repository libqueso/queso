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

#include <queso/TensorProductQuadrature.h>
#include <queso/asserts.h>
#include <queso/1DQuadrature.h>
#include <queso/MultiDimensionalIndexing.h>
#include <queso/VectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO
{
  template<class V,class M>
  TensorProductQuadrature<V,M>::TensorProductQuadrature( const VectorSubset<V,M> & domain,
                                                         const std::vector<QUESO::SharedPtr<Base1DQuadrature>::Type> & q_rules )
    : MultiDQuadratureBase<V,M>(domain)
  {
    const unsigned int dim = domain.vectorSpace().dimGlobal();

    queso_require_equal_to_msg(dim, q_rules.size(), "Mismatched quadrature rule size and vector space dimension!");

    // We figure the number of q_points in each dimension and the total
    unsigned int total_n_q_points = 1;
    std::vector<unsigned int> n_q_points(dim);
    for( unsigned int i = 0; i < dim; i++ )
      {
        n_q_points[i] = q_rules[i]->positions().size();
        total_n_q_points *= q_rules[i]->positions().size();
      }

    this->m_positions.resize(total_n_q_points,SharedPtr<GslVector>::Type());
    this->m_weights.resize(total_n_q_points);

    // The positions are just a tensor product of the positions for each of the 1-D rules.
    // Thus, we leverage the MultiDimensionalIndexing to move from "global" quadrature point
    // number to the index in each of the dimensions. Each of these dimensional indices tells
    // us which component of the corresponding 1-D rule to extract.
    // The final quadrature weight will be the product of each of the corresponding weights
    // for the 1-D rule, with the correspondence dictated by the "local" index in that dimension
    // (as with the position).
    std::vector<unsigned int> indices(dim);
    for( unsigned int q = 0; q < total_n_q_points; q++ )
      {
        MultiDimensionalIndexing::globalToCoord( q, n_q_points, indices );

        // Make sure indices size isn't something we don't expect after population
        queso_assert_equal_to(indices.size(),dim);

        // Create new vector from our vector space
        // We are taking ownership of this pointer.
        typename QUESO::SharedPtr<V>::Type domain_vec(domain.vectorSpace().newVector());

        // Now populate the weights and positions
        double weight = 1.0;
        for( unsigned int i = 0; i < dim; i++ )
          {
            unsigned int idx = indices[i];
            (*domain_vec)[i] = q_rules[i]->positions()[idx];
            weight *= q_rules[i]->weights()[idx];
          }

        // m_positions is taking ownership of the V*
        this->m_positions[q] = domain_vec;
        this->m_weights[q] = weight;
      }
  }

  // Instantiate
  template class TensorProductQuadrature<GslVector,GslMatrix>;

} // end namespace QUESO
