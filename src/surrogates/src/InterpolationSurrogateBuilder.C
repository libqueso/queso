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

// This class
#include <queso/InterpolationSurrogateBuilder.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateBuilder<V,M>::InterpolationSurrogateBuilder( const BoxSubset<V,M> & domain,
                                                                     const std::vector<unsigned int>& n_points )
    : SurrogateBuilderBase<V>(),
    m_domain(domain),
    m_n_points(n_points)
  {
    this->init_values(this->m_n_points);
  }

  template<class V, class M>
  void InterpolationSurrogateBuilder<V,M>::init_values( const std::vector<unsigned int>& n_points )
  {
    unsigned int n_total_points = 0;
    for( std::vector<unsigned int>::const_iterator it = n_points.begin();
         it != n_points.end(); ++it )
      {
        n_total_points += *it;
      }

    this->m_values.resize(n_total_points);
  }

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateBuilder<QUESO::GslVector,QUESO::GslMatrix>;
