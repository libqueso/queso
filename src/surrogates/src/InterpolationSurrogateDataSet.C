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

// This class
#include <queso/InterpolationSurrogateDataSet.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateDataSet<V,M>::InterpolationSurrogateDataSet(const BoxSubset<V,M> & domain,
                                                                    const std::vector<unsigned int>& n_points,
                                                                    unsigned int n_datasets )
    : m_datasets(n_datasets,NULL)
  {
    for( unsigned int s = 0; s < n_datasets; s++ )
      m_datasets[s] = new InterpolationSurrogateData<V,M>( domain, n_points );
  }

  template<class V, class M>
  InterpolationSurrogateDataSet<V,M>::~InterpolationSurrogateDataSet()
  {
    for( typename std::vector<InterpolationSurrogateData<V,M>*>::iterator it = m_datasets.begin();
         it != m_datasets.end(); ++it )
      delete *it;
  }

} // end namespace QUESO

// Instantiate
template class QUESO::InterpolationSurrogateDataSet<QUESO::GslVector,QUESO::GslMatrix>;
