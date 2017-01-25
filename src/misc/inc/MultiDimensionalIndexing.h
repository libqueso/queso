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

#ifndef UQ_MULTI_DIMENSIONAL_INDEXING_H
#define UQ_MULTI_DIMENSIONAL_INDEXING_H

// C++
#include <vector>

namespace QUESO
{
  class MultiDimensionalIndexing
  {
  public:

    MultiDimensionalIndexing(){};

    ~MultiDimensionalIndexing(){};

    //! Map coordinate indices to a singal global index.
    /*! e.g. in 3-D, you pass in i,j,k and the n_points in each of the 3 directions.
        This function will return the corresponding global index to which the value
        at i,j,k can be indexed. Ordering must be consistent between coord_indices
        and n_points. Must be ordered in increasing dimension. e.g. x,y,z,...
        The user shouldn't need to call this method, this is public mainly to
        facilitate testing. */
    static unsigned int coordToGlobal( const std::vector<unsigned int>& coord_indices,
                                       const std::vector<unsigned int>& n_points );

    //! Inverse of coordToGlobal map
    /*! Given the global index and the n_points in each direction, we
        back out the coordinate indices and return them in coord_indices. */
    static void globalToCoord( unsigned int global,
                               const std::vector<unsigned int>& n_points,
                               std::vector<unsigned int>& coord_indices );

  private:

    //! Helper function
    /*! coordToGlobal computes: i + j*n_points[0] + k*n_points[0]*n_points[1] + ...
        To compute inverse we need the product of n_points terms, depending on which
        term we are computing. So, given n_points and the current term, we
        compute the correct number of products of n_points. */
    static unsigned int compute_npoints_factor( const std::vector<unsigned int>& n_points,
                                                unsigned int term );
  };

} // end namespace QUESO

#endif // UQ_MULTI_DIMENSIONAL_INDEXING_H
