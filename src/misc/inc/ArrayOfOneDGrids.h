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

#ifndef UQ_ARRAY_OF_ONE_D_GRIDS_H
#define UQ_ARRAY_OF_ONE_D_GRIDS_H

#include <queso/OneDGrid.h>
#include <queso/VectorSpace.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\file ArrayOfOneDGrids.h
 * \brief Class to accommodate arrays of one-dimensional grid.
 *
 * \class ArrayOfOneDGrids
 * \brief Class to accommodate arrays of one-dimensional grid.
 *
 * Arrays of one-dimensional grids are necessary in the calculation, for instance, of CDFs
 * and MDF of vector functions (refer to BaseVectorCdf, BaseVectorMdf, and
 * derived classes).
 */
template<class V = GslVector, class M = GslMatrix>
class ArrayOfOneDGrids
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  ArrayOfOneDGrids(const char* prefix, const VectorSpace<V,M>& rowSpace);

  //! Destructor.
  ~ArrayOfOneDGrids();
  //@}

  //! @name Property methods
  //@{
  //! Returns the (vector) space to which the row belongs to.
  const VectorSpace<V,M>&     rowSpace       () const;

  //! Returns an array with the sizes of the grids.
  const V&  sizes          () const;

  //! Returns an array with the minimum position of each grid.
  const V&  minPositions   () const;

  //! Returns an array with the maximum position of each grid.
  const V&  maxPositions   () const;
  //@}

  //! @name Math methods
  //@{
//  void      setGrid        (unsigned int                 rowId,
//			    BaseOneDGrid<double>& oneDGrid);

  //! Sets an array of uniform grids.
  void      setUniformGrids(const V& sizesVec,
			    const V& minPositionsVec,
			    const V& maxPositionsVec);
  //@}

  //! @name Accessor methods
  //@{
  //! Returns the grid stored in the <c>rowId</c>-th position of the array of grids.
  const BaseOneDGrid<double>& grid           (unsigned int rowId) const;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the values of the array of grids (points).
  void      print          (std::ostream& os) const;
  friend std::ostream& operator<< (std::ostream& os,
      const ArrayOfOneDGrids<V,M>& obj)
  {
    obj.print(os);
    return os;
  }
  //@}

private:
  const BaseEnvironment&                  m_env;
        std::string                              m_prefix;
  const VectorSpace<V,M>&                 m_rowSpace;
  DistArray<BaseOneDGrid<double>*> m_oneDGrids;

  V* m_sizes;
  V* m_minPositions;
  V* m_maxPositions;

};

}  // End namespace QUESO

#endif // UQ_ARRAY_OF_ONE_D_GRIDS_H
