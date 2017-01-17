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

#ifndef UQ_UNIFORM_ONE_D_GRID_FUNCTION_H
#define UQ_UNIFORM_ONE_D_GRID_FUNCTION_H

#include <queso/Environment.h>
#include <queso/OneDGrid.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Uniform grid class
//*****************************************************
/*!\class UniformOneDGrid
 * \brief Class for accommodating uniform one-dimensional grids.*/

template<class T>
class UniformOneDGrid : public BaseOneDGrid<T> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Constructs a uniform 1D grid between \c minPosition and \c maxPosition, with \c size points.*/
  UniformOneDGrid(const BaseEnvironment& env,
                         const char*               prefix,
                               unsigned int        size,
                               T                   minPosition,
                               T                   maxPosition);
 //! Destructor
  ~UniformOneDGrid();
  //@}

  //! @name Accessor methods
  //@{
  //! Returns the position of the i-th point in the grid.
  T    operator[]    (unsigned int i)      const;
  //@}

  //! @name Mathematical methods
  //@{
  //! Grid size; the amount of points that defines the grid.
  unsigned int size          ()                    const;

  //! Finds the ID of an interval. See template specialization.
  /*! This function finds to which interval the parameter value belongs to.*/
  unsigned int findIntervalId(const T& paramValue) const;
  //@}

protected:
  using BaseOneDGrid<T>::m_env;
  using BaseOneDGrid<T>::m_prefix;

  unsigned int m_size;
  T            m_minPosition;
  T            m_maxPosition;
};

}  // End namespace QUESO

#endif // UQ_UNIFORM_ONE_D_GRID_FUNCTION_H
