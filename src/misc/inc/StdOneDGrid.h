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

#ifndef UQ_STD_ONE_D_GRID_FUNCTION_H
#define UQ_STD_ONE_D_GRID_FUNCTION_H

#include <queso/Environment.h>
#include <queso/OneDGrid.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Std grid class
//*****************************************************
/*!\class StdOneDGrid
 * \brief Class for accommodating standard one-dimensional grids.
 *
 * This class implements a standard one-dimensional grid, which is required, for instance,
 * in the evaluation of the cumulative distribution function (CDF) of a random variable.
 */

template<class T>
class StdOneDGrid : public BaseOneDGrid<T> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  StdOneDGrid(const BaseEnvironment& env,
                     const char*                   prefix,
                     const std::vector<T>&         points);
 //! Destructor.
  ~StdOneDGrid();
  //@}

  //! @name Accessor methods
  //@{
  //! Returns the position of the i-th point in the grid.
  T            operator[]    (unsigned int i)      const;
  //@}

  //! @name Mathematical methods
  //@{
  //! Grid size; the amount of points which defines the grid.
  unsigned int size          ()                    const;

  //! Finds the ID of an interval. See template specialization.
  unsigned int findIntervalId(const T& paramValue) const;
  //@}

protected:
  using BaseOneDGrid<T>::m_env;
  using BaseOneDGrid<T>::m_prefix;

  std::vector<T> m_points;
};

}  // End namespace QUESO

#endif // UQ_STD_ONE_D_GRID_FUNCTION_H
