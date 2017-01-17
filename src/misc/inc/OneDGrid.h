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

#ifndef UQ_ONE_D_GRID_FUNCTION_H
#define UQ_ONE_D_GRID_FUNCTION_H

#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Classes to accommodate a one dimensional grid
//*****************************************************
/*!\file OneDGrid.h
 * \brief Classes to accommodate a one dimensional grid.
 *
 * \class BaseOneDGrid
 * \brief Base class for accommodating one-dimensional grids.*/

//*****************************************************
// Base class
//*****************************************************
template<class T>
class BaseOneDGrid {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  BaseOneDGrid(const BaseEnvironment& env,
		      const char* prefix);
  //! Virtual destructor.
  virtual ~BaseOneDGrid();
  //@}
  //! @name Accessor methods
  //@{
  //! Returns the position of the i-th point in the grid. See template specialization.
  virtual T            operator[]    (unsigned int i)       const = 0;
  //@}
  //! @name Mathematical methods
  //@{
  //! Grid size; the amount of points which defines the grid. See template specialization.
  virtual unsigned int size          ()                     const = 0;

  //! Finds the ID of an interval. See template specialization.
  virtual unsigned int findIntervalId(const T& paramValue)  const = 0;
  //@}
  //! @name I/O methods
  //@{
  //! Prints the values of the grid points.
  void         print         (std::ostream& ofsvar) const;
  friend std::ostream& operator<< (std::ostream& os,
      const BaseOneDGrid<T>& obj)
  {
    obj.print(os);
    return os;
  }
  //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
};

}  // End namespace QUESO

#endif // UQ_ONE_D_GRID_FUNCTION_H
