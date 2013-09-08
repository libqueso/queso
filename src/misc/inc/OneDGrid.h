//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_ONE_D_GRID_FUNCTION_H__
#define __UQ_ONE_D_GRID_FUNCTION_H__

#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Classes to accommodate a one dimensional grid
//*****************************************************
/*!\file uqOneDGrid.h
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
  //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
};

template<class T>
BaseOneDGrid<T>::BaseOneDGrid(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix)+"grid")
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering BaseOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving BaseOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
BaseOneDGrid<T>::~BaseOneDGrid()
{
}

template <class T>
void
BaseOneDGrid<T>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_prefix << "_sub" << m_env.subIdString() << " = zeros(" << this->size()
     << ","                                                      << 1
     << ");"
     << std::endl;
  os << m_prefix << "_sub" << m_env.subIdString() << " = [";
  for (unsigned int j = 0; j < this->size(); ++j) {
    os << (*this)[j] << " ";
  }
  os << "];"
     << std::endl;

  return;
}

template <class T>
std::ostream& operator<< (std::ostream& os, const BaseOneDGrid<T>& obj)
{
  obj.print(os);
  return os;
}

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

// Constructor-------------------------------------------
template<class T>
UniformOneDGrid<T>::UniformOneDGrid(
  const BaseEnvironment& env,
  const char*               prefix,
        unsigned int        size,
        T                   minPosition,
        T                   maxPosition)
  :
  BaseOneDGrid<T>(env,prefix),
  m_size       (size),
  m_minPosition(minPosition),
  m_maxPosition(maxPosition)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering UniformOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving UniformOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor--------------------------------------------
template<class T>
UniformOneDGrid<T>::~UniformOneDGrid()
{
}
// Math methods------------------------------------------
template<class T>
unsigned int
UniformOneDGrid<T>::size() const
{
  return m_size;
}
//-------------------------------------------------------
template<class T>
T
UniformOneDGrid<T>::operator[](unsigned int i) const
{
  UQ_FATAL_TEST_MACRO(i >= m_size,
                      m_env.worldRank(),
                      "UniformOneDGrid<V,M>::operator[]",
                      "too large i");

  T ratio = ((T) i)/(((T)m_size)-1.); // IMPORTANT: Yes, '-1.'
  T position = (1.-ratio)*m_minPosition + ratio*m_maxPosition;
  return position;
}
//-------------------------------------------------------
template<class T>
unsigned int
UniformOneDGrid<T>::findIntervalId(const T& paramValue) const
{
  UQ_FATAL_TEST_MACRO((paramValue < m_minPosition) || (m_maxPosition < paramValue),
                      m_env.worldRank(),
                      "UniformOneDGrid<V,M>::findIntervalId[]",
                      "paramValue is out of domain");

  T ratio = (paramValue - m_minPosition)/(m_maxPosition - m_minPosition);
  unsigned int i = (unsigned int) (ratio*(m_size-1.));
  if ((i > 0                  ) && 
      ((*this)[i] > paramValue)) {
    i--;
  }

  return i;
}

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

template<class T>
StdOneDGrid<T>::StdOneDGrid(
  const BaseEnvironment& env,
  const char*                   prefix,
  const std::vector<T>&         points)
  :
  BaseOneDGrid<T>(env,prefix),
  m_points              (points)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering StdOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving StdOneDGrid<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
StdOneDGrid<T>::~StdOneDGrid()
{
}

template<class T>
unsigned int
StdOneDGrid<T>::size() const
{
  return m_points.size();
}

template<class T>
T
StdOneDGrid<T>::operator[](unsigned int i) const
{
  UQ_FATAL_TEST_MACRO(i >= m_points.size(),
                      m_env.worldRank(),
                      "StdOneDGrid<V,M>::operator[]",
                      "too large i");

  return m_points[i];
}

template<class T>
unsigned int
StdOneDGrid<T>::findIntervalId(const T& paramValue) const
{
  UQ_FATAL_TEST_MACRO((paramValue < m_points[0]) || (m_points[m_points.size()-1] < paramValue),
                      m_env.worldRank(),
                      "StdOneDGrid<V,M>::findIntervalId[]",
                      "paramValue is out of domain");

  unsigned int iMax = m_points.size();
  unsigned int i = 1; // Yes, '1'
  for (i = 1; i < iMax; ++i) { // Yes, '1'
    if (paramValue < m_points[i]) {
      i--;
      break;
    }
  }

  return i;
}

}  // End namespace QUESO

#endif // __UQ_ONE_D_GRID_FUNCTION_H__
