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

#ifndef UQ_2D_ARRAY_OF_STUFF_H
#define UQ_2D_ARRAY_OF_STUFF_H

namespace QUESO {

/*! \file 2dArrayOfStuff.h
 * \brief A templated class for handling arrays of data
 *
 * \class TwoDArray
 * \brief Class for handling arrays of generic data.
 *
 * This class handles array of generic data (doubles, ints, strings, structs, etc.).*/

template <class T>
class TwoDArray
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  TwoDArray(unsigned int numRows, unsigned int numCols);

  //! Destructor.
 ~TwoDArray();
  //@}

  //! @name Attribute methods
  //@{
  //! Number of rows in the array.
  unsigned int numRows    ()                                 const;

  //! Number of columns in the array.
  unsigned int numCols    ()                                 const;

  //! Sets the data in a specific location.
  /*! This method sets the generic (templated) data \c info, in the position <c>(i,j)</c>
   * of the array. */
  void         setLocation(unsigned int i, unsigned int j, T* info);
  //@}

  //! @name Accessor methods
  //@{
  //! Returns data stored in a specific location (non-const).
  /*! This non-const method returns the generic (templated) data that is stored the
   * position  <c>(i,j)</c> of the array. */
  T&           operator   ()(unsigned int i, unsigned int j);

  //! Returns data stored in a specific location (const).
  /*! This const method returns the generic (templated) data that is stored the
   * position  <c>(i,j)</c> of the array. */
  const T&     operator   ()(unsigned int i, unsigned int j) const;
  //@}
private:
  unsigned int m_numRows;
  unsigned int m_numCols;

  std::vector<std::vector<T*>* > m_data;
};

}  // End namespace QUESO

#endif // UQ_2D_ARRAY_OF_STUFF_H
