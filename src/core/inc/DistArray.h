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

#ifndef UQ_DIST_ARRAY_H
#define UQ_DIST_ARRAY_H

#include <queso/Defines.h>
#include <queso/Map.h>
#ifdef QUESO_HAS_TRILINOS
#include <EpetraExt_DistArray.h>
#endif
#include <ostream>

namespace QUESO {

/*! \file DistArray.h
    \brief A class to store row-oriented multi-vectors of type T.
*/

/*! \class DistArray
    \brief A class for partitioning vectors and matrices.

    Class DistArray allows the construction and usage of multi-vectors.
    These vectors contain element of type T, and the storage is row-oriented
    (instead of and not column-oriented; thus his class should be used as a container
    for data, on which no BLAS-like operations are performed).

    DistArray objects are identified by an Map and a RowSize. The map specifies
    the distribution of the elements across the processors and therefore the number
    of local elements, while the RowSize gives the total number of data assigned to
    each node. RowSize is constant for all elements.
*/

template <typename T>
class DistArray
{
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Constructor for a given inputMap and inputRowSize.
  DistArray(const Map& inputMap,
                   const int         inputRowSize);

  //! Destructor
 ~DistArray();
 //@}

  //! @name Query methods
  //@{

  //! Returns a reference to the colId column component of the localElementId local element.
        T&   operator    ()(int localElementId, int colId);

  //! Returns a reference to the colId column component of the localElementId local element.(const)
  const T&   operator    ()(int localElementId, int colId) const;

  //! Returns the global length of the array.
  int  GlobalLength() const;

  //! Returns the length of the locally owned array.
  int  MyLength    () const;

  //! Returns the row size, that is, the amount of data associated with each element.
  int  RowSize     () const;

  //@}

  //! @name I/O methods
  //@{
  void print       (std::ostream& os) const;
  friend std::ostream & operator<<(std::ostream& os, const DistArray<T>& obj)
  {
    obj.print(os);

    return os;
  }
  //@}

private:
  //! Default constructor. Do not call this directly.
  DistArray();

  //! Copy constructor.  Private.
  DistArray(const DistArray<T>& src) : m_Map(src.m_Map) { }

  //! Assignment operator.
  DistArray<T>& operator=(const DistArray<T>& rhs);

  Map                   m_Map;
#ifdef QUESO_HAS_TRILINOS
  EpetraExt::DistArray<T>*     m_epetraDistArray;
#else
  unsigned int                 m_rowSize;
  std::vector<std::vector<T> > m_elements;
#endif
};

}  // End namespace QUESO

#endif // UQ_DIST_ARRAY_H
