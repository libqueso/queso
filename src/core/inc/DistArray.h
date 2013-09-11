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

#ifndef UQ_DIST_ARRAY_H
#define UQ_DIST_ARRAY_H

#include <queso/Defines.h>
#include <queso/Map.h>
#ifdef QUESO_HAS_TRILINOS
#include <EpetraExt_DistArray.h>
#endif
#include <ostream>

namespace QUESO {

/*! \file uqDistArray.h
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

template<typename T>
class DistArray
{
public:
  
    //! @name Constructor/Destructor methods
  //@{

  //! Default constructor. Do not call this directly.
  DistArray();
  
  //! Constructor for a given inputMap and inputRowSize. 
  DistArray(const Map& inputMap,
                   const int         inputRowSize);
  
  //! Copy constructor
  DistArray(const DistArray<T>& src);
  
  //! Destructor
 ~DistArray();
 //@}
 
  //! @name Set methods
  //@{
  //! Assignment operator.
  DistArray<T>& operator= (const DistArray<T>& rhs);
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
  //@}

private:
//! Copies the array.

        void copy        (const DistArray<T>& src);

  Map                   m_Map;
#ifdef QUESO_HAS_TRILINOS
  EpetraExt::DistArray<T>*     m_epetraDistArray;
#else
  unsigned int                 m_rowSize;
  std::vector<std::vector<T> > m_elements;
#endif
};


// --------------------------------------------------
// Constructor/Destructor methods -------------------

// Default constructor ------------------------------
template<typename T>
DistArray<T>::DistArray()
  :
  m_Map()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "DistArray<T>::constructor()",
                      "should not be called");
}
// Constructor for a given inputMap and inputRowSize. 
template<typename T>
DistArray<T>::DistArray(
  const Map& inputMap,
  const int         inputRowSize)
  :
  m_Map  (inputMap),
#ifdef QUESO_HAS_TRILINOS
  m_epetraDistArray(new EpetraExt::DistArray<T>(inputMap.epetraMap(),inputRowSize))
#else
  m_rowSize(inputRowSize)
#endif
{
  //std::cout << "Entering DistArray<T>::constructor(1)" << std::endl;
#ifdef QUESO_HAS_TRILINOS
#else
  m_elements.resize(m_Map.NumGlobalElements());
  for (int i = 0; i < m_Map.NumGlobalElements(); ++i) {
    m_elements[i].resize(m_rowSize);
  }
#endif
  //std::cout << "Leaving DistArray<T>::constructor(1)" << std::endl;
}
// Copy constructor ---------------------------------
template<typename T>
DistArray<T>::DistArray(const DistArray<T>& src)
  :
  m_Map(src.m_Map)
#ifdef QUESO_HAS_TRILINOS
  ,
  m_epetraDistArray(NULL)
#endif
{
  //std::cout << "Entering DistArray<T>::constructor(2)" << std::endl;
#ifdef QUESO_HAS_TRILINOS
#else
  m_elements.clear();
#endif
  this->copy(src);
  //std::cout << "Leaving DistArray<T>::constructor(2)" << std::endl;
}

// Destructor ---------------------------------------
template<typename T>
DistArray<T>::~DistArray()
{
  //std::cout << "Entering DistArray<T>::destructor()" << std::endl;
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraDistArray;
  m_epetraDistArray = NULL;
#else
  for (int i = 0; i < m_Map.NumGlobalElements(); ++i) {
    m_elements[i].clear();
  }
  m_elements.clear();
#endif
  //std::cout << "Leaving DistArray<T>::destructor()" << std::endl;
}

// ---------------------------------------------------
// Set methods----------------------------------------
template<typename T>
DistArray<T>&
DistArray<T>::operator=(const DistArray<T>& rhs)
{
#ifdef QUESO_HAS_TRILINOS
#else
  for (int i = 0; i < m_Map.NumGlobalElements(); ++i) {
    m_elements[i].clear();
  }
  m_elements.clear();
#endif
  this->copy(rhs);
  return *this;
}

// ---------------------------------------------------
// Query methods -------------------------------------
template<typename T>
T&
DistArray<T>::operator()(int localElementId, int colId)
{
#ifdef QUESO_HAS_TRILINOS
  return (*m_epetraDistArray)(localElementId,colId);
#else
  return m_elements[localElementId][colId];
#endif
}
// ---------------------------------------------------
template<typename T>
const T&
DistArray<T>::operator()(int localElementId, int colId) const
{
#ifdef QUESO_HAS_TRILINOS
  return (*m_epetraDistArray)(localElementId,colId);
#else
  return m_elements[localElementId][colId];
#endif
}
// ---------------------------------------------------
template<typename T>
int
DistArray<T>::GlobalLength() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->GlobalLength();
#else
  return m_Map.NumGlobalElements();
#endif
}
// ---------------------------------------------------
template<typename T>
int
DistArray<T>::MyLength() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->MyLength();
#else
  return m_Map.NumMyElements();
#endif
}
// ---------------------------------------------------
template<typename T>
int
DistArray<T>::RowSize() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->RowSize();
#else
  return m_rowSize;
#endif
}
// ---------------------------------------------------
// I/O methods----------------------------------------
template<typename T>
void
DistArray<T>::print(std::ostream& os) const
{
#ifdef QUESO_HAS_TRILINOS
  os << *m_epetraDistArray;
#else
  os << "m_rowSize = "           << m_rowSize
     << ", m_elements.size() = " << m_elements.size()
     << std::endl;
#endif

  return;
}
// ---------------------------------------------------
template<typename T>
std::ostream& operator<<(std::ostream& os, const DistArray<T>& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // UQ_DIST_ARRAY_H
