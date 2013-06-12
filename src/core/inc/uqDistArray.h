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

#ifndef __UQ_DIST_ARRAY_H__
#define __UQ_DIST_ARRAY_H__

#include <uqDefines.h>
#include <uqMap.h>
#ifdef QUESO_HAS_TRILINOS
#include <EpetraExt_DistArray.h>
#endif
#include <ostream>


/*! \file uqDistArray.h
    \brief A class to store row-oriented multi-vectors of type T.
*/

/*! \class uqDistArrayClass
    \brief A class for partitioning vectors and matrices.
    
    Class DistArray allows the construction and usage of multi-vectors. 
    These vectors contain element of type T, and the storage is row-oriented 
    (instead of and not column-oriented; thus his class should be used as a container 
    for data, on which no BLAS-like operations are performed). 
      
    DistArray objects are identified by an uqMap and a RowSize. The map specifies 
    the distribution of the elements across the processors and therefore the number 
    of local elements, while the RowSize gives the total number of data assigned to
    each node. RowSize is constant for all elements.
*/

template<typename T>
class uqDistArrayClass
{
public:
  
    //! @name Constructor/Destructor methods
  //@{

  //! Default constructor. Do not call this directly.
  uqDistArrayClass();
  
  //! Constructor for a given inputMap and inputRowSize. 
  uqDistArrayClass(const uqMapClass& inputMap, 
                   const int         inputRowSize);
  
  //! Copy constructor
  uqDistArrayClass(const uqDistArrayClass<T>& src);
  
  //! Destructor
 ~uqDistArrayClass();
 //@}
 
  //! @name Set methods
  //@{
  //! Assignment operator.
  uqDistArrayClass<T>& operator= (const uqDistArrayClass<T>& rhs);
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

        void copy        (const uqDistArrayClass<T>& src);

  uqMapClass                   m_uqMap;
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
uqDistArrayClass<T>::uqDistArrayClass()
  :
  m_uqMap()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqDistArrayClass<T>::constructor()",
                      "should not be called");
}
// Constructor for a given inputMap and inputRowSize. 
template<typename T>
uqDistArrayClass<T>::uqDistArrayClass(
  const uqMapClass& inputMap, 
  const int         inputRowSize)
  :
  m_uqMap  (inputMap),
#ifdef QUESO_HAS_TRILINOS
  m_epetraDistArray(new EpetraExt::DistArray<T>(inputMap.epetraMap(),inputRowSize))
#else
  m_rowSize(inputRowSize)
#endif
{
  //std::cout << "Entering uqDistArrayClass<T>::constructor(1)" << std::endl;
#ifdef QUESO_HAS_TRILINOS
#else
  m_elements.resize(m_uqMap.NumGlobalElements());
  for (int i = 0; i < m_uqMap.NumGlobalElements(); ++i) {
    m_elements[i].resize(m_rowSize);
  }
#endif
  //std::cout << "Leaving uqDistArrayClass<T>::constructor(1)" << std::endl;
}
// Copy constructor ---------------------------------
template<typename T>
uqDistArrayClass<T>::uqDistArrayClass(const uqDistArrayClass<T>& src)
  :
  m_uqMap(src.m_uqMap)
#ifdef QUESO_HAS_TRILINOS
  ,
  m_epetraDistArray(NULL)
#endif
{
  //std::cout << "Entering uqDistArrayClass<T>::constructor(2)" << std::endl;
#ifdef QUESO_HAS_TRILINOS
#else
  m_elements.clear();
#endif
  this->copy(src);
  //std::cout << "Leaving uqDistArrayClass<T>::constructor(2)" << std::endl;
}

// Destructor ---------------------------------------
template<typename T>
uqDistArrayClass<T>::~uqDistArrayClass()
{
  //std::cout << "Entering uqDistArrayClass<T>::destructor()" << std::endl;
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraDistArray;
  m_epetraDistArray = NULL;
#else
  for (int i = 0; i < m_uqMap.NumGlobalElements(); ++i) {
    m_elements[i].clear();
  }
  m_elements.clear();
#endif
  //std::cout << "Leaving uqDistArrayClass<T>::destructor()" << std::endl;
}

// ---------------------------------------------------
// Set methods----------------------------------------
template<typename T>
uqDistArrayClass<T>&
uqDistArrayClass<T>::operator=(const uqDistArrayClass<T>& rhs)
{
#ifdef QUESO_HAS_TRILINOS
#else
  for (int i = 0; i < m_uqMap.NumGlobalElements(); ++i) {
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
uqDistArrayClass<T>::operator()(int localElementId, int colId)
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
uqDistArrayClass<T>::operator()(int localElementId, int colId) const
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
uqDistArrayClass<T>::GlobalLength() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->GlobalLength();
#else
  return m_uqMap.NumGlobalElements();
#endif
}
// ---------------------------------------------------
template<typename T>
int
uqDistArrayClass<T>::MyLength() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->MyLength();
#else
  return m_uqMap.NumMyElements();
#endif
}
// ---------------------------------------------------
template<typename T>
int
uqDistArrayClass<T>::RowSize() const
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
uqDistArrayClass<T>::print(std::ostream& os) const
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
std::ostream& operator<<(std::ostream& os, const uqDistArrayClass<T>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_DIST_ARRAY_H__