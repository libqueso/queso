//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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
#ifdef QUESO_HAS_TRILINOS

// 'uqDistArrayClass<T>::type' is just an alias to the 'EpetraExt::DistArray<T>' class of Trilinos
//#include <uqMpiComm.h>
#include <uqMap.h>
#include <EpetraExt_DistArray.h>

template<typename T>
struct uqDistArrayClass
{
  typedef EpetraExt::DistArray<T> type;
};

#else // QUESO_HAS_TRILINOS

//#include <uqMpiComm.h>
#include <uqMap.h>
#include <ostream>

template<typename T>
class uqDistArrayCoreClass
{
public:
  uqDistArrayCoreClass();
  uqDistArrayCoreClass(const uqMapClass& inputMap, 
                       const int         inputRowSize);
  uqDistArrayCoreClass(const uqDistArrayCoreClass<T>& src);
 ~uqDistArrayCoreClass();

  uqDistArrayCoreClass<T>& operator= (const uqDistArrayCoreClass<T>& rhs);

        T&  operator    ()(int localElementId, int rowId);
  const T&  operator    ()(int localElementId, int rowId) const;
        int GlobalLength() const;
        int MyLength    () const;
        int RowSize     () const;

private:
  void copy             (const uqDistArrayCoreClass<T>& src);

  uqMapClass                   m_uqMap;
  unsigned int                 m_rowSize;
  std::vector<std::vector<T> > m_elements;
};

template<typename T>
uqDistArrayCoreClass<T>::uqDistArrayCoreClass()
  :
  m_uqMap()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqDistArrayCoreClass<T>::constructor()",
                      "should not be called");
}

template<typename T>
uqDistArrayCoreClass<T>::uqDistArrayCoreClass(
  const uqMapClass& inputMap, 
  const int         inputRowSize)
  :
  m_uqMap  (inputMap),
  m_rowSize(inputRowSize)
{
}

template<typename T>
uqDistArrayCoreClass<T>::uqDistArrayCoreClass(const uqDistArrayCoreClass<T>& src)
  :
  m_uqMap(src.m_uqMap)
{
  this->copy(src);
}

template<typename T>
uqDistArrayCoreClass<T>&
uqDistArrayCoreClass<T>::operator=(const uqDistArrayCoreClass<T>& rhs)
{
  this->copy(rhs);
  return *this;
}

template<typename T>
uqDistArrayCoreClass<T>::~uqDistArrayCoreClass()
{
  // Nothing to do
}

template<typename T>
void
uqDistArrayCoreClass<T>::copy(const uqDistArrayCoreClass<T>& src)
{
  m_uqMap   = src.m_uqMap;
  m_rowSize = src.m_rowSize;

  return;
}

template<typename T>
T&
uqDistArrayCoreClass<T>::operator()(int localElementId, int rowId)
{
  return m_elements[localElementId][rowId];
}

template<typename T>
const T&
uqDistArrayCoreClass<T>::operator()(int localElementId, int rowId) const
{
  return m_elements[localElementId][rowId];
}

template<typename T>
int
uqDistArrayCoreClass<T>::GlobalLength() const
{
  return m_uqMap.NumGlobalElements();
}

template<typename T>
int
uqDistArrayCoreClass<T>::MyLength() const
{
  return m_uqMap.NumMyElements();
}

template<typename T>
int
uqDistArrayCoreClass<T>::RowSize() const
{
  return m_rowSize;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const uqDistArrayCoreClass<T>& obj)
{
  return os;
}


template<typename T>
struct uqDistArrayClass
{
  typedef uqDistArrayCoreClass<T> type;
};

#endif // QUESO_HAS_TRILINOS

#endif // __UQ_DIST_ARRAY_H__
