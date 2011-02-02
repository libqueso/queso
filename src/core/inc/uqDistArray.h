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
#include <uqMpiComm.h>
#include <uqMap.h>
#include <EpetraExt_DistArray.h>
template<typename T>
struct uqDistArrayClass
{
  typedef EpetraExt::DistArray<T> type;
};

#else // QUESO_HAS_TRILINOS

#include <uqMpiComm.h>
class uqDistArrayClass
{
public:
  uqDistArrayClass();
  uqDistArrayClass(unsigned int          numGlobalElements,
                   unsigned int          numNotUsed,
                   const uqMpiCommClass& comm);
  uqDistArrayClass(const uqDistArrayClass& src);
 ~uqDistArrayClass();

  uqDistArrayClass& operator= (const uqDistArrayClass& rhs);

  const uqMpiCommClass& Comm()              const;
  unsigned int          NumGlobalElements() const;
  unsigned int          NumMyElements()     const;

private:
  void copy             (const uqDistArrayClass& src);

  uqMpiCommClass m_comm;
  unsigned int   m_numGlobalElements;
  unsigned int   m_numMyElements;
};

#endif // QUESO_HAS_TRILINOS

#endif // __UQ_DIST_ARRAY_H__
