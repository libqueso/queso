//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#ifndef __UQ_MAP_H__
#define __UQ_MAP_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_TRILINOS
#include <Epetra_Map.h>
#endif

#include <uqMpiComm.h>
class uqMapClass
{
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Default constructor. Do not call this directly.
  uqMapClass();
  uqMapClass(int                   numGlobalElements,
             int                   indexBase,
             const uqMpiCommClass& comm);
  uqMapClass(const uqMapClass& src);

  //! Destructor
 ~uqMapClass();
 //@}

  uqMapClass& operator= (const uqMapClass& rhs);

  const uqMpiCommClass& Comm             () const;
  int                   NumGlobalElements() const;
  int                   IndexBase        () const;//1st position in the global processor in my processors
  int                   NumMyElements    () const;
  int                   MinMyGID         () const;
#ifdef QUESO_HAS_TRILINOS
  const Epetra_Map&     epetraMap        () const;
#endif

private:
  void                  copy             (const uqMapClass& src);

  uqMpiCommClass m_uqMpiComm;
#ifdef QUESO_HAS_TRILINOS
  Epetra_Map*    m_epetraMap;
#else
  int            m_numGlobalElements;
  int            m_indexBase;
  int            m_numMyElements;
#endif
};

#endif // __UQ_MAP_H__
