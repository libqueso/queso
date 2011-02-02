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

#ifndef __UQ_MAP_H__
#define __UQ_MAP_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_TRILINOS

// 'uqMapClass' is just an alias to the 'Epetra_Map' class of Trilinos
#include <Epetra_Map.h>
typedef Epetra_Map uqMapClass ;

#else // QUESO_HAS_TRILINOS

#include <uqMpiComm.h>
class uqMapClass
{
public:
  uqMapClass();
  uqMapClass(unsigned int          numGlobalElements,
             unsigned int          numNotUsed,
             const uqMpiCommClass& comm);
  uqMapClass(const uqMapClass& src);
 ~uqMapClass();

  uqMapClass& operator= (const uqMapClass& rhs);

  const uqMpiCommClass& Comm()              const;
  unsigned int          NumGlobalElements() const;
  unsigned int          NumMyElements()     const;

private:
  void copy             (const uqMapClass& src);

  const uqMpiCommClass& m_comm;
  unsigned int          m_numGlobalElements;
  unsigned int          m_numMyElements;
};

#endif // QUESO_HAS_TRILINOS

#endif // __UQ_MAP_H__
