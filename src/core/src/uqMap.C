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

#include <uqMap.h>

#include <uqMpiComm.h>
uqMapClass::uqMapClass()
  :
  m_uqMpiComm() // *(new uqMpiCommClass()) )
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqMapClass::constructor()",
                      "should not be called");
}

uqMapClass::uqMapClass(
  int                   numGlobalElements,
  int                   indexBase,
  const uqMpiCommClass& comm)
  :
  m_uqMpiComm        (comm),
#ifdef QUESO_HAS_TRILINOS
  m_epetraMap        ( new Epetra_Map(numGlobalElements,0,comm.epetraMpiComm()) )
#else
  m_numGlobalElements(numGlobalElements),
  m_indexBase        (indexBase),
  m_numMyElements    (numGlobalElements)
#endif
{
}

uqMapClass::uqMapClass(const uqMapClass& src)
  :
  m_uqMpiComm(src.m_uqMpiComm)
#ifdef QUESO_HAS_TRILINOS
  ,
  m_epetraMap(NULL)
#endif
{
  this->copy(src);
}

uqMapClass&
uqMapClass::operator=(const uqMapClass& rhs)
{
  this->copy(rhs);
  return *this;
}

uqMapClass::~uqMapClass()
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMap;
  m_epetraMap = NULL;
#else
  // Nothing to do
#endif
}

void
uqMapClass::copy(const uqMapClass& src)
{
  m_uqMpiComm         = src.m_uqMpiComm;
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMap;
  m_epetraMap         = new Epetra_Map(*src.m_epetraMap);
#else
  m_numGlobalElements = src.m_numGlobalElements;
  m_indexBase         = src.m_indexBase;
  m_numMyElements     = src.m_numMyElements;
#endif

  return;
}

const uqMpiCommClass&
uqMapClass::Comm() const
{
  return m_uqMpiComm;
}

#ifdef QUESO_HAS_TRILINOS
const Epetra_Map&
uqMapClass::epetraMap() const
{
  return *m_epetraMap;
}
#endif

int
uqMapClass::NumGlobalElements() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->NumGlobalElements();
#else
  return m_numGlobalElements;
#endif
}

int
uqMapClass::IndexBase() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->IndexBase();
#else
  return m_indexBase;
#endif
}

int
uqMapClass::NumMyElements() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->NumMyElements();
#else
  return m_numMyElements;
#endif
}

int
uqMapClass::MinMyGID() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->MinMyGID();
#else
  return 0;
#endif
}
