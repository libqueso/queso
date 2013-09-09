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

#include <uqMap.h>
#include <uqMpiComm.h>

namespace QUESO {

// --------------------------------------------------
// Constructor/Destructor methods -------------------

// Default constructor ------------------------------
MapClass::MapClass()
  :
  m_MpiComm() // *(new MpiCommClass()) )
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "MapClass::constructor()",
                      "should not be called");
}

// Usable constructor ------------------------------
MapClass::MapClass(
  int                   numGlobalElements,
  int                   indexBase,
  const MpiCommClass& comm)
  :
  m_MpiComm        (comm),
#ifdef QUESO_HAS_TRILINOS
  m_epetraMap        ( new Epetra_Map(numGlobalElements,0,comm.epetraMpiComm()) )
#else
  m_numGlobalElements(numGlobalElements),
  m_indexBase        (indexBase),
  m_numMyElements    (numGlobalElements)
#endif
{
}

// Copy constructor ---------------------------------
MapClass::MapClass(const MapClass& src)
  :
  m_MpiComm(src.m_MpiComm)
#ifdef QUESO_HAS_TRILINOS
  ,
  m_epetraMap(NULL)
#endif
{
  this->copy(src);
}

// Destructor ---------------------------------------
MapClass::~MapClass()
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMap;
  m_epetraMap = NULL;
#else
  // Nothing to do
#endif
}

// --------------------------------------------------
// Set methodos -------------------------------------
MapClass&
MapClass::operator=(const MapClass& rhs)
{
  this->copy(rhs);
  return *this;
}

// --------------------------------------------------
// Size, dimension and local ID accessor methods ----
int
MapClass::NumGlobalElements() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->NumGlobalElements();
#else
  return m_numGlobalElements;
#endif
}

// --------------------------------------------------
int
MapClass::IndexBase() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->IndexBase();
#else
  return m_indexBase;
#endif
}

// --------------------------------------------------
int
MapClass::NumMyElements() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->NumMyElements();
#else
  return m_numMyElements;
#endif
}

// --------------------------------------------------
int
MapClass::MinMyGID() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMap->MinMyGID();
#else
  return 0;
#endif
}


// --------------------------------------------------
// Misc methods -------------------------------------
const MpiCommClass&
MapClass::Comm() const
{
  return m_MpiComm;
}

#ifdef QUESO_HAS_TRILINOS
const Epetra_Map&
MapClass::epetraMap() const
{
  return *m_epetraMap;
}
#endif

// --------------------------------------------------
void
MapClass::copy(const MapClass& src)
{
  m_MpiComm         = src.m_MpiComm;
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

}  // End namespace QUESO
