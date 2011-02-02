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

#include <uqMpiComm.h>
#ifdef QUESO_HAS_TRILINOS

// 'uqMpiCommClass' is just an alias to the 'Epetra_Map' class of Trilinos

#else // QUESO_HAS_TRILINOS

uqMpiCommClass::uqMpiCommClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqMpiCommClass::constructor()"
                      "should not be called");
}

uqMpiCommClass::uqMpiCommClass(
  unsigned int          numGlobalElements,
  unsigned int          numNotUsed,
  const uqMpiCommClass& comm)
  :
  m_comm             (comm),
  m_numGlobalElements(numGlobalElements),
  m_numMyElements    (numGlobalElements)
{
  if (numNotUsed > 1) { // just to avoid compiler warning
    // Do nothing
  }
}

uqMpiCommClass::uqMpiCommClass(const uqMpiCommClass& src)
{
  this->copy(src);
}

uqMpiCommClass&
uqMpiCommClass::operator=(const uqMpiCommClass& rhs)
{
  this->copy(rhs);
  return *this;
}

uqMpiCommClass::~uqMpiCommClass()
{
  // Nothing to do
}

void
uqMpiCommClass::copy(const uqMpiCommClass& src)
{
  m_comm              = src.m_comm;
  m_numGlobalElements = src.m_numGlobalElements;
  m_numMyElements     = src.m_numMyElements;

  return;
}

const uqMpiCommClass&
uqMpiCommClass::Comm() const
{
  return m_comm;
}

unsigned int
uqMpiCommClass::NumGlobalElements() const
{
  return m_numGlobalElements;
}

unsigned int
uqMpiCommClass::NumMyElements() const
{
  return m_numMyElements;
}

#endif // QUESO_HAS_TRILINOS

int UQ_MPI_Barrier(MPI_Comm comm)
{
  int mpiRC = MPI_Barrier(comm);

  return mpiRC;
}

int UQ_MPI_Barrier(const uqMpiCommClass& comm)
{
  comm.Barrier();

  return 0;
}
