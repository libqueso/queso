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

// 'uqMpiCommClass' is just an alias to the 'Epetra_MpiComm' class of Trilinos

#else // QUESO_HAS_TRILINOS

uqMpiCommClass::uqMpiCommClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqMpiCommClass::constructor()",
                      "should not be called");
}

uqMpiCommClass::uqMpiCommClass(MPI_Comm inputRawComm)
  :
  m_rawComm  (inputRawComm),
  m_worldRank(-1),
  m_myPid    (-1),
  m_numProc  (-1)
{
  int mpiRC = MPI_Comm_rank(MPI_COMM_WORLD,&m_worldRank);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      "uqMpiCommClass::constructor()",
                      "failed MPI_Comm_rank() on MPI_COMM_WORLD");

  mpiRC = MPI_Comm_rank(inputRawComm,&m_myPid);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqMpiCommClass::constructor()",
                      "failed MPI_Comm_rank() on inputRawComm");

  mpiRC = MPI_Comm_size(inputRawComm,&m_numProc);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqMpiCommClass::constructor()",
                      "failed MPI_Comm_size() on inputRawComm");
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
  m_rawComm   = src.m_rawComm;
  m_worldRank = src.m_worldRank;
  m_myPid     = src.m_myPid;
  m_numProc   = src.m_numProc;

  return;
}

MPI_Comm
uqMpiCommClass::Comm() const
{
  return m_rawComm;
}

int
uqMpiCommClass::MyPID() const
{
  return m_myPid;
}

int
uqMpiCommClass::NumProc() const
{
  return m_numProc;
}

void
uqMpiCommClass::Barrier() const
{
  int mpiRC = MPI_Barrier(m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqMpiCommClass::Barrier()",
                      "mpiRC indicates no success");

  return;
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
