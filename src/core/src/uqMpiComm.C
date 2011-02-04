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

uqMpiCommClass::uqMpiCommClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqMpiCommClass::constructor()",
                      "should not be called");
}

uqMpiCommClass::uqMpiCommClass(MPI_Comm inputRawComm)
  :
#ifdef QUESO_HAS_TRILINOS
  m_epetraMpiComm( new Epetra_MpiComm(inputRawComm) )
#else
  m_rawComm  (inputRawComm),
  m_worldRank(-1),
  m_myPid    (-1),
  m_numProc  (-1)
#endif
{
#ifdef QUESO_HAS_TRILINOS
#else
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
#endif
}

uqMpiCommClass::uqMpiCommClass(const uqMpiCommClass& src)
#ifdef QUESO_HAS_TRILINOS
  :
  m_epetraMpiComm(NULL)
#endif
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
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMpiComm;
  m_epetraMpiComm = NULL;
#else
  // Nothing to do
#endif
}

void
uqMpiCommClass::copy(const uqMpiCommClass& src)
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMpiComm;
  m_epetraMpiComm = new Epetra_MpiComm(*src.m_epetraMpiComm);
#else
  m_rawComm   = src.m_rawComm;
  m_worldRank = src.m_worldRank;
  m_myPid     = src.m_myPid;
  m_numProc   = src.m_numProc;
#endif

  return;
}

MPI_Comm
uqMpiCommClass::Comm() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->Comm();
#else
  return m_rawComm;
#endif
}

#ifdef QUESO_HAS_TRILINOS
const Epetra_MpiComm&
uqMpiCommClass::epetraMpiComm() const
{
  return *m_epetraMpiComm;
}
#endif

int
uqMpiCommClass::MyPID() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->MyPID();
#else
  return m_myPid;
#endif
}

int
uqMpiCommClass::NumProc() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->NumProc();
#else
  return m_numProc;
#endif
}

void
uqMpiCommClass::Barrier() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->Barrier();
#else
  int mpiRC = MPI_Barrier(m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqMpiCommClass::Barrier()",
                      "mpiRC indicates no success");
#endif
  return;
}

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
