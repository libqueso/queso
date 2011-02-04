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

#ifndef __UQ_MPI_COMM_H__
#define __UQ_MPI_COMM_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_TRILINOS
#include <Epetra_MpiComm.h>
#endif

#include <mpi.h>

class uqMpiCommClass
{
public:
  uqMpiCommClass();
  uqMpiCommClass(MPI_Comm inputRawComm);
  uqMpiCommClass(const uqMpiCommClass& src);
 ~uqMpiCommClass();

  uqMpiCommClass& operator= (const uqMpiCommClass& rhs);

  MPI_Comm Comm     () const;
  int      MyPID    () const;
  int      NumProc  () const;

  void     Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                     const char* whereMsg, const char* whatMsg) const;
  void     Barrier  () const; // const char* whereMsg, const char* whatMsg) const;
  void     Bcast    (void* buffer, int count, MPI_Datatype datatype, int root,
                     const char* whereMsg, const char* whatMsg) const;
  void     Gather   (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                     void *recvbuf, int recvcount, MPI_Datatype recvtype, 
                     int root,
                     const char* whereMsg, const char* whatMsg) const;
  void     Gatherv  (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, 
                     int root,
                     const char* whereMsg, const char* whatMsg) const;
  void     Recv     (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Status *status,
                     const char* whereMsg, const char* whatMsg) const;
  void     Send     (void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                     const char* whereMsg, const char* whatMsg) const;

#ifdef QUESO_HAS_TRILINOS
  const Epetra_MpiComm& epetraMpiComm() const;
#endif

private:
  void     copy   (const uqMpiCommClass& src);

#ifdef QUESO_HAS_TRILINOS
  Epetra_MpiComm* m_epetraMpiComm;
#endif
  MPI_Comm m_rawComm;
  int      m_worldRank;
  int      m_myPid;
  int      m_numProc;
};

#endif // __UQ_MPI_COMM_H__
