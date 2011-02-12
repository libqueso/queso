//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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

#ifdef QUESO_HAS_MPI

#include <mpi.h>
typedef MPI_Comm     uqRawType_MPI_Comm ;
typedef MPI_Group    uqRawType_MPI_Group ;
typedef MPI_Datatype uqRawType_MPI_Datatype ;
typedef MPI_Op       uqRawType_MPI_Op ;
typedef MPI_Status   uqRawType_MPI_Status ;
#define uqRawValue_MPI_COMM_SELF  MPI_COMM_SELF 
#define uqRawValue_MPI_IN_PLACE   MPI_IN_PLACE
#define uqRawValue_MPI_ANY_SOURCE MPI_ANY_SOURCE
#define uqRawValue_MPI_CHAR       MPI_CHAR
#define uqRawValue_MPI_INT        MPI_INT
#define uqRawValue_MPI_DOUBLE     MPI_DOUBLE
#define uqRawValue_MPI_UNSIGNED   MPI_UNSIGNED
#define uqRawValue_MPI_MIN        MPI_MIN
#define uqRawValue_MPI_MAX        MPI_MAX
#define uqRawValue_MPI_SUM        MPI_SUM

#else // QUESO_HAS_MPI

typedef int uqRawType_MPI_Comm ;
typedef int uqRawType_MPI_Group ;
typedef int uqRawType_MPI_Datatype ;
typedef int uqRawType_MPI_Op ;
typedef int uqRawType_MPI_Status ;
#define uqRawValue_MPI_COMM_SELF  0
#define uqRawValue_MPI_IN_PLACE   0
#define uqRawValue_MPI_ANY_SOURCE -1
#define uqRawValue_MPI_CHAR       0
#define uqRawValue_MPI_INT        1
#define uqRawValue_MPI_DOUBLE     2
#define uqRawValue_MPI_UNSIGNED   3
#define uqRawValue_MPI_MIN        0
#define uqRawValue_MPI_MAX        1
#define uqRawValue_MPI_SUM        2

#endif // QUESO_HAS_MPI

class uqBaseEnvironmentClass;

class uqMpiCommClass
{
public:
  uqMpiCommClass();
  uqMpiCommClass(const uqBaseEnvironmentClass& env, uqRawType_MPI_Comm inputRawComm);
  uqMpiCommClass(const uqMpiCommClass& src);
 ~uqMpiCommClass();

  uqMpiCommClass& operator= (const uqMpiCommClass& rhs);

  uqRawType_MPI_Comm Comm     () const;
  int                MyPID    () const;
  int                NumProc  () const;

  void               Allreduce(void* sendbuf, void* recvbuf, int count, uqRawType_MPI_Datatype datatype, uqRawType_MPI_Op op,
                               const char* whereMsg, const char* whatMsg) const;
  void               Barrier  () const; // const char* whereMsg, const char* whatMsg) const;
  void               Bcast    (void* buffer, int count, uqRawType_MPI_Datatype datatype, int root,
                               const char* whereMsg, const char* whatMsg) const;
  void               Gather   (void *sendbuf, int sendcnt, uqRawType_MPI_Datatype sendtype, 
                               void *recvbuf, int recvcount, uqRawType_MPI_Datatype recvtype, 
                               int root,
                               const char* whereMsg, const char* whatMsg) const;
  void               Gatherv  (void *sendbuf, int sendcnt, uqRawType_MPI_Datatype sendtype, 
                               void *recvbuf, int *recvcnts, int *displs, uqRawType_MPI_Datatype recvtype, 
                               int root,
                               const char* whereMsg, const char* whatMsg) const;
  void               Recv     (void *buf, int count, uqRawType_MPI_Datatype datatype, int source, int tag, uqRawType_MPI_Status *status,
                               const char* whereMsg, const char* whatMsg) const;
  void               Send     (void *buf, int count, uqRawType_MPI_Datatype datatype, int dest, int tag,
                               const char* whereMsg, const char* whatMsg) const;

  void               syncPrintDebugMsg(const char* msg, unsigned int msgVerbosity, unsigned int numUSecs) const;

#ifdef QUESO_HAS_TRILINOS
  const Epetra_MpiComm& epetraMpiComm() const;
#endif

private:
  void               copy          (const uqMpiCommClass& src);
#ifdef QUESO_HAS_MPI
#else
  size_t             sizeOfDataType(uqRawType_MPI_Datatype datatype, const char* whereMsg, const char* whatMsg) const;
#endif

  const uqBaseEnvironmentClass& m_env;
#ifdef QUESO_HAS_TRILINOS
  Epetra_MpiComm*               m_epetraMpiComm;
#endif
  uqRawType_MPI_Comm            m_rawComm;
  int                           m_worldRank;
  int                           m_myPid;
  int                           m_numProc;
};

#endif // __UQ_MPI_COMM_H__
