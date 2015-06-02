//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <unistd.h>
#include <queso/MpiComm.h>
#include <queso/Environment.h>

namespace QUESO {

// QUESO MpiComm MPI Constructor ------------------
MpiComm::MpiComm(const BaseEnvironment& env, RawType_MPI_Comm inputRawComm)
  :
  m_env          (env),
#ifdef QUESO_HAS_TRILINOS
  m_epetraMpiComm( new Epetra_MpiComm(inputRawComm) ),
#endif
  m_rawComm      (inputRawComm),
  m_worldRank    (-1),
  m_myPid        (-1),
  m_numProc      (-1)
{
  int mpiRC = MPI_Comm_rank(inputRawComm,&m_worldRank);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_rank() on full rank");

  mpiRC = MPI_Comm_rank(inputRawComm,&m_myPid);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_rank() on inputRawComm");

  mpiRC = MPI_Comm_size(inputRawComm,&m_numProc);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_size() on inputRawComm");
}

// Copy constructor ---------------------------------
MpiComm::MpiComm(const MpiComm& src)
  :
  m_env          (src.m_env)
#ifdef QUESO_HAS_TRILINOS
  ,
  m_epetraMpiComm(NULL)
#endif
{
  this->copy(src);
}

// Destructor ---------------------------------------
MpiComm::~MpiComm()
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMpiComm;
  m_epetraMpiComm = NULL;
#endif
}

// --------------------------------------------------
// Set methodos -------------------------------------
MpiComm&
MpiComm::operator=(const MpiComm& rhs)
{
  this->copy(rhs);
  return *this;
}

// Attribute access methods -------------------------
RawType_MPI_Comm
MpiComm::Comm() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->Comm();
#endif
  return m_rawComm;
}
// --------------------------------------------------
int
MpiComm::MyPID() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->MyPID();
#endif
  return m_myPid;
}
// --------------------------------------------------
int
MpiComm::NumProc() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->NumProc();
#endif
  return m_numProc;
}
// Methods overridden from Comm ---------------------

void
MpiComm::Allreduce(void* sendbuf, void* recvbuf, int count, RawType_MPI_Datatype datatype, RawType_MPI_Op op, const char* whereMsg, const char* whatMsg) const
{
  int mpiRC = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, m_rawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);

  return;
}
//--------------------------------------------------
void
MpiComm::Barrier() const // const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->Barrier();
#endif
  int mpiRC = MPI_Barrier(m_rawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "mpiRC indicates failure");  // whatMsg);
  return;
}
//--------------------------------------------------
void
MpiComm::Bcast(void* buffer, int count, RawType_MPI_Datatype datatype, int root, const char* whereMsg, const char* whatMsg) const
{
  int mpiRC = MPI_Bcast(buffer, count, datatype, root, m_rawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
  return;
}
//--------------------------------------------------
void
MpiComm::Gather(
  void* sendbuf, int sendcnt, RawType_MPI_Datatype sendtype,
  void* recvbuf, int recvcount, RawType_MPI_Datatype recvtype,
  int root,
  const char* whereMsg, const char* whatMsg) const
{
  //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype,
  //                void *recvbuf, int recvcount, MPI_Datatype recvtype,
  //                int root, MPI_Comm comm )
  int mpiRC = MPI_Gather(sendbuf, sendcnt, sendtype,
                         recvbuf, recvcount, recvtype,
                         root, m_rawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
  return;
}
//--------------------------------------------------
void
MpiComm::Gatherv(
  void* sendbuf, int sendcnt, RawType_MPI_Datatype sendtype,
  void* recvbuf, int* recvcnts, int* displs, RawType_MPI_Datatype recvtype,
  int root,
  const char* whereMsg, const char* whatMsg) const
{
  //int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
  //                void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype,
  //                int root, MPI_Comm comm )
  int mpiRC = MPI_Gatherv(sendbuf, sendcnt, sendtype,
                          recvbuf, recvcnts, displs, recvtype,
                          root, m_rawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
  return;
}
//--------------------------------------------------
void
MpiComm::Recv(
  void* buf, int count, RawType_MPI_Datatype datatype, int source, int tag, RawType_MPI_Status* status,
  const char* whereMsg, const char* whatMsg) const
{
  int mpiRC = MPI_Recv(buf, count, datatype, source, tag, m_rawComm, status);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
  return;
}
//--------------------------------------------------
void
MpiComm::Send(
  void* buf, int count, RawType_MPI_Datatype datatype, int dest, int tag,
  const char* whereMsg, const char* whatMsg) const
{
  int mpiRC = MPI_Send(buf, count, datatype, dest, tag, m_rawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
  return;
}
// Misc methods ------------------------------------
void
MpiComm::syncPrintDebugMsg(const char* msg, unsigned int msgVerbosity, unsigned int numUSecs) const
{
  if (m_env.syncVerbosity() >= msgVerbosity) {
    this->Barrier();
    for (int i = 0; i < this->NumProc(); ++i) {
      if (i == this->MyPID()) {
        std::cout << msg
                  << ": fullRank "       << m_env.fullRank()
                  << ", subEnvironment " << m_env.subId()
                  << ", subRank "        << m_env.subRank()
                  << ", inter0Rank "     << m_env.inter0Rank()
                  << std::endl;
      }
      usleep(numUSecs);
      this->Barrier();
    }
    //if (this->fullRank() == 0) std::cout << "Sleeping " << numUSecs << " microseconds..."
    //                                     << std::endl;
    //usleep(numUSecs);
    this->Barrier();
  }

  return;
}
// -------------------------------------------------
#ifdef QUESO_HAS_TRILINOS
const Epetra_MpiComm&
MpiComm::epetraMpiComm() const
{
  return *m_epetraMpiComm;
}
#endif

// Private methods----------------------------------
void
MpiComm::copy(const MpiComm& src)
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraMpiComm;
  m_epetraMpiComm = new Epetra_MpiComm(*src.m_epetraMpiComm);
#endif
  m_rawComm   = src.m_rawComm;
  m_worldRank = src.m_worldRank;
  m_myPid     = src.m_myPid;
  m_numProc   = src.m_numProc;

  return;
}
// -------------------------------------------------

}  // End namespace QUESO
