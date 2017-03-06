//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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
#include <cstring>
#include <queso/MpiComm.h>
#include <queso/Environment.h>

#ifdef QUESO_HAS_TRILINOS
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#endif

namespace QUESO {

// QUESO MpiComm MPI Constructor ------------------
MpiComm::MpiComm(const BaseEnvironment& env, RawType_MPI_Comm inputRawComm)
  :
  m_env          (env),
#ifdef QUESO_HAS_TRILINOS
  m_epetraComm( new Epetra_MpiComm(inputRawComm) ),
#endif
  m_rawComm      (inputRawComm),
  m_worldRank    (-1),
  m_myPid        (-1),
  m_numProc      (-1)
{
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Comm_rank(inputRawComm,&m_worldRank);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_rank() on full rank");

  mpiRC = MPI_Comm_rank(inputRawComm,&m_myPid);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_rank() on inputRawComm");

  mpiRC = MPI_Comm_size(inputRawComm,&m_numProc);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_size() on inputRawComm");
#else
  m_worldRank = 0;
  m_myPid     = 0;
  m_numProc   = 1;
#endif
}

MpiComm::MpiComm(const BaseEnvironment& env)
  :
  m_env          (env),
#ifdef QUESO_HAS_TRILINOS
  m_epetraComm( new Epetra_SerialComm() ),
#endif
  m_rawComm      (RawValue_MPI_COMM_SELF),
  m_worldRank    (0),
  m_myPid        (0),
  m_numProc      (1)
{
}

// Copy constructor ---------------------------------
MpiComm::MpiComm(const MpiComm& src)
  :
  m_env          (src.m_env)
#ifdef QUESO_HAS_TRILINOS
  ,
  m_epetraComm(NULL)
#endif
{
  this->copy(src);
}

// Destructor ---------------------------------------
MpiComm::~MpiComm()
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraComm;
  m_epetraComm = NULL;
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
#ifdef QUESO_HAS_MPI
RawType_MPI_Comm
MpiComm::Comm() const
{
#ifdef QUESO_HAS_TRILINOS
#ifdef QUESO_HAS_MPI
  return dynamic_cast<Epetra_MpiComm *>(m_epetraComm)->Comm();
#endif
#endif
  return m_rawComm;
}
#endif  // QUESO_HAS_MPI

// --------------------------------------------------
int
MpiComm::MyPID() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraComm->MyPID();
#endif
  return m_myPid;
}
// --------------------------------------------------
int
MpiComm::NumProc() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraComm->NumProc();
#endif
  return m_numProc;
}
// Methods overridden from Comm ---------------------

void
MpiComm::Allreduce(void* sendbuf, void* recvbuf, int count, RawType_MPI_Datatype datatype, RawType_MPI_Op op, const char* whereMsg, const char* whatMsg) const
{
  queso_deprecated();
  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    int mpiRC = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
}

template <typename T>
void
MpiComm::Allreduce(const T* sendbuf, T* recvbuf, int count, RawType_MPI_Op op,
                   const char* whereMsg, const char* whatMsg) const
{
  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    T * sendbuf_noconst = const_cast<T *>(sendbuf);
    int mpiRC = MPI_Allreduce(sendbuf_noconst, recvbuf, count, StandardType<T>(sendbuf), op, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
  else {
    size_t dataTypeSize = sizeof(T);
    size_t dataTotal = dataTypeSize*count;
    std::memcpy(recvbuf, sendbuf, dataTotal);
  }
}
//--------------------------------------------------
void
MpiComm::Barrier() const // const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraComm->Barrier();
#endif

  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    int mpiRC = MPI_Barrier(m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "mpiRC indicates failure");  // whatMsg);
#endif
  }

  return;
}
//--------------------------------------------------
void
MpiComm::Bcast(void* buffer, int count, RawType_MPI_Datatype datatype, int root, const char* whereMsg, const char* whatMsg) const
{
  if (NumProc() > 1) {  // Necesarrily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    int mpiRC = MPI_Bcast(buffer, count, datatype, root, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
}
//--------------------------------------------------
void
MpiComm::Gather(
  void* sendbuf, int sendcnt, RawType_MPI_Datatype sendtype,
  void* recvbuf, int recvcount, RawType_MPI_Datatype recvtype,
  int root,
  const char* whereMsg, const char* whatMsg) const
{
  queso_deprecated();
  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype,
    //                void *recvbuf, int recvcount, MPI_Datatype recvtype,
    //                int root, MPI_Comm comm )
    int mpiRC = MPI_Gather(sendbuf, sendcnt, sendtype,
                           recvbuf, recvcount, recvtype,
                           root, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
}

template <typename T>
void
MpiComm::Gather(const T * sendbuf, int sendcnt, T * recvbuf, int recvcount,
                int root, const char* whereMsg, const char* whatMsg) const
{
  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype,
    //                void *recvbuf, int recvcount, MPI_Datatype recvtype,
    //                int root, MPI_Comm comm )
    T * sendbuf_noconst = const_cast<T *>(sendbuf);
    int mpiRC = MPI_Gather(sendbuf_noconst, sendcnt, StandardType<T>(sendbuf),
                           recvbuf, recvcount, StandardType<T>(sendbuf),
                           root, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
  else {
    size_t dataTypeSize = sizeof(T);
    size_t sendTotal = dataTypeSize*sendcnt;
    size_t recvTotal = dataTypeSize*recvcount;
    if (sendTotal != recvTotal) {
      std::cerr << "MpiCommClass::Gather()"
                << ": sendTotal != recvTotal"
                << std::endl;
    }
    queso_require_equal_to_msg(sendTotal, recvTotal, whatMsg);
    std::memcpy(recvbuf, sendbuf, sendTotal);
  }
}
//--------------------------------------------------
void
MpiComm::Gatherv(
  void* sendbuf, int sendcnt, RawType_MPI_Datatype sendtype,
  void* recvbuf, int* recvcnts, int* displs, RawType_MPI_Datatype recvtype,
  int root,
  const char* whereMsg, const char* whatMsg) const
{
  queso_deprecated();
  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    //int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
    //                void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype,
    //                int root, MPI_Comm comm )
    int mpiRC = MPI_Gatherv(sendbuf, sendcnt, sendtype,
                            recvbuf, recvcnts, displs, recvtype,
                            root, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
}

template <typename T>
void
MpiComm::Gatherv(const T * sendbuf, int sendcnt, T * recvbuf, int * recvcnts,
                 int * displs, int root, const char * whereMsg,
                 const char * whatMsg) const
{
  if (NumProc() > 1) {  // Necessarily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    //int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
    //                void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype,
    //                int root, MPI_Comm comm )
    T * sendbuf_noconst = const_cast<T *>(sendbuf);
    int mpiRC = MPI_Gatherv(sendbuf_noconst, sendcnt, StandardType<T>(sendbuf),
                            recvbuf, recvcnts, displs,
                            StandardType<T>(recvbuf), root, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
  else {
    size_t dataTypeSize = sizeof(T);
    size_t sendTotal = dataTypeSize*sendcnt;
    size_t recvTotal = dataTypeSize*recvcnts[0];
    if (sendTotal != recvTotal) {
      std::cerr << "MpiCommClass::Gatherv()"
                << ": sendTotal != recvTotal"
                << std::endl;
    }
    queso_require_equal_to_msg(sendTotal, recvTotal, whatMsg);
    std::memcpy(recvbuf, sendbuf, sendTotal);
  }
}
//--------------------------------------------------
void
MpiComm::Recv(
  void* buf, int count, RawType_MPI_Datatype datatype, int source, int tag, RawType_MPI_Status* status,
  const char* whereMsg, const char* whatMsg) const
{
  if (NumProc() > 1) {  // Necesarrily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    int mpiRC = MPI_Recv(buf, count, datatype, source, tag, m_rawComm, status);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
}
//--------------------------------------------------
void
MpiComm::Send(
  void* buf, int count, RawType_MPI_Datatype datatype, int dest, int tag,
  const char* whereMsg, const char* whatMsg) const
{
  if (NumProc() > 1) {  // Necesarrily true if QUESO_HAS_MPI
#ifdef QUESO_HAS_MPI
    int mpiRC = MPI_Send(buf, count, datatype, dest, tag, m_rawComm);
    queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, whatMsg);
#endif
  }
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
const Epetra_Comm&
MpiComm::epetraMpiComm() const
{
  return *m_epetraComm;
}
#endif

// Private methods----------------------------------
void
MpiComm::copy(const MpiComm& src)
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraComm;
  m_epetraComm = src.m_epetraComm->Clone();
#endif
  m_rawComm   = src.m_rawComm;
  m_worldRank = src.m_worldRank;
  m_myPid     = src.m_myPid;
  m_numProc   = src.m_numProc;

  return;
}
// -------------------------------------------------

// Explicit template function instantiations
template void MpiComm::Allreduce<int>(const int *,
                                      int *,
                                      int,
                                      RawType_MPI_Op,
                                      const char *,
                                      const char*) const;
template void MpiComm::Allreduce<char>(const char *,
                                       char *,
                                       int,
                                       RawType_MPI_Op,
                                       const char *,
                                       const char*) const;
template void MpiComm::Allreduce<unsigned int>(const unsigned int *,
                                               unsigned int *,
                                               int,
                                               RawType_MPI_Op,
                                               const char *,
                                               const char *) const;
template void MpiComm::Allreduce<double>(const double *,
                                         double *,
                                         int,
                                         RawType_MPI_Op,
                                         const char *,
                                         const char *) const;

template void MpiComm::Gather<int>(const int * sendbuf,
                                   int sendcnt,
                                   int * recvbuf,
                                   int recvcount,
                                   int root,
                                   const char * whereMsg,
                                   const char * whatMsg) const;

template void MpiComm::Gather<char>(const char * sendbuf,
                                    int sendcnt,
                                    char * recvbuf,
                                    int recvcount,
                                    int root,
                                    const char * whereMsg,
                                    const char * whatMsg) const;

template void MpiComm::Gather<unsigned int>(const unsigned int * sendbuf,
                                            int sendcnt,
                                            unsigned int * recvbuf,
                                            int recvcount,
                                            int root,
                                            const char * whereMsg,
                                            const char * whatMsg) const;

template void MpiComm::Gather<double>(const double * sendbuf,
                                      int sendcnt,
                                      double * recvbuf,
                                      int recvcount,
                                      int root,
                                      const char * whereMsg,
                                      const char * whatMsg) const;

template void MpiComm::Gatherv<int>(const int * sendbuf,
                                    int sendcnt,
                                    int * recvbuf,
                                    int * recvcnts,
                                    int * displs,
                                    int root,
                                    const char * whereMsg,
                                    const char * whatMsg) const;

template void MpiComm::Gatherv<char>(const char * sendbuf,
                                     int sendcnt,
                                     char * recvbuf,
                                     int * recvcnts,
                                     int * displs,
                                     int root,
                                     const char * whereMsg,
                                     const char * whatMsg) const;

template void MpiComm::Gatherv<unsigned int>(const unsigned int * sendbuf,
                                             int sendcnt,
                                             unsigned int * recvbuf,
                                             int * recvcnts,
                                             int * displs,
                                             int root,
                                             const char * whereMsg,
                                             const char * whatMsg) const;

template void MpiComm::Gatherv<double>(const double * sendbuf,
                                       int sendcnt,
                                       double * recvbuf,
                                       int * recvcnts,
                                       int * displs,
                                       int root,
                                       const char * whereMsg,
                                       const char * whatMsg) const;

}  // End namespace QUESO
