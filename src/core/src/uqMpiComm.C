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

#include <uqMpiComm.h>
#include <uqEnvironment.h>

uqMpiCommClass::uqMpiCommClass()
  :
  m_env( *(new uqEmptyEnvironmentClass()) )
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqMpiCommClass::constructor()",
                      "should not be called");
}

uqMpiCommClass::uqMpiCommClass(const uqBaseEnvironmentClass& env, uqRawType_MPI_Comm inputRawComm)
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
#ifdef QUESO_HAS_MPI
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
#else
  m_worldRank = 0;
  m_myPid     = 0;
  m_numProc   = 1;
#endif
}

uqMpiCommClass::uqMpiCommClass(const uqMpiCommClass& src)
  :
  m_env          (src.m_env)
#ifdef QUESO_HAS_TRILINOS
  ,
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
#endif
}

void
uqMpiCommClass::copy(const uqMpiCommClass& src)
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

uqRawType_MPI_Comm
uqMpiCommClass::Comm() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->Comm();
#endif
  return m_rawComm;
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
#endif
  return m_myPid;
}

int
uqMpiCommClass::NumProc() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->NumProc();
#endif
  return m_numProc;
}

void
uqMpiCommClass::Allreduce(void* sendbuf, void* recvbuf, int count, uqRawType_MPI_Datatype datatype, uqRawType_MPI_Op op, const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#else
  size_t dataTypeSize = sizeOfDataType(datatype, whereMsg, whatMsg);
  size_t dataTotal = dataTypeSize*count;
  memcpy(recvbuf, sendbuf, dataTotal);
#endif

  return;
}

void
uqMpiCommClass::Barrier() const // const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraMpiComm->Barrier();
#endif
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Barrier(m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqMPICommClass::Barrier()", // whereMsg,
                      "mpiRC indicates failure");  // whatMsg);
#else
  // Nothing needs to be done
#endif
  return;
}

void
uqMpiCommClass::Bcast(void* buffer, int count, uqRawType_MPI_Datatype datatype, int root, const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Bcast(buffer, count, datatype, root, m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#else
  // Nothing needs to be done
#endif
  return;
}

void
uqMpiCommClass::Gather(
  void* sendbuf, int sendcnt, uqRawType_MPI_Datatype sendtype, 
  void* recvbuf, int recvcount, uqRawType_MPI_Datatype recvtype, 
  int root,
  const char* whereMsg, const char* whatMsg) const
{
  //int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
  //                void *recvbuf, int recvcount, MPI_Datatype recvtype, 
  //                int root, MPI_Comm comm )
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Gather(sendbuf, sendcnt, sendtype,
                         recvbuf, recvcount, recvtype,
                         root, m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#else
  size_t sendDataTypeSize = sizeOfDataType(sendtype, whereMsg, whatMsg);
  size_t recvDataTypeSize = sizeOfDataType(recvtype, whereMsg, whatMsg);
  size_t sendTotal = sendDataTypeSize*sendcnt;
  size_t recvTotal = recvDataTypeSize*recvcount;
  if (sendTotal != recvTotal) {
    std::cerr << "uqMpiCommClass::Gather()"
              << ": sendTotal != recvTotal"
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO(sendTotal != recvTotal,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
  memcpy(recvbuf, sendbuf, sendTotal);
#endif
  return;
}
 
void
uqMpiCommClass::Gatherv(
  void* sendbuf, int sendcnt, uqRawType_MPI_Datatype sendtype, 
  void* recvbuf, int* recvcnts, int* displs, uqRawType_MPI_Datatype recvtype, 
  int root,
  const char* whereMsg, const char* whatMsg) const
{
  //int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
  //                void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, 
  //                int root, MPI_Comm comm )
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Gatherv(sendbuf, sendcnt, sendtype,
                          recvbuf, recvcnts, displs, recvtype,
                          root, m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#else
  size_t sendDataTypeSize = sizeOfDataType(sendtype, whereMsg, whatMsg);
  size_t recvDataTypeSize = sizeOfDataType(recvtype, whereMsg, whatMsg);
  size_t sendTotal = sendDataTypeSize*sendcnt;
  size_t recvTotal = recvDataTypeSize*recvcnts[0];
  if (sendTotal != recvTotal) {
    std::cerr << "uqMpiCommClass::Gatherv()"
              << ": sendTotal != recvTotal"
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO(sendTotal != recvTotal,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
  memcpy(recvbuf, sendbuf, sendTotal);
#endif
  return;
}

void
uqMpiCommClass::Recv(
  void* buf, int count, uqRawType_MPI_Datatype datatype, int source, int tag, uqRawType_MPI_Status* status,
  const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Recv(buf, count, datatype, source, tag, m_rawComm, status);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#else
  std::cerr << "uqMpiCommClass::Recv()"
            << ": should note be used if there is no 'mpi'"
            << std::endl;
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#endif
  return;
}
 
void
uqMpiCommClass::Send(
  void* buf, int count, uqRawType_MPI_Datatype datatype, int dest, int tag,
  const char* whereMsg, const char* whatMsg) const
{
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Send(buf, count, datatype, dest, tag, m_rawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#else
  std::cerr << "uqMpiCommClass::Send()"
            << ": should note be used if there is no 'mpi'"
            << std::endl;
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      whereMsg,
                      whatMsg);
#endif
  return;
}

void
uqMpiCommClass::syncPrintDebugMsg(const char* msg, unsigned int msgVerbosity, unsigned int numUSecs) const
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

#ifdef QUESO_HAS_MPI
#else
size_t 
uqMpiCommClass::sizeOfDataType(uqRawType_MPI_Datatype datatype, const char* whereMsg, const char* whatMsg) const
{
  size_t dataTypeSize = 0;

  switch (datatype) {
    case uqRawValue_MPI_CHAR:
      dataTypeSize = sizeof(char);
    break;

    case uqRawValue_MPI_INT:
      dataTypeSize = sizeof(int);
    break;

    case uqRawValue_MPI_DOUBLE:
      dataTypeSize = sizeof(double);
    break;

    case uqRawValue_MPI_UNSIGNED:
      dataTypeSize = sizeof(unsigned int);
    break;

    default: 
      std::cerr << "uqMpiCommClass::Allreduce()"
                << ": datatype not supported yet"
                << std::endl;
      UQ_FATAL_TEST_MACRO(true,
                          m_worldRank,
                          whereMsg,
                          whatMsg);
    break;
  }

  return dataTypeSize;
}
#endif
