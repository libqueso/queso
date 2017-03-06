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

#ifndef UQ_MPI_COMM_H
#define UQ_MPI_COMM_H

#include <queso/Defines.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

#ifdef QUESO_HAS_TRILINOS
class Epetra_Comm;
#endif

namespace QUESO {

#ifdef QUESO_HAS_MPI
typedef MPI_Comm     RawType_MPI_Comm ;
typedef MPI_Group    RawType_MPI_Group ;
typedef MPI_Datatype RawType_MPI_Datatype ;
typedef MPI_Datatype data_type ;
typedef MPI_Op       RawType_MPI_Op ;
typedef MPI_Status   RawType_MPI_Status ;
#define RawValue_MPI_COMM_SELF  MPI_COMM_SELF
#define RawValue_MPI_ANY_SOURCE MPI_ANY_SOURCE
#define RawValue_MPI_CHAR       MPI_CHAR
#define RawValue_MPI_INT        MPI_INT
#define RawValue_MPI_DOUBLE     MPI_DOUBLE
#define RawValue_MPI_UNSIGNED   MPI_UNSIGNED
#define RawValue_MPI_MIN        MPI_MIN
#define RawValue_MPI_MAX        MPI_MAX
#define RawValue_MPI_SUM        MPI_SUM
#else
typedef int RawType_MPI_Comm;
typedef int RawType_MPI_Group;
typedef int RawType_MPI_Datatype;
struct data_type { };
typedef int RawType_MPI_Op;
typedef int RawType_MPI_Status;
#define RawValue_MPI_COMM_SELF   0
#define RawValue_MPI_ANY_SOURCE -1
#define RawValue_MPI_CHAR        0
#define RawValue_MPI_INT         1
#define RawValue_MPI_DOUBLE      2
#define RawValue_MPI_UNSIGNED    3
#define RawValue_MPI_MIN         0
#define RawValue_MPI_MAX         1
#define RawValue_MPI_SUM         2
#endif

/**
 * Encapsulates the MPI_Datatype.  Taken from libmesh.
 */
class DataType
{
public:
  DataType () : _datatype() {}

  DataType (const DataType &other) :
    _datatype(other._datatype)
  {}

  DataType (const RawType_MPI_Datatype &type) :
    _datatype(type)
  {}

#ifdef QUESO_HAS_MPI
  DataType (const DataType &other, unsigned int count)
  {
    // FIXME - if we nest an inner type here will we run into bug
    // https://github.com/libMesh/libmesh/issues/631 again?
    MPI_Type_contiguous(count, other._datatype, &_datatype);
    this->commit();
  }
#else
  DataType (const DataType &, unsigned int)
  {
  }
#endif

  DataType & operator = (const DataType &other)
  { _datatype = other._datatype; return *this; }

  DataType & operator = (const RawType_MPI_Datatype &type)
  { _datatype = type; return *this; }

  operator const RawType_MPI_Datatype & () const
  { return _datatype; }

  operator RawType_MPI_Datatype & ()
  { return _datatype; }

  void commit ()
  {
#ifdef QUESO_HAS_MPI
    MPI_Type_commit (&_datatype);
#endif
  }

  void free ()
  {
#ifdef QUESO_HAS_MPI
    MPI_Type_free (&_datatype);
#endif
  }

protected:
  RawType_MPI_Datatype _datatype;
};

/**
 * Templated class to provide the appropriate MPI datatype
 * for use with built-in C types or simple C++ constructions.
 *
 * More complicated data types may need to provide a pointer-to-T so
 * that we can use MPI_Address without constructing a new T.
 */
template <typename T>
class StandardType : public DataType
{
#ifdef QUESO_HAS_CXX11  // This macro isn't defined (yet)
    // Get a slightly better compiler diagnostic if we have C++11
  static_assert(dependent_false<T>::value,
                "Only specializations of StandardType may be used, did you forget to include a header file (e.g. parallel_algebra.h)?");
#endif

 /*
  * The unspecialized class is useless, so we make its constructor
  * private to catch mistakes at compile-time rather than link-time.
  * Specializations should have a public constructor of the same
  * form.
  */
private:
  StandardType(const T* example = NULL);
};

#ifdef QUESO_HAS_MPI

#define QUESO_STANDARD_TYPE(cxxtype,mpitype)                                  \
  template<>                                                            \
  class StandardType<cxxtype> : public DataType                         \
  {                                                                     \
  public:                                                               \
    explicit                                                            \
      StandardType(const cxxtype* = NULL) : DataType(mpitype) {}        \
  }

#else

#define QUESO_STANDARD_TYPE(cxxtype,mpitype)                          \
  template<>                                                    \
  class StandardType<cxxtype> : public DataType                 \
  {                                                             \
  public:                                                       \
    explicit                                                    \
      StandardType(const cxxtype* = NULL) : DataType() {}       \
  }

#endif

QUESO_STANDARD_TYPE(char,MPI_CHAR);
QUESO_STANDARD_TYPE(int,MPI_INT);
QUESO_STANDARD_TYPE(unsigned int,MPI_UNSIGNED);
QUESO_STANDARD_TYPE(double,MPI_DOUBLE);

class BaseEnvironment;

/*! \file MpiComm.h
    \brief MPI Communicator Class.
*/

/*! \class MpiComm
    \brief The QUESO MPI Communicator Class.

    This class uses MPI (the Message Passing Interface) for distributed-memory
    communication between one or more parallel processes. It is meant to insulate
    the user from the specifics of communication that are not required for normal
    manipulation of linear algebra objects.
*/
class MpiComm
{
public:
   //! @name Constructor/Destructor methods
  //@{

  //! QUESO MpiComm MPI parallel constructor.
  /*!
   * This constructs an MpiComm that uses the given "raw" MPI communicator
   * underneath.  MPI_Init *must* have been called before instantiating an
   * object of this type.
   *
   * The MPI_Comm must be valid for the lifetime of this MpiComm.
   */
  MpiComm(const BaseEnvironment& env, RawType_MPI_Comm inputRawComm);

  //! QUESO MpiComm MPI serial constructor.
  /*!
   * This constructs an MpiComm that defaults to MPI_COMM_SELF
   * underneath.  MPI_Init need not be called before using this constructor.
   *
   * The MPI_Comm must be valid for the lifetime of this MpiComm.
   */
  MpiComm(const BaseEnvironment& env);

  //! Copy Constructor.
  /** Makes an exact copy of an existing MpiComm instance.*/
  MpiComm(const MpiComm& src);

  //! Destructor
 ~MpiComm();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator.
  MpiComm& operator= (const MpiComm& rhs);
  //@}


  //! @name Attribute Accessor Methods
  //@{
#ifdef QUESO_HAS_MPI
  //! Extract MPI Communicator from a MpiComm object.
  RawType_MPI_Comm Comm     () const;
#endif  // QUESO_HAS_MPI

  //! Return my process ID.
  int                MyPID    () const;

  //! Returns total number of processes.
  int                NumProc  () const;
  //@}

  //! @name Methods Overridden from Comm
  //@{
  //! Combines values from all processes and distributes the result back to all processes
  /*! \param sendbuf starting address of send buffer
   * \param count number of elements in send buffer
   * \param datatype data type of elements of send buffer
   * \param op operation
   * \param recvbuf (output) starting address of receive buffer
   *
   * This method is deprecated.  Use the templated Allreduce method instead.
   */
  void               Allreduce(void* sendbuf, void* recvbuf, int count, RawType_MPI_Datatype datatype,
			       RawType_MPI_Op op, const char* whereMsg, const char* whatMsg) const;

  //! Combines values from all processes and distributes the result back to all processes
  /*!
   * \param sendbuf starting address of send buffer containing elements of type T
   * \param count number of elements in send buffer
   * \param op operation
   * \param recvbuf (output) starting address of receive buffer containing elements of type T
   */
  template <typename T>
  void Allreduce(const T * sendbuf, T * recvbuf, int count, RawType_MPI_Op op,
                 const char* whereMsg, const char* whatMsg) const;

  //! Pause every process in *this communicator until all the processes reach this point.
  /*! Blocks the caller until all processes in the communicator have called it; that is,
   * the call returns at any process only after all members of the communicator have entered the call.*/
  void               Barrier  () const; // const char* whereMsg, const char* whatMsg) const;

  //! Broadcast values from the root process to the slave processes.
  /*! Broadcasts a message from the process with rank "root" to all other processes of the communicator.
   * \param buffer (input/output) starting address of buffer
   * \param count number of entries in buffer
   * \param datatype data type of buffer
   * \param root rank of broadcast root */
  void               Bcast    (void* buffer, int count, RawType_MPI_Datatype datatype, int root,
                               const char* whereMsg, const char* whatMsg) const;

  //! Gather values from each process to collect on all processes.
  /*!\param sendbuf starting address of send buffer
   * \param sendcnt number of elements in send buffer
   * \param sendtype data type of send buffer elements
   * \param recvcount number of elements for any single receive
   * \param recvtype data type of recv buffer elements
   * \param root rank of receiving process
   * \param recvbuf (output) address of receive buffer
   *
   * This method is deprecated.  Use the templated Gather method instead.
   */
  void               Gather   (void *sendbuf, int sendcnt, RawType_MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, RawType_MPI_Datatype recvtype,
                               int root,
                               const char* whereMsg, const char* whatMsg) const;

  //! Gather values from each process to collect on all processes.
  /*!
   * \param sendbuf starting address of send buffer containing elements of type T
   * \param sendcnt number of elements in send buffer
   * \param recvcount number of elements for any single receive
   * \param root rank of receiving process
   * \param recvbuf (output) address of receive buffer containing elements of type T
   */
  template <typename T>
  void Gather(const T * sendbuf, int sendcnt, T * recvbuf, int recvcount, int root,
              const char * whereMsg, const char * whatMsg) const;

 //! Gathers into specified locations from all processes in a group
 /*! \param sendbuf starting address of send buffer
  * \param sendcount number of elements in send buffer
  * \param sendtype data type of send buffer elements
  * \param recvcounts integer array (of length group size) containing the number of elements
  * that are received from each process
  * \param displs integer array (of length group size). Entry i specifies the displacement
  * relative to recvbuf at which to place the incoming data from process i
  * \param recvtype data type of recv buffer elements
  * \param root rank of receiving process*/
  void               Gatherv  (void *sendbuf, int sendcnt, RawType_MPI_Datatype sendtype,
                               void *recvbuf, int *recvcnts, int *displs, RawType_MPI_Datatype recvtype,
                               int root,
                               const char* whereMsg, const char* whatMsg) const;

  //! Gathers into specified locations from all processes in a group
  /*!
   * \param sendbuf starting address of send buffer containing elements of type T
   * \param sendcnt number of elements in send buffer
   * \param recvcnts integer array (of length group size) containing the number
   * of elements that are received from each process
   * \param displs integer array (of length group size). Entry i specifies the
   * displacement relative to recvbuf at which to place the incoming data from
   * process i
   * \param root rank of receiving process
   */
  template <typename T>
  void Gatherv(const T * sendbuf, int sendcnt, T * recvbuf, int * recvcnts,
               int * displs, int root, const char * whereMsg,
               const char * whatMsg) const;

  //! Blocking receive of data from this process to another process.
  /*!\param buf (output) initial address of receive buffer
   * \param status (output) status object
   * \param count maximum number of elements in receive buffer
   * \param datatype datatype of each receive buffer element
   * \param source rank of source
   * \param tag message tag */
  void               Recv     (void *buf, int count, RawType_MPI_Datatype datatype, int source, int tag, RawType_MPI_Status *status,
                               const char* whereMsg, const char* whatMsg) const;

  //! Possibly blocking send of data from this process to another process.
  /*!\param buf initial address of send buffer
   * \param count number of elements in send buffer
   * \param datatype datatype of each send buffer element
   * \param dest rank of destination
   * \param tag message tag*/
  void               Send     (void *buf, int count, RawType_MPI_Datatype datatype, int dest, int tag,
                               const char* whereMsg, const char* whatMsg) const;
 //@}

//! @name Miscellaneous Methods
  //@{
  //! Synchronizes all the processes and print debug message.
  void               syncPrintDebugMsg(const char* msg, unsigned int msgVerbosity, unsigned int numUSecs) const;

#ifdef QUESO_HAS_TRILINOS
  //! Extract MPI Communicator from a Epetra_MpiComm object.
  const Epetra_Comm& epetraMpiComm() const;
#endif

 //@}
private:
  //! Default Constructor
  /*! It should not be used by user.*/
  MpiComm();

  //! Copies from an existing MpiComm instance.
  void               copy          (const MpiComm& src);

  // QUESO environment
  const BaseEnvironment& m_env;

#ifdef QUESO_HAS_TRILINOS
  // Epetra communicator
  Epetra_Comm*               m_epetraComm;
#endif

  //! Embedded wrapped opaque MPI_Comm object.
  RawType_MPI_Comm            m_rawComm;

  //! World rank
  int                           m_worldRank;

  //! Process ID of this process
  int                           m_myPid;

  // Total number of processes
  int                           m_numProc;
};

}  // End namespace QUESO

#endif // UQ_MPI_COMM_H
