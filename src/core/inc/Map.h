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

#ifndef UQ_MAP_H
#define UQ_MAP_H

#include <queso/Defines.h>
#ifdef QUESO_HAS_TRILINOS
#include <Epetra_Map.h>
#endif

#include <queso/MpiComm.h>

namespace QUESO {

/*! \file Map.h
    \brief A class for partitioning vectors and matrices.
*/

/*! \class Map
    \brief A class for partitioning vectors and matrices.

    It is often the case that multiple matrix and vector objects have an identical distribution
    of elements on a parallel machine. The Map keeps information that describes this
    distribution for matrices and vectors. Inspired by Trilinos Epetra_Map class.
*/

class Map
{
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Constructor for a uniform linear distribution of elements.
  Map(int                   numGlobalElements,
             int                   indexBase,
             const MpiComm& comm);

  //! Copy constructor.
  Map(const Map& src);

  //! Destructor
 ~Map();
 //@}

 //! @name Set methods
  //@{
  //! Assignment operator.
  Map& operator= (const Map& rhs);
  //@}


  //! @name Size, dimension and local ID accessor methods
  //@{
  //! Returns the total number of elements across all processors.
  int                   NumGlobalElements() const;

  //! Returns the base integer value for indexed array references.
  //! The first position in the global processor in my processors.
  int                   IndexBase        () const;//1st position in the global processor in my processors

  //! Returns the number of elements owned by the calling processor.
  int                   NumMyElements    () const;

  //!The minimum global index value on the calling processor.
  int                   MinMyGID         () const;
  //@}

  //! @name Miscellaneous  methods
  //@{
  //! Access function for MpiComm communicator.
  const MpiComm& Comm             () const;

#ifdef QUESO_HAS_TRILINOS
  //! Trilinos Epetra_Map: A class for partitioning vectors and matrices.
  const Epetra_Map&     epetraMap        () const;
#endif

  //@}
private:
  //! Default constructor. Do not call this directly.
  Map();

  //! Copies the map.
  void                  copy             (const Map& src);

  //! This communicator can be queried for processor rank and size information.
  MpiComm m_MpiComm;

#ifdef QUESO_HAS_TRILINOS
  //! Epetra_Map
  Epetra_Map*    m_epetraMap;
#else
  //! Total number of elements across all processors.
  int            m_numGlobalElements;

  //!  Base integer value for indexed array references.
  int            m_indexBase;

  //! Number of elements owned by the calling processor.
  int            m_numMyElements;
#endif
};

}  // End namespace QUESO

#endif // UQ_MAP_H
