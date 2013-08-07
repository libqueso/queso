//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef __UQ_MAP_H__
#define __UQ_MAP_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_TRILINOS
#include <Epetra_Map.h>
#endif

#include <uqMpiComm.h>


/*! \file uqMap.h
    \brief A class for partitioning vectors and matrices.
*/

/*! \class uqMapClass
    \brief A class for partitioning vectors and matrices.
    
    It is often the case that multiple matrix and vector objects have an identical distribution 
    of elements on a parallel machine. The uqMapClass keeps information that describes this 
    distribution for matrices and vectors. Inspired by Trilinos Epetra_Map class.
*/

class uqMapClass
{
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Default constructor. Do not call this directly.
  uqMapClass();
  
  //! Constructor for a uniform linear distribution of elements. 
  uqMapClass(int                   numGlobalElements,
             int                   indexBase,
             const uqMpiCommClass& comm);
  
  //! Copy constructor.
  uqMapClass(const uqMapClass& src);

  //! Destructor
 ~uqMapClass();
 //@}

 //! @name Set methods
  //@{
  //! Assignment operator. 
  uqMapClass& operator= (const uqMapClass& rhs);
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
  //! Access function for uqMpiComm communicator. 
  const uqMpiCommClass& Comm             () const;
  
#ifdef QUESO_HAS_TRILINOS
  //! Trilinos Epetra_Map: A class for partitioning vectors and matrices.
  const Epetra_Map&     epetraMap        () const;
#endif

  //@}
private:
  //! Copies the map.
  void                  copy             (const uqMapClass& src);

  //! This communicator can be queried for processor rank and size information. 
  uqMpiCommClass m_uqMpiComm;
  
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

#endif // __UQ_MAP_H__
