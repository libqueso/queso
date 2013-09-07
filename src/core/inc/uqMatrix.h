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

#ifndef __UQ_MATRIX_H__
#define __UQ_MATRIX_H__

#include <uqEnvironment.h>
#include <uqVector.h>
#include <iostream>

namespace QUESO {

/*! \file uqMatrix.h
    \brief Matrix class.
*/

/*! \class MatrixClass
    \brief Class for matrix operations (virtual). 
    
    Base matrix class. The matrix class is an abstract class designed to be used as a base class 
    for different matrix implementations (all actual matrix classes in QUESO).
*/


class MatrixClass
{
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  MatrixClass();
	   
  //! Shaped constructor.
  MatrixClass(const BaseEnvironmentClass& env, const MapClass& map);
	   
  //! Copy constructor.
  MatrixClass(const MatrixClass& rhs);
	   
  //! Virtual Destructor
  virtual ~MatrixClass();
  //@}

  //! @name Set methods
  //@{    
  //! Operator for copying a matrix.
  MatrixClass& operator= (const MatrixClass& rhs);
  
  //! Operator for multiplication of the matrix by a scalar.
  MatrixClass& operator*=(double a);
  
  //! Operator for addition (element-wise) of two matrices.
  MatrixClass& operator+=(const MatrixClass& rhs);
  
  //! Operator for subtraction (element-wise) of two matrices.
  MatrixClass& operator-=(const MatrixClass& rhs);
  //@}

  //! @name Environment and Map methods
  //@{ 
  const BaseEnvironmentClass& env                 ()           const;
  const MapClass&             map                 ()           const;
        unsigned int            numOfProcsForStorage()           const;
  //@}       
  
  //! @name Attribute methods
  //@{ 
  //! Returns the number of rows local of the matrix. 	  
  virtual unsigned int            numRowsLocal        () const = 0;
  
  //! Returns the number of rows global of the matrix. 
  virtual unsigned int            numRowsGlobal       () const = 0;
  
  //! Returns the number of columns local of the matrix. 	  
  virtual unsigned int            numCols             () const = 0;
  //@}
    
  //! @name Mathematical  methods
  //@{ 
  //! Performs Cholesky factorization of the matrix.
  virtual int                     chol                () = 0;
  //@}
  
  
  //! @name Get/Set methods
  //@{ 
  //! Sets to zero all the elements bellow (including or not) the diagonal of the matrix.
  virtual void                    zeroLower           (bool includeDiagonal = false) = 0;
  
  //! Sets to zero all the elements above (including or not) the diagonal of the matrix.
  virtual void                    zeroUpper           (bool includeDiagonal = false) = 0;
  //@}
  
  //! @name I/O and Miscellaneous methods
  //@{ 
  //! Print this matrix.
  virtual void                    print               (std::ostream& os) const = 0;

  //! Determines whether the matrix should be printed horizontally.  
          void                    setPrintHorizontally(bool value) const; // Yes, 'const'
          
  //! Checks if matrix will be is printed horizontally.          
          bool                    getPrintHorizontally()           const;
          
  //! Determines whether QUESO will run through this class in debug mode.        
          void                    setInDebugMode      (bool value) const; // Yes, 'const'
          
  //! Checks if QUESO will run through this class in debug mode.        
          bool                    getInDebugMode      ()           const;
  //@}
  
protected:
  //! Copies matrix \c src to \c this matrix.
  virtual void                    copy                (const MatrixClass& src);

  //! QUESO environment variable.
  const   BaseEnvironmentClass& m_env;
  
#ifdef QUESO_CLASSES_INSTANTIATE_NEW_MAPS

  //! Mapping variable.
  const   MapClass              m_map;

#else

  //! Mapping variable.
  const   MapClass&             m_map;

#endif

  //! Flag for either or not print this matrix.
  mutable bool                    m_printHorizontally;
  
  //! Flag for either or not QUESO is in debug mode.
  mutable bool                    m_inDebugMode;
};

}  // End namespace QUESO

#endif // __UQ_MATRIX_H__
