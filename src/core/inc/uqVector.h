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

#ifndef __UQ_VECTOR_H__
#define __UQ_VECTOR_H__

#include <uqEnvironment.h>
#include <uqMap.h>
#include <iostream>
#include <uqDefines.h>
/*! \file uqVector.h
    \brief Vector class.
*/

/*! \class uqVectorClass
    \brief Class for vector operations (virtual). 
    
    Base vector class. The vector class is an abstract class designed to be used as a base class 
    for different vector implementations (all actual vector classes in QUESO). 
*/


class uqVectorClass
{
public:
     
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default Constructor
  uqVectorClass();
	   
  //! Shaped Constructor
  uqVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
          	   
  //! Copy Constructor 
  uqVectorClass(const uqVectorClass& rhs);
	   
  //! Virtual Destructor
  virtual ~uqVectorClass();
  //@}
  
  //! @name Set methods
  //@{ 
    
  //! Operator for copying a vector.
  uqVectorClass& operator=(const uqVectorClass& rhs);
  
  //! Operator for multiplication of the vector by a scalar.
  uqVectorClass& operator*=(double a);
  
  //! Operator for division of the vector by a scalar.
  uqVectorClass& operator/=(double a);
  
  //! Operator for addition (element-wise) of two vectors.
  uqVectorClass& operator+=(const uqVectorClass& rhs);
  
  //! Operator for subtraction (element-wise) of two vectors.
  uqVectorClass& operator-=(const uqVectorClass& rhs);
  //@}
  
  //! @name Environment and Map methods
  //@{ 
    const uqBaseEnvironmentClass& env                 ()           const;
    const uqMapClass&             map                 ()           const;
          unsigned int            numOfProcsForStorage()           const;
  //@}
  
  //! @name Attribute methods
  //@{ 
  //! Returns the local size of the vector.
  virtual unsigned int            sizeLocal           () const = 0;
  
  //! Returns the global size of the vector.
  virtual unsigned int            sizeGlobal          () const = 0;
  //@}
  
    //! @name Get/Set methods
  //@{ 
  //! Component-wise set the values of the vector to \c value.
  virtual void                    cwSet               (double value) = 0;
  
  //! Component-wise set the values of the vector to a random number from a Gaussian distribution.
  virtual void                    cwSetGaussian       (double mean, double stdDev) = 0;
  
  //! Component-wise invert the values of the vector.
  virtual void                    cwInvert            () = 0;
  //@}
  
  
  //! @name I/O and Miscellaneous methods
  //@{ 
  //! Sort the elements.
  virtual void                    sort                () = 0;
  
  //! Print this vector.
  virtual void                    print               (std::ostream& os) const = 0;
 
  //! Determines whether vector should be printed horizontally.
          void                    setPrintHorizontally(bool value) const; // Yes, 'const'
  
  //! Checks if vector is printed horizontally.
          bool                    getPrintHorizontally()           const;
	  
  //! Determines whether vector should be printed in Scientific Notation.
          void                    setPrintScientific  (bool value) const; // Yes, 'const'
          
  //! Checks if the vector should be printed in Scientific Notation.
          bool                    getPrintScientific  ()           const;
  
  //@}
  
protected:
  
  //! Copies vector \c src to \c this matrix.
  virtual void                    copy                (const uqVectorClass& src);

  //! Environment variable.
  const uqBaseEnvironmentClass& m_env;
  
#ifdef QUESO_CLASSES_INSTANTIATE_NEW_MAPS
  
   //! Mapping variable.
  const uqMapClass              m_map;
#else
  
  //! Mapping variable.
  const uqMapClass&             m_map;
#endif
  
  //! Flag for either or not print this matrix horizontally.
  mutable bool                  m_printHorizontally;
  
  //! Flag for either or not print this matrix in scientific notation.
  mutable bool                  m_printScientific;
};

#endif // __UQ_VECTOR_H__
