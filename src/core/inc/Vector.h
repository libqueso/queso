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

#ifndef UQ_VECTOR_H
#define UQ_VECTOR_H

#include <queso/Environment.h>
#include <queso/Map.h>
#include <iostream>
#include <queso/Defines.h>

namespace QUESO {

/*! \file Vector.h
    \brief Vector class.
*/

/*! \class Vector
    \brief Class for vector operations (virtual).

    Base vector class. The vector class is an abstract class designed to be used as a base class
    for different vector implementations (all actual vector classes in QUESO).
*/


class Vector
{
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Shaped Constructor
  Vector(const BaseEnvironment& env, const Map& map);

  //! Virtual Destructor
  virtual ~Vector();
  //@}

  //! @name Environment and Map methods
  //@{
    const BaseEnvironment& env                 ()           const;
    const Map&             map                 ()           const;
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

  //! Copies base data from vector \c src to \c this vector
  virtual void                    base_copy           (const Vector& src);

  //! Environment variable.
  const BaseEnvironment& m_env;

#ifdef QUESO_CLASSES_INSTANTIATE_NEW_MAPS

   //! Mapping variable.
  const Map              m_map;
#else

  //! Mapping variable.
  const Map&             m_map;
#endif

  //! Flag for either or not print this matrix horizontally.
  mutable bool                  m_printHorizontally;

  //! Flag for either or not print this matrix in scientific notation.
  mutable bool                  m_printScientific;

private:
  //! Default Constructor
  Vector();
};

}  // End namespace QUESO

#endif // UQ_VECTOR_H
