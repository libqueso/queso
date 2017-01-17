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

#ifndef QUESO_FUNCTION_BASE_H
#define QUESO_FUNCTION_BASE_H

#include <string>

#include <queso/SharedPtr.h>

namespace QUESO {

class FunctionOperatorBuilder;

/*!
 * \file FunctionBase.h
 * \brief Abstract base class for function objects
 * \class FunctionBase
 * \brief Abstract base class for function objects
 *
 * One needs to sublclass this abstract class to implement their own functions
 * using a backend not supported by QUESO.  If QUESO is linked against libMesh
 * (and libMesh was compiled with SLEPc), you may use LibMeshFunction.
 */

class FunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  FunctionBase();

  //! Destructor
  virtual ~FunctionBase();
  //@}

  //! Execute \c this += \c scale * \c rhs
  virtual void add(double scale, const FunctionBase & rhs) = 0;

  //! Pointwise multiply \c f1 and \c f2 and store the result in \c *this
  virtual void pointwise_mult(const FunctionBase & f1,
      const FunctionBase & f2) = 0;

  //! Execute \c this *= \c scale
  virtual void scale(double scale) = 0;

  //! Set \c this to zero everywhere
  virtual void zero() = 0;

  //! Return the L2-norm of \c this
  virtual double L2_norm() const = 0;

  //! Create a zero function copy of \c this and return pointer to it
  /*!
   * Create a new instance of FunctionBase representing the function that is
   * identically zero (by copying \c this) everywhere and return a boost shared
   * pointer to it
   */
  virtual SharedPtr<FunctionBase>::Type zero_clone() const = 0;

  //! Save the current function to an Exodus file called \c filename.  \c time is the time to attach to the function and is usually the iteration number
  virtual void save_function(const std::string & filename, double time) const = 0;
};

}  // End namespace QUESO

#endif // QUESO_FUNCTION_BASE_H
