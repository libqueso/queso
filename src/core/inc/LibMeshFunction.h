//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef QUESO_LIBMESHFUNCTION_H
#define QUESO_LIBMESHFUNCTION_H

#include <string>
#include <boost/shared_ptr.hpp>
#include <queso/FunctionBase.h>

namespace libMesh {
  class MeshBase;
  class EquationSystems;
}

namespace QUESO {

class FunctionOperatorBuilder;

/*!
 * \file LibMeshFunction.h
 * \brief Function objects using libMesh for the backend
 *
 * \class LibMeshFunction
 * \brief Function objects using libMesh for the backend
 */

class LibMeshFunction : public FunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Construct a function with a user-provided mesh \c m and builder
  /*!
   * It is expected the lifetime of \c m will outlive \c this
   *
   * A \c builder object is just one that a FEM library backend can use to set
   * up various options. Polynomial type, polynomial order, and the number of
   * eigenpairs to request are good examples.
   */
  LibMeshFunction(const FunctionOperatorBuilder & builder,
      libMesh::MeshBase & m);

  //! Destructor
  ~LibMeshFunction();
  //@}

  //! Will print mesh-related libMesh foo to \c std::cerr
  void print_info() const;

  //! Save the current function to an Exodus file called \c filename
  virtual void save_function(const std::string & filename) const;

  //! Execute \c this += \c scale * \c rhs
  virtual void add(double scale, const FunctionBase & rhs);

  //! Pointwise multiply \c f1 and \c f2 and store the result in \c *this
  virtual void pointwise_mult(const FunctionBase & f1,
      const FunctionBase & f2);

  //! Execute \c this *= \c scale
  virtual void scale(double scale);

  //! Set \c this to zero everywhere
  virtual void zero();

  //! Return the L2-norm of \c this
  virtual double L2_norm() const;

  //! Create a zero function copy of \c this and return pointer to it
  /*!
   * Create a new instance of FunctionBase representing the function that is
   * identically zero (by copying \c this) everywhere and return a boost shared
   * pointer to it
   */
  virtual boost::shared_ptr<FunctionBase> zero_clone() const;

  //! This is public for now, but it should be encapsulated. Don't touch it.
  boost::shared_ptr<libMesh::EquationSystems> equation_systems;
};

}  // End namespace QUESO

#endif // QUESO_LIBMESHFUNCTION_H
