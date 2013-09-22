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

#ifndef __QUESO_LIBMESHFUNCTION__
#define __QUESO_LIBMESHFUNCTION__

#include <string>
#include <boost/shared_ptr.hpp>
#include <queso/uqFunctionBase.h>

namespace libMesh {
  class MeshBase;
  class EquationSystems;
}

namespace QUESO {

class uqFunctionOperatorBuilder;

/*!
 * \file uqLibMeshFunction.h
 * \brief Function objects using libMesh for the backend
 */

class uqLibMeshFunction : public uqFunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Construct a function with a user-provided mesh \c m and builder
  /*!
   * It is expected the lifetime of \c m will outlive \c this
   */
  uqLibMeshFunction(const uqFunctionOperatorBuilder & builder,
      libMesh::MeshBase & m);

  //! Destructor
  ~uqLibMeshFunction();
  //@}

  //! Will print mesh-related libMesh foo to \c std::cerr
  void print_info() const;

  //! Save the current function to an Exodus file called \c filename
  virtual void save_function(const std::string & filename) const;

  virtual void add(double scale, const uqFunctionBase & rhs);

  virtual void pointwise_mult(const uqFunctionBase & f1,
      const uqFunctionBase & f2);

  virtual void scale(double scale);

  virtual void zero();

  virtual double L2_norm() const;

  virtual boost::shared_ptr<uqFunctionBase> zero_clone() const;

  //! This is public for now, but it should be encapsulated. Don't touch it.
  boost::shared_ptr<libMesh::EquationSystems> equation_systems;
};

}  // End namespace QUESO

#endif // __QUESO_LIBMESHFUNCTION__
