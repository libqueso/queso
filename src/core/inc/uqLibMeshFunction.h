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
#include <uqFunctionBase.h>

namespace libMesh {
  class Mesh;
  class EquationSystems;
}

/*!
 * \file uqLibMeshFunction.h
 * \brief Function objects using lib mesh for the backend
 */

class uqLibMeshFunction : public uqFunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. Zero everywhere.
  uqLibMeshFunction();

  //! Create a function that is equal to \c value everywhere
  uqLibMeshFunction(double value);

  //! Destructor
  ~uqLibMeshFunction();

  //! Will print mesh-related libMesh foo to the screen
  void print_info();

  //@}

  virtual void save_function(const std::string & filename) const;

  libMesh::Mesh *mesh;
  libMesh::EquationSystems *equation_systems;
};

#endif // __QUESO_LIBMESHFUNCTION__
