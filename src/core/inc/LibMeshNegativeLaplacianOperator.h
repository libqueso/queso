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

#include <queso/Defines.h>

#ifdef QUESO_HAVE_LIBMESH_SLEPC

#ifndef QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H
#define QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H

#include <string>
#include <queso/LibMeshOperatorBase.h>

namespace libMesh {
  class MeshBase;
  class EquationSystems;
}


namespace QUESO {

class FunctionOperatorBuilder;

/*!
 * \file LibMeshNegativeLaplacianOperator.h
 * \brief Class describing negative Laplacian operator using libmesh backend
 *
 * \class LibMeshNegativeLaplacianOperator
 * \brief Class describing negative Laplacian operator using libmesh backend
 */

class LibMeshNegativeLaplacianOperator : public LibMeshOperatorBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Construct the negative Laplacian operator on the libmesh mesh \c m with builder \c builder
  /*!
   * The negative Laplacian operator is ( - \Delta )
   *
   * A \c builder object is just one that a FEM library backend can use to set
   * up various options. Polynomial type, polynomial order, and the number of
   * eigenpairs to request are good examples.
   */
  LibMeshNegativeLaplacianOperator(const FunctionOperatorBuilder & builder,
      libMesh::MeshBase & m);

  //! Destructor
  ~LibMeshNegativeLaplacianOperator();

  //! Method to assemble the mass and stiffness matrices
  virtual void assemble();

  //! Print libmesh related foo
  virtual void print_info() const;
};

}  // End namespace QUESO

#endif // QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H

#endif  // QUESO_HAVE_LIBMESH_SLEPC
