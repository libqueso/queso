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

#ifndef QUESO_LIBMESHOPERATOR_BASE_H
#define QUESO_LIBMESHOPERATOR_BASE_H

#include <string>
#include <set>
#include <vector>
#include <queso/SharedPtr.h>
#include <queso/FunctionBase.h>
#include <queso/OperatorBase.h>
#include <libmesh/system.h>

namespace libMesh {
  class Mesh;
  class EquationSystems;
}

namespace QUESO {

class FunctionOperatorBuilder;

/*!
 * \file LibMeshOperatorBase.h
 * \brief Abstract base class for operator objects using libmesh in the backend
 *
 * \class LibMeshOperatorBase
 * \brief Abstract base class for operator objects using libmesh in the backend
 */

class LibMeshOperatorBase : public OperatorBase,
                              public libMesh::System::Assembly {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constuct an operator on the mesh \c m using a builder \c builder
  /*!
   * A FunctionOperatorBuilder object is just one that a FEM library backend can use to set
   * up various options. Polynomial type, polynomial order, and the number of
   * eigenpairs to request are good examples.
   */
  LibMeshOperatorBase(const FunctionOperatorBuilder & builder,
      libMesh::MeshBase & m);

  //! Destructor
  ~LibMeshOperatorBase();
  //@}

  //! Must implement this for the solve to work
  /*!
   * This gets called by libMesh to do the assembly for the eigenvalue problem
   */
  virtual void assemble() = 0;

  //! Print libmesh related information
  virtual void print_info() const = 0;

  //! Save the eigenvalues to file \c filename
  virtual void save_converged_evals(const std::string &filename) const;

  //! Save converged eigenfunction \c i to \c filename
  virtual void
  save_converged_evec(const std::string &filename, unsigned int i) const;

  //! Return the number of converged eigenpairs
  virtual unsigned int get_num_converged() const;

  //! Return eigenvalue \c i.
  /*!
   * You can store them however you want, but having some kind of order to them
   * is useful for \c InfiniteDimensionalMeasure
   */
  virtual double get_eigenvalue(unsigned int i) const;

  //! Return the reciprocal of eigenvalue \c i.
  virtual double get_inverted_eigenvalue(unsigned int i) const;

  //! Return the internal libmesh equation systems object
  virtual libMesh::EquationSystems & get_equation_systems() const;

  //! Given coefficients \c xi, computes the Karhunen-Loeve transform
  /*!
   *  This transform goes from coefficient space to physical space using
   *  \c this as the precision operator:
   *  \sum_k \xi_k / pow(\lambda_k, \c alpha / 2.0) \phi_k(x)
   *  where the lambda are eigenvalues of the precision operator, \c this, and
   *  the \phi(x) are eigenfunctions of the precision operator, \c this
   */
  virtual SharedPtr<FunctionBase>::Type
  inverse_kl_transform(std::vector<double> & xi, double alpha) const;

protected:
  SharedPtr<libMesh::EquationSystems>::Type equation_systems;

  const FunctionOperatorBuilder & builder;

  //! The number of requested eigenpairs
  unsigned int num_req_pairs;

  //! The number of converged eigenpairs
  unsigned int nconv;
};

}  // End namespace QUESO

#endif // QUESO_LIBMESHOPERATOR_BASE_H

#endif  // QUESO_HAVE_LIBMESH_SLEPC
