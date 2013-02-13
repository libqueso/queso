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

#ifndef __QUESO_LIBMESHOPERATOR_BASE__
#define __QUESO_LIBMESHOPERATOR_BASE__

/*!
 * \file uqLibMeshOperatorBase.h
 * \brief Abstract base class for operator objects using libmesh in the
 *        backend
 */

#include <string>
#include <set>
#include <vector>
#include <memory>
#include <uqFunctionBase.h>
#include <uqOperatorBase.h>
#include <libmesh/system.h>

namespace libMesh {
  class Mesh;
  class EquationSystems;
}

class uqLibMeshOperatorBase : public uqOperatorBase, public libMesh::System::Assembly {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constuct an operator on the mesh \c m
  uqLibMeshOperatorBase(libMesh::MeshBase & m);

  //! Construct an operator on the meshed domain in \c filename
  uqLibMeshOperatorBase(const std::string& filename);

  //! Destructor
  ~uqLibMeshOperatorBase();
  //@}

  //! Must implement this for the solve to work
  virtual void assemble() = 0;

  //! Print libmesh related information
  virtual void print_info() const = 0;

  //! Save the eigenvalues to file \c filename
  virtual void save_converged_evals(const std::string &filename) const;

  //! Save converged eigenfunction \c i to \c filename
  virtual void save_converged_evec(const std::string &filename, unsigned int i) const;

  //! Return the number of converged eigenpairs
  virtual unsigned int get_num_converged() const;

  //! Return eigenvalue \c i. You can store them however you want, but
  //! having some kind of order to them is useful for \c uqInfiniteDimensionalMeasure
  virtual double get_eigenvalue(unsigned int i) const;

  //! Return the reciprocal of eigenvalue \c i.
  virtual double get_inverted_eigenvalue(unsigned int i) const;

  //! Given coefficients \c xi, compute the inverse Karhunen-Loeve transform using \c num_terms terms.
  /*!
   *  This transform goes from coefficient space to physical space:
   *  \sum_k \xi_k \phi_k(x) / sqrt(\lambda_k)
   *  where the lambda are eigenvalues of \c this and the \phi(x) are
   *  eigenfunctions of \c this
   */
  virtual std::auto_ptr<uqFunctionBase> inverse_kl_transform(const std::vector<double> & xi) const;

protected:
  libMesh::EquationSystems * equation_systems;

  //! The number of converged eigenvalue/eigenvector pairs
  unsigned int nconv;
};

#endif // __QUESO_LIBMESHOPERATOR_BASE__
