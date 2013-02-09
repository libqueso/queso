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
  //! Default constructor
  uqLibMeshOperatorBase();

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
  unsigned int get_num_converged() const;

protected:
  libMesh::Mesh *mesh;
  libMesh::EquationSystems *equation_systems;

  //! The number of converged eigenvalue/eigenvector pairs
  unsigned int nconv;
};

#endif // __QUESO_LIBMESHOPERATOR_BASE__
