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

#include <iostream>
#include <fstream>
#include <vector>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshFunction.h>
#include <queso/LibMeshOperatorBase.h>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/eigen_system.h>
#include <libmesh/condensed_eigen_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/parallel.h>
#include <libmesh/parallel_implementation.h>

namespace QUESO {

LibMeshOperatorBase::LibMeshOperatorBase(
    const FunctionOperatorBuilder & builder, libMesh::MeshBase & m)
  : OperatorBase(),
    equation_systems(new libMesh::EquationSystems(m)),
    builder(builder)
{
#ifndef LIBMESH_HAVE_SLEPC
  if (m.processor_id() == 0)
    std::cerr << "ERROR: This example requires libMesh to be\n"
              << "compiled with SLEPc eigen solvers support!"
              << std::endl;
#else

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_assert(false, "--disable-singleprecision");
#endif

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  // SLEPc currently gives us an "inner product not well defined" with
  // Number==complex
  libmesh_example_assert(false, "--disable-complex");
#endif

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  this->equation_systems->add_system<libMesh::CondensedEigenSystem>("Eigensystem");

#endif // LIBMESH_HAVE_SLEPC
}

LibMeshOperatorBase::~LibMeshOperatorBase()
{
}

void LibMeshOperatorBase::save_converged_evals(const std::string & filename) const
{
  unsigned int i;
  std::ofstream evals_file(filename.c_str());

  std::pair<libMesh::Real, libMesh::Real> eval;
  for (i = 0; i < this->nconv; i++) {
    eval = this->equation_systems
               ->get_system<libMesh::EigenSystem>("Eigensystem").get_eigenpair(i);
    if (equation_systems->processor_id() == 0) {
      evals_file << eval.first << " " << eval.second << std::endl;
    }
  }
  if (this->equation_systems->processor_id() == 0) {
    evals_file.close();
  }
}

void LibMeshOperatorBase::save_converged_evec(const std::string & filename,
    unsigned int i) const
{
  if (i < this->nconv) {
    SharedPtr<libMesh::EquationSystems>::Type es(this->equation_systems);
    es->get_system<libMesh::EigenSystem>("Eigensystem").get_eigenpair(i);
    libMesh::ExodusII_IO(es->get_mesh()).write_equation_systems(filename, *es);
  }
  else {
    std::cerr << "Warning: eigenpair " << i
              << " did not converge. Not saving."
              << std::endl;
  }
}

unsigned int LibMeshOperatorBase::get_num_converged() const {
  return this->nconv;
}

double LibMeshOperatorBase::get_eigenvalue(unsigned int i) const
{
  if (i < this->nconv) {
    std::pair<libMesh::Real, libMesh::Real> eval;
    SharedPtr<libMesh::EquationSystems>::Type es(this->equation_systems);
    eval = es->get_system<libMesh::EigenSystem>("Eigensystem").get_eigenpair(i);
    return eval.first;
  }
  else {
    return -1;
  }
}

double LibMeshOperatorBase::get_inverted_eigenvalue(unsigned int i) const
{
  return 1.0 / this->get_eigenvalue(i);
}

libMesh::EquationSystems & LibMeshOperatorBase::get_equation_systems() const
{
  return *this->equation_systems;
}

SharedPtr<FunctionBase>::Type
LibMeshOperatorBase::inverse_kl_transform(std::vector<double> & xi,
    double alpha) const
{
  unsigned int i;
  SharedPtr<libMesh::EquationSystems>::Type es(this->equation_systems);
  LibMeshFunction *kl = new LibMeshFunction(this->builder, es->get_mesh());

  // Make sure all procs in libmesh mpi communicator all have the same xi.  No,
  // I can't set the seed in QUESO.  That would mess with the QUESO
  // communicator.
  this->equation_systems->comm().broadcast(xi);

  SharedPtr<libMesh::EquationSystems>::Type kl_eq_sys(kl->get_equation_systems());

  std::pair<libMesh::Real, libMesh::Real> eval;
  for (i = 0; i < this->get_num_converged(); i++) {
    eval = es->get_system<libMesh::EigenSystem>("Eigensystem").get_eigenpair(i);
    kl_eq_sys->get_system<libMesh::ExplicitSystem>("Function").solution->add(
        xi[i] / pow(eval.first, alpha / 2.0),
        *es->get_system<libMesh::EigenSystem>("Eigensystem").solution);
  }

  SharedPtr<FunctionBase>::Type ap(kl);
  return ap;
}

}  // End namespace QUESO

#endif  // QUESO_HAVE_LIBMESH_SLEPC
