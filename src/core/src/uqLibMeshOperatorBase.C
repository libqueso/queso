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

#include <iostream>
#include <uqFunctionOperatorBuilder.h>
#include <uqLibMeshFunction.h>
#include <uqLibMeshOperatorBase.h>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/condensed_eigen_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/parallel.h>
#include <libmesh/parallel_implementation.h>

using namespace std;
using namespace libMesh;

uqLibMeshOperatorBase::uqLibMeshOperatorBase(
    const uqFunctionOperatorBuilder & builder, MeshBase & m)
  : uqOperatorBase(builder)
{
#ifndef LIBMESH_HAVE_SLEPC
  if (processor_id() == 0)
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

  // Create an equation systems object.
  this->equation_systems = new EquationSystems(m);

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  this->equation_systems->add_system<CondensedEigenSystem>("Eigensystem");

#endif // LIBMESH_HAVE_SLEPC
}

uqLibMeshOperatorBase::~uqLibMeshOperatorBase()
{
  delete this->equation_systems;
}

void uqLibMeshOperatorBase::save_converged_evals(const string & filename) const
{
  unsigned int i;
  ofstream evals_file(filename.c_str());

  pair<Real, Real> eval;
  for (i = 0; i < this->nconv; i++) {
    eval = this->equation_systems
               ->get_system<EigenSystem>("Eigensystem").get_eigenpair(i);
    if (processor_id() == 0) {
      evals_file << eval.first << " " << eval.second << endl;
    }
  }
  if (processor_id() == 0) {
    evals_file.close();
  }
}

void uqLibMeshOperatorBase::save_converged_evec(const string & filename,
    unsigned int i) const
{
  if (i < this->nconv) {
    EquationSystems * es = this->equation_systems;
    es->get_system<EigenSystem>("Eigensystem").get_eigenpair(i);
    ExodusII_IO(es->get_mesh()).write_equation_systems(filename, *es);
  }
  else {
    std::cerr << "Warning: eigenpair " << i
              << " did not converge. Not saving."
              << std::endl;
  }
}

unsigned int uqLibMeshOperatorBase::get_num_converged() const {
  return this->nconv;
}

double uqLibMeshOperatorBase::get_eigenvalue(unsigned int i) const
{
  if (i < this->nconv) {
    pair<Real, Real> eval;
    EquationSystems * es = this->equation_systems;
    eval = es->get_system<EigenSystem>("Eigensystem").get_eigenpair(i);
    return eval.first;
  }
  else {
    return -1;
  }
}

double uqLibMeshOperatorBase::get_inverted_eigenvalue(unsigned int i) const
{
  return 1.0 / this->get_eigenvalue(i);
}

EquationSystems & uqLibMeshOperatorBase::get_equation_systems() const
{
  return *this->equation_systems;
}

boost::shared_ptr<uqFunctionBase>
uqLibMeshOperatorBase::inverse_kl_transform(vector<double> & xi,
    double alpha) const
{
  unsigned int i;
  EquationSystems * es = this->equation_systems;
  uqLibMeshFunction *kl = new uqLibMeshFunction(this->builder, es->get_mesh());

  // Make sure all procs in libmesh mpi communicator all have the same xi.  No,
  // I can't set the seed in QUESO.  That would mess with the QUESO
  // communicator.
  CommWorld.broadcast(xi);

  EquationSystems *kl_eq_sys = kl->equation_systems;

  pair<Real, Real> eval;
  for (i = 0; i < this->get_num_converged(); i++) {
    eval = es->get_system<EigenSystem>("Eigensystem").get_eigenpair(i);
    kl_eq_sys->get_system<ExplicitSystem>("Function").solution->add(
        xi[i] / pow(eval.first, alpha / 2.0),
        *es->get_system<EigenSystem>("Eigensystem").solution);
  }

  boost::shared_ptr<uqFunctionBase> ap(kl);
  return ap;
}
