//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#include <iostream>
#include <uqLibMeshOperatorBase.h>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/condensed_eigen_system.h>
#include <libmesh/exodusII_io.h>

uqLibMeshOperatorBase::uqLibMeshOperatorBase()
  : uqOperatorBase()
{
#ifndef LIBMESH_HAVE_SLEPC
  if (libMesh::processor_id() == 0)
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

  this->mesh = new libMesh::Mesh;
  libMesh::MeshTools::Generation::build_square(*this->mesh, 20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);

  // Create an equation systems object.
  this->equation_systems = new libMesh::EquationSystems(*this->mesh);

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  this->equation_systems->add_system<libMesh::CondensedEigenSystem>("Eigensystem");

#endif // LIBMESH_HAVE_SLEPC
}

uqLibMeshOperatorBase::uqLibMeshOperatorBase(const std::string& filename)
  : uqOperatorBase(filename)
{
#ifndef LIBMESH_HAVE_SLEPC
  if (libMesh::processor_id() == 0)
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

  this->mesh = new libMesh::Mesh;
  this->mesh->read(filename);

  // Create an equation systems object.
  this->equation_systems = new libMesh::EquationSystems(*this->mesh);

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  this->equation_systems->add_system<libMesh::EigenSystem> ("Eigensystem");

#endif // LIBMESH_HAVE_SLEPC
}

uqLibMeshOperatorBase::~uqLibMeshOperatorBase()
{
  delete this->mesh;
  delete this->equation_systems;
}

void uqLibMeshOperatorBase::save_converged_evals(const std::string &filename) const
{
  unsigned int i;
  std::ofstream evals_file(filename.c_str());

  for (i = 0; i < this->nconv; i++) {
    std::pair<libMesh::Real, libMesh::Real> eval =
      this->equation_systems->get_system<libMesh::EigenSystem>("Eigensystem").get_eigenpair(i);
    evals_file << eval.first << " " << eval.second << std::endl;
  }
  evals_file.close();
}

void uqLibMeshOperatorBase::save_converged_evec(const std::string &filename, unsigned int i) const
{
  if (i < this->nconv) {
    this->equation_systems->get_system<libMesh::EigenSystem>("Eigensystem").get_eigenpair(i);
    libMesh::ExodusII_IO(*this->mesh).write_equation_systems(filename, *this->equation_systems);
  }
  else {
    std::cerr << "Warning: eigenpair" << i
              << "did not converge. Not saving."
              << std::endl;
  }
}

unsigned int uqLibMeshOperatorBase::get_num_converged() const {
  return this->nconv;
}
