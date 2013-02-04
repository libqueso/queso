//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#include <iostream>
#include <uqLibMeshOperatorBase.h>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/eigen_system.h>

uqLibMeshOperatorBase::uqLibMeshOperatorBase()
  : uqOperatorBase()
{
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

  // Finally, read in the number of eigenpairs we want to compute!
  // Refactor this out
  unsigned int n_evals = 10;

  this->mesh = new libMesh::Mesh;
  this->mesh->read(filename);

  // Create an equation systems object.
  this->equation_systems = new libMesh::EquationSystems(*this->mesh);

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  libMesh::EigenSystem & eigen_system =
    this->equation_systems->add_system<libMesh::EigenSystem> ("Eigensystem");

  // Declare the system variables.
  // Adds the variable "u" to "Eigensystem".   "u"
  // will be approximated using second-order approximation.
  eigen_system.add_variable("u", SECOND);

#endif // LIBMESH_HAVE_SLEPC
}

uqLibMeshOperatorBase::~uqLibMeshOperatorBase()
{
  delete this->mesh;
  delete this->equation_systems;
}
