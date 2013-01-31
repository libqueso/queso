//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#include <iostream>
#include <uqFunctionBase.h>
#include <uqLibMeshFunction.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>

// Define the Finite Element object.
#include <libmesh/fe.h>

// Define Gauss quadrature rules.
#include <libmesh/quadrature_gauss.h>

// Define useful datatypes for finite element
// matrix and vector components.
#include <libmesh/sparse_matrix.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/elem.h>

// Define the DofMap, which handles degree of freedom
// indexing.
#include <libmesh/dof_map.h>

uqLibMeshFunction::uqLibMeshFunction()
  : uqFunctionBase()
{
  mesh = new libMesh::Mesh;
  
  // Use the MeshTools::Generation mesh generator to create a uniform 1D grid
  // on the line [0,1].  We instruct the mesh generator to build a mesh of 15
  // QUAD9 elements.  Building QUAD9 elements instead of the default QUAD4s
  // allow us to use higher-order approximation.
  libMesh::MeshTools::Generation::build_square(*mesh, 15, 0.0, 1.0, QUAD9);

  // Create an equation systems object.
  equations_systems = new libMesh::EquationSystems(*mesh);
  
  // Declare the zero funtion equations system
  equations_systems->add_system<libMesh::ExplicitSystem>("Zero function");

  // Adds the variable "u".  "u" will be approximated using second-order
  // approximation.
  equations_systems->get_system("Zero function").add_variable("u", SECOND);

  // Initialize the data structures for the equation system.
  equations_systems->init();
}

uqLibMeshFunction::~uqLibMeshFunction()
{
  delete mesh;
  delete equations_systems;
}

void uqLibMeshFunction::print_info()
{
  // Print information about the mesh to the screen.
  mesh->print_info(std::cerr);
}
