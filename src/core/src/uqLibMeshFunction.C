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
  // QUAD9 elements.  Building QUAD4 elements instead of the default QUAD4s
  // allow us to use higher-order approximation.
  libMesh::MeshTools::Generation::build_square(*mesh, 20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);

  // Create an equation systems object.
  this->equation_systems = new libMesh::EquationSystems(*mesh);
  
  // Declare the zero funtion equations system
  this->equation_systems->add_system<libMesh::ExplicitSystem>("Function");

  // Adds the variable "u".  "u" will be approximated using second-order
  // approximation.
  this->equation_systems->get_system("Function").add_variable("u", FIRST);

  // Initialize the data structures for the equation system.
  this->equation_systems->init();
}

uqLibMeshFunction::~uqLibMeshFunction()
{
  delete this->mesh;
  delete this->equation_systems;
}

void uqLibMeshFunction::print_info()
{
  // Print information about the mesh to the screen.
  this->mesh->print_info(std::cerr);
}
