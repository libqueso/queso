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
#include <queso/uqLibMeshFunction.h>
#include <queso/uqFunctionOperatorBuilder.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/system_norm.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/vtk_io.h>

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
#include <libmesh/utility.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_fe_family.h>

namespace QUESO {

uqLibMeshFunction::uqLibMeshFunction(
    const uqFunctionOperatorBuilder & builder, libMesh::MeshBase & m)
  : uqFunctionBase(builder),
    equation_systems(new libMesh::EquationSystems(m))
{
  this->equation_systems->add_system<libMesh::ExplicitSystem>("Function");
  this->equation_systems->get_system("Function").add_variable("u",
      libMesh::Utility::string_to_enum<libMeshEnums::Order>(this->builder.order),
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>(this->builder.family));
  this->equation_systems->init();
}

uqLibMeshFunction::~uqLibMeshFunction()
{
}

void uqLibMeshFunction::print_info() const
{
  // Print information about the mesh to the screen.
  this->equation_systems->get_mesh().print_info(std::cerr);
}

void uqLibMeshFunction::save_function(const std::string & filename) const
{
  libMesh::VTKIO(this->equation_systems->get_mesh()).write_equation_systems(
      filename, *this->equation_systems);
}

void uqLibMeshFunction::add(double scale, const uqFunctionBase & rhs) {
  // We know we're dealing with a derived class type, so cast
  const uqLibMeshFunction & rhs_derived = dynamic_cast<
    const uqLibMeshFunction &>(rhs);

  this->equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution->add(
      scale, *(rhs_derived.equation_systems->get_system<libMesh::ExplicitSystem>(
          "Function").solution));
}

void uqLibMeshFunction::pointwise_mult(const uqFunctionBase & f1,
    const uqFunctionBase & f2)
{
  const uqLibMeshFunction & f1_derived = static_cast<
    const uqLibMeshFunction &>(f1);
  const uqLibMeshFunction & f2_derived = static_cast<
    const uqLibMeshFunction &>(f2);

  this->equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution->pointwise_mult(
      *(f1_derived.equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution),
      *(f2_derived.equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution));
}

void uqLibMeshFunction::scale(double scale) {
  this->equation_systems->get_system<libMesh::ExplicitSystem>(
      "Function").solution->scale(scale);
}

void uqLibMeshFunction::zero() {
  this->equation_systems->get_system<libMesh::ExplicitSystem>(
      "Function").solution->zero();
}

double uqLibMeshFunction::L2_norm() const {
  libMesh::ExplicitSystem & system =
    this->equation_systems->get_system<libMesh::ExplicitSystem>("Function");

  double norm = system.calculate_norm(*system.solution, libMesh::SystemNorm(L2));
  return norm;
}

boost::shared_ptr<uqFunctionBase> uqLibMeshFunction::zero_clone() const
{
  uqLibMeshFunction * clone = new uqLibMeshFunction(this->builder,
      this->equation_systems->get_mesh());
  clone->equation_systems->get_system<libMesh::ExplicitSystem>(
      "Function").solution->zero();

  boost::shared_ptr<uqFunctionBase> ptr(clone);
  return ptr;
}

}  // End namespace QUESO
