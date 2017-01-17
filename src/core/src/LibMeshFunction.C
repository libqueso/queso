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
#include <queso/LibMeshFunction.h>
#include <queso/FunctionOperatorBuilder.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/system_norm.h>
#include <libmesh/exodusII_io.h>

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

LibMeshFunction::LibMeshFunction(
    const FunctionOperatorBuilder & builder, libMesh::MeshBase & m)
  : FunctionBase(),
    builder(builder),
    equation_systems(new libMesh::EquationSystems(m))
{
  this->equation_systems->add_system<libMesh::ExplicitSystem>("Function");
  this->equation_systems->get_system("Function").add_variable("u",
      libMesh::Utility::string_to_enum<libMeshEnums::Order>(this->builder.order),
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>(this->builder.family));
  this->equation_systems->init();
}

LibMeshFunction::~LibMeshFunction()
{
}

void LibMeshFunction::print_info() const
{
  // Print information about the mesh to the screen.
  this->equation_systems->get_mesh().print_info(std::cerr);
}

void LibMeshFunction::save_function(const std::string & filename, double time) const
{
  // The "1" is hardcoded for the number of time steps because the ExodusII
  // manual states that it should be the number of timesteps within the file.
  // Here, we are explicitly only doing one timestep per file.
  libMesh::ExodusII_IO(this->equation_systems->get_mesh()).write_timestep(
      filename, *this->equation_systems, 1, time);
}

void LibMeshFunction::add(double scale, const FunctionBase & rhs) {
  // We know we're dealing with a derived class type, so cast
  const LibMeshFunction & rhs_derived = dynamic_cast<
    const LibMeshFunction &>(rhs);

  this->equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution->add(
      scale, *(rhs_derived.equation_systems->get_system<libMesh::ExplicitSystem>(
          "Function").solution));
}

void LibMeshFunction::pointwise_mult(const FunctionBase & f1,
    const FunctionBase & f2)
{
  const LibMeshFunction & f1_derived = static_cast<
    const LibMeshFunction &>(f1);
  const LibMeshFunction & f2_derived = static_cast<
    const LibMeshFunction &>(f2);

  this->equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution->pointwise_mult(
      *(f1_derived.equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution),
      *(f2_derived.equation_systems->get_system<libMesh::ExplicitSystem>("Function").solution));
}

void LibMeshFunction::scale(double scale) {
  this->equation_systems->get_system<libMesh::ExplicitSystem>(
      "Function").solution->scale(scale);
}

void LibMeshFunction::zero() {
  this->equation_systems->get_system<libMesh::ExplicitSystem>(
      "Function").solution->zero();
}

double LibMeshFunction::L2_norm() const {
  libMesh::ExplicitSystem & system =
    this->equation_systems->get_system<libMesh::ExplicitSystem>("Function");

  double norm = system.calculate_norm(*system.solution,
                                      libMesh::SystemNorm(libMeshEnums::L2));
  return norm;
}

SharedPtr<FunctionBase>::Type LibMeshFunction::zero_clone() const
{
  LibMeshFunction * clone = new LibMeshFunction(this->builder,
      this->equation_systems->get_mesh());
  clone->equation_systems->get_system<libMesh::ExplicitSystem>(
      "Function").solution->zero();

  SharedPtr<FunctionBase>::Type ptr(clone);
  return ptr;
}

SharedPtr<libMesh::EquationSystems>::Type
LibMeshFunction::get_equation_systems() const {
  return this->equation_systems;
}

}  // End namespace QUESO

#endif  // QUESO_HAVE_LIBMESH_SLEPC
