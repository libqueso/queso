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

#include <set>
#include <vector>

#include <queso/SharedPtr.h>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshNegativeLaplacianOperator.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/condensed_eigen_system.h>
#include <libmesh/auto_ptr.h>
#include <libmesh/fe_base.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>
#include <libmesh/id_types.h>
#include <libmesh/elem.h>
#include <libmesh/zero_function.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/utility.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_fe_family.h>

namespace QUESO {

LibMeshNegativeLaplacianOperator::LibMeshNegativeLaplacianOperator(
    const FunctionOperatorBuilder & builder, libMesh::MeshBase & m)
  : LibMeshOperatorBase(builder, m)
{
  SharedPtr<libMesh::EquationSystems>::Type es(this->equation_systems);

  // Give the system a pointer to the matrix assembly
  // function defined below.
  libMesh::CondensedEigenSystem & eigen_system =
    es->get_system<libMesh::CondensedEigenSystem>("Eigensystem");

  // Declare the system variables.
  // Adds the variable "u" to "Eigensystem".   "u"
  // will be approximated using second-order approximation.
  unsigned int u_var = eigen_system.add_variable("u",
      libMesh::Utility::string_to_enum<libMeshEnums::Order>(this->builder.order),
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>(this->builder.family));

  // This works because *this is a subclass of System::Assembly
  // and requires the class to implement an 'assemble'
  eigen_system.attach_assemble_object(*this);

  // Set the number of requested eigenpairs \p n_evals and the number
  // of basis vectors used in the solution algorithm.
  es->parameters.set<unsigned int>("eigenpairs") =
    this->builder.num_req_eigenpairs;
  es->parameters.set<unsigned int>("basis vectors") =
    this->builder.num_req_eigenpairs * 3;

  // Set the solver tolerance and the maximum number of iterations.
  es->parameters.set<libMesh::Real>("linear solver tolerance") =
    pow(libMesh::TOLERANCE, 5./3.);
  es->parameters.set<unsigned int>("linear solver maximum iterations") = 1000;

  // Set the type of the problem, here we deal with
  // a generalized Hermitian problem.
  eigen_system.set_eigenproblem_type(libMeshEnums::GHEP);

  // Order the eigenvalues "smallest first"
  // This hoses performance?
  eigen_system.eigen_solver->set_position_of_spectrum
    (libMeshEnums::SMALLEST_MAGNITUDE);

  // Set up the boundary (only works if this->m is a square)
  // We'll just the whole boundary to be Dirichlet, because why not
  std::set<libMesh::boundary_id_type> boundary_ids;
  boundary_ids.insert(0);
  boundary_ids.insert(1);
  boundary_ids.insert(2);
  boundary_ids.insert(3);

  // Assign which variables must satisfy the boundary constraints
  std::vector<unsigned int> vars;
  vars.push_back(u_var);

  // Add the boundary object to the eigensystem
  libMesh::ZeroFunction<> zf;
  libMesh::DirichletBoundary dirichlet_bc(boundary_ids, vars, &zf);
  eigen_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Initialize the data structures for the equation system.
  es->init();

  // Here we obtain the ids of the Dirichlet dofs
  std::set<unsigned int> global_dirichlet_dof_ids;
  libMesh::DofConstraints::const_iterator it =
    eigen_system.get_dof_map().constraint_rows_begin();

  for ( ; it != eigen_system.get_dof_map().constraint_rows_end(); ++it) {
    global_dirichlet_dof_ids.insert(it->first);
  }

  eigen_system.initialize_condensed_dofs(global_dirichlet_dof_ids);

  // Solve the system "Eigensystem".
  eigen_system.solve();

  // Get the number of converged eigen pairs.
  this->nconv = eigen_system.get_n_converged();
}

LibMeshNegativeLaplacianOperator::~LibMeshNegativeLaplacianOperator()
{
}

void LibMeshNegativeLaplacianOperator::assemble()
{
#ifdef LIBMESH_HAVE_SLEPC

  SharedPtr<libMesh::EquationSystems>::Type es(this->equation_systems);

  // Get a constant reference to the mesh object.
  const libMesh::MeshBase& mesh = es->get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system.
  libMesh::EigenSystem & eigen_system = es->get_system<libMesh::EigenSystem>("Eigensystem");

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  libMesh::FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  // A reference to the two system matrices
  libMesh::SparseMatrix<libMesh::Number>&  matrix_A = *eigen_system.matrix_A;
  libMesh::SparseMatrix<libMesh::Number>&  matrix_B = *eigen_system.matrix_B;

  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  libMesh::QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<libMesh::Real>& JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const libMesh::DofMap& dof_map = eigen_system.get_dof_map();

  // The element stiffness matrix.
  libMesh::DenseMatrix<libMesh::Number> Me;
  libMesh::DenseMatrix<libMesh::Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<libMesh::dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the \p active_elem_iterator.
  libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el) {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const libMesh::Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices(elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit(elem);

      // Zero the element matrices before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize(dof_indices.size(), dof_indices.size());
      Me.resize(dof_indices.size(), dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //
      // We will build the element matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp < qrule.n_points(); qp++)
        for (unsigned int i=0; i < phi.size(); i++)
          for (unsigned int j=0; j < phi.size(); j++)
            {
              Me(i, j) += JxW[qp] * phi[i][qp] * phi[j][qp];
              Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
            }

      // On an unrefined mesh, constrain_element_matrix does
      // nothing.  If this assembly function is ever repurposed to
      // run on a refined mesh, getting the hanging node constraints
      // right will be important.  Note that, even with
      // asymmetric_constraint_rows = false, the constrained dof
      // diagonals still exist in the matrix, with diagonal entries
      // that are there to ensure non-singular matrices for linear
      // solves but which would generate positive non-physical
      // eigenvalues for eigensolves.
      // dof_map.constrain_element_matrix(Ke, dof_indices, false);
      // dof_map.constrain_element_matrix(Me, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrices A and B.
      matrix_A.add_matrix (Ke, dof_indices);
      matrix_B.add_matrix (Me, dof_indices);

    } // end of element loop

#endif // LIBMESH_HAVE_SLEPC
}

void LibMeshNegativeLaplacianOperator::print_info() const
{
  // Prints information about the system to the screen.
  this->equation_systems->get_mesh().print_info();
  this->equation_systems->print_info();
}

}  // End namespace QUESO

#endif  // QUESO_HAVE_LIBMESH_SLEPC
