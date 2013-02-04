//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#include <uqLibMeshNegativeLaplacianOperator.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/equation_systems.h>
#include <libmesh/eigen_system.h>
#include <libmesh/auto_ptr.h>
#include <libmesh/fe_base.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>
#include <libmesh/id_types.h>
#include <libmesh/elem.h>

uqLibMeshNegativeLaplacianOperator::uqLibMeshNegativeLaplacianOperator()
  : uqLibMeshOperatorBase()
{
}

uqLibMeshNegativeLaplacianOperator::uqLibMeshNegativeLaplacianOperator(const std::string& filename)
  : uqLibMeshOperatorBase(filename)
{
  // Refactor this out
  unsigned int n_evals = 10;
  // Give the system a pointer to the matrix assembly
  // function defined below.
  libMesh::EigenSystem & eigen_system =
    this->equation_systems->get_system<libMesh::EigenSystem>("Eigensystem");
  eigen_system.attach_assemble_object(*this);

  // Set the number of requested eigenpairs \p n_evals and the number
  // of basis vectors used in the solution algorithm.
  equation_systems->parameters.set<unsigned int>("eigenpairs")    = n_evals;
  equation_systems->parameters.set<unsigned int>("basis vectors") = n_evals*3;

  // Set the solver tolerance and the maximum number of iterations. 
  equation_systems->parameters.set<libMesh::Real>("linear solver tolerance") = pow(libMesh::TOLERANCE, 5./3.);
  equation_systems->parameters.set<unsigned int>
    ("linear solver maximum iterations") = 1000;

  // Set the type of the problem, here we deal with
  // a generalized Hermitian problem.
  eigen_system.set_eigenproblem_type(GHEP);

  // Order the eigenvalues "smallest first"
  eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);

  // Initialize the data structures for the equation system.
  equation_systems->init();

  // Prints information about the system to the screen.
  equation_systems->print_info();

  // Pass the Dirichlet dof IDs to the CondensedEigenSystem
  // std::set<unsigned int> dirichlet_dof_ids;
  // get_boundary_dofs(*equation_systems, "Eigensystem", dirichlet_dof_ids);
  // eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);

  // Solve the system "Eigensystem".
  eigen_system.solve();

  // Get the number of converged eigen pairs.
  unsigned int nconv = eigen_system.get_n_converged();

  std::cout << "Number of converged eigenpairs: " << nconv
            << "\n" << std::endl;

  // if (plotting_index > n_evals)
  //   {
  //     std::cout << "WARNING: Solver did not converge for the requested eigenvector!" << std::endl;
  //   }

  // // write out all of the computed eigenvalues and plot the specified eigenvector
  // std::ostringstream eigenvalue_output_name;
  // eigenvalue_output_name << mesh_name << "_evals.txt";
  // std::ofstream evals_file(eigenvalue_output_name.str().c_str());

  // for(unsigned int i=0; i<nconv; i++)
  // {
  //   std::pair<Real,Real> eval = eigen_system.get_eigenpair(i);
  //   
  //   // The eigenvalues should be real!
  //   libmesh_assert_less (eval.second, TOLERANCE);
  //   evals_file << eval.first << std::endl;
  //   
  // }
  // 
  // evals_file.close();
}

uqLibMeshNegativeLaplacianOperator::~uqLibMeshNegativeLaplacianOperator()
{
  // delete this->mesh;
  // delete this->equations_systems;
}

void uqLibMeshNegativeLaplacianOperator::assemble()
{
#ifdef LIBMESH_HAVE_SLEPC
  // Get a constant reference to the mesh object.
  const libMesh::MeshBase& mesh = this->equation_systems->get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system.
  libMesh::EigenSystem & eigen_system = this->equation_systems->get_system<libMesh::EigenSystem>("Eigensystem");

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  libMesh::FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  // A reference to the two system matrices
  libMesh::SparseMatrix<libMesh::Number>&  matrix_A = *eigen_system.matrix_A;
  libMesh::SparseMatrix<libMesh::Number>&  matrix_B = *eigen_system.matrix_B;

  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p AutoPtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  libMesh::AutoPtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
  
  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  libMesh::QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<libMesh::Real>& JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const libMesh::DofMap& dof_map = eigen_system.get_dof_map();

  // The element mass and stiffness matrices.
  libMesh::DenseMatrix<libMesh::Number>   Me;
  libMesh::DenseMatrix<libMesh::Number>   Ke;

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
  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
 
  for ( ; el != end_el; ++el)
    {
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
              Me(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp];
              Ke(i,j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
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

void uqLibMeshNegativeLaplacianOperator::get_boundary_dofs(
    std::set<unsigned int>& global_boundary_dofs_set)
{
//   // This function gets the Dirichlet degrees of freedom
// #ifdef LIBMESH_HAVE_SLEPC
//   dirichlet_dof_ids.clear();
// 
//   // It is a good idea to make sure we are assembling
//   // the proper system.
//   // libmesh_assert_equal_to(system_name, "Eigensystem");
// 
//   // Get a constant reference to the mesh object.
//   const MeshBase& mesh = this->equations_systems->get_mesh();
// 
//   // The dimension that we are running.
//   const unsigned int dim = mesh.mesh_dimension();
// 
//   // Get a reference to our system.
//   EigenSystem & eigen_system = this->equations_systems->get_system<EigenSystem>(system_name);
// 
//   // Get a constant reference to the Finite Element type
//   // for the first (and only) variable in the system.
//   FEType fe_type = eigen_system.get_dof_map().variable_type(0);
//   
//   const DofMap& dof_map = eigen_system.get_dof_map();
// 
//   // This vector will hold the degree of freedom indices for
//   // the element.  These define where in the global system
//   // the element degrees of freedom get mapped.
//   std::vector<unsigned int> dof_indices;
// 
//   // Now we will loop over all the elements in the mesh that
//   // live on the local processor. We will compute the element
//   // matrix and right-hand-side contribution.  In case users
//   // later modify this program to include refinement, we will
//   // be safe and will only consider the active elements;
//   // hence we use a variant of the \p active_elem_iterator.
//   MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//  
//   for ( ; el != end_el; ++el)
//     {
//       // Store a pointer to the element we are currently
//       // working on.  This allows for nicer syntax later.
//       const Elem* elem = *el;
// 
//       // Get the degree of freedom indices for the
//       // current element.  These define where in the global
//       // matrix and right-hand-side this element will
//       // contribute to.
//       dof_map.dof_indices(elem, dof_indices);
// 
//      {
//         // All boundary dofs are Dirichlet dofs in this case
//         for (unsigned int s = 0; s < elem->n_sides(); s++)
//           if (elem->neighbor(s) == NULL)
//             {
//               std::vector<unsigned int> side_dofs;
//               FEInterface::dofs_on_side(elem, dim, fe_type,
//                                         s, side_dofs);
// 
//               for(unsigned int ii = 0; ii < side_dofs.size(); ii++)
//                 dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
//             }
//       }
// 
//     } // end of element loop
// 
// #endif // LIBMESH_HAVE_SLEPC
}
