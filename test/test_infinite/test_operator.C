#include <memory>
#include <vector>
#include <cmath>

#include <queso/Environment.h>

#ifdef QUESO_HAVE_LIBMESH
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/condensed_eigen_system.h>
#include <libmesh/fe.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/elem.h>
#include <libmesh/dof_map.h>
#include <libmesh/system_norm.h>
#include <queso/EnvironmentOptions.h>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshFunction.h>
#include <queso/LibMeshNegativeLaplacianOperator.h>
#include <queso/InfiniteDimensionalGaussian.h>
#endif  // QUESO_HAVE_LIBMESH

#define TEST_TOL 1e-8
#define INTEGRATE_TOL 1e-2

int main(int argc, char **argv)
{
#ifdef QUESO_HAVE_LIBMESH
  unsigned int i, j;
  QUESO::EnvOptionsValues opts;
  opts.m_seed = -1;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &opts);
#else
  QUESO::FullEnvironment env("", "", &opts);
#endif

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_assert(false, "--disable-singleprecision");
#endif

// Need an artificial block here because libmesh needs to
// call PetscFinalize before we call MPI_Finalize
#ifdef LIBMESH_HAVE_SLEPC
{
  libMesh::LibMeshInit init(argc, argv);

  libMesh::Mesh mesh(init.comm());
  libMesh::MeshTools::Generation::build_square(mesh,
      20, 20, 0.0, 1.0, 0.0, 1.0, libMeshEnums::QUAD4);

  QUESO::FunctionOperatorBuilder builder;

  builder.order = "FIRST";
  builder.family = "LAGRANGE";
  builder.num_req_eigenpairs = 10;

  QUESO::LibMeshNegativeLaplacianOperator precision(builder, mesh);

  libMesh::EquationSystems & es = precision.get_equation_systems();
  libMesh::CondensedEigenSystem & eig_sys = es.get_system<libMesh::CondensedEigenSystem>(
      "Eigensystem");

  // Check all eigenfunctions have unit L2 norm
  std::vector<double> norms(builder.num_req_eigenpairs, 0);
  for (i = 0; i < builder.num_req_eigenpairs; i++) {
    eig_sys.get_eigenpair(i);
    norms[i] = eig_sys.calculate_norm(*eig_sys.solution,
                                      libMesh::SystemNorm(libMeshEnums::L2));
    if (abs(norms[i] - 1.0) > TEST_TOL) {
      return 1;
    }
  }

  const unsigned int dim = mesh.mesh_dimension();
  const libMesh::DofMap & dof_map = eig_sys.get_dof_map();
  libMesh::FEType fe_type = dof_map.variable_type(0);
  libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
  libMesh::QGauss qrule(dim, libMeshEnums::FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<libMesh::Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();

  libMesh::UniquePtr<libMesh::NumericVector<libMesh::Real> > u, v;
  double ui = 0.0;
  double vj = 0.0;
  double ip = 0.0;

  for (i = 0; i < builder.num_req_eigenpairs - 1; i++) {
    eig_sys.get_eigenpair(i);
    u = eig_sys.solution->clone();
    for (j = i + 1; j < builder.num_req_eigenpairs; j++) {
      libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
      libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      eig_sys.get_eigenpair(j);
      v = eig_sys.solution->clone();
      for ( ; el != end_el; ++el) {
        const libMesh::Elem * elem = *el;
        fe->reinit(elem);
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
          for (unsigned int dof = 0; dof < phi.size(); dof++) {
            ui += (*u)(dof) * phi[dof][qp];
            vj += (*v)(dof) * phi[dof][qp];
          }
          ip += ui * vj * JxW[qp];
          ui = 0.0;
          vj = 0.0;
        }
      }
      std::cerr << "INTEGRAL of " << i << " against " << j << " is: " << ip << std::endl;
      if (abs(ip) > INTEGRATE_TOL) {
        return 1;
      }
      ip = 0.0;
    }
  }
}
#endif  // LIBMESH_HAVE_SLEPC
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
#else
  return 77;
#endif
}
