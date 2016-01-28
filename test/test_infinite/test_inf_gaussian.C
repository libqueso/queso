#include <memory>
#include <sstream>
#include <queso/Environment.h>

#ifdef QUESO_HAVE_LIBMESH
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <queso/EnvironmentOptions.h>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshFunction.h>
#include <queso/LibMeshNegativeLaplacianOperator.h>
#include <queso/InfiniteDimensionalGaussian.h>
#endif  // QUESO_HAVE_LIBMESH

int main(int argc, char **argv)
{
#ifdef QUESO_HAVE_LIBMESH
  unsigned int i;
  unsigned int j;
  const unsigned int num_pairs = 5;
  const unsigned int num_samples = 1e4;
  const double alpha = 3.0;
  const double beta = 1.0;
  QUESO::EnvOptionsValues opts;
  opts.m_seed = -1;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &opts);
#else
  QUESO::FullEnvironment env("", "", &opts);
#endif

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc farts with libMesh::Real==float
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

  QUESO::FunctionOperatorBuilder fobuilder;

  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = num_pairs;

  QUESO::LibMeshFunction mean(fobuilder, mesh);
  QUESO::LibMeshNegativeLaplacianOperator precision(fobuilder, mesh);
  QUESO::InfiniteDimensionalGaussian mu(env, mean, precision, alpha, beta);

  // Vector to hold all KL coeffs
  std::vector<double> means(num_pairs, 0.0);
  std::vector<double> sumsqs(num_pairs, 0.0);
  std::vector<double> deltas(num_pairs, 0.0);
  double draw;

  for (i = 1; i < num_samples + 1; i++) {
    mu.draw();
    for (j = 0; j < num_pairs; j++) {
      draw = mu.get_kl_coefficient(j);
      deltas[j] = draw - means[j];
      means[j] += (double) deltas[j] / i;
      sumsqs[j] += deltas[j] * (draw - means[j]);
    }
    // std::cerr << "MEAN IS: " << means[0] << std::endl;
  }

  std::vector<double> vars(num_pairs, 0.0);
  for (j = 0; j < num_pairs; j++) {
    vars[j] = sumsqs[j] / (num_samples - 1);
  }

  double sigma = beta / std::pow(precision.get_eigenvalue(j), alpha / 2.0);
  double sigmasq = sigma * sigma;
  double mean_min;
  double mean_max;

  for (j = 0; j < num_pairs; j++) {
    // Mean is N(0, (lambda_j^{- alpha / 2} * beta)^2 / n)
    mean_min = -3.0 * sigma / std::sqrt(num_samples);
    mean_max =  3.0 * sigma / std::sqrt(num_samples);
    if (means[j] < mean_min || means[j] > mean_max) {
      std::cerr << "mean kl test failed" << std::endl;
      return 1;
    }
  }

  double var_min;
  double var_max;

  // var[j] should be approximately ~ N(sigma^2, 2 sigma^4 / (num_samples - 1))
  for (j = 0; j < num_pairs; j++) {
    var_min = sigmasq - 3.0 * sigmasq * std::sqrt(2.0 / (num_samples - 1));
    var_max = sigmasq + 3.0 * sigmasq * std::sqrt(2.0 / (num_samples - 1));
    if (vars[j] < var_min || vars[j] > var_max) {
      std::cerr << "variance kl test failed" << std::endl;
      return 1;
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
