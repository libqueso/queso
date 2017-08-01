#include <queso/Environment.h>

#ifdef QUESO_HAVE_LIBMESH
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshFunction.h>
#include <queso/LibMeshNegativeLaplacianOperator.h>
#include <queso/InfiniteDimensionalGaussian.h>
#include <queso/InfiniteDimensionalLikelihoodBase.h>
#include <queso/InfiniteDimensionalMCMCSampler.h>
#include <queso/InfiniteDimensionalMCMCSamplerOptions.h>

#include <cstdlib>
#include <memory>
#include <sstream>

class Likelihood : public QUESO::InfiniteDimensionalLikelihoodBase {
public:
  Likelihood();
  ~Likelihood();
  virtual double evaluate(QUESO::FunctionBase &flow);
};

Likelihood::Likelihood()
  : QUESO::InfiniteDimensionalLikelihoodBase(1.0)
{
}

Likelihood::~Likelihood() {}

double Likelihood::evaluate(QUESO::FunctionBase &flow)
{
  return 1.0;
}
#endif

#ifdef QUESO_HAVE_LIBMESH
int main(int argc, char **argv)
{
  std::string in_file_name = "test_infinite/inf_options";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    in_file_name = test_srcdir + ('/' + in_file_name);

  const char * prefix = "";
  const unsigned int num_pairs = 5;
  const double alpha = 3.0;
  const double beta = 1.0;
  // EnvOptionsValuesClass opts;
  // opts.m_seed = -1;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, in_file_name, prefix, NULL);
#else
  QUESO::FullEnvironment env(in_file_name, prefix, NULL);
#endif

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc farts with libMesh::Real==float
  libmesh_example_assert(false, "--disable-singleprecision");
#endif


// Need an artificial block here because libmesh needs to
// call PetscFinalize before we call MPI_Finalize
#if defined(LIBMESH_HAVE_SLEPC) && defined(QUESO_HAVE_HDF5)
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

  Likelihood llhd;

  QUESO::InfiniteDimensionalMCMCSamplerOptions opts(env, "");
  QUESO::InfiniteDimensionalMCMCSampler sampler(env, mu, llhd, &opts);

  queso_require_equal_to(opts.m_num_iters, 10);
  queso_require_equal_to(opts.m_save_freq, 5);
  queso_require_equal_to(opts.m_rwmh_step, 0.1);

  return 0;
}
#endif  // LIBMESH_HAVE_SLEPC

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
#else
int main()
{
  std::cout << "QUESO was configured without libMesh support." << std::endl;
  return 77;
}
#endif
