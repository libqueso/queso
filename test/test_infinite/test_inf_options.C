#include <memory>
#include <sstream>
#include <uqEnvironment.h>

#ifdef QUESO_HAVE_LIBMESH
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <uqFunctionOperatorBuilder.h>
#include <uqLibMeshFunction.h>
#include <uqLibMeshNegativeLaplacianOperator.h>
#include <uqInfiniteDimensionalGaussian.h>
#include <uqInfiniteDimensionalLikelihoodBase.h>
#include <uqInfiniteDimensionalMCMCSampler.h>
#include <uqInfiniteDimensionalMCMCSamplerOptions.h>
#endif

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

#ifdef QUESO_HAVE_LIBMESH
using namespace libMesh;
#endif

class Likelihood : public uqInfiniteDimensionalLikelihoodBase {
public:
  Likelihood();
  ~Likelihood();
  virtual double evaluate(uqFunctionBase &flow);
};

Likelihood::Likelihood()
  : uqInfiniteDimensionalLikelihoodBase(1.0)
{
}

Likelihood::~Likelihood() {}

double Likelihood::evaluate(uqFunctionBase &flow)
{
  return 1.0;
}

int main(int argc, char **argv)
{
  const char * in_file_name = "test_infinite/inf_options_in.txt";
  const char * prefix = "";
  const unsigned int num_pairs = 5;
  const double alpha = 3.0;
  const double beta = 1.0;
  // uqEnvOptionsValuesClass opts;
  // opts.m_seed = -1;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef QUESO_HAS_MPI
  uqFullEnvironmentClass env(MPI_COMM_WORLD, in_file_name, prefix, NULL);
#else
  uqFullEnvironmentClass env(0, in_file_name, prefix, NULL);
#endif

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc farts with libMesh::Real==float
  libmesh_example_assert(false, "--disable-singleprecision");
#endif


// Need an artificial block here because libmesh needs to
// call PetscFinalize before we call MPI_Finalize
#ifdef LIBMESH_HAVE_SLEPC
{
  LibMeshInit init(argc, argv);

  Mesh mesh;
  MeshTools::Generation::build_square(mesh,
      20, 20, 0.0, 1.0, 0.0, 1.0, QUAD4);

  uqFunctionOperatorBuilder fobuilder;

  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = num_pairs;

  uqLibMeshFunction mean(fobuilder, mesh);
  uqLibMeshNegativeLaplacianOperator precision(fobuilder, mesh);
  uqInfiniteDimensionalGaussian mu(env, mean, precision, alpha, beta);

  Likelihood llhd;

  uqInfiniteDimensionalMCMCSamplerOptions opts(env, "");
  uqInfiniteDimensionalMCMCSampler sampler(env, mu, llhd, &opts);

  if (opts.m_num_iters != 10 ||  // Input file value is 10
      opts.m_save_freq != 5 ||  // Ditto 5
      opts.m_rwmh_step != 0.1) {
    return 1;
  }

  return 0;
}
#endif  // LIBMESH_HAVE_SLEPC

  MPI_Finalize();
  return 0;
}
