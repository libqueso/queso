#include <memory>
#include <sstream>
#include <queso/Environment.h>

#ifdef QUESO_HAVE_LIBMESH
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <queso/uqFunctionOperatorBuilder.h>
#include <queso/uqLibMeshFunction.h>
#include <queso/uqLibMeshNegativeLaplacianOperator.h>
#include <queso/uqInfiniteDimensionalGaussian.h>
#include <queso/uqInfiniteDimensionalLikelihoodBase.h>
#include <queso/uqInfiniteDimensionalMCMCSampler.h>
#include <queso/uqInfiniteDimensionalMCMCSamplerOptions.h>
#endif

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

#ifdef QUESO_HAVE_LIBMESH
class Likelihood : public QUESO::uqInfiniteDimensionalLikelihoodBase {
public:
  Likelihood();
  ~Likelihood();
  virtual double evaluate(QUESO::uqFunctionBase &flow);
};

Likelihood::Likelihood()
  : QUESO::uqInfiniteDimensionalLikelihoodBase(1.0)
{
}

Likelihood::~Likelihood() {}

double Likelihood::evaluate(QUESO::uqFunctionBase &flow)
{
  return 1.0;
}
#endif

int main(int argc, char **argv)
{
#ifdef QUESO_HAVE_LIBMESH
  const char * in_file_name = "test_infinite/inf_options";
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
  QUESO::FullEnvironment env(MPI_COMM_WORLD, in_file_name, prefix, NULL);
#else
  QUESO::FullEnvironment env(0, in_file_name, prefix, NULL);
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

  libMesh::Mesh mesh;
  libMesh::MeshTools::Generation::build_square(mesh,
      20, 20, 0.0, 1.0, 0.0, 1.0, QUAD4);

  QUESO::uqFunctionOperatorBuilder fobuilder;

  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = num_pairs;

  QUESO::uqLibMeshFunction mean(fobuilder, mesh);
  QUESO::uqLibMeshNegativeLaplacianOperator precision(fobuilder, mesh);
  QUESO::uqInfiniteDimensionalGaussian mu(env, mean, precision, alpha, beta);

  Likelihood llhd;

  QUESO::uqInfiniteDimensionalMCMCSamplerOptions opts(env, "");
  QUESO::uqInfiniteDimensionalMCMCSampler sampler(env, mu, llhd, &opts);

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
#else
  return 77;
#endif
}
