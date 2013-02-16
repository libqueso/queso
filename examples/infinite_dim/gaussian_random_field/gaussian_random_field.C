#include <memory>
#include <sstream>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <uqEnvironment.h>
#include <uqFunctionOperatorBuilder.h>
#include <uqLibMeshFunction.h>
#include <uqLibMeshNegativeLaplacianOperator.h>
#include <uqInfiniteDimensionalGaussian.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

using namespace libMesh;

int main(int argc, char **argv)
{
  unsigned int i;
  uqEnvOptionsValuesClass opts;
  opts.m_seed = -1;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef QUESO_HAS_MPI
  uqFullEnvironmentClass env(MPI_COMM_WORLD, "", "", &opts);
#else
  uqFullEnvironmentClass env(0, "", "", &opts);
#endif


// Need an artificial block here because libmesh needs to
// call PetscFinalize before we call MPI_Finalize
{
  LibMeshInit init(argc, argv);

  Mesh mesh;
  MeshTools::Generation::build_square(mesh,
      20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);

  uqFunctionOperatorBuilder fobuilder;

  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = 360;

  uqLibMeshFunction mean(fobuilder, mesh);
  uqLibMeshNegativeLaplacianOperator precision(fobuilder, mesh);

  precision.print_info();
  precision.save_converged_evals("evals.txt");

  uqInfiniteDimensionalGaussian mu(env, mean, precision, 3.0, 1.0);

  for (i = 0; i < 5; i++) {
    std::ostringstream number;
    number << i;
    mu.draw()->save_function("rand_draw" + number.str() + ".e");
  }
}

  MPI_Finalize();

  return 0;
}
