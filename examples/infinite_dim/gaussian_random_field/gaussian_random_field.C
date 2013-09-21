#include <memory>
#include <sstream>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <queso/Environment.h>
#include <queso/uqFunctionOperatorBuilder.h>
#include <queso/uqLibMeshFunction.h>
#include <queso/uqLibMeshNegativeLaplacianOperator.h>
#include <queso/uqInfiniteDimensionalGaussian.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

using namespace libMesh;

int main(int argc, char **argv)
{
  unsigned int i;
  QUESO::EnvOptionsValues opts;
  opts.m_seed = -1;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &opts);
#else
  QUESO::FullEnvironment env(0, "", "", &opts);
#endif


// Need an artificial block here because libmesh needs to
// call PetscFinalize before we call MPI_Finalize
{
  LibMeshInit init(argc, argv);

  Mesh mesh;
  MeshTools::Generation::build_square(mesh,
      20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);

  QUESO::uqFunctionOperatorBuilder fobuilder;

  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = 360;

  QUESO::uqLibMeshFunction mean(fobuilder, mesh);
  QUESO::uqLibMeshNegativeLaplacianOperator precision(fobuilder, mesh);

  precision.print_info();
  precision.save_converged_evals("evals.txt");

  QUESO::uqInfiniteDimensionalGaussian mu(env, mean, precision, 3.0, 1.0);

  for (i = 0; i < 5; i++) {
    std::ostringstream number;
    number << i;
    mu.draw()->save_function("rand_draw" + number.str() + ".e");
  }
}

  MPI_Finalize();

  return 0;
}
