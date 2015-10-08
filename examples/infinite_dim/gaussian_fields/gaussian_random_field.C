#include <memory>
#include <sstream>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshFunction.h>
#include <queso/LibMeshNegativeLaplacianOperator.h>
#include <queso/InfiniteDimensionalGaussian.h>

#include <mpi.h>

using namespace libMesh;

int main(int argc, char **argv)
{
  unsigned int i;
  QUESO::EnvOptionsValues opts;
  opts.m_seed = -1;

  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &opts);

// Need an artificial block here because libmesh needs to
// call PetscFinalize before we call MPI_Finalize
{
  LibMeshInit init(argc, argv);

  Mesh mesh(init.comm());
  MeshTools::Generation::build_square(mesh,
      20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);

  QUESO::FunctionOperatorBuilder fobuilder;

  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = 360;

  QUESO::LibMeshFunction mean(fobuilder, mesh);
  QUESO::LibMeshNegativeLaplacianOperator precision(fobuilder, mesh);

  precision.print_info();
  precision.save_converged_evals("evals.txt");

  QUESO::InfiniteDimensionalGaussian mu(env, mean, precision, 3.0, 1.0);

  for (i = 0; i < 5; i++) {
    std::ostringstream number;
    number << i;
    mu.draw()->save_function("rand_draw" + number.str() + ".e", 0);
  }
}

  MPI_Finalize();

  return 0;
}
