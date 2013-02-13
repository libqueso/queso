#include <memory>
#include <sstream>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <uqEnvironment.h>
#include <uqLibMeshFunction.h>
#include <uqLibMeshNegativeLaplacianOperator.h>
#include <uqInfiniteDimensionalGaussian.h>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

using namespace libMesh;

int main(int argc, char **argv)
{
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
  MeshTools::Generation::build_square(mesh, 20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);
  uqLibMeshFunction mean(mesh);
  uqLibMeshNegativeLaplacianOperator precision(mesh);

  precision.print_info();
  precision.save_converged_evals("evals.txt");

  std::vector<double> xi(precision.get_num_converged(), 0.5);
  uqInfiniteDimensionalGaussian mu(env, mean, precision);
  mu.draw()->save_function("rand_draw.e");
  }

  MPI_Finalize();

  return 0;
}
