#include <sstream>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <uqLibMeshNegativeLaplacianOperator.h>
#include <mpi.h>

using namespace libMesh;

int main(int argc, char **argv)
{
  unsigned int i;

  MPI_Init(&argc, &argv);

  // Need an artificial block here because libmesh needs to
  // call PetscFinalize before we call MPI_Finalize
  {
  LibMeshInit init(argc, argv);

  Mesh mesh;
  MeshTools::Generation::build_square(mesh, 20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);
  uqLibMeshNegativeLaplacianOperator C(mesh);

  C.print_info();
  C.save_converged_evals("evals.txt");

  std::vector<double> xi(C.get_num_converged(), 0.5);
  C.inverse_kl_transform(xi)->save_function("draw.e");
  }

  MPI_Finalize();

  return 0;
}
