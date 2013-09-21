#include <sstream>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <queso/uqFunctionOperatorBuilder.h>
#include <queso/uqLibMeshNegativeLaplacianOperator.h>
#include <mpi.h>

using namespace libMesh;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  // Need an artificial block here because libmesh needs to
  // call PetscFinalize before we call MPI_Finalize
  {
  LibMeshInit init(argc, argv);

  Mesh mesh;
  MeshTools::Generation::build_square(mesh, 20, 20, -1.0, 1.0, -1.0, 1.0, QUAD4);

  QUESO::uqFunctionOperatorBuilder fobuilder;
  fobuilder.num_req_eigenpairs = 5;

  QUESO::uqLibMeshNegativeLaplacianOperator C(fobuilder, mesh);
  C.print_info();
  C.save_converged_evals("evals.txt");

  for (unsigned int i = 0; i < 5; i++) {
    std::ostringstream base_name("evec");
    std::ostringstream number;
    number << i;
    C.save_converged_evec("evec" + number.str() + ".e", i);
  }

  }

  MPI_Finalize();

  return 0;
}
