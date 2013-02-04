#include <libmesh/libmesh.h>
#include <uqLibMeshNegativeLaplacianOperator.h>
#include <mpi.h>

using namespace libMesh;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  // Need an artificial block here because libmesh needs to
  // call PetscFinalize before we call MPI_Finalize
  {
  LibMeshInit init(argc, argv);

  uqLibMeshNegativeLaplacianOperator *C = new uqLibMeshNegativeLaplacianOperator("mesh.e");
  C->print_info();
  C->save_converged_evals("evals.txt");
  C->save_converged_evec("evec.e", 1);

  delete C;
  }

  MPI_Finalize();

  return 0;
}
