#include <sstream>
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

  uqLibMeshNegativeLaplacianOperator *C = new uqLibMeshNegativeLaplacianOperator();
  C->print_info();
  C->save_converged_evals("evals.txt");


  for (unsigned int i = 0; i < 5; i++) {
    std::ostringstream base_name("evec");
    std::ostringstream number;
    number << i;
    C->save_converged_evec("evec" + number.str() + ".e", i);
  }

  delete C;
  }

  MPI_Finalize();

  return 0;
}
