#include <slope_compute.h>

int main(int argc, char** argv)
{
	// Initialize QUESO env
	MPI_Init(&argc, &argv);
	QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);
	
	// Call application
	infer_slope(env);

	// Finalize QUESO environment
	// delete env;
	MPI_Finalize();

	return 0;
}

