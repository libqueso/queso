#include <slope_compute.h>

int main(int argc, char** argv)
{
	// Initialize QUESO env
	MPI_Init(&argc, &argv);
	QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

	// Check point for input files
	std::ifstream chk_xs("xs.txt");
	if(chk_xs.fail()){
	     std::cout << "ERROR: Input truth data files are missing (See README). Run: gen_truth.py to generate them." << std::endl;
	     return 0;
	}
	
	std::ifstream chk_obs("obs.txt");
	if(chk_obs.fail()){
	     std::cout << "ERROR: Input truth data files are missing (See README). Run: gen_truth.py to generate them." << std::endl;
	     return 0;
	}

	std::ifstream chk_sigmas("sigmas.txt");
	if(chk_sigmas.fail()){
	     std::cout << "ERROR: Input truth data files are missing (See README). Run: gen_truth.py to generate them." << std::endl;
	     return 0;
	}
	// Call application
	infer_slope(env);

	// Finalize QUESO environment
	// delete env;
	MPI_Finalize();

	return 0;
}

