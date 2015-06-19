#include <iostream>
#include <cmath>
#include <string>

// Boost includes
#include <boost/math/constants/constants.hpp>

// LibMesh includes
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/vector_value.h>
#include <libmesh/exact_solution.h>
#include <libmesh/exact_error_estimator.h>
#include <libmesh/explicit_system.h>

// QUESO includes
#include <queso/Defines.h>
#include <queso/FunctionOperatorBuilder.h>
#include <queso/LibMeshFunction.h>
#include <queso/LibMeshNegativeLaplacianOperator.h>
#include <queso/InfiniteDimensionalGaussian.h>
#include <queso/InfiniteDimensionalLikelihoodBase.h>
#include <queso/InfiniteDimensionalMCMCSamplerOptions.h>
#include <queso/InfiniteDimensionalMCMCSampler.h>

// GSL includes
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// The likelihood object subclass
class Likelihood : public QUESO::InfiniteDimensionalLikelihoodBase
{
public:
  Likelihood(double obs_stddev);
  ~Likelihood();
  virtual double evaluate(QUESO::FunctionBase & state);
  gsl_rng *r;
};

Likelihood::Likelihood(double obs_stddev)
  : QUESO::InfiniteDimensionalLikelihoodBase(obs_stddev)
{
  this->r = gsl_rng_alloc(gsl_rng_default);
}

Likelihood::~Likelihood()
{
  gsl_rng_free(this->r);
}

double Likelihood::evaluate(QUESO::FunctionBase & flow)
{
  // This is where you call your forward problem, make observations and compute
  // the likelihood norm.  Here we'll just pretend we observe a zero with
  // error.
  const double obs_stddev = this->obs_stddev();
  const double obs = gsl_ran_gaussian(this->r, obs_stddev);
  return obs * obs / (2.0 * obs_stddev * obs_stddev);
}

int main(int argc, char **argv) {
  unsigned int num_pts = 100;
  unsigned int num_pairs = num_pts / 2.0;
  const double alpha = 3.0;
  const double beta = 10.0;
  const double obs_stddev = 1e-3;

  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

#ifndef QUESO_HAVE_HDF5
  std::cerr << "Cannot run infinite dimensional inverse problems\n" <<
               "without HDF5 enabled." << std::endl;
  return 0;
#else

  libMesh::LibMeshInit init(argc, argv);

  // Generate the mesh on which the samples live.
  libMesh::Mesh mesh(init.comm());
  libMesh::MeshTools::Generation::build_line(mesh, num_pts, 0.0, 1.0,
                                             libMeshEnums::EDGE2);

  // Use a helper object to define some of the properties of our samples
  QUESO::FunctionOperatorBuilder fobuilder;
  fobuilder.order = "FIRST";
  fobuilder.family = "LAGRANGE";
  fobuilder.num_req_eigenpairs = num_pairs;

  // Define the mean of the prior
  QUESO::LibMeshFunction mean(fobuilder, mesh);

  // Define the precision operator of the prior
  QUESO::LibMeshNegativeLaplacianOperator precision(fobuilder, mesh);

  // Define the prior measure
  QUESO::InfiniteDimensionalGaussian mu(env, mean, precision, alpha, beta);

  // Create likelihood object
  Likelihood llhd(obs_stddev);

  // The worst hack in the world ever
  // boost::shared_ptr<LibMeshFunction> lm_draw(boost::static_pointer_cast<LibMeshFunction>(mu.draw()));
  // std::cout << "before zero" << std::endl;
  // lm_draw->equation_systems->get_system<ExplicitSystem>("Function").solution->zero();
  // std::cout << "after zero" << std::endl;

  // Create the options helper object that determines what options to pass to
  // the sampler
  QUESO::InfiniteDimensionalMCMCSamplerOptions opts(env, "");

  // Set the number of iterations to do
  opts.m_num_iters = 1000;

  // Set the frequency with which we save samples
  opts.m_save_freq = 10;

  // Set the RWMH step size
  opts.m_rwmh_step = 0.1;

  // Construct the sampler, and set the name of the output file (will only
  // write HDF5 files)
  QUESO::InfiniteDimensionalMCMCSampler s(env, mu, llhd, &opts);

  for (unsigned int i = 0; i < opts.m_num_iters; i++) {
    s.step();
    if (i % 100 == 0) {
      std::cout << "sampler iteration: " << i << std::endl;
      std::cout << "avg acc prob is: " << s.avg_acc_prob() << std::endl;
      std::cout << "l2 norm is: " << s.llhd_val() << std::endl;
    }
  }

#endif // QUESO_HAVE_HDF5

  return 0;
}
