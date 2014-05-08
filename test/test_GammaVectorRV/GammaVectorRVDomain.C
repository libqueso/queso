#include <queso/exceptions.h>
#include <queso/asserts.h>
#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>
#include <queso/GammaVectorRV.h>

int main(int argc, char **argv) {
  int num_samples = 1e5;

  QUESO::EnvOptionsValues *opts = new QUESO::EnvOptionsValues();
  opts->m_seed = 1;

  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", opts);

  // One-dimensional domain
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> param_space(env,
      "param_", 1, NULL);

  // Full domain
  QUESO::GslVector mins(param_space.zeroVector());
  QUESO::GslVector maxs(param_space.zeroVector());
  mins.cwSet(0);
  maxs.cwSet(INFINITY);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> param_domain("param_",
      param_space, mins, maxs);

  // Truncated domain
  QUESO::GslVector mins2(param_space.zeroVector());
  QUESO::GslVector maxs2(param_space.zeroVector());
  mins2.cwSet(0.3);
  maxs2.cwSet(INFINITY);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> param_domain2("param_",
      param_space, mins2, maxs2);

  // Compact domain
  QUESO::GslVector mins3(param_space.zeroVector());
  QUESO::GslVector maxs3(param_space.zeroVector());
  mins3.cwSet(0.5);
  maxs3.cwSet(1.0);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> param_domain3("param_",
      param_space, mins3, maxs3);

  // Gamma shape parameter
  QUESO::GslVector shape_param(param_space.zeroVector());
  shape_param.cwSet(1.0);

  // Gamma scale parameter
  QUESO::GslVector scale_param(param_space.zeroVector());
  scale_param.cwSet(0.1);

  // Full Gamma
  QUESO::GammaVectorRV<QUESO::GslVector, QUESO::GslMatrix> gamma_rv("prior_",
      param_domain, shape_param, scale_param);

  // Truncated Gamma
  QUESO::GammaVectorRV<QUESO::GslVector, QUESO::GslMatrix> gamma_rv2("prior_",
      param_domain2, shape_param, scale_param);

  // Compact Gamma
  QUESO::GammaVectorRV<QUESO::GslVector, QUESO::GslMatrix> gamma_rv3("prior_",
      param_domain3, shape_param, scale_param);

  QUESO::GslVector draw(param_space.zeroVector());

  std::cout << "Full domain is:" << param_domain << std::endl;
  std::cout << "Truncated domain is:" << param_domain2 << std::endl;
  std::cout << "Compact domain is:" << param_domain3 << std::endl;

  for (unsigned int i = 0; i < num_samples; i++) {
    // Only one-dimensional, so only testing draw[0] is ok

    // Full domain
    gamma_rv.realizer().realization(draw);
    if (draw[0] < mins[0]) {
      std::cerr << "Test failing: generated Gamma samples outside of domain.  "
                << "Domain is:"
                << std::endl
                << param_domain
                << std::endl;
      queso_error();
    }

    // Truncated domain
    gamma_rv2.realizer().realization(draw);
    if (draw[0] < mins2[0]) {
      std::cerr << "Test failing: generated Gamma samples outside of domain.  "
                << "Domain is:"
                << std::endl
                << param_domain2
                << std::endl;
      queso_error();
    }

    // Compact domain
    gamma_rv3.realizer().realization(draw);
    if (draw[0] < mins3[0] || draw[0] > maxs3[0]) {
      std::cerr << "Test failing: generated Gamma samples outside of domain.  "
                << "Domain is:"
                << std::endl
                << param_domain3
                << std::endl;
      queso_error();
    }
  }

  MPI_Finalize();

  return 0;
}
