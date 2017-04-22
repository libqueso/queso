//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/BoxSubset.h>
#include <queso/GaussianLikelihoodDiagonalCovariance.h>

#include <cstdlib>

#define TOL 1e-8

template<class V, class M>
class Likelihood : public QUESO::GaussianLikelihoodDiagonalCovariance<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain,
      const V & observations, const V & covariance)
    : QUESO::GaussianLikelihoodDiagonalCovariance<V, M>(prefix, domain,
        observations, covariance)
  {
  }

  virtual ~Likelihood()
  {
  }

  virtual void evaluateModel(const V & domainVector, V & modelOutput) const
  {
    // Model is a map from R to R^2

    // Evaluate model and fill up the m_modelOutput member variable
    for (unsigned int i = 0; i < modelOutput.sizeLocal(); i++) {
      modelOutput[i] = domainVector[0] + 3.0;
    }
  }
};

int main(int argc, char ** argv) {
  std::string inputFileName = "test_gaussian_likelihoods/queso_input.txt";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);
#else
  QUESO::FullEnvironment env(inputFileName, "", NULL);
#endif

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 1, NULL);

  double min_val = -INFINITY;
  double max_val = INFINITY;

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(min_val);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  // Set up observation space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> obsSpace(env,
      "obs_", 2, NULL);

  // Fill up observation vector
  QUESO::GslVector observations(obsSpace.zeroVector());
  observations[0] = 1.0;
  observations[1] = 1.0;

  // Fill up covariance 'matrix'
  QUESO::GslVector covariance(obsSpace.zeroVector());
  covariance[0] = 1.0;
  covariance[1] = 2.0;

  // Pass in observations to Gaussian likelihood object
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood("llhd_", paramDomain,
      observations, covariance);

  double lhood_value;
  double truth_value;
  QUESO::GslVector point(paramSpace.zeroVector());
  point[0] = 0.0;
  lhood_value = lhood.actualValue(point, NULL, NULL, NULL, NULL);
  truth_value = std::exp(-3.0);

  if (std::abs(lhood_value - truth_value) > TOL) {
    std::cerr << "Scalar Gaussian test case failure." << std::endl;
    std::cerr << "Computed likelihood value is: " << lhood_value << std::endl;
    std::cerr << "Likelihood value should be: " << truth_value << std::endl;
    queso_error();
  }

  point[0] = -2.0;
  lhood_value = lhood.actualValue(point, NULL, NULL, NULL, NULL);
  truth_value = 1.0;

  if (std::abs(lhood_value - truth_value) > TOL) {
    std::cerr << "Scalar Gaussian test case failure." << std::endl;
    std::cerr << "Computed likelihood value is: " << lhood_value << std::endl;
    std::cerr << "Likelihood value should be: " << truth_value << std::endl;
    queso_error();
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
