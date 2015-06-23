//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/GaussianLikelihoodFullCovarianceRandomCoefficient.h>
#include <queso/StatisticalInverseProblem.h>

template<class V, class M>
class Likelihood : public QUESO::GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain,
      const V & observations, const M & covariance)
    : QUESO::GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>(prefix,
        domain, observations, covariance)
  {
    // Default covariance coefficient is 1.0
  }

  virtual ~Likelihood()
  {
  }

  virtual void evaluateModel(const V & domainVector, const V * domainDirection,
      V & modelOutput, V * gradVector, M * hessianMatrix,
      V * hessianEffect) const
  {
    // domainVector is of size 2 because the first element is the model
    // parameter and the second (last) element is the multiplicative
    // coefficient of the observational error covariance matrix.  Therefore,
    // for calling the model code, we need only concern ourselves with the
    // first element of domainVector.

    // Evaluate model and fill up the m_modelOutput member variable
    for (unsigned int i = 0; i < modelOutput.sizeLocal(); i++) {
      modelOutput[i] = 1.0;
    }
  }
};

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 2, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  // Model parameter between 0 and 1
  paramMins[0] = 0.0;
  paramMaxs[0] = 1.0;

  // Hyperparameter (multiplicative coefficient of observational error
  // covariance matrix) between 0.0 and \infty
  paramMins[1] = 0.0;
  paramMaxs[1] = INFINITY;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorRv("prior_",
      paramDomain);

  // Set up observation space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> obsSpace(env,
      "obs_", 2, NULL);

  // Fill up observation vector
  QUESO::GslVector observations(obsSpace.zeroVector());
  observations[0] = 1.0;
  observations[1] = 1.0;

  // Fill up covariance 'matrix'
  QUESO::GslMatrix covariance(obsSpace.zeroVector());
  covariance(0, 0) = 1.0;
  covariance(0, 1) = 0.0;
  covariance(1, 0) = 0.0;
  covariance(1, 1) = 1.0;

  // Pass in observations to Gaussian likelihood object
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood("llhd_", paramDomain,
      observations, covariance);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix>
    postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix>
    ip("", NULL, priorRv, lhood, postRv);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  paramInitials[0] = 0.0;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  for (unsigned int i = 0; i < 1; i++) {
    proposalCovMatrix(i, i) = 0.1;
  }

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  MPI_Finalize();

  return 0;
}
