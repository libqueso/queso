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
#include <queso/GslBlockMatrix.h>
#include <queso/VectorSet.h>
#include <queso/BoxSubset.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/ScalarFunction.h>
#include <queso/GaussianLikelihoodBlockDiagonalCovariance.h>
#include <queso/StatisticalInverseProblem.h>

#include <cstdlib>

#define TOL 1e-15

// A Gaussian likelihood
template<class V, class M>
class Likelihood : public QUESO::GaussianLikelihoodBlockDiagonalCovariance<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain,
      const V & observations, const QUESO::GslBlockMatrix & covariance)
    : QUESO::GaussianLikelihoodBlockDiagonalCovariance<V, M>(prefix, domain,
        observations, covariance)
  {
  }

  virtual ~Likelihood()
  {
  }

  virtual void evaluateModel(const V & domainVector, V & modelOutput) const
  {
    for (unsigned int i = 0; i < modelOutput.sizeLocal(); i++) {
      modelOutput[i] = domainVector[i];
    }
  }
};

// A custom likelihood (that just so happens to also be Gaussian)
template<class V, class M>
class CustomLikelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  CustomLikelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain)
    : QUESO::BaseScalarFunction<V, M>(prefix, domain)
  {
  }

  virtual ~CustomLikelihood()
  {
  }

  virtual double lnValue(const V & domainVector) const
  {
    double d1 = domainVector[0] - 1.0;
    double d2 = domainVector[1] - 1.0;
    double d3 = domainVector[2] - 1.0;

    // Sigma inverse multiplied by d
    double r1 = 2.0 * d1 - 0.5 * d2;
    double r2 = -0.5 * d1 + 0.25 * d2;
    double r3 = 0.5 * d3;

    r1 *= d1;
    r2 *= d2;
    r3 *= d3;

    // The square of the 2-norm of the misfit
    double misfit = r1 + r2 + r3;
    misfit /= -2.0;
    return misfit;
  }

  virtual double actualValue(const V & domainVector, const V * /*domainDirection*/,
      V * /*gradVector*/, M * /*hessianMatrix*/, V * /*hessianEffect*/) const
  {
    return this->lnValue(domainVector);
  }
};

template<class V, class M>
class BayesianInverseProblem
{
public:
  BayesianInverseProblem(unsigned int likelihoodFlag)
  {
    std::string inputFileName =
      "test_gaussian_likelihoods/gaussian_consistency_input.txt";
    const char * test_srcdir = std::getenv("srcdir");
    if (test_srcdir)
      inputFileName = test_srcdir + ('/' + inputFileName);

#ifdef QUESO_HAS_MPI
    this->env = new QUESO::FullEnvironment
      (MPI_COMM_WORLD, inputFileName, "", NULL);
#else
    this->env = new QUESO::FullEnvironment
      (inputFileName, "", NULL);
#endif

    this->paramSpace =
      new QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>(*env,
          "param_", 3, NULL);

    double min_val = -INFINITY;
    double max_val = INFINITY;

    this->paramMins = new QUESO::GslVector(this->paramSpace->zeroVector());
    this->paramMins->cwSet(min_val);
    this->paramMaxs = new QUESO::GslVector(this->paramSpace->zeroVector());
    this->paramMaxs->cwSet(max_val);

    this->paramDomain =
      new QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>("param_",
          *(this->paramSpace), *(this->paramMins), *(this->paramMaxs));

    this->priorRv =
      new QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix>("prior_",
          *(this->paramDomain));

    // Set up observation space
    this->obsSpace =
      new QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>(*env, "obs_",
          3, NULL);

    // Fill up observation vector
    this->observations = new QUESO::GslVector(this->obsSpace->zeroVector());
    (*(this->observations))[0] = 1.0;
    (*(this->observations))[1] = 1.0;
    (*(this->observations))[2] = 1.0;

    // Fill up covariance 'matrix'
    std::vector<unsigned int> blockSizes(2);
    blockSizes[0] = 2;  // First block is 2x2
    blockSizes[1] = 1;  // Second block is 1x1 (scalar)

    this->covariance = new QUESO::GslBlockMatrix(blockSizes,
        this->obsSpace->zeroVector(), 1.0);
    (this->covariance->getBlock(0))(0, 0) = 1.0;
    (this->covariance->getBlock(0))(0, 1) = 2.0;
    (this->covariance->getBlock(0))(1, 0) = 2.0;
    (this->covariance->getBlock(0))(1, 1) = 8.0;
    (this->covariance->getBlock(1))(0, 0) = 2.0;

    // Construct whatever likelihood we need

    if (likelihoodFlag == 1) {
      // Do canned Gaussian likelihood
      this->lhood = new Likelihood<QUESO::GslVector, QUESO::GslMatrix>("llhd_",
          *(this->paramDomain), *(this->observations), *(this->covariance));
    }
    else if (likelihoodFlag == 2) {
      // Do the other, custom, likelihood (which just so happens to be
      // Gaussian)
      this->lhood = new CustomLikelihood<QUESO::GslVector, QUESO::GslMatrix>(
          "llhd_", *(this->paramDomain));
    }
    else {
      // Wat
      queso_error_msg("DIE!");
    }

    this->postRv =
      new QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix>("post_",
        *(this->paramSpace));

    this->ip =
      new QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix>(
          "", NULL, *(this->priorRv), *(this->lhood), *(this->postRv));

    this->paramInitials = new QUESO::GslVector(this->paramSpace->zeroVector());
    (*(this->paramInitials))[0] = 1.0;
    (*(this->paramInitials))[1] = 1.0;

    this->proposalCovMatrix = new QUESO::GslMatrix(
        this->paramSpace->zeroVector());
    (*(this->proposalCovMatrix))(0, 0) = 1.0;
    (*(this->proposalCovMatrix))(0, 1) = 2.0;
    (*(this->proposalCovMatrix))(0, 2) = 0.0;
    (*(this->proposalCovMatrix))(1, 0) = 2.0;
    (*(this->proposalCovMatrix))(1, 1) = 8.0;
    (*(this->proposalCovMatrix))(1, 2) = 0.0;
    (*(this->proposalCovMatrix))(2, 0) = 0.0;
    (*(this->proposalCovMatrix))(2, 1) = 0.0;
    (*(this->proposalCovMatrix))(2, 2) = 2.0;

    this->ip->solveWithBayesMetropolisHastings(NULL, *(this->paramInitials),
        this->proposalCovMatrix);

    this->draw = new QUESO::GslVector(this->paramSpace->zeroVector());
  }

  virtual ~BayesianInverseProblem()
  {
  }

  QUESO::FullEnvironment * env;
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> * paramSpace;
  QUESO::GslVector * paramMins;
  QUESO::GslVector * paramMaxs;
  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> * paramDomain;
  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> * priorRv;
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> * obsSpace;
  QUESO::GslVector * observations;
  QUESO::GslBlockMatrix * covariance;
  QUESO::BaseScalarFunction<QUESO::GslVector, QUESO::GslMatrix> * lhood;
  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> * ip;
  QUESO::GslVector * paramInitials;
  QUESO::GslMatrix * proposalCovMatrix;
  QUESO::GslVector * draw;
  QUESO::GenericVectorRV<V, M> * postRv;
};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  // Instantiate each inverse problem
  BayesianInverseProblem<QUESO::GslVector, QUESO::GslMatrix> b1(1);
  BayesianInverseProblem<QUESO::GslVector, QUESO::GslMatrix> b2(2);

  // Compare each draw
  for (unsigned int i = 0; i < 1000; i++) {
    b1.postRv->realizer().realization(*(b1.draw));
    b2.postRv->realizer().realization(*(b2.draw));

    // Test draws
    QUESO::GslVector diff(*(b1.draw));
    diff -= *(b2.draw);
    if (diff.norm2() > TOL) {
      queso_error_msg("test_fullCovarianceChain: Samples from the canned likelihood do not match samples from the equivalent custom likelihood");
    }
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
