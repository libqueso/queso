#include <vector>

#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/SharedPtr.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/RandomWalkSG.h>
#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalInverseProblemOptions.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/VectorSequence.h>
#include <queso/TKGroup.h>

#define MEAN 10.0
#define STDDEV 0.5
#define DISPLAY_FREQ 100000

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain)
    : QUESO::BaseScalarFunction<V, M>(prefix, domain)
  {
  }

  virtual ~Likelihood()
  {
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    double x1 = domainVector[0] - MEAN;
    double x2 = domainVector[1] - MEAN;

    return -0.5 * (x1 * x1 + x2 * x2) / (STDDEV * STDDEV);
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }
};

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class InformedRandomWalkSG : public QUESO::RandomWalkSG<V, M>
{
public:
  InformedRandomWalkSG(const char * prefix,
                       const QUESO::MhOptionsValues * options,
                       const QUESO::BaseVectorRV<V, M> & sourceRv,
                       const V & initialPosition,
                       const M * inputProposalCovMatrix)
    : QUESO::RandomWalkSG<V, M>(prefix, options, sourceRv,
                                initialPosition, inputProposalCovMatrix)
  {
  }

  virtual ~InformedRandomWalkSG()
  {
  }

  virtual void updateProposalCovariance(unsigned int positionId,
      const QUESO::BaseVectorSequence<V, M> & workingChain,
      M & covarianceMatrix)
  {
    const double var = STDDEV * STDDEV;

    for (unsigned int i = 0; i < this->m_vectorSpace.dimLocal(); i++) {
      covarianceMatrix(i, i) = var + (100.0 / (double) (positionId + 1));
    }

    if (positionId % DISPLAY_FREQ == 0) {
      std::cout << "iter "
                << positionId
                << " updating prop cov matrix to: "
                << covarianceMatrix
                << std::endl;
    }
  }
};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues envOptions;
  envOptions.m_numSubEnvironments = 1;
  envOptions.m_subDisplayFileName = "output_test_randomWalkTK/display";
  envOptions.m_subDisplayAllowAll = 1;
  envOptions.m_displayVerbosity = 2;
  envOptions.m_syncVerbosity = 0;
  envOptions.m_seed = 0;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &envOptions);
#else
  QUESO::FullEnvironment env("", "", &envOptions);
#endif

  unsigned int dim = 2;
  QUESO::VectorSpace<> paramSpace(env, "param_", dim, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  double min_val = -INFINITY;
  double max_val = INFINITY;
  paramMins.cwSet(min_val);
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  Likelihood<> lhood("llhd_", paramDomain);

  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  QUESO::SipOptionsValues sipOptions;
  sipOptions.m_computeSolution = 1;
  sipOptions.m_dataOutputFileName = "output_test_randomWalkTK/sipOutput";
  sipOptions.m_dataOutputAllowedSet.clear();
  sipOptions.m_dataOutputAllowedSet.insert(0);

  QUESO::StatisticalInverseProblem<> ip("", &sipOptions, priorRv, lhood,
      postRv);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 10.0;
  paramInitials[1] = 10.0;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  proposalCovMatrix(0, 0) = 1.0;
  proposalCovMatrix(0, 1) = 0.0;
  proposalCovMatrix(1, 0) = 0.0;
  proposalCovMatrix(1, 1) = 1.0;

  QUESO::MhOptionsValues mhOptions;
  mhOptions.m_dataOutputFileName = "output_test_randomWalkTK/sipOutput";
  mhOptions.m_dataOutputAllowAll = 1;
  mhOptions.m_rawChainGenerateExtra = 0;
  mhOptions.m_rawChainDisplayPeriod = 50000;
  mhOptions.m_rawChainMeasureRunTimes = 1;
  mhOptions.m_rawChainDataOutputFileName = "output_test_randomWalkTK/ip_raw_chain";
  mhOptions.m_rawChainDataOutputFileType = "txt";
  mhOptions.m_rawChainDataOutputAllowAll = 1;
  mhOptions.m_displayCandidates = 0;
  mhOptions.m_tkUseLocalHessian = 0;
  mhOptions.m_tkUseNewtonComponent = 1;
  mhOptions.m_filteredChainGenerate = 0;
  mhOptions.m_rawChainSize = 1000000;
  mhOptions.m_putOutOfBoundsInChain = false;
  mhOptions.m_drMaxNumExtraStages = 1;
  mhOptions.m_drScalesForExtraStages.resize(1);
  mhOptions.m_drScalesForExtraStages[0] = 5.0;
  mhOptions.m_amInitialNonAdaptInterval = 100;
  mhOptions.m_amAdaptInterval = 100;
  mhOptions.m_amEta = (double) 2.4 * 2.4 / dim;  // From Gelman 95
  mhOptions.m_amEpsilon = 1.e-8;
  mhOptions.m_doLogitTransform = false;

  QUESO::SharedPtr<InformedRandomWalkSG<> >::Type
    sequenceGenerator(new InformedRandomWalkSG<>("",
                                                &mhOptions,
                                                ip.postRv(),
                                                paramInitials,
                                                &proposalCovMatrix));

  ip.setSequenceGenerator(sequenceGenerator);
  ip.solveWithBayesMetropolisHastings(&mhOptions, paramInitials,
      &proposalCovMatrix);

  QUESO::GslVector postSample(paramSpace.zeroVector());
  QUESO::GslVector mean(paramSpace.zeroVector());
  QUESO::GslVector delta(paramSpace.zeroVector());
  for (unsigned int i = 1; i < mhOptions.m_rawChainSize + 1; i++) {
    ip.postRv().realizer().realization(postSample);
    for (unsigned int j = 0; j < dim; j++) {
      delta[j] = postSample[j] - mean[j];
      mean[j] += (double) delta[j] / i;
    }
  }

  double mean_lowerbound = MEAN - (3.0 * STDDEV / std::sqrt(mhOptions.m_rawChainSize));
  double mean_upperbound = MEAN + (3.0 * STDDEV / std::sqrt(mhOptions.m_rawChainSize));

  int return_val = 0;
  for (unsigned int j = 0; j < dim; j++) {
    if (mean[0] < mean_lowerbound || mean[0] > mean_upperbound) {
      return_val = 1;
      break;
    }
  }

  std::cout << "Mean lower bound: " << mean_lowerbound << std::endl;
  std::cout << "Mean: " << mean[0] << std::endl;
  std::cout << "Mean upper bound: " << mean_upperbound << std::endl;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return return_val;
}
