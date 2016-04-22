#include <queso/Environment.h>
#include <queso/EnvironmentOptions.h>
#include <queso/SharedPtr.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/SequenceGenerator.h>
#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalInverseProblemOptions.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/VectorSequence.h>

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
    double x1 = domainVector[0];
    double x2 = domainVector[1];

    return -0.5 * (x1 * x1 + x2 * x2);
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }
};

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class MySequenceGenerator : public QUESO::SequenceGenerator<V, M>
{
public:
  MySequenceGenerator(const char * prefix,
                       const QUESO::MhOptionsValues * options,
                       const QUESO::BaseVectorRV<V, M> & sourceRv,
                       const V & initialPosition,
                       const M * inputProposalCovMatrix)
    : QUESO::SequenceGenerator<V, M>(prefix, options, sourceRv, initialPosition, inputProposalCovMatrix)
  {
  }

  virtual ~MySequenceGenerator()
  {
  }

  virtual void generateSequence(
      QUESO::BaseVectorSequence<V, M> & workingChain,
      QUESO::ScalarSequence<double> * workingLogLikelihoodValues,
      QUESO::ScalarSequence<double> * workingLogTargetValues)
  {
    unsigned int chainSize = 10;

    QUESO::GslVector v(workingChain.vectorSpace().zeroVector());
    v.cwSet(1.0);

    workingChain.resizeSequence(chainSize);
    workingLogLikelihoodValues->resizeSequence(chainSize);
    workingLogTargetValues->resizeSequence(chainSize);

    for (unsigned int i = 0; i < chainSize; i++) {
      workingChain.setPositionValues(i, v);
    }

    for (unsigned int i = 0; i < chainSize; i++) {
      (*workingLogLikelihoodValues)[i] = 1.0;
    }

    for (unsigned int i = 0; i < chainSize; i++) {
      (*workingLogTargetValues)[i] = 1.0;
    }
  }

  virtual void propose(unsigned int positionId,
      const QUESO::BaseVectorSequence<V, M> & workingChain,
      V & proposedState)
  {
    proposedState.cwSet(1.0);
  }

};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  QUESO::EnvOptionsValues envOptions;
  envOptions.m_numSubEnvironments = 1;
  envOptions.m_subDisplayFileName = "output_test_setSequenceGenerator/display";
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

  double min_val = -10.0;
  double max_val = 10.0;
  paramMins.cwSet(min_val);
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

  Likelihood<> lhood("llhd_", paramDomain);

  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  QUESO::SipOptionsValues sipOptions;
  sipOptions.m_computeSolution = 1;
  sipOptions.m_dataOutputFileName = "output_test_setSequenceGenerator/sipOutput";
  sipOptions.m_dataOutputAllowedSet.clear();
  sipOptions.m_dataOutputAllowedSet.insert(0);

  QUESO::StatisticalInverseProblem<> ip("", &sipOptions, priorRv, lhood,
      postRv);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.0;
  paramInitials[1] = 0.0;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  proposalCovMatrix(0, 0) = 1.0;
  proposalCovMatrix(0, 1) = 0.0;
  proposalCovMatrix(1, 0) = 0.0;
  proposalCovMatrix(1, 1) = 1.0;

  QUESO::MhOptionsValues mhOptions;
  mhOptions.m_dataOutputFileName = "output_test_setSequenceGenerator/sipOutput";
  mhOptions.m_dataOutputAllowAll = 1;
  mhOptions.m_rawChainGenerateExtra = 0;
  mhOptions.m_rawChainDisplayPeriod = 50000;
  mhOptions.m_rawChainMeasureRunTimes = 1;
  mhOptions.m_rawChainDataOutputFileName = "output_test_setSequenceGenerator/ip_raw_chain";
  mhOptions.m_rawChainDataOutputAllowAll = 1;
  mhOptions.m_displayCandidates = 0;
  mhOptions.m_tkUseLocalHessian = 0;
  mhOptions.m_tkUseNewtonComponent = 1;
  mhOptions.m_filteredChainGenerate = 0;
  mhOptions.m_rawChainSize = 1000;
  mhOptions.m_putOutOfBoundsInChain = false;
  mhOptions.m_drMaxNumExtraStages = 1;
  mhOptions.m_drScalesForExtraStages.resize(1);
  mhOptions.m_drScalesForExtraStages[0] = 5.0;
  mhOptions.m_amInitialNonAdaptInterval = 100;
  mhOptions.m_amAdaptInterval = 100;
  mhOptions.m_amEta = (double) 2.4 * 2.4 / dim;  // From Gelman 95
  mhOptions.m_amEpsilon = 1.e-8;
  mhOptions.m_doLogitTransform = false;

  QUESO::SharedPtr<MySequenceGenerator<> >::Type
    sequenceGenerator(new MySequenceGenerator<>("",
                                                NULL,
                                                ip.postRv(),
                                                paramInitials,
                                                &proposalCovMatrix));

  ip.setSequenceGenerator(sequenceGenerator);
  ip.solveWithBayesMetropolisHastings(&mhOptions, paramInitials,
      &proposalCovMatrix);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
