//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <queso/GaussianProcessEmulator.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class V, class M>
GaussianProcessEmulator<V, M>::GaussianProcessEmulator(
    const VectorSet<V, M> & domain,
    const VectorSpace<V, M> & m_scenarioSpace,
    const VectorSpace<V, M> & m_parameterSpace,
    const VectorSpace<V, M> & m_simulationOutputSpace,
    const VectorSpace<V, M> & m_experimentOutputSpace,
    const unsigned int m_numSimulations,
    const unsigned int m_numExperiments,
    const std::vector<V *> & m_simulationScenarios,
    const std::vector<V *> & m_simulationParameters,
    const std::vector<V *> & m_simulationOutputs,
    const std::vector<V *> & m_experimentScenarios,
    const std::vector<V *> & m_experimentOutputs,
    const M & m_experimentErrors,
    const ConcatenatedVectorRV<V, M> & m_totalPrior)
  :
  BaseScalarFunction<V, M>("", m_totalPrior.imageSet()),
  m_scenarioSpace(m_scenarioSpace),
  m_parameterSpace(m_parameterSpace),
  m_simulationOutputSpace(m_simulationOutputSpace),
  m_experimentOutputSpace(m_experimentOutputSpace),
  m_numSimulations(m_numSimulations),
  m_numExperiments(m_numExperiments),
  m_simulationScenarios(m_simulationScenarios),
  m_simulationParameters(m_simulationParameters),
  m_simulationOutputs(m_simulationOutputs),
  m_experimentScenarios(m_experimentScenarios),
  m_experimentOutputs(m_experimentOutputs),
  m_experimentErrors(m_experimentErrors),
  m_totalPrior(m_totalPrior)
{
}

template <class V, class M>
GaussianProcessEmulator<V, M>::~GaussianProcessEmulator()
{
  // Do nothing?
}

template <class V, class M>
double
GaussianProcessEmulator<V, M>::lnValue(const V & domainVector,
                                       const V * domainDirection,
                                       V * gradVector,
                                       M * hessianMatrix,
                                       V * hessianEffect) const
{
  // Components of domainVector:
  // theta(1)
  // theta(2)
  // ...
  // theta(dimParameterSpace)
  // emulator_mean
  // emulator_precision
  // emulator_corr_strength(1)
  // ...
  // emulator_corr_strength(dimScenario + dimParameter)
  // discrepancy_precision
  // discrepancy_corr_strength(1)
  // ...
  // discrepancy_corr_strength(dimScenario)

  // Construct covariance matrix
  unsigned int totalDim = this->m_numExperiments + this->m_numSimulations;
  double prodScenario = 1.0;
  double prodParameter = 1.0;
  double prodDiscrepancy = 1.0;
  unsigned int dimScenario = (this->m_scenarioSpace).dimLocal();
  unsigned int dimParameter = (this->m_parameterSpace).dimLocal();

  // This is cumbersome.  All I want is a matrix.
  VectorSpace<V, M> gpSpace(this->m_scenarioSpace.env(), "", totalDim, NULL);
  V residual(gpSpace.zeroVector());
  M covMatrix(residual);

  V *scenario1;
  V *scenario2;
  V *parameter1;
  V *parameter2;

  // This for loop is a disaster and could do with a *lot* of optimisation
  for (unsigned int i = 0; i < totalDim; i++) {
    for (unsigned int j = 0; j < totalDim; j++) {
      // Get i-th and j-th simulation and parameter
      if (i < this->m_numExperiments) {
        scenario1 = new V(*((this->m_experimentScenarios)[i]));
        parameter1 = new V(*((this->m_simulationParameters)[0]));
        for (unsigned int k = 0; k < dimParameter; k++) {
          (*parameter1)[k] = domainVector[k];
        }
      }
      else {
        scenario1 =
          new V(*((this->m_simulationScenarios)[i-this->m_numExperiments]));
        parameter1 =
          new V(*((this->m_simulationParameters)[i-this->m_numExperiments]));
      }

      if (j < this->m_numExperiments) {
        scenario2 = new V(*((this->m_experimentScenarios)[j]));
        parameter2 = new V(*((this->m_simulationParameters)[0]));
        for (unsigned int k = 0; k < dimParameter; k++) {
          (*parameter2)[k] = domainVector[k];
        }
      }
      else {
        scenario2 =
          new V(*((this->m_simulationScenarios)[j-this->m_numExperiments]));
        parameter2 =
          new V(*((this->m_simulationParameters)[j-this->m_numExperiments]));
      }

      // Emulator component
      prodScenario = 1.0;
      prodParameter = 1.0;
      unsigned int emulatorCorrStrStart = dimParameter + 2;
      for (unsigned int k = 0; k < dimScenario; k++) {
        prodScenario *= std::pow(domainVector[emulatorCorrStrStart+k],
                                 4.0 * ((*scenario1)[k] - (*scenario2)[k]) *
                                       ((*scenario1)[k] - (*scenario2)[k]));
      }

      for (unsigned int k = 0; k < dimParameter; k++) {
        prodParameter *= std::pow(
            domainVector[emulatorCorrStrStart+dimScenario+k],
            4.0 * ((*parameter1)[k] - (*parameter2)[k]) *
                  ((*parameter1)[k] - (*parameter2)[k]));
      }

      covMatrix(i, j) = prodScenario * prodParameter /
                        domainVector[dimParameter+1];  // emulator precision

      delete scenario1;
      delete scenario2;
      delete parameter1;
      delete parameter2;

      // If we're in the experiment cross correlation part, need extra foo
      if (i < this->m_numExperiments && j < this->m_numExperiments) {
        scenario1 = new V(*((this->m_simulationScenarios)[i]));
        scenario2 = new V(*((this->m_simulationScenarios)[j]));
        prodDiscrepancy = 1.0;
        unsigned int discrepancyCorrStrStart = dimParameter +
                                               dimParameter +
                                               dimScenario + 3;
        for (unsigned int k = 0; k < dimScenario; k++) {
          prodDiscrepancy *= std::pow(domainVector[discrepancyCorrStrStart+k],
                                      4.0 * ((*scenario1)[k] - (*scenario2)[k]) *
                                            ((*scenario1)[k] - (*scenario2)[k]));
        }

        covMatrix(i, j) += prodDiscrepancy /
                           domainVector[discrepancyCorrStrStart-1];
        covMatrix(i, j) += (this->m_experimentErrors)(i, j);

        delete scenario1;
        delete scenario2;
      }
    }

    // Add small component to diagonal to make stuff +ve def
    covMatrix(i, i) += 0.0001;
  }

  // Form residual = D - mean
  for (unsigned int i = 0; i < this->m_numExperiments; i++) {
    // Scalar so ok -- will need updating for nonscalar case
    residual[i] = (*((this->m_experimentOutputs)[i]))[0] -
      domainVector[dimParameter];
  }
  for (unsigned int i = 0; i < this->m_numSimulations; i++) {
    // Scalar so ok -- will need updating for nonscalar case
    residual[i+this->m_numExperiments] =
      (*((this->m_simulationOutputs)[i]))[0] - domainVector[dimParameter];
  }

  // Solve covMatrix * sol = residual
  V sol(covMatrix.invertMultiply(residual));

  // There's no dot product function in GslVector.
  double minus_2_log_lhd = 0.0;
  for (unsigned int i = 0; i < totalDim; i++) {
    minus_2_log_lhd += sol[i] * residual[i];
  }

  if (minus_2_log_lhd < 0) {
    std::cout << " oh noes!  ln value negative: " << minus_2_log_lhd << std::endl;
  }

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  return -0.5 * minus_2_log_lhd;
#else
  return minus_2_log_lhd;
#endif
}

template <class V, class M>
double
GaussianProcessEmulator<V, M>::actualValue(const V & domainVector,
                                           const V * domainDirection,
                                           V * gradVector,
                                           M * hessianMatrix,
                                           V * hessianEffect) const
{
  // Do nothing?
  return 1.0;
}

template <class V, class M>
GaussianProcessFactory<V, M>::GaussianProcessFactory(
    const char * prefix,
    const BaseVectorRV<V, M> & parameterPrior,
    const VectorSpace<V, M> & scenarioSpace,
    const VectorSpace<V, M> & parameterSpace,
    const VectorSpace<V, M> & simulationOutputSpace,
    const VectorSpace<V, M> & experimentOutputSpace,
    unsigned int numSimulations,
    unsigned int numExperiments)
  :
    m_prefix(prefix),
    m_parameterPrior(parameterPrior),
    m_env(scenarioSpace.env()),
    m_scenarioSpace(scenarioSpace),
    m_parameterSpace(parameterSpace),
    m_simulationOutputSpace(simulationOutputSpace),
    m_experimentOutputSpace(experimentOutputSpace),
    m_numSimulations(numSimulations),
    m_numExperiments(numExperiments),
    m_simulationScenarios(numSimulations, (V *)NULL),
    m_simulationParameters(numSimulations, (V *)NULL),
    m_simulationOutputs(numSimulations, (V *)NULL),
    m_experimentScenarios(numExperiments, (V *)NULL),
    m_experimentOutputs(numExperiments, (V *)NULL),
    m_numSimulationAdds(0),
    m_numExperimentAdds(0),
    priors(6, NULL)
{
  this->setUpHyperpriors();
  this->m_constructedGP = false;
}

template <class V, class M>
GaussianProcessFactory<V, M>::~GaussianProcessFactory()
{
  // Do nothing
}

template <class V, class M>
unsigned int
GaussianProcessFactory<V, M>::numSimulations() const
{
  return this->m_numSimulations;
}

template <class V, class M>
unsigned int
GaussianProcessFactory<V, M>::numExperiments() const
{
  return this->m_numExperiments;
}

template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessFactory<V, M>::scenarioSpace() const
{
  return this->m_scenarioSpace;
}

template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessFactory<V, M>::parameterSpace() const
{
  return this->m_parameterSpace;
}
 
template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessFactory<V, M>::simulationOutputSpace() const
{
  return this->m_simulationOutputSpace;
}

template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessFactory<V, M>::experimentOutputSpace() const
{
  return this->m_experimentOutputSpace;
}

template <class V, class M>
const V &
GaussianProcessFactory<V, M>::simulationScenario(
    unsigned int simulationId) const
{
  UQ_FATAL_TEST_MACRO(simulationId >= m_simulationScenarios.size(),
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::scenarioVec_original()",
                      "simulationId is too large");

  UQ_FATAL_TEST_MACRO(m_simulationScenarios[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::scenarioVec_original()",
                      "vector is NULL");

  return *(this->m_simulationScenarios[simulationId]);
}

template <class V, class M>
const std::vector<V *> &
GaussianProcessFactory<V, M>::simulationScenarios() const
{
  return this->m_simulationScenarios;
}

template <class V, class M>
const V &
GaussianProcessFactory<V, M>::simulationParameter(
    unsigned int simulationId) const
{
  UQ_FATAL_TEST_MACRO(simulationId >= m_simulationParameters.size(),
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::parameterVec_original()",
                      "simulationId is too large");

  UQ_FATAL_TEST_MACRO(m_simulationParameters[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::parameterVec_original()",
                      "vector is NULL");

  return *(this->m_simulationParameters[simulationId]);
}

template <class V, class M>
const std::vector<V *> &
GaussianProcessFactory<V, M>::simulationParameters() const
{
  return this->m_simulationParameters;
}

template <class V, class M>
const V &
GaussianProcessFactory<V, M>::simulationOutput(
    unsigned int simulationId) const
{
  UQ_FATAL_TEST_MACRO(simulationId >= m_simulationOutputs.size(),
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::outputVec_original()",
                      "simulationId is too large");

  UQ_FATAL_TEST_MACRO(m_simulationOutputs[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::outputVec_original()",
                      "vector is NULL");

  return *(this->m_simulationOutputs[simulationId]);
}

template <class V, class M>
const std::vector<V *> &
GaussianProcessFactory<V, M>::simulationOutputs() const
{
  return this->m_simulationOutputs;
}

template <class V, class M>
const V &
GaussianProcessFactory<V, M>::experimentScenario(
    unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= (this->m_experimentScenarios).size(),
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::experimentScenario()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(this->m_experimentScenarios[experimentId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::experimentScenario()",
                      "vector is NULL");

  return *(this->m_experimentScenarios[experimentId]);
}

template <class V, class M>
const std::vector<V *> &
GaussianProcessFactory<V, M>::experimentScenarios() const
{
  return this->m_experimentScenarios;
}

template <class V, class M>
const V &
GaussianProcessFactory<V, M>::experimentOutput(
    unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= (this->m_experimentOutputs).size(),
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::experimentOutput()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(this->m_experimentOutputs[experimentId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::experimentOutput()",
                      "vector is NULL");

  return *(this->m_experimentOutputs[experimentId]);
}

template <class V, class M>
const std::vector<V *> &
GaussianProcessFactory<V, M>::experimentOutputs() const
{
  return this->m_experimentOutputs;
}

template <class V, class M>
const M &
GaussianProcessFactory<V, M>::experimentErrors() const
{
  return *(this->m_experimentErrors);
}

template <class V, class M>
const BaseEnvironment &
GaussianProcessFactory<V, M>::env() const
{
  return this->m_env;
}

template <class V, class M>
const GaussianProcessEmulator<V, M> &
GaussianProcessFactory<V, M>::getGaussianProcess() const
{
  return *(this->gaussianProcess);
}

template <class V, class M>
void
GaussianProcessFactory<V, M>::addSimulation(V & simulationScenario,
                                            V & simulationParameter,
                                            V & simulationOutput)
{
  UQ_FATAL_TEST_MACRO(this->m_numSimulationAdds >= this->m_numSimulations,
                      this->m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::addSimulation()",
                      "too many simulation adds...");

  this->m_simulationScenarios[this->m_numSimulationAdds] = &simulationScenario;
  this->m_simulationParameters[this->m_numSimulationAdds] = &simulationParameter;
  this->m_simulationOutputs[this->m_numSimulationAdds] = &simulationOutput;
  this->m_numSimulationAdds++;

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->m_constructedGP = true;
    this->gaussianProcess = new GaussianProcessEmulator<V, M>(
        this->prior().imageSet(),
        this->m_scenarioSpace,
        this->m_parameterSpace,
        this->m_simulationOutputSpace,
        this->m_experimentOutputSpace,
        this->m_numSimulations,
        this->m_numExperiments,
        this->m_simulationScenarios,
        this->m_simulationParameters,
        this->m_simulationOutputs,
        this->m_experimentScenarios,
        this->m_experimentOutputs,
        *(this->m_experimentErrors),
        *(this->m_totalPrior));
  }
}

template <class V, class M>
void
GaussianProcessFactory<V, M>::addSimulations(
    const std::vector<V *> & simulationScenarios,
    const std::vector<V *> & simulationParameters,
    const std::vector<V *> & simulationOutputs)
{
  for (unsigned int i = 0; i < this->m_numSimulations; i++) {
    this->addSimulation(*(simulationScenarios[i]), *(simulationParameters[i]),
        *(simulationOutputs[i]));
  }
}

template <class V, class M>
void
GaussianProcessFactory<V, M>::addExperiments(
    const std::vector<V *> & experimentScenarios,
    const std::vector<V *> & experimentOutputs,
    const M * experimentErrors)
{
  UQ_FATAL_TEST_MACRO(experimentScenarios.size() > this->m_numExperiments,
                      this->m_env.worldRank(),
                      "GaussianProcessFactory<V, M>::addExperiment()",
                      "too many experiments...");

  for (unsigned int i = 0; i < this->m_experimentScenarios.size(); i++) {
    this->m_experimentScenarios[i] = experimentScenarios[i];
    this->m_experimentOutputs[i] = experimentOutputs[i];
  }
  this->m_experimentErrors = experimentErrors;
  this->m_numExperimentAdds += experimentScenarios.size();

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->m_constructedGP = true;
    this->gaussianProcess = new GaussianProcessEmulator<V, M>(
        this->prior().imageSet(),
        this->m_scenarioSpace,
        this->m_parameterSpace,
        this->m_simulationOutputSpace,
        this->m_experimentOutputSpace,
        this->m_numSimulations,
        this->m_numExperiments,
        this->m_simulationScenarios,
        this->m_simulationParameters,
        this->m_simulationOutputs,
        this->m_experimentScenarios,
        this->m_experimentOutputs,
        *(this->m_experimentErrors),
        *(this->m_totalPrior));
  }
}

template <class V, class M>
const ConcatenatedVectorRV<V, M> &
GaussianProcessFactory<V, M>::prior() const
{
  return *(this->m_totalPrior);
}

template <class V, class M>
void
GaussianProcessFactory<V, M>::print(std::ostream& os) const
{
  // Do nothing
}

// Private methods follow
template <class V, class M>
void
GaussianProcessFactory<V, M>::setUpHyperpriors()
{
  // Default value for hyperprior parameters
  this->m_emulatorPrecisionShape = 5.0;
  this->m_emulatorPrecisionScale = 1.0 / 5.0;
  this->m_emulatorCorrelationStrengthAlpha = 1.0;
  this->m_emulatorCorrelationStrengthBeta = 0.1;
  this->m_discrepancyPrecisionShape = 1.0;
  this->m_discrepancyPrecisionScale = 1.0 / 0.0001;
  this->m_discrepancyCorrelationStrengthAlpha = 1.0;
  this->m_discrepancyCorrelationStrengthBeta = 0.1;

  this->oneDSpace = new VectorSpace<V, M>(this->m_env, "", 1, NULL);

  // Emulator mean
  this->emulatorMeanMin = new V(this->oneDSpace->zeroVector());
  this->emulatorMeanMax = new V(this->oneDSpace->zeroVector());
  this->emulatorMeanMin->cwSet(-INFINITY);
  this->emulatorMeanMax->cwSet(INFINITY);

  this->emulatorMeanDomain = new BoxSubset<V, M>(
      "",
      *(this->oneDSpace),
      *(this->emulatorMeanMin),
      *(this->emulatorMeanMax));

  this->m_emulatorMean = new UniformVectorRV<V, M>(
      "",
      *(this->emulatorMeanDomain));

  // Emulator precision
  this->emulatorPrecisionMin = new V(this->oneDSpace->zeroVector());
  this->emulatorPrecisionMax = new V(this->oneDSpace->zeroVector());
  this->m_emulatorPrecisionShapeVec = new V(this->oneDSpace->zeroVector());
  this->m_emulatorPrecisionScaleVec = new V(this->oneDSpace->zeroVector());
  this->emulatorPrecisionMin->cwSet(0);
  this->emulatorPrecisionMax->cwSet(INFINITY);
  this->m_emulatorPrecisionShapeVec->cwSet(this->m_emulatorPrecisionShape);
  this->m_emulatorPrecisionScaleVec->cwSet(this->m_emulatorPrecisionScale);

  this->emulatorPrecisionDomain = new BoxSubset<V, M>(
      "",
      *(this->oneDSpace),
      *(this->emulatorPrecisionMin),
      *(this->emulatorPrecisionMax));

  this->m_emulatorPrecision = new GammaVectorRV<V, M>("",
      *(this->emulatorPrecisionDomain),
      *(this->m_emulatorPrecisionShapeVec),
      *(this->m_emulatorPrecisionScaleVec));

  // Emulator correlation strength
  unsigned int dimScenario = (this->scenarioSpace()).dimLocal();
  unsigned int dimParameter = (this->parameterSpace()).dimLocal();
  this->emulatorCorrelationSpace = new VectorSpace<V, M>(
      this->m_env,
      "",
      dimScenario + dimParameter,
      NULL);

  this->emulatorCorrelationMin = new V(
      this->emulatorCorrelationSpace->zeroVector());
  this->emulatorCorrelationMax = new V(
      this->emulatorCorrelationSpace->zeroVector());
  this->m_emulatorCorrelationStrengthAlphaVec = new V(
      this->emulatorCorrelationSpace->zeroVector());
  this->m_emulatorCorrelationStrengthBetaVec = new V(
      this->emulatorCorrelationSpace->zeroVector());
  this->emulatorCorrelationMin->cwSet(0);
  this->emulatorCorrelationMax->cwSet(1);
  this->m_emulatorCorrelationStrengthAlphaVec->cwSet(this->m_emulatorCorrelationStrengthAlpha);
  this->m_emulatorCorrelationStrengthBetaVec->cwSet(this->m_emulatorCorrelationStrengthBeta);

  this->emulatorCorrelationDomain = new BoxSubset<V, M>(
      "",
      *(this->emulatorCorrelationSpace),
      *(this->emulatorCorrelationMin),
      *(this->emulatorCorrelationMax));

  this->m_emulatorCorrelationStrength = new BetaVectorRV<V, M>("",
      *(this->emulatorCorrelationDomain),
      *(this->m_emulatorCorrelationStrengthAlphaVec),
      *(this->m_emulatorCorrelationStrengthBetaVec));

  // Discrepancy precision
  this->discrepancyPrecisionMin = new V(this->oneDSpace->zeroVector());
  this->discrepancyPrecisionMax = new V(this->oneDSpace->zeroVector());
  this->m_discrepancyPrecisionShapeVec = new V(this->oneDSpace->zeroVector());
  this->m_discrepancyPrecisionScaleVec = new V(this->oneDSpace->zeroVector());
  this->discrepancyPrecisionMin->cwSet(0);
  this->discrepancyPrecisionMax->cwSet(INFINITY);
  this->m_discrepancyPrecisionShapeVec->cwSet(this->m_discrepancyPrecisionShape);
  this->m_discrepancyPrecisionScaleVec->cwSet(this->m_discrepancyPrecisionScale);

  this->discrepancyPrecisionDomain = new BoxSubset<V, M>(
      "",
      *(this->oneDSpace),
      *(this->discrepancyPrecisionMin),
      *(this->emulatorPrecisionMax));

  this->m_discrepancyPrecision = new GammaVectorRV<V, M>("",
      *(this->discrepancyPrecisionDomain),
      *(this->m_discrepancyPrecisionShapeVec),
      *(this->m_discrepancyPrecisionScaleVec));

  // Discrepancy correlation strength
  this->discrepancyCorrelationSpace = new VectorSpace<V, M>(
      this->m_env,
      "",
      dimScenario,
      NULL);

  this->discrepancyCorrelationMin = new V(
      this->discrepancyCorrelationSpace->zeroVector());
  this->discrepancyCorrelationMax = new V(
      this->discrepancyCorrelationSpace->zeroVector());
  this->m_discrepancyCorrelationStrengthAlphaVec = new V(
      this->discrepancyCorrelationSpace->zeroVector());
  this->m_discrepancyCorrelationStrengthBetaVec = new V(
      this->discrepancyCorrelationSpace->zeroVector());
  this->discrepancyCorrelationMin->cwSet(0);
  this->discrepancyCorrelationMax->cwSet(1);
  this->m_discrepancyCorrelationStrengthAlphaVec->cwSet(this->m_discrepancyCorrelationStrengthAlpha);
  this->m_discrepancyCorrelationStrengthBetaVec->cwSet(this->m_discrepancyCorrelationStrengthBeta);

  this->discrepancyCorrelationDomain = new BoxSubset<V, M>(
      "",
      *(this->discrepancyCorrelationSpace),
      *(this->discrepancyCorrelationMin),
      *(this->discrepancyCorrelationMax));

  this->m_discrepancyCorrelationStrength = new BetaVectorRV<V, M>("",
      *(this->discrepancyCorrelationDomain),
      *(this->m_discrepancyCorrelationStrengthAlphaVec),
      *(this->m_discrepancyCorrelationStrengthBetaVec));

  // Now form full prior
  unsigned int dimSum = 3 +
                        dimParameter +
                        dimParameter +
                        dimScenario +
                        dimScenario;  // yum

  this->totalSpace = new VectorSpace<V, M>(
      this->m_env,
      "",
      dimSum,
      NULL);
  this->totalMins = new V(this->totalSpace->zeroVector());
  this->totalMaxs = new V(this->totalSpace->zeroVector());

  // Hackety hack McHackington.  There's no better way to do this unfortunately
  this->totalMins->cwSet(0);
  this->totalMaxs->cwSet(1);

  const BoxSubset<GslVector, GslMatrix> * box_prior =
    dynamic_cast<const BoxSubset<GslVector, GslMatrix>* >(&(this->m_parameterPrior.imageSet()));
  std::cout << "CTOR prior min[0]: " << box_prior->minValues()[0] << std::endl;
  for (unsigned int i = 0; i < dimParameter; i++) {
    (*(this->totalMins))[i] = box_prior->minValues()[i];
    (*(this->totalMaxs))[i] = box_prior->maxValues()[i];
  }

  (*(this->totalMins))[dimParameter] = -INFINITY;  // Min mean
  (*(this->totalMaxs))[dimParameter] = INFINITY;  // Max mean
  (*(this->totalMins))[dimParameter+1] = 0;  // Min emulator precision
  (*(this->totalMaxs))[dimParameter+1] = INFINITY;  // Max emulator precision

  // Min discrepancy precision
  (*(this->totalMins))[dimScenario+dimParameter+dimParameter+2] = 0;
  // Max discrepancy precision
  (*(this->totalMaxs))[dimScenario+dimParameter+dimParameter+2] = INFINITY;

  this->totalDomain = new BoxSubset<V, M>(
      "",
      *(this->totalSpace),
      *(this->totalMins),
      *(this->totalMaxs));

  this->priors[0] = &(this->m_parameterPrior);
  this->priors[1] = this->m_emulatorMean;
  this->priors[2] = this->m_emulatorPrecision;
  this->priors[3] = this->m_emulatorCorrelationStrength;
  this->priors[4] = this->m_discrepancyPrecision;
  this->priors[5] = this->m_discrepancyCorrelationStrength;

  // Finally
  this->m_totalPrior = new ConcatenatedVectorRV<V, M>(
      "",
      this->priors,
      *(this->totalDomain));
}

}  // End namespace QUESO

template class QUESO::GaussianProcessFactory<QUESO::GslVector, QUESO::GslMatrix>;
