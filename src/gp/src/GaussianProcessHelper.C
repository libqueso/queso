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

#include <queso/GaussianProcessHelper.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class V, class M>
GaussianProcessHelper<V, M>::GaussianProcessHelper(
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
    m_simulationScenarios(m_paper_m, (V*) NULL),
    m_simulationParameters(m_paper_m, (V*) NULL),
    m_simulationOutputs(m_paper_m, (V*) NULL),
    m_experimentScenarios(m_paper_m, (V*) NULL),
    m_experimentOutputs(m_paper_m, (V*) NULL),
    m_experimentErrors(m_paper_m, (M*) NULL),
    m_numSimulationAdds(0),
    m_numExperimentAdds(0),
    m_emulatorSpace(m_env, "emulator_",
        m_numSimulations * (m_simulationOutputSpace).dimLocal(), NULL),
    m_emulator(m_emulatorSpace.zeroVector())
{
  // Do nothing
  this->setUpHyperpriors();
}

template <class V, class M>
GaussianProcessHelper<V, M>::~GaussianProcessHelper()
{
  // Do nothing
}

template <class V, class M>
unsigned int
GaussianProcessHelper<V, M>::numSimulations() const
{
  return this->m_numSimulations;
}

template <class V, class M>
unsigned int
GaussianProcessHelper<V, M>::numExperiments() const
{
  return this->m_numExperiments;
}

template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessHelper<V, M>::scenarioSpace() const
{
  return this->m_scenarioSpace;
}

template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessHelper<V, M>::parameterSpace() const
{
  return this->m_parameterSpace;
}
 
template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessHelper<V, M>::simulationOutputSpace() const
{
  return this->m_simulationOutputSpace;
}

template <class V, class M>
const VectorSpace<V, M> &
GaussianProcessHelper<V, M>::experimentOutputSpace() const
{
  return this->m_experimentOutputSpace;
}

template <class V, class M>
const V &
GaussianProcessHelper<V, M>::simulationScenario(unsigned int simulationId) const
{
  UQ_FATAL_TEST_MACRO(simulationId >= m_scenarioVecs_original.size(),
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::scenarioVec_original()",
                      "simulationId is too large");

  UQ_FATAL_TEST_MACRO(m_scenarioVecs_original[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::scenarioVec_original()",
                      "vector is NULL");

  return *(this->m_simulationScenarios[simulationId]);
}

template <class V, class M>
const std::vector<const V *> &
GaussianProcessHelper<V, M>::simulationScenarios() const
{
  return this->m_simulationScenarios;
}

template <class V, class M>
const V &
GaussianProcessHelper<V, M>::simulationParameter(unsigned int simulationId) const
{
  UQ_FATAL_TEST_MACRO(simulationId >= m_parameterVecs_original.size(),
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::parameterVec_original()",
                      "simulationId is too large");

  UQ_FATAL_TEST_MACRO(m_parameterVecs_original[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::parameterVec_original()",
                      "vector is NULL");

  return *(this->m_simulationParameters[simulationId]);
}

template <class V, class M>
const std::vector<const V *> &
GaussianProcessHelper<V, M>::simulationParameters() const
{
  return this->m_simulationParameters;
}

template <class V, class M>
const V &
GaussianProcessHelper<V, M>::simulationOutput(unsigned int simulationId) const
{
  UQ_FATAL_TEST_MACRO(simulationId >= m_outputVecs_original.size(),
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::outputVec_original()",
                      "simulationId is too large");

  UQ_FATAL_TEST_MACRO(m_outputVecs_original[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::outputVec_original()",
                      "vector is NULL");

  return *(this->m_simulationOutputs[simulationId]);
}

template <class V, class M>
const std::vector<const V *> &
GaussianProcessHelper<V, M>::simulationOutputs() const
{
  return this->m_simulationOutputs;
}

template <class V, class M>
const V &
GaussianProcessHelper<V, M>::experimentScenario(unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= (this->m_experimentScenarios).size(),
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::experimentScenario()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(this->m_experimentScenarios[simulationId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::experimentScenario()",
                      "vector is NULL");

  return *(this->m_experimentScenarios[simulationId]);
}

template <class V, class M>
const std::vector<const V *> &
GaussianProcessHelper<V, M>::experimentScenarios() const
{
  return this->m_experimentScenarios;
}

template <class V, class M>
const V &
GaussianProcessHelper<V, M>::experimentOutput(unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= (this->m_experimentOutputs).size(),
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::experimentOutput()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(this->m_experimentOutputs[experimentId] == NULL,
                      m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::experimentOutput()",
                      "vector is NULL");

  return *(this->m_experimentOutputs[experimentId]);
}

template <class V, class M>
const std::vector<const V *> &
GaussianProcessHelper<V, M>::experimentOutputs() const
{
  return this->m_experimentOutputs;
}

template <class V, class M>
const M &
GaussianProcessHelper<V, M>::experimentError(unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= (this->m_experimentErrors).size(),
                      this->m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::experimentError()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(this->m_experimentErrors[experimentId] == NULL,
                      this->m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::experimentError()",
                      "matrix is NULL");

  return *(this->m_experimentErrors[experimentId]);
}

template <class V, class M>
const std::vector<const M *> &
GaussianProcessHelper<V, M>::experimentErrors() const
{
  return this->m_experimentErrors;
}

template <class V, class M>
const BaseEnvironment &
GaussianProcessHelper<V, M>::env() const
{
  return this->m_env;
}

template <class V, class M>
void
GaussianProcessHelper<V, M>::addSimulation(const V & simulationScenario,
                                           const V & simulationParameter,
                                           const V & simulationOutput)
{
  UQ_FATAL_TEST_MACRO(this->m_numSimulationAdds >= this->m_numSimulations,
                      this->m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::addSimulation()",
                      "too many simulation adds...");

  this->m_simulationScenarios[this->m_numSimulationAdds] = &simulationScenario;
  this->m_simulationParameters[this->m_numSimulationAdds] = &simulationParameter;
  this->m_simulationOutputs[this->m_numSimulationAdds] = &simulationOutput;
  this->m_numSimulationAdds++;

  // Done adding so form the emulator state
  if (this->m_numSimulationAdds == this->m_numSimulations) {
    this->m_emulator.cwSetConcatenated(this->m_simulationOutputs);
  }
}

template <class V, class M>
void
GaussianProcessHelper<V, M>::addExperiment(const V & experimentScenario,
                                           const V & experimentParameter,
                                           const M & experimentError)
{
  UQ_FATAL_TEST_MACRO(this->m_numExperimentAdds >= this->m_numExperiments,
                      this->m_env.worldRank(),
                      "GaussianProcessHelper<V, M>::addExperiment()",
                      "too many experiment adds...");

  this->m_experimentScenarios[this->m_numExperimentAdds] = &experimentScenario;
  this->m_experimentParameters[this->m_numExperimentAdds] = &experimentParameter;
  this->m_experimentErrors[this->m_numExperimentAdds] = &experimentError;
  this->m_numExperimentAdds++;
}

template <class V, class M>
void
GaussianProcessHelper<V, M>::addSimulations(
    const std::vector<const V *> & simulationScenarios,
    const std::vector<const V *> & simulationParameters,
    const std::vector<const V *> & simulationOutputs)
{
  for (unsigned int i = 0; i < this->m_numSimulations; i++) {
    this->addSimulation(simulationScenarios[i], simulationParameters[i],
        simulationOutputs[i]);
  }
}

template <class V, class M>
void
GaussianProcessHelper<V, M>::addExperiments(
    const std::vector<const V *> & experimentScenarios,
    const std::vector<const V *> & experimentOutputs,
    const std::vector<const M *> & experimentErrors)
{
  for (unsigned int i = 0; i < this->m_numExperiments; i++) {
    this->addExperiment(experimentScenarios[i], experimentOutputs[i],
        experimentErrors[i]);
  }
}

template <class V, class M>
void
GaussianProcessHelper<V, M>::print(std::ostream& os) const
{
  // Do nothing
}

}  // End namespace QUESO

// Private methods follow
template <class V, class M>
void
GaussianProcessHelper<V, M>::setUpHyperpriors()
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

  VectorSpace<V, M> oneDSpace(this->m_env, "", 1, NULL);

  // Emulator mean
  V emulatorMeanMin(oneDSpace.zeroVector());
  V emulatorMeanMax(oneDSpace.zeroVector());
  emulatorMeanMin.cwSet(-INFINITY);
  emulatorMeanMax.cwSet(INFINITY);

  BoxSubset<V, M> emulatorMeanDomain("", oneDSpace, emulatorMeanMin,
      emulatorMeanMax);

  this->m_emulatorMean = new UniformVectorRV<V, M>("", emulatorMeanDomain);

  // Emulator precision
  V emulatorPrecisionMin(oneDSpace.zeroVector());
  V emulatorPrecisionMax(oneDSpace.zeroVector());
  emulatorPrecisionMin.cwSet(0);
  emulatorPrecisionMax.cwSet(INFINITY);

  BoxSubset<V, M> emulatorPrecisionDomain("", oneDSpace, emulatorPrecisionMin,
      emulatorPrecisionMax);

  this->m_emulatorPrecision = new GammaVectorRV<V, M>("",
      emulatorPrecisionDomain,
      this->m_emulatorPrecisionShape,
      this->m_emulatorPrecisionScale);

  // Emulator correlation strength
  unsigned int dimScenario = (this->scenarioSpace).dimLocal();
  unsigned int dimParameter = (this->parameterSpace).dimLocal();
  VectorSpace<V, M> emulatorCorrelationSpace(this->m_env, "",
      dimScenario + dimParameter,
      NULL);

  V emulatorCorrelationMin(emulatorCorrelationSpace.zeroVector());
  V emulatorCorrelationMax(emulatorCorrelationSpace.zeroVector());
  emulatorCorrelationMin.cwSet(0);
  emulatorCorrelationMax.cwSet(1);

  BoxSubset<V, M> emulatorCorrelationDomain("", emulatorCorrelationSpace,
      emulatorCorrelationMin,
      emulatorCorrelationMax);

  this->m_emulatorCorrelationStrength = new BetaVectorRV<V, M>("",
      emulatorCorrelationDomain,
      this->m_emulatorCorrelationStrengthAlpha,
      this->m_emulatorCorrelationStrengthBeta);

  // Discrepancy precision
  V discrepancyPrecisionMin(oneDSpace.zeroVector());
  V discrepancyPrecisionMax(oneDSpace.zeroVector());
  discrepancyPrecisionMin.cwSet(0);
  descrepancyPrecisionMax.cwSet(INFINITY);

  BoxSubset<V, M> discrepancyPrecisionDomain("", oneDSpace,
      discrepancyPrecisionMin,
      emulatorPrecisionMax);

  this->m_discrepancyPrecision = new GammaVectorRV<V, M>("",
      discrepancyPrecisionDomain,
      this->m_discrepancyPrecisionShape,
      this->m_discrepancyPrecisionScale);

  // Discrepancy correlation strength
  unsigned int dimScenario = (this->scenarioSpace).dimLocal();
  VectorSpace<V, M> discrepancyCorrelationSpace(this->m_env, "",
      dimScenario,
      NULL);

  V discrepancyCorrelationMin(discrepancyCorrelationSpace.zeroVector());
  V discrepancyCorrelationMax(discrepancyCorrelationSpace.zeroVector());
  discrepancyCorrelationMin.cwSet(0);
  discrepancyCorrelationMax.cwSet(1);

  BoxSubset<V, M> discrepancyCorrelationDomain("", discrepancyCorrelationSpace,
      discrepancyCorrelationMin,
      discrepancyCorrelationMax);

  this->m_discrepancyCorrelationStrength = new BetaVectorRV<V, M>("",
      discrepancyCorrelationDomain,
      this->m_discrepancyCorrelationStrengthAlpha,
      this->m_discrepancyCorrelationStrengthBeta);

  // Now form full prior
  unsigned int dimSum = 3 + dimScenario + dimParameter + dimScenario;  // yum
  VectorSpace<V, M> totalSpace(this->m_env, "", dimSum, NULL);
  V totalMins(totalSpace.zeroVector());
  V totalMaxs(totalSpace.zeroVector());

  // Hackety hack McHackington.  There's no better way to do this unfortunately
  totalMins.cwSet(0);
  totalMaxs.cwSet(1);
  totalMins[0] = -INFINITY;
  totalMaxs[0] = INFINITY;
  totalMins[1] = 0;
  totalMaxs[1] = INFINITY;
  totalMins[dimScenario + dimParameter + 2] = 0;
  totalMaxs[dimScenario + dimParameter + 2] = INFINITY;

  BoxSubset<V, M> totalDomain("", totalSpace, totalMins, totalMaxs);

  std::vector<BaseVectorRV *> priors;
  priors.push_back(this->m_parameterPrior);
  priors.push_back(this->m_emulatorMean);
  priors.push_back(this->m_emulatorPrecision);
  priors.push_back(this->m_emulatorCorrelationStrength);
  priors.push_back(this->m_discrepancyPrecision);
  priors.push_back(this->m_discrepancyCorrelationStrength);

  // Finally
  this->m_totalPrior = new ConcatenatedVectorRV<V, M>("", priors, totalDomain);
}

template class QUESO::GaussianProcessHelper<QUESO::GslVector, QUESO::GslMatrix>;
