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

#ifndef UQ_GP_HELPER_H
#define UQ_GP_HELPER_H

#include <vector>

#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSpace.h>
#include <queso/VectorRV.h>
#include <queso/ConcatenatedVectorRV.h>
#include <queso/GammaVectorRV.h>
#include <queso/BetaVectorRV.h>
#include <queso/UniformVectorRV.h>

namespace QUESO {

template <class V, class M>
class GaussianProcessEmulator : public BaseScalarFunction<V, M>
{
public:
  GaussianProcessEmulator(const VectorSet<V, M> & domain,
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
      const ConcatenatedVectorRV<V, M> & m_totalPrior);

  virtual ~GaussianProcessEmulator();

  virtual double lnValue(const V & domainVector,
                         const V * domainDirection,
                         V * gradVector,
                         M * hessianMatrix,
                         V * hessianEffect) const;

  virtual double actualValue(const V & domainVector,
                             const V * domainDirection,
                             V * gradVector,
                             M * hessianMatrix,
                             V * hessianEffect) const;

  const VectorSpace<V, M> & m_scenarioSpace;
  const VectorSpace<V, M> & m_parameterSpace;
  const VectorSpace<V, M> & m_simulationOutputSpace;
  const VectorSpace<V, M> & m_experimentOutputSpace;

  const unsigned int m_numSimulations;
  const unsigned int m_numExperiments;

  const std::vector<V *> & m_simulationScenarios;
  const std::vector<V *> & m_simulationParameters;
  const std::vector<V *> & m_simulationOutputs;
  const std::vector<V *> & m_experimentScenarios;
  const std::vector<V *> & m_experimentOutputs;

  // Total observation error covriance matrix
  const M & m_experimentErrors;

  const ConcatenatedVectorRV<V, M> & m_totalPrior;
};

template <class V, class M>
class GaussianProcessFactory
{
public:
  //! Constructor
  GaussianProcessFactory(const char * prefix,
                         const BaseVectorRV<V, M> & parameterPrior,
                         const VectorSpace<V, M> & scenarioSpace,
                         const VectorSpace<V, M> & parameterSpace,
                         const VectorSpace<V, M> & simulationOutputSpace,
                         const VectorSpace<V, M> & experimentOutputSpace,
                         unsigned int numSimulations,
                         unsigned int numExperiments);

  //! Destructor
  ~GaussianProcessFactory();

  //! @name Getters
  //@{

  //! Return number of simulations
  unsigned int numSimulations() const;

  //! Return number of experiments
  unsigned int numExperiments() const;

  //! Return the vector space in which scenarios live
  const VectorSpace<V, M> & scenarioSpace() const;

  //! Return the vector space in which parameters live
  const VectorSpace<V, M> & parameterSpace() const;

  //! Return the vector space in which simulations live
  const VectorSpace<V, M> & simulationOutputSpace() const;

  //! Return the vector space in which experiments live
  const VectorSpace<V, M> & experimentOutputSpace() const;

  //! Return the point in \c scenarioSpace for simulation \c simulationId
  /*!
   * This returns the point in scenario space at which simulation
   * \c simulationId was executed
   */
  const V & simulationScenario(unsigned int simulationId) const;

  //! Return all points in \c scenarioSpace for all simulations
  /*!
   * This returns all points in scenario space at which simulations were
   * executed
   */
  const std::vector<V *> & simulationScenarios() const;

  //! Return the point in \c parameterSpace for simulation \c simulationId
  /*!
   * This returns the point in parameter space at which simulation
   * \c simulationId was executed
   */
  const V & simulationParameter(unsigned int simulationId) const;

  //! Return all points in \c parameterSpace for all simulations
  /*!
   * This returns all points in parameter space at which simulations were
   * executed
   */
  const std::vector<V *> & simulationParameters() const;

  //! Return the simulation output for simulation \c simulationId
  /*!
   * The returned vector is a point in \c simulationOutputSpace
   */
  const V & simulationOutput(unsigned int simulationId) const;

  //! Return all points in \c simulationOutputSpace for all simulations
  /*!
   * This returns all points in simulation output space at which simulations
   * were executed
   */
  const std::vector<V *> & simulationOutputs() const;

  //! Return the point in \c scenarioSpace for experiment \c experimentId
  /*!
   * This returns the point in scenario space at which experiment
   * \c experimentId was executed
   */
  const V & experimentScenario(unsigned int experimentId) const;

  //! Return all points in \c scenarioSpace for all experiments
  /*!
   * This returns all points in scenario space at which experiments were
   * executed
   */
  const std::vector<V *> & experimentScenarios() const;

  //! Return the experiment output for experiment \c experimentId
  /*!
   * The returned vector is a point in \c experimentOutputSpace
   */
  const V & experimentOutput(unsigned int experimentId) const;

  //! Return all points in \c experimentOutputSpace for all experiments
  /*!
   * This returns all points in experiment output space at which experiments
   * were executed
   */
  const std::vector<V *> & experimentOutputs() const;

  //! Return all observation error covarince matrices for all experiments
  const M & experimentErrors() const;

  //! Return the QUESO environment
  const BaseEnvironment & env() const;

  //! Return the Gaussian process likelihood object
  const GaussianProcessEmulator<V, M> & getGaussianProcess() const;

  //@}

  //! Add a simulation to \c this
  /*!
   * The simulation added to \c this is assumed to correspond to the point
   * \c simulationScenario in scenario space and \c simulationParameter in
   * parameter space.  The simulation output is assumed to be stored in
   * \c simulationOutput.
   */
  void addSimulation(V & simulationScenario,
                     V & simulationParameter,
                     V & simulationOutput);

  //! Adds multiple simulations to \c this
  /*!
   * This method takes a vector of simulations and calls \c addSimulation on
   * each element
   */
  void addSimulations(const std::vector<V *> & simulationScenarios,
                      const std::vector<V *> & simulationParameters,
                      const std::vector<V *> & simulationOutputs);

  //! Add all experiments to \c this
  /*!
   * This method takes a vector of *all* the experimental data and associated
   * observation errors/correlations and stores them.  This cannot be done
   * piecemeal like the simulation data.
   *
   * Each experiment (\experimentOutputs[i]) is assumed to correspond to the
   * point \c expermientScenarios[i] in scenario space.  The observation error
   * covariance matrix is assumed to be stored in \c experimentErrors.
   */
  void addExperiments(const std::vector<V *> & experimentScenarios,
                      const std::vector<V *> & experimentOutputs,
                      const M * experimentErrors);

  const ConcatenatedVectorRV<V, M> & prior() const;

  void print(std::ostream& os) const;
  friend std::ostream & operator<<(std::ostream& os,
                                   const GaussianProcessFactory<V, M> & obj)
  {
    obj.print(os);
    return os;
  }

  const char * m_prefix;

  const BaseVectorRV<V, M> & m_parameterPrior;

  const BaseEnvironment & m_env;

  const VectorSpace<V, M> & m_scenarioSpace;
  const VectorSpace<V, M> & m_parameterSpace;
  const VectorSpace<V, M> & m_simulationOutputSpace;
  const VectorSpace<V, M> & m_experimentOutputSpace;

  unsigned int m_numSimulations;
  unsigned int m_numExperiments;

  std::vector<V *> m_simulationScenarios;
  std::vector<V *> m_simulationParameters;
  std::vector<V *> m_simulationOutputs;
  std::vector<V *> m_experimentScenarios;
  std::vector<V *> m_experimentOutputs;

  // Total observation error covriance matrix
  const M * m_experimentErrors;

  // Counter for the number of adds that happen
  unsigned int m_numSimulationAdds;
  unsigned int m_numExperimentAdds;

  // The space in which the emulator (simulator) lives
  // const VectorSpace<V, M> & m_emulatorSpace;
  
  // The emulator state
  // const V & m_emulator;

  // All the GP priors information for a scalar GP follows:
  void setUpHyperpriors();

  // Domains for all the hyperpriors
  VectorSpace<V, M> * oneDSpace;

  // Emulator mean
  V * emulatorMeanMin;
  V * emulatorMeanMax;
  BoxSubset<V, M> * emulatorMeanDomain;

  // Emulator precision
  V * emulatorPrecisionMin;
  V * emulatorPrecisionMax;
  BoxSubset<V, M> * emulatorPrecisionDomain;

  // Emulator correlation strength
  VectorSpace<V, M> * emulatorCorrelationSpace;

  V * emulatorCorrelationMin;
  V * emulatorCorrelationMax;

  BoxSubset<V, M> * emulatorCorrelationDomain;

  // Discrepancy precision
  V * discrepancyPrecisionMin;
  V * discrepancyPrecisionMax;

  BoxSubset<V, M> * discrepancyPrecisionDomain;

  // Discrepancy correlation strength
  VectorSpace<V, M> * discrepancyCorrelationSpace;

  V * discrepancyCorrelationMin;
  V * discrepancyCorrelationMax;

  BoxSubset<V, M> * discrepancyCorrelationDomain;

  // Now form full prior
  VectorSpace<V, M> * totalSpace;
  V * totalMins;
  V * totalMaxs;

  BoxSubset<V, M> * totalDomain;

  std::vector<const BaseVectorRV<V, M> *> * priors;

  // The hyperpriors
  UniformVectorRV<V, M> * m_emulatorMean;  // scalar
  GammaVectorRV<V, M> * m_emulatorPrecision;  // (scalar) gamma(a, b) shape-rate
  BetaVectorRV<V, M> * m_emulatorCorrelationStrength;  // (dim scenariosspace + dim parameterspace)
  GammaVectorRV<V, M> * m_discrepancyPrecision;  // (scalar) shape-rate
  BetaVectorRV<V, M> * m_discrepancyCorrelationStrength;  // (dim scenariospace)
  ConcatenatedVectorRV<V, M> * m_totalPrior;  // prior for joint parameters and hyperparameters

  // Parameters for hyperpriors.  Hyper-hyperparameters, if you will.  Yo dawg.
  double m_emulatorPrecisionShape;
  double m_emulatorPrecisionScale;
  double m_emulatorCorrelationStrengthAlpha;
  double m_emulatorCorrelationStrengthBeta;
  double m_discrepancyPrecisionShape;
  double m_discrepancyPrecisionScale;
  double m_discrepancyCorrelationStrengthAlpha;
  double m_discrepancyCorrelationStrengthBeta;

  V * m_emulatorPrecisionShapeVec;
  V * m_emulatorPrecisionScaleVec;
  V * m_emulatorCorrelationStrengthAlphaVec;
  V * m_emulatorCorrelationStrengthBetaVec;
  V * m_discrepancyPrecisionShapeVec;
  V * m_discrepancyPrecisionScaleVec;
  V * m_discrepancyCorrelationStrengthAlphaVec;
  V * m_discrepancyCorrelationStrengthBetaVec;

  // The gaussian process object to build
  GaussianProcessEmulator<V, M> * gaussianProcess;
};

}  // End namespace QUESO

#endif // UQ_GP_HELPER_H
