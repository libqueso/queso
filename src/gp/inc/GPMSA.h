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

#ifndef UQ_GPMSA_HELPER_H
#define UQ_GPMSA_HELPER_H

#include <vector>

#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSpace.h>
#include <queso/VectorRV.h>
#include <queso/ConcatenatedVectorRV.h>
#include <queso/GammaVectorRV.h>
#include <queso/BetaVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/GPMSAOptions.h>
#include <queso/ScopedPtr.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class V = GslVector, class M = GslMatrix>
class GPMSAEmulator : public BaseScalarFunction<V, M>
{
public:
  GPMSAEmulator(const VectorSet<V, M> & domain,
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
                const std::vector<V>   & m_discrepancyBases,
                const M & m_experimentErrors,
                const ConcatenatedVectorRV<V, M> & m_totalPrior);

  virtual ~GPMSAEmulator();

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

        std::vector<V>     m_discrepancyBases;

  unsigned int num_svd_terms;
  typename ScopedPtr<M>::Type m_TruncatedSVD_simulationOutputs;

  // Matrix of svd basis vectors
  typename ScopedPtr<M>::Type K;

  // Saved calculation of K*K^T, to be used for LU invert-multiplies
  typename ScopedPtr<M>::Type KKT;

  // Total observation error covriance matrix
  const M & m_experimentErrors;

  const ConcatenatedVectorRV<V, M> & m_totalPrior;
};

template <class V = GslVector, class M = GslMatrix>
class GPMSAFactory
{
public:
  //! Constructor
  GPMSAFactory(const BaseEnvironment & env,
               GPMSAOptions * opts,
               const BaseVectorRV<V, M> & parameterPrior,
               const VectorSpace<V, M> & scenarioSpace,
               const VectorSpace<V, M> & parameterSpace,
               const VectorSpace<V, M> & simulationOutputSpace,
               const VectorSpace<V, M> & experimentOutputSpace,
               unsigned int numSimulations,
               unsigned int numExperiments);

  //! Destructor
  ~GPMSAFactory();

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

  //! Return the GPMSAEmulator likelihood object
  const GPMSAEmulator<V, M> & getGPMSAEmulator() const;

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

  //! Add all discrepancy bases to \c this
  /*!
   * This method takes a vector of *all* the bases to use in the
   * discrepancy model and stores a copy.
   *
   * The user is responsible for normalizing each basis vector to be
   * consistent with the discrepancy precision coefficients which will
   * multiply them.
   *
   * If no discrepancy basis is provided, a single "1 for each output"
   * vector will be used.
   */
  void setDiscrepancyBases(const std::vector<V *> & discrepancyBases);

  const ConcatenatedVectorRV<V, M> & prior() const;

  void print(std::ostream& os) const;
  friend std::ostream & operator<<(std::ostream& os,
                                   const GPMSAFactory<V, M> & obj)
  {
    obj.print(os);
    return os;
  }

  const BaseEnvironment & m_env;

  const BaseVectorRV<V, M> & m_parameterPrior;

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

  std::vector<V> m_discrepancyBases;

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
  typename ScopedPtr<VectorSpace<V, M> >::Type oneDSpace;

  // Emulator mean
  typename ScopedPtr<V>::Type emulatorMeanMin;
  typename ScopedPtr<V>::Type emulatorMeanMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type emulatorMeanDomain;

  // Emulator precision
  typename ScopedPtr<V>::Type emulatorPrecisionMin;
  typename ScopedPtr<V>::Type emulatorPrecisionMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type emulatorPrecisionDomain;

  // Emulator correlation strength
  typename ScopedPtr<VectorSpace<V, M> >::Type emulatorCorrelationSpace;
  typename ScopedPtr<V>::Type emulatorCorrelationMin;
  typename ScopedPtr<V>::Type emulatorCorrelationMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type emulatorCorrelationDomain;

  // Discrepancy precision
  typename ScopedPtr<V>::Type discrepancyPrecisionMin;
  typename ScopedPtr<V>::Type discrepancyPrecisionMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type discrepancyPrecisionDomain;

  // Discrepancy correlation strength
  typename ScopedPtr<VectorSpace<V, M> >::Type discrepancyCorrelationSpace;
  typename ScopedPtr<V>::Type discrepancyCorrelationMin;
  typename ScopedPtr<V>::Type discrepancyCorrelationMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type discrepancyCorrelationDomain;

  // Emulator data precision
  typename ScopedPtr<V>::Type emulatorDataPrecisionMin;
  typename ScopedPtr<V>::Type emulatorDataPrecisionMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type emulatorDataPrecisionDomain;

  // Now form full prior
  typename ScopedPtr<VectorSpace<V, M> >::Type totalSpace;
  typename ScopedPtr<V>::Type totalMins;
  typename ScopedPtr<V>::Type totalMaxs;

  typename ScopedPtr<BoxSubset<V, M> >::Type totalDomain;

  std::vector<const BaseVectorRV<V, M> *> priors;

  // The hyperpriors
  typename ScopedPtr<UniformVectorRV<V, M> >::Type m_emulatorMean;  // scalar
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_emulatorPrecision;  // (scalar) gamma(a, b) shape-rate
  typename ScopedPtr<BetaVectorRV<V, M> >::Type m_emulatorCorrelationStrength;  // (dim scenariosspace + dim parameterspace)
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_discrepancyPrecision;  // (scalar) shape-rate
  typename ScopedPtr<BetaVectorRV<V, M> >::Type m_discrepancyCorrelationStrength;  // (dim scenariospace)
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_emulatorDataPrecision;  // (scalar) shape-rate
  typename ScopedPtr<ConcatenatedVectorRV<V, M> >::Type m_totalPrior;  // prior for joint parameters and hyperparameters

  typename ScopedPtr<V>::Type m_emulatorPrecisionShapeVec;
  typename ScopedPtr<V>::Type m_emulatorPrecisionScaleVec;
  typename ScopedPtr<V>::Type m_emulatorCorrelationStrengthAlphaVec;
  typename ScopedPtr<V>::Type m_emulatorCorrelationStrengthBetaVec;
  typename ScopedPtr<V>::Type m_discrepancyPrecisionShapeVec;
  typename ScopedPtr<V>::Type m_discrepancyPrecisionScaleVec;
  typename ScopedPtr<V>::Type m_discrepancyCorrelationStrengthAlphaVec;
  typename ScopedPtr<V>::Type m_discrepancyCorrelationStrengthBetaVec;
  typename ScopedPtr<V>::Type m_emulatorDataPrecisionShapeVec;
  typename ScopedPtr<V>::Type m_emulatorDataPrecisionScaleVec;

  // The gaussian process object to build
  typename ScopedPtr<GPMSAEmulator<V, M> >::Type gpmsaEmulator;
  bool m_constructedGP;

private:
  bool allocated_m_opts;
  GPMSAOptions * m_opts;

};

}  // End namespace QUESO

#endif // UQ_GPMSA_HELPER_H
