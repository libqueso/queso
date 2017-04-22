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
#include <queso/SharedPtr.h>

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
                const std::vector<typename SharedPtr<V>::Type> & m_simulationScenarios,
                const std::vector<typename SharedPtr<V>::Type> & m_simulationParameters,
                const std::vector<typename SharedPtr<V>::Type> & m_simulationOutputs,
                const std::vector<typename SharedPtr<V>::Type> & m_experimentScenarios,
                const std::vector<typename SharedPtr<V>::Type> & m_experimentOutputs,
                const std::vector<typename SharedPtr<V>::Type> & m_discrepancyBases,
                const std::vector<typename SharedPtr<M>::Type> & m_observationErrorMatrices,
                const typename SharedPtr<M>::Type & m_observationErrorMatrix,
                const ConcatenatedVectorRV<V, M> & m_totalPrior,
                const V & residual_in,
                const M & BT_Wy_B_inv_in,
                const M & KT_K_inv_in,
                const GPMSAOptions & opts);

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

  const std::vector<typename SharedPtr<V>::Type> & m_simulationScenarios;
  const std::vector<typename SharedPtr<V>::Type> & m_simulationParameters;
  const std::vector<typename SharedPtr<V>::Type> & m_simulationOutputs;
  const std::vector<typename SharedPtr<V>::Type> & m_experimentScenarios;
  const std::vector<typename SharedPtr<V>::Type> & m_experimentOutputs;

        std::vector<typename SharedPtr<V>::Type>   m_discrepancyBases;

  const std::vector<typename SharedPtr<M>::Type> & m_observationErrorMatrices;

  // Block diagonal matrix; sacrificing efficiency for clarity
  typename SharedPtr<M>::Type m_observationErrorMatrix;

  //
  // Intermediate calculations we can cache
  //
  unsigned int num_svd_terms;

  const ConcatenatedVectorRV<V, M> & m_totalPrior;

  //
  // Intermediate calculations cached by factory
  //
  const V & residual;

  const M & BT_Wy_B_inv;

  const M & KT_K_inv;

  const GPMSAOptions & m_opts;
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

  //! Return GPMSAOptions structure
  const GPMSAOptions & options() const;

  GPMSAOptions & options();

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
  const std::vector<typename SharedPtr<V>::Type> & simulationScenarios() const;

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
  const std::vector<typename SharedPtr<V>::Type> & simulationParameters() const;

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
  const std::vector<typename SharedPtr<V>::Type> & simulationOutputs() const;

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
  const std::vector<typename SharedPtr<V>::Type> & experimentScenarios() const;

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
  const std::vector<typename SharedPtr<V>::Type> & experimentOutputs() const;

  //! Return all observation error covarince matrices for all experiments
  const typename SharedPtr<M>::Type experimentErrors() const;

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
  void addSimulation(typename SharedPtr<V>::Type simulationScenario,
                     typename SharedPtr<V>::Type simulationParameter,
                     typename SharedPtr<V>::Type simulationOutput);

  //! Adds multiple simulations to \c this
  /*!
   * This method takes a vector of simulations and calls \c addSimulation on
   * each element
   */
  void addSimulations(const std::vector<typename SharedPtr<V>::Type> & simulationScenarios,
                      const std::vector<typename SharedPtr<V>::Type> & simulationParameters,
                      const std::vector<typename SharedPtr<V>::Type> & simulationOutputs);

  //! Add all experiments to \c this
  /*!
   * This method takes a vector of *all* the experimental data and associated
   * observation errors/correlations and stores them.  This cannot be done
   * piecemeal like the simulation data.
   *
   * Each experiment (\experimentOutputs[i]) is assumed to correspond to the
   * point \c expermientScenarios[i] in scenario space.  The observation error
   * covariance matrix is assumed to be stored in \c experimentErrors.
   *
   * This method is solely for backward compatibility.  Covariances
   * between errors from different experiments will be ignored.
   */
  void addExperiments(const std::vector<typename SharedPtr<V>::Type> & experimentScenarios,
                      const std::vector<typename SharedPtr<V>::Type> & experimentOutputs,
                      const typename SharedPtr<M>::Type experimentErrors);

  //! Add all experiments to \c this
  /*!
   * This method takes a vector of *all* the experimental data and associated
   * observation errors/correlations and stores them.  This cannot be done
   * piecemeal like the simulation data.
   *
   * Each experiment (\experimentOutputs[i]) is assumed to correspond to the
   * point \c expermientScenarios[i] in scenario space.  Each
   * experiment has a corresponding error covariance matrix stored in
   * \c experimentErrors[i]
   */
  void addExperiments(const std::vector<typename SharedPtr<V>::Type> & experimentScenarios,
                      const std::vector<typename SharedPtr<V>::Type> & experimentOutputs,
                      const std::vector<typename SharedPtr<M>::Type> & experimentErrors);

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
   *
   * For now we will assume that outputs are the same for each
   * experiment (tau_ij and phi_ij are independent of i, in the
   * notation of Higdon et al), and so each discrepancy basis can be
   * expressed as a vector indexed by output index.
   */
  void setDiscrepancyBases(const std::vector<typename SharedPtr<V>::Type> & discrepancyBases);

  M & getObservationErrorCovariance(unsigned int simulationNumber);

  const M & getObservationErrorCovariance(unsigned int simulationNumber) const;

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

  std::vector<typename SharedPtr<V>::Type> m_simulationScenarios;
  std::vector<typename SharedPtr<V>::Type> m_simulationParameters;
  std::vector<typename SharedPtr<V>::Type> m_simulationOutputs;
  std::vector<typename SharedPtr<V>::Type> m_experimentScenarios;
  std::vector<typename SharedPtr<V>::Type> m_experimentOutputs;

  // We will be recentering data around the simulation output mean
  typename ScopedPtr<V>::Type simulationOutputMeans;

  std::vector<typename SharedPtr<V>::Type> m_discrepancyBases;

  std::vector<typename SharedPtr<M>::Type> m_observationErrorMatrices;

  // Counter for the number of adds that happen
  unsigned int m_numSimulationAdds;
  unsigned int m_numExperimentAdds;

  // The space in which the emulator (simulator) lives
  // const VectorSpace<V, M> & m_emulatorSpace;

  // The emulator state
  // const V & m_emulator;

  // Build the emulator once all data has been added
  void setUpEmulator();

  // All the GP priors information for a scalar GP follows:
  void setUpHyperpriors();

  // Number of dimensions preserved from the SVD of simulation outputs
  unsigned int num_svd_terms;

  // Domains for all the hyperpriors
  typename ScopedPtr<VectorSpace<V, M> >::Type oneDSpace;

  // Truncation error precision
  typename ScopedPtr<V>::Type truncationErrorPrecisionMin;
  typename ScopedPtr<V>::Type truncationErrorPrecisionMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type truncationErrorPrecisionDomain;

  // Emulator precision
  typename ScopedPtr<VectorSpace<V, M> >::Type emulatorPrecisionSpace;
  typename ScopedPtr<V>::Type emulatorPrecisionMin;
  typename ScopedPtr<V>::Type emulatorPrecisionMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type emulatorPrecisionDomain;

  // Emulator correlation strength
  typename ScopedPtr<VectorSpace<V, M> >::Type emulatorCorrelationSpace;
  typename ScopedPtr<V>::Type emulatorCorrelationMin;
  typename ScopedPtr<V>::Type emulatorCorrelationMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type emulatorCorrelationDomain;

  // Observational precision
  typename ScopedPtr<VectorSpace<V, M> >::Type observationalPrecisionSpace;
  typename ScopedPtr<V>::Type observationalPrecisionMin;
  typename ScopedPtr<V>::Type observationalPrecisionMax;
  typename ScopedPtr<BoxSubset<V, M> >::Type observationalPrecisionDomain;

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
  typename ScopedPtr<UniformVectorRV<V, M> >::Type m_truncationErrorPrecision;  // scalar
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_emulatorPrecision;  // (dim num_svd_terms) gamma(a, b) shape-rate
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_observationalPrecision;  // scalar gamma(a, b) shape-rate
  typename ScopedPtr<BetaVectorRV<V, M> >::Type m_emulatorCorrelationStrength;  // (dim scenariosspace + dim parameterspace)
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_discrepancyPrecision;  // (scalar) shape-rate
  typename ScopedPtr<BetaVectorRV<V, M> >::Type m_discrepancyCorrelationStrength;  // (dim scenariospace)
  typename ScopedPtr<GammaVectorRV<V, M> >::Type m_emulatorDataPrecision;  // (scalar) shape-rate
  typename ScopedPtr<ConcatenatedVectorRV<V, M> >::Type m_totalPrior;  // prior for joint parameters and hyperparameters

  typename ScopedPtr<V>::Type m_emulatorPrecisionShapeVec;
  typename ScopedPtr<V>::Type m_emulatorPrecisionScaleVec;
  typename ScopedPtr<V>::Type m_observationalPrecisionShapeVec;
  typename ScopedPtr<V>::Type m_observationalPrecisionScaleVec;
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

  // Block diagonal matrix; sacrificing efficiency for clarity
  typename SharedPtr<M>::Type m_observationErrorMatrix;

  //
  // Intermediate calculations we can cache
  //
  typename ScopedPtr<M>::Type m_TruncatedSVD_simulationOutputs;

  std::vector<M> m_discrepancyMatrices;

  typename ScopedPtr<M>::Type m_BMatrix;

  // Matrix of svd basis vectors
  typename ScopedPtr<M>::Type K;

  typename ScopedPtr<V>::Type residual;

  // Cached calculation of (B^T*W_y*B)^-1
  typename ScopedPtr<M>::Type BT_Wy_B_inv;

  // Cached calculation of (K^T*K)^-1
  typename ScopedPtr<M>::Type KT_K_inv;


private:
  bool allocated_m_opts;
  GPMSAOptions * m_opts;
};

}  // End namespace QUESO

#endif // UQ_GPMSA_HELPER_H
