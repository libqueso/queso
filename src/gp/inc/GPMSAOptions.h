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

#include <queso/Environment.h>
#include <queso/SharedPtr.h>

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#include <queso/ScopedPtr.h>
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

#ifndef UQ_GPMSA_OPTIONS_H
#define UQ_GPMSA_OPTIONS_H

namespace QUESO {

template <typename V>
class SimulationOutputMesh;

/*!
 * \file GPMSAOptions.h
 * \brief This class defines the options that specify the behaviour of the Gaussian process emulator
 *
 * \class GPMSAOptions
 * \brief This class defines the options that specify the behaviour of the Gaussian process emulator
 */

class GPMSAOptions
{
public:
  //! Given prefix, read the input file for parameters named "prefix"+*
  GPMSAOptions(const BaseEnvironment& env, const char* prefix);

  //! Construct with default parameters
  GPMSAOptions();

  //! Set parameter option names to begin with prefix
  void set_prefix(const char* prefix);

  //! Set default values for parameter options
  void set_defaults();

  //! Given prefix, read the input file for parameters named "prefix"+*
  void parse(const BaseEnvironment& env, const char* prefix);

  //! Destructor
  virtual ~GPMSAOptions();

  //! Prints \c this to \c os
  void print(std::ostream& os) const;

  //! The prefix to look for in the input file
  std::string m_prefix;

  //! If this string is non-empty, print the options object to the output file
  std::string m_help;

  //! The maximum number of basis vectors to use for approximating
  //  emulator output.  If this number is set to be zero (as it is by
  //  default), then the number of basis vectors will be determined at
  //  run time to preserve some minimum fraction of the variance in
  //  simulation output.
  int m_maxEmulatorBasisVectors;

  //! The minimum fraction of the variance in simulation output to
  //  capture with the emulator basis.  By default this is 1.0, i.e.
  //  100%, in which case there will be one basis vector for each
  //  dimension in the simulation output space.
  //
  //  Currently the only supported value is 1.0
  double m_emulatorBasisVarianceToCapture;

  //! Do automatic normalization, using minimum and maximum values in
  //  the supplied simulator data, for uncertain parameter i.
  //
  //  Normalized values will range from 0 to 1.
  void set_autoscale_minmax_uncertain_parameter(unsigned int i);

  //! Do automatic normalization, using minimum and maximum values in
  //  the supplied simulator data, for scenario parameter i.
  //
  //  Normalized values will range from 0 to 1.
  void set_autoscale_minmax_scenario_parameter(unsigned int i);

  //! Do automatic normalization, using minimum and maximum values in
  //  the supplied simulator data, for output variable i.  Note that,
  //  if functional data exists, variable i will *not* be the same as
  //  index i.
  //
  //  Normalized values will range from 0 to 1.
  void set_autoscale_minmax_output(unsigned int i);

  //! Do automatic normalization, using minimum and maximum values in
  //  the supplied simulator data, for all input uncertain parameters,
  //  all input scenario parameters, and all output values.
  //
  //  Normalized values will range from 0 to 1.
  void set_autoscale_minmax();

  //! Do automatic normalization, using mean and variance of the
  //  supplied simulator data, for uncertain parameter i.
  //
  //  Normalized values will have mean 0 and standard deviation 1.
  void set_autoscale_meanvar_uncertain_parameter(unsigned int i);

  //! Do automatic normalization, using mean and variance of the
  //  supplied simulator data, for scenario parameter i.
  //
  //  Normalized values will have mean 0 and standard deviation 1.
  void set_autoscale_meanvar_scenario_parameter(unsigned int i);

  //! Do automatic normalization, using mean and variance of the
  //  supplied simulator data, for output variable i.  Note that, when
  //  functional data exists, variable i will *not* be the same as
  //  index i.
  //
  //  Normalized values will have mean 0 and standard deviation 1.
  void set_autoscale_meanvar_output(unsigned int i);

  //! Do automatic normalization, using mean and variance of the
  //  supplied simulator data, for all input uncertain parameters, all
  //  input scenario parameters, and all output values.
  //
  //  Normalized values will have mean 0 and standard deviation 1.
  void set_autoscale_meanvar();

  //! Set a value, for uncertain parameter i in simulation inputs, of
  //  the physical parameter range (range_min, range_max) which should
  //  be rescaled to (0,1)
  //
  //  If no value is set and no automatic normalization is specified,
  //  a default of (0,1) will be used; i.e. no normalization.
  //
  //  This option is mutually exclusive with the above automatic
  //  scaling options; trying to use both is a runtime error.
  void set_uncertain_parameter_scaling(unsigned int i,
                                       double range_min,
                                       double range_max);

  //! Set a value, for scenario parameter i in simulation and
  //  experimental inputs, of the physical parameter range (range_min,
  //  range_max) which should be rescaled to (0,1)
  //
  //  If no value is set and no automatic normalization is specified,
  //  a default of (0,1) will be used; i.e. no normalization.
  //
  //  This option is mutually exclusive with the above automatic
  //  scaling options; trying to use both is a runtime error.
  void set_scenario_parameter_scaling(unsigned int i,
                                      double range_min,
                                      double range_max);

  //! Set a value, for output variable i in simulation and
  //  experimental outputs, of the physical output range (range_min,
  //  range_max) which should be rescaled to (0,1)
  //
  //  Note that, when functional data exists, variable i will *not* be
  //  the same as output vector index i.
  //
  //  If no value is set and no automatic normalization is specified,
  //  a default of (0,1) will be used; i.e. no normalization.
  //
  //  This option is mutually exclusive with the above automatic
  //  scaling options; trying to use both is a runtime error.
  void set_output_scaling(unsigned int i,
                          double range_min,
                          double range_max);

  //! Determine the physical parameter ranges (range_min, range_max)
  //  which should be rescaled to (0,1), based on our autoscaling option
  //  settings.
  template <typename V>
  void set_final_scaling
    (const std::vector<typename SharedPtr<V>::Type> & m_simulationScenarios,
     const std::vector<typename SharedPtr<V>::Type> & m_simulationParameters,
     const std::vector<typename SharedPtr<V>::Type> & m_simulationOutputs,
     const std::vector<typename SharedPtr<V>::Type> & m_experimentScenarios,
     const std::vector<typename SharedPtr<V>::Type> & m_experimentOutputs,
     const std::vector<typename SharedPtr<SimulationOutputMesh<V> >::Type> & m_simulationMeshes);

  //! Calculate a normalized value from a physical value for the
  //  specified scenario parameter.
  double normalized_scenario_parameter(unsigned int i,
                                       double physical_param) const;

  //! Calculate a normalized value from a physical value for the
  //  specified uncertain parameter.
  double normalized_uncertain_parameter(unsigned int i,
                                        double physical_param) const;

  //! Calculate a normalized value from a physical value for the
  //  output variable at vector index i.
  double normalized_output(unsigned int i,
                           double output_data) const;

  //! Calculate a normalized value from a physical value for the
  //  output variable i.
  double normalized_output_variable(unsigned int i,
                                    double output_data) const;

  //! Returns the scale, in physical units, corresponding to a single
  //  nondimensionalized unit for the output at vector index i.
  double output_scale(unsigned int i) const;

  //! Returns the scale, in physical units, corresponding to a single
  //  nondimensionalized unit for output variable i.
  double output_scale_variable(unsigned int i) const;

  //! The shape parameter for the Gamma hyperprior for the truncation error precision
  double m_truncationErrorPrecisionShape;

  //! The scale parameter for the Gamma hyperprior for the truncation error precision
  double m_truncationErrorPrecisionScale;

  //! The shape parameter for the Gamma hyperprior for the emulator precision
  double m_emulatorPrecisionShape;

  //! The scale parameter for the Gamma hyperprior for the emulator precision
  double m_emulatorPrecisionScale;

  //! Whether to use an observational error precision hyperparameter
  //  lambda_y rather than a fixed multiplier of 1.
  //
  // False by default.
  bool m_calibrateObservationalPrecision;

  //! The shape parameter for the Gamma hyperprior for the observational precision
  double m_observationalPrecisionShape;

  //! The scale parameter for the Gamma hyperprior for the observational precision
  double m_observationalPrecisionScale;


  //! The alpha paramter for the Beta hyperprior for the emulator correlation strength
  double m_emulatorCorrelationStrengthAlpha;

  //! The beta paramter for the Beta hyperprior for the emulator correlation strength
  double m_emulatorCorrelationStrengthBeta;

  //! The shape parameter for the Gamma hyperprior for the discrepancy precision
  double m_discrepancyPrecisionShape;

  //! The scale parameter for the Gamma hyperprior for the discrepancy precision
  double m_discrepancyPrecisionScale;

  //! The alpha paramter for the Beta hyperprior for the discrepancy correlation strength
  double m_discrepancyCorrelationStrengthAlpha;

  //! The beta paramter for the Beta hyperprior for the discrepancy correlation strength
  double m_discrepancyCorrelationStrengthBeta;

  //! Returns the QUESO environment
  const BaseEnvironment& env() const;

  //! The shape parameter for the Gamma hyperprior for the emulator data precision
  double m_emulatorDataPrecisionShape;

  //! The scale parameter for the Gamma hyperprior for the emulator data precision
  double m_emulatorDataPrecisionScale;

  //! The ridge to add to B^T*W_y*B before inverting it
  double m_observationalPrecisionRidge;

  //! The ridge to add to (B^T*W_y*B)^-1 before using it
  double m_observationalCovarianceRidge;

  //! The distance between gaussian discrepancy kernels, on each
  //  provided simulation mesh, in each direction.
  std::vector<double> m_gaussianDiscrepancyDistanceX,
                      m_gaussianDiscrepancyDistanceY,
                      m_gaussianDiscrepancyDistanceZ,
                      m_gaussianDiscrepancyDistanceT;

  //! Whether or not gaussian discrepancy kernels should be periodic,
  //  on each simulation mesh, in each direction.
  std::vector<bool> m_gaussianDiscrepancyPeriodicX,
                    m_gaussianDiscrepancyPeriodicY,
                    m_gaussianDiscrepancyPeriodicZ,
                    m_gaussianDiscrepancyPeriodicT;

  //! Gaussian discrepancy kernels have infinite support, and we want
  //  to include some kernels which are not centered within a bounded
  //  domain, so we have a cutoff (default 5% == 0.05) for what the
  //  value of a gaussian at the domain boundary must be, relative to
  //  the value of its out-of-boundary peak, for us to include it.
  //
  //  Currently this is a single scalar value, because I can't imagine
  //  users wanting it to differ from mesh to mesh (or from direction
  //  to direction) within a single problem.
  double m_gaussianDiscrepancySupportThreshold;

  friend std::ostream & operator<<(std::ostream& os, const GPMSAOptions & obj);

private:
  const BaseEnvironment * m_env;

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  QUESO::ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  // For fast lookups when doing normalization
  std::vector<unsigned int> m_output_index_to_variable_index;

  // True if the specified autoscaling should be done for all inputs
  bool m_autoscaleMinMaxAll;
  bool m_autoscaleMeanVarAll;

  // True if the specified autoscaling should be done for the specific
  // uncertain input parameter index
  std::set<unsigned int> m_autoscaleMinMaxUncertain;
  std::set<unsigned int> m_autoscaleMeanVarUncertain;

  // True if the specified autoscaling should be done for the specific
  // scenario input parameter index
  std::set<unsigned int> m_autoscaleMinMaxScenario;
  std::set<unsigned int> m_autoscaleMeanVarScenario;

  // True if the specified autoscaling should be done for the specific
  // output variable index; this index could apply to many vector
  // indices in the case of functional output.
  std::set<unsigned int> m_autoscaleMinMaxOutput;
  std::set<unsigned int> m_autoscaleMeanVarOutput;

  // The point in each uncertain input parameter range (typically min
  // or mean) corresponding to a normalized parameter of 0
  std::vector<double> m_uncertainScaleMin;

  // The width of the physical uncertain input parameter range
  // (typically max minus min or standard deviation) corresponding to
  // a normalized width of 1
  std::vector<double> m_uncertainScaleRange;

  // The point in each input scenario parameter range (typically min
  // or mean) corresponding to a normalized parameter of 0
  std::vector<double> m_scenarioScaleMin;

  // The width of the physical input scenario parameter range
  // (typically max minus min or standard deviation) corresponding to
  // a normalized width of 1
  std::vector<double> m_scenarioScaleRange;

  // The point in each output range (typically min
  // or mean) corresponding to a normalized parameter of 0
  std::vector<double> m_outputScaleMin;

  // The width of the physical output range (typically max minus min
  // or standard deviation) corresponding to a normalized width of 1
  std::vector<double> m_outputScaleRange;

  std::string m_option_help;
  std::string m_option_maxEmulatorBasisVectors;
  std::string m_option_truncationErrorPrecisionShape;
  std::string m_option_truncationErrorPrecisionScale;
  std::string m_option_emulatorBasisVarianceToCapture;
  std::string m_option_emulatorPrecisionShape;
  std::string m_option_emulatorPrecisionScale;
  std::string m_option_calibrateObservationalPrecision;
  std::string m_option_observationalPrecisionShape;
  std::string m_option_observationalPrecisionScale;
  std::string m_option_emulatorCorrelationStrengthAlpha;
  std::string m_option_emulatorCorrelationStrengthBeta;
  std::string m_option_discrepancyPrecisionShape;
  std::string m_option_discrepancyPrecisionScale;
  std::string m_option_discrepancyCorrelationStrengthAlpha;
  std::string m_option_discrepancyCorrelationStrengthBeta;
  std::string m_option_emulatorDataPrecisionShape;
  std::string m_option_emulatorDataPrecisionScale;
  std::string m_option_observationalPrecisionRidge;
  std::string m_option_observationalCovarianceRidge;

  std::string m_option_autoscaleMinMaxAll;
  std::string m_option_autoscaleMeanVarAll;

  std::string m_option_gaussianDiscrepancyDistanceX,
              m_option_gaussianDiscrepancyDistanceY,
              m_option_gaussianDiscrepancyDistanceZ,
              m_option_gaussianDiscrepancyDistanceT;

  std::string m_option_gaussianDiscrepancyPeriodicX,
              m_option_gaussianDiscrepancyPeriodicY,
              m_option_gaussianDiscrepancyPeriodicZ,
              m_option_gaussianDiscrepancyPeriodicT;

  std::string m_option_gaussianDiscrepancySupportThreshold;

  bool options_have_been_used;

  void checkOptions();
};

}  // End namespace QUESO

#endif // UQ_GPMSA_OPTIONS_H
