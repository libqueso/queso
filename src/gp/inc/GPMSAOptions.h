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

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#include <queso/ScopedPtr.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#ifndef UQ_GPMSA_OPTIONS_H
#define UQ_GPMSA_OPTIONS_H

namespace QUESO {

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

  friend std::ostream & operator<<(std::ostream& os, const GPMSAOptions & obj);

private:
  const BaseEnvironment * m_env;

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  QUESO::ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  std::string m_option_help;
  std::string m_option_maxEmulatorBasisVectors;
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

  void checkOptions();
};

}  // End namespace QUESO

#endif // UQ_GPMSA_OPTIONS_H
