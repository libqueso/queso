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
#include <queso/BoostInputOptionsParser.h>

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
  //! Given prefix, read the input file for parameters named prefix_*
  GPMSAOptions(const BaseEnvironment& env, const char* prefix);

  //! Destructor
  virtual ~GPMSAOptions();

  //! Prints \c this to \c os
  void print(std::ostream& os) const;

  //! The prefix to look for in the input file
  std::string m_prefix;

  //! If this string is non-empty, print the options object to the output file
  std::string m_help;

  //! The shape parameter for the Gamma hyperprior for the emulator precision
  double m_emulatorPrecisionShape;

  //! The scale parameter for the Gamma hyperprior for the emulator precision
  double m_emulatorPrecisionScale;

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
  const BaseEnvironment& m_env;

  BoostInputOptionsParser * m_parser;

  std::string m_option_help;
  std::string m_option_emulatorPrecisionShape;
  std::string m_option_emulatorPrecisionScale;
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
