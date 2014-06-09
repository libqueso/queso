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

#include <queso/GPMSAOptions.h>

// ODV = option default value
#define UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV 5.0
#define UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV 0.2
#define UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_ALPHA_ODV 1.0
#define UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_BETA_ODV 0.1
#define UQ_GPMSA_DISCREPANCY_PRECISION_SHAPE_ODV 1.0
#define UQ_GPMSA_DISCREPANCY_PRECISION_SCALE_ODV 1e4
#define UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_ALPHA_ODV 1.0
#define UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_BETA_ODV 0.1
#define UQ_GPMSA_EMULATOR_DATA_PRECISION_SHAPE_ODV 3.0
#define UQ_GPMSA_EMULATOR_DATA_PRECISION_SCALE_ODV 333.333

namespace QUESO {

GPMSAOptions::GPMSAOptions(
  const BaseEnvironment & env,
  const char * prefix)
  :
  m_prefix((std::string)(prefix) + "gpmsa_"),
  m_env(env),
  m_optionsDesc(new po::options_description("Gaussian process emulator options")),
  m_option_help(m_prefix + "help"),
  m_option_emulatorPrecisionShape(m_prefix + "emulator_precision_shape"),
  m_option_emulatorPrecisionScale(m_prefix + "emulator_precision_scale"),
  m_option_emulatorCorrelationStrengthAlpha(m_prefix + "emulator_correlation_strength_alpha"),
  m_option_emulatorCorrelationStrengthBeta(m_prefix + "emulator_correlation_strength_beta"),
  m_option_discrepancyPrecisionShape(m_prefix + "discrepancy_precision_shape"),
  m_option_discrepancyPrecisionScale(m_prefix + "discrepancy_precision_scale"),
  m_option_discrepancyCorrelationStrengthAlpha(m_prefix + "discrepancy_correlation_strength_alpha"),
  m_option_discrepancyCorrelationStrengthBeta(m_prefix + "discrepancy_correlation_strength_beta"),
  m_option_emulatorDataPrecisionShape(m_prefix + "emulator_data_precision_shape"),
  m_option_emulatorDataPrecisionScale(m_prefix + "emulator_data_precision_scale")
{
  if (m_env.optionsInputFileName() == "") {
    queso_error_msg("Missing input file is required");
  }
}

GPMSAOptions::~GPMSAOptions()
{
  if (m_optionsDesc) {
    delete m_optionsDesc;
  }
}

void
GPMSAOptions::scanOptionsValues()
{
  if (m_optionsDesc == NULL) {
    queso_error_msg("m_optionsDesc variable is NULL");
  }

  defineMyOptions(*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues(*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In GPMSAOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

void
GPMSAOptions::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()
    (m_option_help.c_str(), "produce help message Gaussian process emulator")
    (m_option_emulatorPrecisionShape.c_str(), po::value<double>()->default_value(UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV), "shape hyperprior (Gamma) parameter for emulator precision")
    (m_option_emulatorPrecisionScale.c_str(), po::value<double>()->default_value(UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV), "scale hyperprior (Gamma) parameter for emulator precision")
    (m_option_emulatorCorrelationStrengthAlpha.c_str(), po::value<double>()->default_value(UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_ALPHA_ODV), "alpha hyperprior (Beta) parameter for emulator correlation strength")
    (m_option_emulatorCorrelationStrengthBeta.c_str(), po::value<double>()->default_value(UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_BETA_ODV), "beta hyperprior (Beta) parameter for emulator correlation strength")
    (m_option_discrepancyPrecisionShape.c_str(), po::value<double>()->default_value(UQ_GPMSA_DISCREPANCY_PRECISION_SHAPE_ODV), "shape hyperprior (Gamma) parameter for discrepancy precision")
    (m_option_discrepancyPrecisionScale.c_str(), po::value<double>()->default_value(UQ_GPMSA_DISCREPANCY_PRECISION_SCALE_ODV), "scale hyperprior (Gamma) parameter for discrepancy precision")
    (m_option_discrepancyCorrelationStrengthAlpha.c_str(), po::value<double>()->default_value(UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_ALPHA_ODV), "alpha hyperprior (Beta) parameter for discrepancy correlation strength")
    (m_option_discrepancyCorrelationStrengthBeta.c_str(), po::value<double>()->default_value(UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_BETA_ODV), "beta hyperprior (Beta) parameter for discrepancy correlation strength")
    (m_option_emulatorDataPrecisionShape.c_str(), po::value<double>()->default_value(UQ_GPMSA_EMULATOR_DATA_PRECISION_SHAPE_ODV), "shape hyperprior (Gamma) parameter for emulator data precision")
    (m_option_emulatorDataPrecisionScale.c_str(), po::value<double>()->default_value(UQ_GPMSA_EMULATOR_DATA_PRECISION_SCALE_ODV), "scale hyperprior (Gamma) parameter for emulator data precision");
}

void
GPMSAOptions::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_emulatorPrecisionShape)) {
    this->m_emulatorPrecisionShape = ((const po::variable_value &) m_env.allOptionsMap()[m_option_emulatorPrecisionShape]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_emulatorPrecisionScale)) {
    this->m_emulatorPrecisionScale = ((const po::variable_value &) m_env.allOptionsMap()[m_option_emulatorPrecisionScale]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_emulatorCorrelationStrengthAlpha)) {
    this->m_emulatorCorrelationStrengthAlpha = ((const po::variable_value &) m_env.allOptionsMap()[m_option_emulatorCorrelationStrengthAlpha]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_emulatorCorrelationStrengthBeta)) {
    this->m_emulatorCorrelationStrengthBeta = ((const po::variable_value &) m_env.allOptionsMap()[m_option_emulatorCorrelationStrengthBeta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_discrepancyPrecisionShape)) {
    this->m_discrepancyPrecisionShape = ((const po::variable_value &) m_env.allOptionsMap()[m_option_discrepancyPrecisionShape]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_discrepancyPrecisionScale)) {
    this->m_discrepancyPrecisionScale = ((const po::variable_value &) m_env.allOptionsMap()[m_option_discrepancyPrecisionScale]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_discrepancyCorrelationStrengthAlpha)) {
    this->m_discrepancyCorrelationStrengthAlpha = ((const po::variable_value &) m_env.allOptionsMap()[m_option_discrepancyCorrelationStrengthAlpha]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_discrepancyCorrelationStrengthBeta)) {
    this->m_discrepancyCorrelationStrengthBeta = ((const po::variable_value &) m_env.allOptionsMap()[m_option_discrepancyCorrelationStrengthBeta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_emulatorDataPrecisionShape)) {
    this->m_emulatorDataPrecisionShape = ((const po::variable_value &) m_env.allOptionsMap()[m_option_emulatorDataPrecisionShape]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_emulatorDataPrecisionScale)) {
    this->m_emulatorDataPrecisionScale = ((const po::variable_value &) m_env.allOptionsMap()[m_option_emulatorDataPrecisionScale]).as<double>();
  }
}

void
GPMSAOptions::print(std::ostream& os) const
{
  os << "\n" << m_option_emulatorPrecisionShape << " = " << this->m_emulatorPrecisionShape
     << "\n" << m_option_emulatorPrecisionScale << " = " << this->m_emulatorPrecisionScale
     << "\n" << m_option_emulatorCorrelationStrengthAlpha << " = " << this->m_emulatorCorrelationStrengthAlpha
     << "\n" << m_option_emulatorCorrelationStrengthBeta << " = " << this->m_emulatorCorrelationStrengthBeta
     << "\n" << m_option_discrepancyPrecisionShape << " = " << this->m_discrepancyPrecisionShape
     << "\n" << m_option_discrepancyPrecisionScale << " = " << this->m_discrepancyPrecisionScale
     << "\n" << m_option_discrepancyCorrelationStrengthAlpha << " = " << this->m_discrepancyCorrelationStrengthAlpha
     << "\n" << m_option_discrepancyCorrelationStrengthBeta << " = " << this->m_discrepancyCorrelationStrengthBeta
     << "\n" << m_option_emulatorDataPrecisionShape << " = " << this->m_emulatorDataPrecisionShape
     << "\n" << m_option_emulatorDataPrecisionScale << " = " << this->m_emulatorDataPrecisionScale
     << std::endl;
}

}  // End namespace QUESO
