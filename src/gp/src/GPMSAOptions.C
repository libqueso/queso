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

#include <boost/program_options.hpp>

#include <queso/GPMSAOptions.h>

// ODV = option default value
#define UQ_GPMSA_HELP ""
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
  m_help(UQ_GPMSA_HELP),
  m_emulatorPrecisionShape(UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV),
  m_emulatorPrecisionScale(UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV),
  m_emulatorCorrelationStrengthAlpha(UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_ALPHA_ODV),
  m_emulatorCorrelationStrengthBeta(UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_BETA_ODV),
  m_discrepancyPrecisionShape(UQ_GPMSA_DISCREPANCY_PRECISION_SHAPE_ODV),
  m_discrepancyPrecisionScale(UQ_GPMSA_DISCREPANCY_PRECISION_SCALE_ODV),
  m_discrepancyCorrelationStrengthAlpha(UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_ALPHA_ODV),
  m_discrepancyCorrelationStrengthBeta(UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_BETA_ODV),
  m_emulatorDataPrecisionShape(UQ_GPMSA_EMULATOR_DATA_PRECISION_SHAPE_ODV),
  m_emulatorDataPrecisionScale(UQ_GPMSA_EMULATOR_DATA_PRECISION_SCALE_ODV),
  m_env(env),
  m_parser(new BoostInputOptionsParser(env.optionsInputFileName())),
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

  m_parser->registerOption<std::string>(m_option_help, UQ_GPMSA_HELP, "produce help message Gaussian process emulator");
  m_parser->registerOption<double>(m_option_emulatorPrecisionShape, UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV, "shape hyperprior (Gamma) parameter for emulator precision");
  m_parser->registerOption<double>(m_option_emulatorPrecisionScale, UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV, "scale hyperprior (Gamma) parameter for emulator precision");
  m_parser->registerOption<double>(m_option_emulatorCorrelationStrengthAlpha, UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_ALPHA_ODV, "alpha hyperprior (Beta) parameter for emulator correlation strength");
  m_parser->registerOption<double>(m_option_emulatorCorrelationStrengthBeta, UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_BETA_ODV, "beta hyperprior (Beta) parameter for emulator correlation strength");
  m_parser->registerOption<double>(m_option_discrepancyPrecisionShape, UQ_GPMSA_DISCREPANCY_PRECISION_SHAPE_ODV, "shape hyperprior (Gamma) parameter for discrepancy precision");
  m_parser->registerOption<double>(m_option_discrepancyPrecisionScale, UQ_GPMSA_DISCREPANCY_PRECISION_SCALE_ODV, "scale hyperprior (Gamma) parameter for discrepancy precision");
  m_parser->registerOption<double>(m_option_discrepancyCorrelationStrengthAlpha, UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_ALPHA_ODV, "alpha hyperprior (Beta) parameter for discrepancy correlation strength");
  m_parser->registerOption<double>(m_option_discrepancyCorrelationStrengthBeta, UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_BETA_ODV, "beta hyperprior (Beta) parameter for discrepancy correlation strength");
  m_parser->registerOption<double>(m_option_emulatorDataPrecisionShape, UQ_GPMSA_EMULATOR_DATA_PRECISION_SHAPE_ODV, "shape hyperprior (Gamma) parameter for emulator data precision");
  m_parser->registerOption<double>(m_option_emulatorDataPrecisionScale, UQ_GPMSA_EMULATOR_DATA_PRECISION_SCALE_ODV, "scale hyperprior (Gamma) parameter for emulator data precision");

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help,                           m_help);
  m_parser->getOption<double>(m_option_emulatorPrecisionShape,              m_emulatorPrecisionShape);
  m_parser->getOption<double>(m_option_emulatorPrecisionScale,              m_emulatorPrecisionScale);
  m_parser->getOption<double>(m_option_emulatorCorrelationStrengthAlpha,    m_emulatorCorrelationStrengthAlpha);
  m_parser->getOption<double>(m_option_emulatorCorrelationStrengthBeta,     m_emulatorCorrelationStrengthBeta);
  m_parser->getOption<double>(m_option_discrepancyPrecisionShape,           m_discrepancyPrecisionShape);
  m_parser->getOption<double>(m_option_discrepancyPrecisionScale,           m_discrepancyPrecisionScale);
  m_parser->getOption<double>(m_option_discrepancyCorrelationStrengthAlpha, m_discrepancyCorrelationStrengthAlpha);
  m_parser->getOption<double>(m_option_discrepancyCorrelationStrengthBeta,  m_discrepancyCorrelationStrengthBeta);
  m_parser->getOption<double>(m_option_emulatorDataPrecisionShape,          m_emulatorDataPrecisionShape);
  m_parser->getOption<double>(m_option_emulatorDataPrecisionScale,          m_emulatorDataPrecisionScale);

  checkOptions();
}

GPMSAOptions::~GPMSAOptions()
{
}

void
GPMSAOptions::checkOptions()
{
  if (m_help != "") {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*this) << std::endl;
    }
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

std::ostream &
operator<<(std::ostream& os, const GPMSAOptions & obj)
{
  os << (*(obj.m_parser)) << std::endl;
  obj.print(os);
  return os;
}


}  // End namespace QUESO
