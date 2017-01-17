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

#include <queso/Defines.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/GPMSAOptions.h>

// ODV = option default value
#define UQ_GPMSA_HELP ""
#define UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV 5.0
#define UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV 0.2
#define UQ_GPMSA_OBSERVATIONAL_PRECISION_SHAPE_ODV 5.0
#define UQ_GPMSA_OBSERVATIONAL_PRECISION_SCALE_ODV 0.2
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
{
  this->set_defaults();
  this->parse(env, prefix);
}


GPMSAOptions::GPMSAOptions()
  :
  m_env(NULL)
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  ,m_parser(new BoostInputOptionsParser())
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
{
  this->set_defaults();
  this->set_prefix("");
}


void
GPMSAOptions::set_prefix(const char * prefix)
{
  m_prefix = std::string(prefix) + "gpmsa_";

  m_option_help = m_prefix + "help";
  m_option_emulatorPrecisionShape = m_prefix + "emulator_precision_shape";
  m_option_emulatorPrecisionScale = m_prefix + "emulator_precision_scale";
  m_option_observationalPrecisionShape = m_prefix + "observational_precision_shape";
  m_option_observationalPrecisionScale = m_prefix + "observational_precision_scale";
  m_option_emulatorCorrelationStrengthAlpha = m_prefix + "emulator_correlation_strength_alpha";
  m_option_emulatorCorrelationStrengthBeta = m_prefix + "emulator_correlation_strength_beta";
  m_option_discrepancyPrecisionShape = m_prefix + "discrepancy_precision_shape";
  m_option_discrepancyPrecisionScale = m_prefix + "discrepancy_precision_scale";
  m_option_discrepancyCorrelationStrengthAlpha = m_prefix + "discrepancy_correlation_strength_alpha";
  m_option_discrepancyCorrelationStrengthBeta = m_prefix + "discrepancy_correlation_strength_beta";
  m_option_emulatorDataPrecisionShape = m_prefix + "emulator_data_precision_shape";
  m_option_emulatorDataPrecisionScale = m_prefix + "emulator_data_precision_scale";
}



void
GPMSAOptions::set_defaults()
{
  m_help = UQ_GPMSA_HELP;
  m_emulatorPrecisionShape = UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV;
  m_emulatorPrecisionScale = UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV;
  m_observationalPrecisionShape = UQ_GPMSA_OBSERVATIONAL_PRECISION_SHAPE_ODV;
  m_observationalPrecisionScale = UQ_GPMSA_OBSERVATIONAL_PRECISION_SCALE_ODV;
  m_emulatorCorrelationStrengthAlpha = UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_ALPHA_ODV;
  m_emulatorCorrelationStrengthBeta = UQ_GPMSA_EMULATOR_CORRELATION_STRENGTH_BETA_ODV;
  m_discrepancyPrecisionShape = UQ_GPMSA_DISCREPANCY_PRECISION_SHAPE_ODV;
  m_discrepancyPrecisionScale = UQ_GPMSA_DISCREPANCY_PRECISION_SCALE_ODV;
  m_discrepancyCorrelationStrengthAlpha = UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_ALPHA_ODV;
  m_discrepancyCorrelationStrengthBeta = UQ_GPMSA_DISCREPANCY_CORRELATION_STRENGTH_BETA_ODV;
  m_emulatorDataPrecisionShape = UQ_GPMSA_EMULATOR_DATA_PRECISION_SHAPE_ODV;
  m_emulatorDataPrecisionScale = UQ_GPMSA_EMULATOR_DATA_PRECISION_SCALE_ODV;

  checkOptions();
}


void
GPMSAOptions::parse(const BaseEnvironment & env,
                    const char * prefix)
{
  m_env = &env;

  if (m_env->optionsInputFileName() == "") {
    queso_error_msg("Missing input file is required");
  }

  this->set_prefix(prefix);

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser.reset(new BoostInputOptionsParser(env.optionsInputFileName()));

  m_parser->registerOption<std::string>
    (m_option_help,
     m_help,
     "produce help message Gaussian process emulator");

  m_parser->registerOption
    (m_option_emulatorPrecisionShape,
     m_emulatorPrecisionShape,
     "shape hyperprior (Gamma) parameter for emulator precision");
  m_parser->registerOption
    (m_option_emulatorPrecisionScale,
    m_emulatorPrecisionScale,
    "scale hyperprior (Gamma) parameter for emulator precision");

  m_parser->registerOption
    (m_option_observationalPrecisionShape,
    m_observationalPrecisionShape,
    "shape hyperprior (Gamma) parameter for observational precision");
  m_parser->registerOption
    (m_option_observationalPrecisionScale,
    m_observationalPrecisionScale,
    "scale hyperprior (Gamma) parameter for observational precision");

  m_parser->registerOption
    (m_option_emulatorCorrelationStrengthAlpha,
    m_emulatorCorrelationStrengthAlpha,
    "alpha hyperprior (Beta) parameter for emulator correlation strength");
  m_parser->registerOption
    (m_option_emulatorCorrelationStrengthBeta,
    m_emulatorCorrelationStrengthBeta,
    "beta hyperprior (Beta) parameter for emulator correlation strength");

  m_parser->registerOption
    (m_option_discrepancyPrecisionShape,
    m_discrepancyPrecisionShape,
    "shape hyperprior (Gamma) parameter for discrepancy precision");
  m_parser->registerOption
    (m_option_discrepancyPrecisionScale,
    m_discrepancyPrecisionScale,
    "scale hyperprior (Gamma) parameter for discrepancy precision");

  m_parser->registerOption
    (m_option_discrepancyCorrelationStrengthAlpha,
    m_discrepancyCorrelationStrengthAlpha,
    "alpha hyperprior (Beta) parameter for discrepancy correlation strength");
  m_parser->registerOption
    (m_option_discrepancyCorrelationStrengthBeta,
    m_discrepancyCorrelationStrengthBeta,
    "beta hyperprior (Beta) parameter for discrepancy correlation strength");

  m_parser->registerOption
    (m_option_emulatorDataPrecisionShape,
    m_emulatorDataPrecisionShape,
    "shape hyperprior (Gamma) parameter for emulator data precision");
  m_parser->registerOption
    (m_option_emulatorDataPrecisionScale,
    m_emulatorDataPrecisionScale,
    "scale hyperprior (Gamma) parameter for emulator data precision");

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help,                           m_help);
  m_parser->getOption<double>(m_option_emulatorPrecisionShape,              m_emulatorPrecisionShape);
  m_parser->getOption<double>(m_option_emulatorPrecisionScale,              m_emulatorPrecisionScale);
  m_parser->getOption<double>(m_option_observationalPrecisionShape,         m_observationalPrecisionShape);
  m_parser->getOption<double>(m_option_observationalPrecisionScale,         m_observationalPrecisionScale);
  m_parser->getOption<double>(m_option_emulatorCorrelationStrengthAlpha,    m_emulatorCorrelationStrengthAlpha);
  m_parser->getOption<double>(m_option_emulatorCorrelationStrengthBeta,     m_emulatorCorrelationStrengthBeta);
  m_parser->getOption<double>(m_option_discrepancyPrecisionShape,           m_discrepancyPrecisionShape);
  m_parser->getOption<double>(m_option_discrepancyPrecisionScale,           m_discrepancyPrecisionScale);
  m_parser->getOption<double>(m_option_discrepancyCorrelationStrengthAlpha, m_discrepancyCorrelationStrengthAlpha);
  m_parser->getOption<double>(m_option_discrepancyCorrelationStrengthBeta,  m_discrepancyCorrelationStrengthBeta);
  m_parser->getOption<double>(m_option_emulatorDataPrecisionShape,          m_emulatorDataPrecisionShape);
  m_parser->getOption<double>(m_option_emulatorDataPrecisionScale,          m_emulatorDataPrecisionScale);
#else
  m_help = env.input()(m_option_help, UQ_GPMSA_HELP);
  m_emulatorPrecisionShape =
    env.input()(m_option_emulatorPrecisionShape,
                m_emulatorPrecisionShape);
  m_emulatorPrecisionScale =
    env.input()(m_option_emulatorPrecisionScale,
                m_emulatorPrecisionScale);

  m_observationalPrecisionShape =
    env.input()(m_option_observationalPrecisionShape,
                m_observationalPrecisionShape);
  m_observationalPrecisionScale =
    env.input()(m_option_observationalPrecisionScale,
                m_observationalPrecisionScale);

  m_emulatorCorrelationStrengthAlpha =
    env.input()(m_option_emulatorCorrelationStrengthAlpha,
                m_emulatorCorrelationStrengthAlpha);
  m_emulatorCorrelationStrengthBeta =
    env.input()(m_option_emulatorCorrelationStrengthBeta,
                m_emulatorCorrelationStrengthBeta);

  m_discrepancyPrecisionShape =
    env.input()(m_option_discrepancyPrecisionShape,
                m_discrepancyPrecisionShape);
  m_discrepancyPrecisionScale =
    env.input()(m_option_discrepancyPrecisionScale,
                m_discrepancyPrecisionScale);

  m_discrepancyCorrelationStrengthAlpha =
    env.input()(m_option_discrepancyCorrelationStrengthAlpha,
                m_discrepancyCorrelationStrengthAlpha);
  m_discrepancyCorrelationStrengthBeta =
    env.input()(m_option_discrepancyCorrelationStrengthBeta,
                m_discrepancyCorrelationStrengthBeta);

  m_emulatorDataPrecisionShape =
    env.input()(m_option_emulatorDataPrecisionShape,
                m_emulatorDataPrecisionShape);
  m_emulatorDataPrecisionScale =
    env.input()(m_option_emulatorDataPrecisionScale,
                m_emulatorDataPrecisionScale);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

GPMSAOptions::~GPMSAOptions()
{
}

void
GPMSAOptions::checkOptions()
{
  if (m_help != "") {
    if (m_env && m_env->subDisplayFile()) {
      *m_env->subDisplayFile() << (*this) << std::endl;
    }
  }
}

void
GPMSAOptions::print(std::ostream& os) const
{
  os << "\n" << m_option_emulatorPrecisionShape << " = " << this->m_emulatorPrecisionShape
     << "\n" << m_option_emulatorPrecisionScale << " = " << this->m_emulatorPrecisionScale
     << "\n" << m_option_observationalPrecisionShape << " = " << this->m_observationalPrecisionShape
     << "\n" << m_option_observationalPrecisionScale << " = " << this->m_observationalPrecisionScale
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
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  obj.print(os);
  return os;
}


}  // End namespace QUESO
