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

#include <queso/Defines.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#include <queso/getpot.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/GPMSAOptions.h>

#include <queso/GslVector.h>


// ODV = option default value
#define UQ_GPMSA_HELP ""
#define UQ_GPMSA_MAX_SIMULATOR_BASIS_VECTORS_ODV 0
#define UQ_GPMSA_SIMULATOR_BASIS_VARIANCE_TO_CAPTURE 1.0
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

namespace { // Anonymous namespace for helper functions

template <typename V>
void min_max_update(V & min, V & max, const V & new_data)
{
  unsigned int dim = min.sizeGlobal();
  queso_assert_equal_to(dim, max.sizeGlobal());
  queso_assert_equal_to(dim, new_data.sizeGlobal());

  for (unsigned int p=0; p != dim; ++p)
    {
      min[p] = std::min(min[p], new_data[p]);
      max[p] = std::max(max[p], new_data[p]);
    }
}


template <typename V>
void mean_var_update(unsigned int & n, V & mean, V & var, const V & new_data)
{
  ++n;

  const V delta (new_data - mean);

  V delta_n (delta);
  delta_n /= n;

  mean += delta_n;

  V delta2 (new_data - mean);
  delta2 *= delta;

  var += delta2;
}

} // end anonymous namespace



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
  m_option_maxEmulatorBasisVectors = m_prefix + "max_emulator_basis_vectors";
  m_option_emulatorBasisVarianceToCapture = m_prefix + "emulator_basis_variance_to_capture";
  m_option_emulatorPrecisionShape = m_prefix + "emulator_precision_shape";
  m_option_emulatorPrecisionScale = m_prefix + "emulator_precision_scale";
  m_option_calibrateObservationalPrecision = m_prefix + "calibrate_observational_precision";
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
  m_option_autoscaleMinMaxAll = m_prefix + "autoscale_min_max_all";
  m_option_autoscaleMeanVarAll = m_prefix + "autoscale_mean_var_all";
}



void
GPMSAOptions::set_defaults()
{
  m_help = UQ_GPMSA_HELP;
  m_maxEmulatorBasisVectors = UQ_GPMSA_MAX_SIMULATOR_BASIS_VECTORS_ODV;
  m_emulatorBasisVarianceToCapture = UQ_GPMSA_SIMULATOR_BASIS_VARIANCE_TO_CAPTURE;
  m_emulatorPrecisionShape = UQ_GPMSA_EMULATOR_PRECISION_SHAPE_ODV;
  m_emulatorPrecisionScale = UQ_GPMSA_EMULATOR_PRECISION_SCALE_ODV;
  m_calibrateObservationalPrecision = false;
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

  m_autoscaleMinMaxAll = false;
  m_autoscaleMeanVarAll = false;

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
    (m_option_calibrateObservationalPrecision,
    m_calibrateObservationalPrecision,
    "whether to use a calibrated hyperparameter for observational precision");

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

  m_parser->registerOption
    (m_option_autoscaleMinMaxAll,
    m_autoscaleMinMaxAll,
    "option to autoscale all parameters and outputs based on data range");
  m_parser->registerOption
    (m_option_autoscaleMeanVarAll,
    m_autoscaleMeanVarAll,
    "option to autoscale all parameters and outputs based on data statistics");

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help,                           m_help);
  m_parser->getOption<double>(m_option_emulatorPrecisionShape,              m_emulatorPrecisionShape);
  m_parser->getOption<double>(m_option_emulatorPrecisionScale,              m_emulatorPrecisionScale);
  m_parser->getOption<bool>  (m_option_calibrateObservationalPrecision,     m_calibrateObservationalPrecision);
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
  m_parser->getOption<bool>  (m_option_autoscaleMinMaxAll,                  m_autoscaleMinMaxAll);
  m_parser->getOption<bool>  (m_option_autoscaleMeanVarAll,                 m_autoscaleMeanVarAll);
#else
  m_help = env.input()(m_option_help, UQ_GPMSA_HELP);
  m_emulatorPrecisionShape =
    env.input()(m_option_emulatorPrecisionShape,
                m_emulatorPrecisionShape);
  m_emulatorPrecisionScale =
    env.input()(m_option_emulatorPrecisionScale,
                m_emulatorPrecisionScale);

  m_calibrateObservationalPrecision =
    env.input()(m_option_calibrateObservationalPrecision,
                m_calibrateObservationalPrecision);
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

  m_autoscaleMinMaxAll =
    env.input()(m_option_autoscaleMinMaxAll,
                m_autoscaleMinMaxAll);
  m_autoscaleMeanVarAll =
    env.input()(m_option_autoscaleMeanVarAll,
                m_autoscaleMeanVarAll);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

GPMSAOptions::~GPMSAOptions()
{
}


void
GPMSAOptions::set_autoscale_minmax()
{
  this->m_autoscaleMinMaxAll = true;
}


void
GPMSAOptions::set_autoscale_minmax_uncertain_parameter(unsigned int i)
{
  m_autoscaleMinMaxUncertain.insert(i);
}


void
GPMSAOptions::set_autoscale_minmax_scenario_parameter(unsigned int i)
{
  m_autoscaleMinMaxScenario.insert(i);
}


void
GPMSAOptions::set_autoscale_meanvar()
{
  this->m_autoscaleMeanVarAll = true;
}


void
GPMSAOptions::set_autoscale_meanvar_uncertain_parameter(unsigned int i)
{
  m_autoscaleMeanVarUncertain.insert(i);
}


void
GPMSAOptions::set_autoscale_meanvar_scenario_parameter(unsigned int i)
{
  m_autoscaleMeanVarScenario.insert(i);
}


void
GPMSAOptions::set_uncertain_parameter_scaling(unsigned int i,
                                              double range_min,
                                              double range_max)
{
  if (i >= m_uncertainScaleMin.size())
  {
    m_uncertainScaleMin.resize(i+1, 0);
    m_uncertainScaleRange.resize(i+1, 1);
  }
  m_uncertainScaleMin[i] = range_min;
  m_uncertainScaleRange[i] = range_max - range_min;
}


void
GPMSAOptions::set_scenario_parameter_scaling(unsigned int i,
                                             double range_min,
                                             double range_max)
{
  if (i >= m_scenarioScaleMin.size())
  {
    m_scenarioScaleMin.resize(i+1, 0);
    m_scenarioScaleRange.resize(i+1, 1);
  }
  m_scenarioScaleMin[i] = range_min;
  m_scenarioScaleRange[i] = range_max - range_min;
}


template <typename V>
void
GPMSAOptions::set_final_scaling
  (const std::vector<typename SharedPtr<V>::Type> & m_simulationScenarios,
   const std::vector<typename SharedPtr<V>::Type> & m_simulationParameters,
   const std::vector<typename SharedPtr<V>::Type> & m_simulationOutputs,
   const std::vector<typename SharedPtr<V>::Type> & m_experimentScenarios,
   const std::vector<typename SharedPtr<V>::Type> & m_experimentOutputs)
{
  if ((m_autoscaleMinMaxAll && m_autoscaleMeanVarAll) ||
      ((m_autoscaleMinMaxAll || m_autoscaleMeanVarAll) &&
       (!m_autoscaleMinMaxUncertain.empty() ||
        !m_autoscaleMeanVarUncertain.empty() ||
        !m_autoscaleMinMaxScenario.empty() ||
        !m_autoscaleMeanVarScenario.empty())))
    queso_error_msg("Cannot autoscale based on incompatible criteria");

  unsigned int dimScenario = m_simulationScenarios.size() ?
          m_simulationScenarios[0]->sizeGlobal() : 0;

  unsigned int dimParameter = m_simulationParameters.size() ?
          m_simulationParameters[0]->sizeGlobal() : 0;

  if (m_autoscaleMinMaxAll || !m_autoscaleMinMaxScenario.empty())
    {
      V maxScenario(*m_simulationScenarios[0]);

      V minScenario(*m_simulationScenarios[0]);

      for (unsigned int i=1; i < m_simulationScenarios.size(); ++i)
        min_max_update(minScenario, maxScenario,
                       (*m_simulationScenarios[i]));

      for (unsigned int i=0; i < m_experimentScenarios.size(); ++i)
        min_max_update(minScenario, maxScenario,
                       (*m_experimentScenarios[i]));

      for (unsigned int p=0; p != dimScenario; ++p)
        if (m_autoscaleMinMaxAll ||
            m_autoscaleMinMaxScenario.count(p))
          {
            if ((m_scenarioScaleMin.size() > p) &&
                ((m_scenarioScaleMin[p] != 0) ||
                 (m_scenarioScaleRange[p] != 1)))
              queso_error_msg("Cannot autoscale and manually scale the same scenario parameter");

            this->set_scenario_parameter_scaling(p, minScenario[p],
                                                    maxScenario[p]);
          }
    }

  if (m_autoscaleMinMaxAll || !m_autoscaleMinMaxUncertain.empty())
    {
      V maxUncertain(*m_simulationParameters[0]);

      V minUncertain(*m_simulationParameters[0]);

      for (unsigned int i=1; i < m_simulationParameters.size(); ++i)
        min_max_update(minUncertain, maxUncertain,
                       (*m_simulationParameters[i]));

      for (unsigned int p=0; p != dimParameter; ++p)
        if (m_autoscaleMinMaxAll ||
            m_autoscaleMinMaxUncertain.count(p))
          {
            if ((m_uncertainScaleMin.size() > p) &&
                ((m_uncertainScaleMin[p] != 0) ||
                 (m_uncertainScaleRange[p] != 1)))
              queso_error_msg("Cannot autoscale and manually scale the same uncertain parameter");

            this->set_uncertain_parameter_scaling(p, minUncertain[p],
                                                     maxUncertain[p]);
          }
    }


  if (m_autoscaleMeanVarAll || !m_autoscaleMeanVarScenario.empty())
    {
      unsigned int n=0;

      V meanScenario(*m_simulationScenarios[0]);

      V varScenario(m_simulationScenarios[0]->env(),
                    m_simulationScenarios[0]->map());

      for (unsigned int i=0; i < m_simulationScenarios.size(); ++i)
        mean_var_update(n, meanScenario, varScenario,
                        *m_simulationScenarios[i]);

      for (unsigned int i=0; i < m_experimentScenarios.size(); ++i)
        mean_var_update(n, meanScenario, varScenario,
                        *m_experimentScenarios[i]);

      varScenario /= n;

      for (unsigned int p=0; p != dimScenario; ++p)
        if (m_autoscaleMeanVarAll ||
            m_autoscaleMeanVarScenario.count(p))
          {
            if ((m_scenarioScaleMin.size() > p) &&
                ((m_scenarioScaleMin[p] != 0) ||
                 (m_scenarioScaleRange[p] != 1)))
              queso_error_msg("Cannot autoscale and manually scale the same scenario parameter");

            this->set_scenario_parameter_scaling(p, meanScenario[p],
                                                 meanScenario[p] +
                                                 std::sqrt(varScenario[p]));
          }
    }


  if (m_autoscaleMeanVarAll || !m_autoscaleMeanVarUncertain.empty())
    {
      unsigned int n=0;

      V meanUncertain(*m_simulationParameters[0]);

      V varUncertain(m_simulationParameters[0]->env(),
                     m_simulationParameters[0]->map());

      for (unsigned int i=0; i < m_simulationParameters.size(); ++i)
        mean_var_update(n, meanUncertain, varUncertain,
                        *m_simulationParameters[i]);

      varUncertain /= n;

      for (unsigned int p=0; p != dimScenario; ++p)
        if (m_autoscaleMeanVarAll ||
            m_autoscaleMeanVarUncertain.count(p))
          {
            if ((m_uncertainScaleMin.size() > p) &&
                ((m_uncertainScaleMin[p] != 0) ||
                 (m_uncertainScaleRange[p] != 1)))
              queso_error_msg("Cannot autoscale and manually scale the same uncertain parameter");

            this->set_uncertain_parameter_scaling(p, meanUncertain[p],
                                                  meanUncertain[p] +
                                                  std::sqrt(varUncertain[p]));
          }
    }
}


double
GPMSAOptions::normalized_scenario_parameter(unsigned int i,
                                            double physical_param)
const
{
  if (i < m_scenarioScaleMin.size())
    return (physical_param - m_scenarioScaleMin[i]) /
            m_scenarioScaleRange[i];
  return physical_param;
}


double
GPMSAOptions::normalized_uncertain_parameter(unsigned int i,
                                             double physical_param)
const
{
  if (i < m_uncertainScaleMin.size())
    return (physical_param - m_uncertainScaleMin[i]) /
            m_uncertainScaleRange[i];
  return physical_param;
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
     << "\n" << m_option_calibrateObservationalPrecision << " = " << this->m_calibrateObservationalPrecision
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
     << "\n" << m_option_autoscaleMinMaxAll << " = " << this->m_autoscaleMinMaxAll
     << "\n" << m_option_autoscaleMeanVarAll << " = " << this->m_autoscaleMeanVarAll
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



// Template instantiations
template
void
GPMSAOptions::set_final_scaling<GslVector>
  (const std::vector<SharedPtr<GslVector>::Type> &,
   const std::vector<SharedPtr<GslVector>::Type> &,
   const std::vector<SharedPtr<GslVector>::Type> &,
   const std::vector<SharedPtr<GslVector>::Type> &,
   const std::vector<SharedPtr<GslVector>::Type> &);


}  // End namespace QUESO
