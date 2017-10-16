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

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#else
#define GETPOT_NAMESPACE QUESO // So we don't clash with other getpots
#include <queso/getpot.h>
#undef GETPOT_NAMESPACE
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/GPMSAOptions.h>

#include <queso/GslVector.h>
#include <queso/SimulationOutputMesh.h>


// ODV = option default value
#define UQ_GPMSA_HELP ""
#define UQ_GPMSA_MAX_SIMULATOR_BASIS_VECTORS_ODV 0
#define UQ_GPMSA_SIMULATOR_BASIS_VARIANCE_TO_CAPTURE 1.0
#define UQ_GPMSA_TRUNCATION_ERROR_PRECISION_SHAPE_ODV 5.0
#define UQ_GPMSA_TRUNCATION_ERROR_PRECISION_SCALE_ODV 200.0
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
#define UQ_GPMSA_OBSERVATIONAL_PRECISION_RIDGE 1e-4
#define UQ_GPMSA_OBSERVATIONAL_COVARIANCE_RIDGE 0.0
#define UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE 1.0
static const bool UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC = false;
#define UQ_GPMSA_GAUSSIAN_DISCREPANCY_SUPPORT_THRESHOLD 0.05

namespace { // Anonymous namespace for helper functions

template <typename V>
void min_max_update(V & min, V & max, const V & new_data,
                    const char * warn_on_update = NULL)
{
  unsigned int dim = min.sizeGlobal();
  queso_assert_equal_to(dim, max.sizeGlobal());
  queso_assert_equal_to(dim, new_data.sizeGlobal());

  for (unsigned int p=0; p != dim; ++p)
    {
      if (warn_on_update)
        {
          if (min[p] > new_data[p])
            queso_warning("Experimental " << warn_on_update << " " <<
                          new_data[p] << " at index " << p <<
                          " is below minimum simulation " <<
                          warn_on_update << " " << min[p]);
          if (max[p] < new_data[p])
            queso_warning("Experimental " << warn_on_update << " " <<
                          new_data[p] << " at index " << p <<
                          " is above maximum simulation " <<
                          warn_on_update << " " << max[p]);
        }
      else
        {
          min[p] = std::min(min[p], new_data[p]);
          max[p] = std::max(max[p], new_data[p]);
        }
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
  const char * prefix) :
  options_have_been_used(false)
{
  this->set_defaults();
  this->parse(env, prefix);
}


GPMSAOptions::GPMSAOptions()
  :
  m_env(NULL),
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser(new BoostInputOptionsParser()),
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  options_have_been_used(false)
{
  this->set_defaults();
  this->set_prefix("");
}


void
GPMSAOptions::set_prefix(const char * prefix)
{
  queso_require(!options_have_been_used);

  m_prefix = std::string(prefix) + "gpmsa_";

  m_option_help = m_prefix + "help";
  m_option_maxEmulatorBasisVectors = m_prefix + "max_emulator_basis_vectors";
  m_option_emulatorBasisVarianceToCapture = m_prefix + "emulator_basis_variance_to_capture";
  m_option_truncationErrorPrecisionShape = m_prefix + "truncation_error_precision_shape";
  m_option_truncationErrorPrecisionScale = m_prefix + "truncation_error_precision_scale";
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
  m_option_observationalPrecisionRidge = m_prefix + "observational_precision_ridge";
  m_option_observationalCovarianceRidge = m_prefix + "observational_covariance_ridge";
  m_option_autoscaleMinMaxAll = m_prefix + "autoscale_min_max_all";
  m_option_autoscaleMeanVarAll = m_prefix + "autoscale_mean_var_all";
  m_option_gaussianDiscrepancyDistanceX = m_prefix + "gaussian_discrepancy_distance_x";
  m_option_gaussianDiscrepancyDistanceY = m_prefix + "gaussian_discrepancy_distance_y";
  m_option_gaussianDiscrepancyDistanceZ = m_prefix + "gaussian_discrepancy_distance_z";
  m_option_gaussianDiscrepancyDistanceT = m_prefix + "gaussian_discrepancy_distance_t";
  m_option_gaussianDiscrepancyPeriodicX = m_prefix + "gaussian_discrepancy_periodic_x";
  m_option_gaussianDiscrepancyPeriodicY = m_prefix + "gaussian_discrepancy_periodic_y";
  m_option_gaussianDiscrepancyPeriodicZ = m_prefix + "gaussian_discrepancy_periodic_z";
  m_option_gaussianDiscrepancyPeriodicT = m_prefix + "gaussian_discrepancy_periodic_t";
  m_option_gaussianDiscrepancySupportThreshold = m_prefix + "gaussian_discrepancy_support_threshold";
}



void
GPMSAOptions::set_defaults()
{
  queso_require(!options_have_been_used);

  m_help = UQ_GPMSA_HELP;
  m_maxEmulatorBasisVectors = UQ_GPMSA_MAX_SIMULATOR_BASIS_VECTORS_ODV;
  m_emulatorBasisVarianceToCapture = UQ_GPMSA_SIMULATOR_BASIS_VARIANCE_TO_CAPTURE;
  m_truncationErrorPrecisionShape = UQ_GPMSA_TRUNCATION_ERROR_PRECISION_SHAPE_ODV;
  m_truncationErrorPrecisionScale = UQ_GPMSA_TRUNCATION_ERROR_PRECISION_SCALE_ODV;
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
  m_observationalPrecisionRidge = UQ_GPMSA_OBSERVATIONAL_PRECISION_RIDGE;
  m_observationalCovarianceRidge = UQ_GPMSA_OBSERVATIONAL_COVARIANCE_RIDGE;

  m_autoscaleMinMaxAll = false;
  m_autoscaleMeanVarAll = false;

  m_gaussianDiscrepancySupportThreshold = UQ_GPMSA_GAUSSIAN_DISCREPANCY_SUPPORT_THRESHOLD;

  checkOptions();
}


void
GPMSAOptions::parse(const BaseEnvironment & env,
                    const char * prefix)
{
  queso_require(!options_have_been_used);

  m_env = &env;

  if (m_env->optionsInputFileName() == "") {
    queso_error_msg("Missing input file is required");
  }

  this->set_prefix(prefix);

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  m_parser.reset(new BoostInputOptionsParser(env.optionsInputFileName()));

  m_parser->registerOption<std::string>
    (m_option_help,
     m_help,
     "produce help message Gaussian process emulator");

  m_parser->registerOption
    (m_option_truncationErrorPrecisionShape,
     m_truncationErrorPrecisionShape,
     "shape hyperprior (Gamma) parameter for truncation error precision");
  m_parser->registerOption
    (m_option_truncationErrorPrecisionScale,
    m_truncationErrorPrecisionScale,
    "scale hyperprior (Gamma) parameter for truncation error precision");

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
    (m_option_observationalPrecisionRidge,
    m_observationalPrecisionRidge,
    "ridge to add to observational precision matrix");

  m_parser->registerOption
    (m_option_observationalCovarianceRidge,
    m_observationalCovarianceRidge,
    "ridge to add to observational covariance matrix");

  m_parser->registerOption
    (m_option_autoscaleMinMaxAll,
    m_autoscaleMinMaxAll,
    "option to autoscale all parameters and outputs based on data range");
  m_parser->registerOption
    (m_option_autoscaleMeanVarAll,
    m_autoscaleMeanVarAll,
    "option to autoscale all parameters and outputs based on data statistics");

  std::ostringstream defaultDiscrepancyDistanceString;
  defaultDiscrepancyDistanceString << UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE;

  m_parser->registerOption
    (m_option_gaussianDiscrepancyDistanceX,
    defaultDiscrepancyDistanceString.str(),
    "x distance between neighboring gaussian discrepancy kernels");

  m_parser->registerOption
    (m_option_gaussianDiscrepancyDistanceY,
    defaultDiscrepancyDistanceString.str(),
    "y distance between neighboring gaussian discrepancy kernels");

  m_parser->registerOption
    (m_option_gaussianDiscrepancyDistanceZ,
    defaultDiscrepancyDistanceString.str(),
    "z distance between neighboring gaussian discrepancy kernels");

  m_parser->registerOption
    (m_option_gaussianDiscrepancyDistanceT,
    defaultDiscrepancyDistanceString.str(),
    "t distance between neighboring gaussian discrepancy kernels");

  std::ostringstream defaultDiscrepancyPeriodicString;
  defaultDiscrepancyPeriodicString << std::boolalpha << UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC;

  m_parser->registerOption
    (m_option_gaussianDiscrepancyPeriodicX,
    defaultDiscrepancyPeriodicString.str(),
    "whether gaussian discrepancy kernels are periodic in x");

  m_parser->registerOption
    (m_option_gaussianDiscrepancyPeriodicY,
    defaultDiscrepancyPeriodicString.str(),
    "whether gaussian discrepancy kernels are periodic in y");

  m_parser->registerOption
    (m_option_gaussianDiscrepancyPeriodicZ,
    defaultDiscrepancyPeriodicString.str(),
    "whether gaussian discrepancy kernels are periodic in z");

  m_parser->registerOption
    (m_option_gaussianDiscrepancyPeriodicT,
    defaultDiscrepancyPeriodicString.str(),
    "whether gaussian discrepancy kernels are periodic in t");

  m_parser->registerOption
    (m_option_gaussianDiscrepancySupportThreshold,
    m_gaussianDiscrepancySupportThreshold,
    "threshold below which to omit out-of-domain centered gaussian kernels");


  m_parser->registerOption
    (m_option_maxEmulatorBasisVectors,
    m_maxEmulatorBasisVectors,
    "max number of basis vectors to use in SVD of simulation output");

  m_parser->scanInputFile();

  m_parser->getOption<std::string>(m_option_help,                           m_help);
  m_parser->getOption<double>(m_option_truncationErrorPrecisionShape,       m_truncationErrorPrecisionShape);
  m_parser->getOption<double>(m_option_truncationErrorPrecisionScale,       m_truncationErrorPrecisionScale);
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
  m_parser->getOption<double>(m_option_observationalPrecisionRidge,         m_observationalPrecisionRidge);
  m_parser->getOption<double>(m_option_observationalCovarianceRidge,        m_observationalCovarianceRidge);
  m_parser->getOption<bool>  (m_option_autoscaleMinMaxAll,                  m_autoscaleMinMaxAll);
  m_parser->getOption<bool>  (m_option_autoscaleMeanVarAll,                 m_autoscaleMeanVarAll);
  m_parser->getOption<int>   (m_option_maxEmulatorBasisVectors,             m_maxEmulatorBasisVectors);
  m_parser->getOption<std::vector<double> >(m_option_gaussianDiscrepancyDistanceX,        m_gaussianDiscrepancyDistanceX);
  m_parser->getOption<std::vector<double> >(m_option_gaussianDiscrepancyDistanceY,        m_gaussianDiscrepancyDistanceY);
  m_parser->getOption<std::vector<double> >(m_option_gaussianDiscrepancyDistanceZ,        m_gaussianDiscrepancyDistanceZ);
  m_parser->getOption<std::vector<double> >(m_option_gaussianDiscrepancyDistanceT,        m_gaussianDiscrepancyDistanceT);

  // I cannot get boost::program_options to play nicely with
  // vector<bool>, and we don't yet implement periodic discrepancy
  // options, so let's just ignore them here for now.  By the time we
  // implement them we'll have removed the deprecated boost::po code
  // anyway.
/*
  m_parser->getOption<std::vector<bool> > (m_option_gaussianDiscrepancyPeriodicX,        m_gaussianDiscrepancyPeriodicX);
  m_parser->getOption<std::vector<bool> > (m_option_gaussianDiscrepancyPeriodicY,        m_gaussianDiscrepancyPeriodicY);
  m_parser->getOption<std::vector<bool> > (m_option_gaussianDiscrepancyPeriodicZ,        m_gaussianDiscrepancyPeriodicZ);
  m_parser->getOption<std::vector<bool> > (m_option_gaussianDiscrepancyPeriodicT,        m_gaussianDiscrepancyPeriodicT);
*/
  m_parser->getOption<double>(m_option_gaussianDiscrepancySupportThreshold, m_gaussianDiscrepancySupportThreshold);
#else
  m_help = env.input()(m_option_help, UQ_GPMSA_HELP);

  m_truncationErrorPrecisionShape =
    env.input()(m_option_truncationErrorPrecisionShape,
                m_truncationErrorPrecisionShape);
  m_truncationErrorPrecisionScale =
    env.input()(m_option_truncationErrorPrecisionScale,
                m_truncationErrorPrecisionScale);

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

  m_observationalPrecisionRidge =
    env.input()(m_option_observationalPrecisionRidge,
                m_observationalPrecisionRidge);
  m_observationalCovarianceRidge =
    env.input()(m_option_observationalCovarianceRidge,
                m_observationalCovarianceRidge);

  m_autoscaleMinMaxAll =
    env.input()(m_option_autoscaleMinMaxAll,
                m_autoscaleMinMaxAll);
  m_autoscaleMeanVarAll =
    env.input()(m_option_autoscaleMeanVarAll,
                m_autoscaleMeanVarAll);

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyDistanceX);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyDistanceX.size() <= i)
        m_gaussianDiscrepancyDistanceX.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
      m_gaussianDiscrepancyDistanceX[i] =
        env.input()(m_option_gaussianDiscrepancyDistanceX,
                    m_gaussianDiscrepancyDistanceX[i]);
    }

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyDistanceY);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyDistanceY.size() <= i)
        m_gaussianDiscrepancyDistanceY.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
      m_gaussianDiscrepancyDistanceY[i] =
        env.input()(m_option_gaussianDiscrepancyDistanceY,
                    m_gaussianDiscrepancyDistanceY[i]);
    }

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyDistanceZ);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyDistanceZ.size() <= i)
        m_gaussianDiscrepancyDistanceZ.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
      m_gaussianDiscrepancyDistanceZ[i] =
        env.input()(m_option_gaussianDiscrepancyDistanceZ,
                    m_gaussianDiscrepancyDistanceZ[i]);
    }

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyDistanceT);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyDistanceT.size() <= i)
        m_gaussianDiscrepancyDistanceT.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
      m_gaussianDiscrepancyDistanceT[i] =
        env.input()(m_option_gaussianDiscrepancyDistanceT,
                    m_gaussianDiscrepancyDistanceT[i]);
    }


  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyPeriodicX);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyPeriodicX.size() <= i)
        m_gaussianDiscrepancyPeriodicX.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
      m_gaussianDiscrepancyPeriodicX[i] =
        env.input()(m_option_gaussianDiscrepancyPeriodicX,
                    bool(m_gaussianDiscrepancyPeriodicX[i]));
    }

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyPeriodicY);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyPeriodicY.size() <= i)
        m_gaussianDiscrepancyPeriodicY.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
      m_gaussianDiscrepancyPeriodicY[i] =
        env.input()(m_option_gaussianDiscrepancyPeriodicY,
                    bool(m_gaussianDiscrepancyPeriodicY[i]));
    }

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyPeriodicZ);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyPeriodicZ.size() <= i)
        m_gaussianDiscrepancyPeriodicZ.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
      m_gaussianDiscrepancyPeriodicZ[i] =
        env.input()(m_option_gaussianDiscrepancyPeriodicZ,
                    bool(m_gaussianDiscrepancyPeriodicZ[i]));
    }

  for (unsigned int i = 0,
       size = env.input().vector_variable_size(m_option_gaussianDiscrepancyPeriodicT);
       i != size; ++i)
    {
      if (m_gaussianDiscrepancyPeriodicT.size() <= i)
        m_gaussianDiscrepancyPeriodicT.push_back(UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
      m_gaussianDiscrepancyPeriodicT[i] =
        env.input()(m_option_gaussianDiscrepancyPeriodicT,
                    bool(m_gaussianDiscrepancyPeriodicT[i]));
    }

  m_gaussianDiscrepancySupportThreshold =
    env.input()(m_option_gaussianDiscrepancySupportThreshold,
                m_gaussianDiscrepancySupportThreshold);

  m_maxEmulatorBasisVectors =
    env.input()(m_option_maxEmulatorBasisVectors,
                m_maxEmulatorBasisVectors);
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  checkOptions();
}

GPMSAOptions::~GPMSAOptions()
{
}


void
GPMSAOptions::set_autoscale_minmax()
{
  queso_require(!options_have_been_used);

  this->m_autoscaleMinMaxAll = true;
}


void
GPMSAOptions::set_autoscale_minmax_uncertain_parameter(unsigned int i)
{
  queso_require(!options_have_been_used);

  m_autoscaleMinMaxUncertain.insert(i);
}


void
GPMSAOptions::set_autoscale_minmax_scenario_parameter(unsigned int i)
{
  queso_require(!options_have_been_used);

  m_autoscaleMinMaxScenario.insert(i);
}


void
GPMSAOptions::set_autoscale_minmax_output(unsigned int i)
{
  m_autoscaleMinMaxOutput.insert(i);
}


void
GPMSAOptions::set_autoscale_meanvar()
{
  queso_require(!options_have_been_used);

  this->m_autoscaleMeanVarAll = true;
}


void
GPMSAOptions::set_autoscale_meanvar_uncertain_parameter(unsigned int i)
{
  queso_require(!options_have_been_used);

  m_autoscaleMeanVarUncertain.insert(i);
}


void
GPMSAOptions::set_autoscale_meanvar_scenario_parameter(unsigned int i)
{
  queso_require(!options_have_been_used);

  m_autoscaleMeanVarScenario.insert(i);
}


void
GPMSAOptions::set_autoscale_meanvar_output(unsigned int i)
{
  m_autoscaleMeanVarOutput.insert(i);
}


void
GPMSAOptions::set_uncertain_parameter_scaling(unsigned int i,
                                              double range_min,
                                              double range_max)
{
  queso_require(!options_have_been_used);

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
  queso_require(!options_have_been_used);

  if (i >= m_scenarioScaleMin.size())
  {
    m_scenarioScaleMin.resize(i+1, 0);
    m_scenarioScaleRange.resize(i+1, 1);
  }
  m_scenarioScaleMin[i] = range_min;
  m_scenarioScaleRange[i] = range_max - range_min;
}


void
GPMSAOptions::set_output_scaling(unsigned int i,
                                 double range_min,
                                 double range_max)
{
  if (i >= m_outputScaleMin.size())
  {
    m_outputScaleMin.resize(i+1, 0);
    m_outputScaleRange.resize(i+1, 1);
  }
  m_outputScaleMin[i] = range_min;
  m_outputScaleRange[i] = range_max - range_min;
}


template <typename V>
void
GPMSAOptions::set_final_scaling
  (const std::vector<typename SharedPtr<V>::Type> & m_simulationScenarios,
   const std::vector<typename SharedPtr<V>::Type> & m_simulationParameters,
   const std::vector<typename SharedPtr<V>::Type> & m_simulationOutputs,
   const std::vector<typename SharedPtr<V>::Type> & m_experimentScenarios,
   const std::vector<typename SharedPtr<V>::Type> & m_experimentOutputs,
   const std::vector<typename SharedPtr<SimulationOutputMesh<V> >::Type> & m_simulationMeshes)
{
  if ((m_autoscaleMinMaxAll && m_autoscaleMeanVarAll) ||
      ((m_autoscaleMinMaxAll || m_autoscaleMeanVarAll) &&
       (!m_autoscaleMinMaxUncertain.empty() ||
        !m_autoscaleMeanVarUncertain.empty() ||
        !m_autoscaleMinMaxScenario.empty() ||
        !m_autoscaleMeanVarScenario.empty() ||
        !m_autoscaleMinMaxOutput.empty() ||
        !m_autoscaleMeanVarOutput.empty())))
    queso_error_msg("Cannot autoscale based on incompatible criteria");

  unsigned int dimScenario = m_simulationScenarios.size() ?
          m_simulationScenarios[0]->sizeGlobal() : 0;

  unsigned int dimParameter = m_simulationParameters.size() ?
          m_simulationParameters[0]->sizeGlobal() : 0;

  unsigned int dimOutput = m_simulationOutputs.size() ?
          m_simulationOutputs[0]->sizeGlobal() : 0;

  queso_require_msg(dimOutput, "GPMSA needs *some* output data");

  if (m_autoscaleMinMaxAll || !m_autoscaleMinMaxScenario.empty())
    {
      V maxScenario(*m_simulationScenarios[0]);

      V minScenario(*m_simulationScenarios[0]);

      // Only use simulation data for scaling.
      for (unsigned int i=1; i < m_simulationScenarios.size(); ++i)
        min_max_update(minScenario, maxScenario,
                       (*m_simulationScenarios[i]));

      // But check the experimental data, and yell at the user if it
      // falls outside the simulation data range, because they
      // probably don't want to be extrapolating with their GP.
      for (unsigned int i=0; i < m_experimentScenarios.size(); ++i)
        min_max_update(minScenario, maxScenario,
                       (*m_experimentScenarios[i]),
                       /*warn_on_update =*/ "scenario parameter");

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

  if (m_autoscaleMinMaxAll || !m_autoscaleMinMaxOutput.empty())
    {
      V maxOutput(*m_simulationOutputs[0]);

      V minOutput(*m_simulationOutputs[0]);

      // Only use simulation data for scaling.  In this case we won't
      // warn if the experimental data falls outside the simulation
      // data bounds, because the user can't control their outputs,
      // just their inputs.
      for (unsigned int i=1; i < m_simulationOutputs.size(); ++i)
        min_max_update(minOutput, maxOutput,
                       (*m_simulationOutputs[i]));

      // Scale together any functional output data, data that is
      // described by a simulation mesh, then scale multivariate
      // outputs independently.
      unsigned int first_multivariate_index = 0;
      for (unsigned int m=0; m != m_simulationMeshes.size(); ++m)
        {
          const SimulationOutputMesh<V> & mesh = *m_simulationMeshes[m];
          const unsigned int mesh_n_outputs = mesh.n_outputs();
          queso_assert_greater(mesh_n_outputs, 0);
          queso_assert_equal_to(mesh.first_solution_index(),
                                first_multivariate_index);
          first_multivariate_index += mesh_n_outputs;

          if (m_autoscaleMinMaxAll ||
              m_autoscaleMinMaxOutput.count(m))
            {
              if ((m_outputScaleMin.size() > m) &&
                  ((m_outputScaleMin[m] != 0) ||
                   (m_outputScaleRange[m] != 1)))
                queso_error_msg("Cannot autoscale and manually scale the same output data");

              // FIXME - this will break if/when we introduce
              // non-Lagrange solution bases.
              double full_min = minOutput[mesh.first_solution_index()],
                     full_max = maxOutput[mesh.first_solution_index()];
              for (unsigned int i=1; i != mesh_n_outputs; ++i)
                {
                  full_min = std::min(full_min, minOutput[mesh.first_solution_index()+i]);
                  full_max = std::max(full_max, maxOutput[mesh.first_solution_index()+i]);
                }
              this->set_output_scaling(m, full_min, full_max);
            }
        }

      const unsigned int n_multivariate_indices =
        dimOutput - first_multivariate_index;
      const unsigned int n_variables = m_simulationMeshes.size() + n_multivariate_indices;
      const unsigned int variable_index_to_output_index = first_multivariate_index - m_simulationMeshes.size();
      for (unsigned int p=m_simulationMeshes.size(); p != n_variables; ++p)
        {
          if (m_autoscaleMinMaxAll ||
              m_autoscaleMinMaxOutput.count(p))
            {
              if ((m_outputScaleMin.size() > p) &&
                  ((m_outputScaleMin[p] != 0) ||
                   (m_outputScaleRange[p] != 1)))
                queso_error_msg("Cannot autoscale and manually scale the same output data");

              this->set_output_scaling
                (p, minOutput[p+variable_index_to_output_index],
                    maxOutput[p+variable_index_to_output_index]);
            }
        }
    }


  if (m_autoscaleMeanVarAll || !m_autoscaleMeanVarScenario.empty())
    {
      unsigned int n=0;

      V meanScenario(*m_simulationScenarios[0]);

      V varScenario(m_simulationScenarios[0]->env(),
                    m_simulationScenarios[0]->map());

      // For consistency with the min-max behavior, only normalize
      // using simulation data, not experiment data.
      for (unsigned int i=0; i < m_simulationScenarios.size(); ++i)
        mean_var_update(n, meanScenario, varScenario,
                        *m_simulationScenarios[i]);

      varScenario /= n - 1;

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

      varUncertain /= n - 1;

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

  if (m_autoscaleMeanVarAll || !m_autoscaleMeanVarOutput.empty())
    {
      unsigned int n=0;

      V meanOutput(*m_simulationOutputs[0]);

      V varOutput(m_simulationOutputs[0]->env(),
                  m_simulationOutputs[0]->map());

      // This gives us the right mean and var values for any
      // multivariate components of the problem; we'll handle
      // functional components separately shortly.
      for (unsigned int i=0; i < m_simulationOutputs.size(); ++i)
        mean_var_update(n, meanOutput, varOutput,
                        *m_simulationOutputs[i]);

      varOutput /= n - 1;

      // Scale together any functional output data, data that is
      // described by a simulation mesh, then scale multivariate
      // outputs independently.
      unsigned int first_multivariate_index = 0;
      for (unsigned int m=0; m != m_simulationMeshes.size(); ++m)
        {
          const SimulationOutputMesh<V> & mesh = *m_simulationMeshes[m];
          const unsigned int mesh_n_outputs = mesh.n_outputs();
          queso_assert_greater(mesh_n_outputs, 0);
          queso_assert_equal_to(mesh.first_solution_index(),
                                first_multivariate_index);
          first_multivariate_index += mesh_n_outputs;

          if (m_autoscaleMeanVarAll ||
              m_autoscaleMeanVarOutput.count(m))
            {
              if ((m_outputScaleMin.size() > m) &&
                  ((m_outputScaleMin[m] != 0) ||
                   (m_outputScaleRange[m] != 1)))
                queso_error_msg("Cannot autoscale and manually scale the same output data");

              unsigned int n = 0;
              double full_mean = 0, full_var = 0;
              for (unsigned int i=0; i < m_simulationOutputs.size(); ++i)
                for (unsigned int j=0; j != mesh_n_outputs; ++j)
                  mean_var_update
                    (n, full_mean, full_var,
                     (*m_simulationOutputs[i])[mesh.first_solution_index()+i]);

              full_var /= (n - 1);

              this->set_output_scaling
                (m, full_mean, full_mean + std::sqrt(full_var));
            }
        }

      const unsigned int n_multivariate_indices =
        dimOutput - first_multivariate_index;
      const unsigned int n_variables = m_simulationMeshes.size() + n_multivariate_indices;
      const unsigned int variable_index_to_output_index = first_multivariate_index - m_simulationMeshes.size();
      for (unsigned int p=m_simulationMeshes.size(); p != n_variables; ++p)
        if (m_autoscaleMeanVarAll ||
            m_autoscaleMeanVarOutput.count(p))
          {
            if ((m_outputScaleMin.size() > p) &&
                ((m_outputScaleMin[p] != 0) ||
                 (m_outputScaleRange[p] != 1)))
              queso_error_msg("Cannot autoscale and manually scale the same output data");

            this->set_output_scaling
              (p, meanOutput[p+variable_index_to_output_index],
               meanOutput[p+variable_index_to_output_index] +
               std::sqrt(varOutput[p+variable_index_to_output_index]));
          }
    }

  m_output_index_to_variable_index.resize(dimOutput);

  unsigned int next_variable_begin = 0;
  for (unsigned int m=0; m != m_simulationMeshes.size(); ++m)
    {
      const SimulationOutputMesh<V> & mesh = *m_simulationMeshes[m];
      const unsigned int mesh_n_outputs = mesh.n_outputs();
      next_variable_begin += mesh_n_outputs;
      for (unsigned int i = mesh.first_solution_index();
           i != next_variable_begin; ++i)
        m_output_index_to_variable_index[i] = m;
    }

  unsigned int next_mv_index = m_simulationMeshes.size();
  for (unsigned int i = next_variable_begin; i != dimOutput; ++i)
    m_output_index_to_variable_index[i] = next_mv_index++;

  // Make sure we have gaussian discrepancy option values available
  // for every mesh, even if they weren't set by the input file.
  const unsigned int n_meshes = m_simulationMeshes.size();
  m_gaussianDiscrepancyDistanceX.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
  m_gaussianDiscrepancyDistanceY.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
  m_gaussianDiscrepancyDistanceZ.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
  m_gaussianDiscrepancyDistanceT.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_DISTANCE);
  m_gaussianDiscrepancyPeriodicX.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
  m_gaussianDiscrepancyPeriodicY.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
  m_gaussianDiscrepancyPeriodicZ.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);
  m_gaussianDiscrepancyPeriodicT.resize
    (n_meshes, UQ_GPMSA_GAUSSIAN_DISCREPANCY_PERIODIC);

  options_have_been_used = true;
}


double
GPMSAOptions::normalized_scenario_parameter(unsigned int i,
                                            double physical_param)
const
{
  if (i < m_scenarioScaleMin.size())
    return (physical_param - m_scenarioScaleMin[i]) /
            (m_scenarioScaleRange[i] ? m_scenarioScaleRange[i] : 1);
  return physical_param;
}


double
GPMSAOptions::normalized_uncertain_parameter(unsigned int i,
                                             double physical_param)
const
{
  if (i < m_uncertainScaleMin.size())
    return (physical_param - m_uncertainScaleMin[i]) /
            (m_uncertainScaleRange[i] ? m_uncertainScaleRange[i] : 1);
  return physical_param;
}


double
GPMSAOptions::normalized_output(unsigned int i,
                                double output_data)
const
{
  return this->normalized_output_variable
    (m_output_index_to_variable_index[i], output_data);
}



double
GPMSAOptions::normalized_output_variable(unsigned int i,
                                         double output_data)
const
{
  if (i < m_outputScaleMin.size())
    return (output_data - m_outputScaleMin[i]) /
            (m_outputScaleRange[i] ? m_outputScaleRange[i] : 1);
  return output_data;
}



double
GPMSAOptions::output_scale(unsigned int i)
const
{
  return this->output_scale_variable
    (m_output_index_to_variable_index[i]);
}



double
GPMSAOptions::output_scale_variable(unsigned int i)
const
{
  if (i < m_outputScaleRange.size() &&
      m_outputScaleRange[i])
    return m_outputScaleRange[i];
  return 1;
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
  os << "\n" << m_option_truncationErrorPrecisionShape << " = " << this->m_truncationErrorPrecisionShape
     << "\n" << m_option_truncationErrorPrecisionScale << " = " << this->m_truncationErrorPrecisionScale
     << "\n" << m_option_emulatorPrecisionShape << " = " << this->m_emulatorPrecisionShape
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
     << "\n" << m_option_observationalPrecisionRidge << " = " << this->m_observationalPrecisionRidge
     << "\n" << m_option_observationalCovarianceRidge << " = " << this->m_observationalCovarianceRidge
     << "\n" << m_option_autoscaleMinMaxAll << " = " << this->m_autoscaleMinMaxAll
     << "\n" << m_option_autoscaleMeanVarAll << " = " << this->m_autoscaleMeanVarAll
     << "\n" << m_option_maxEmulatorBasisVectors << " = " << this->m_maxEmulatorBasisVectors;

     os << "\n" << m_option_gaussianDiscrepancyDistanceX << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyDistanceX.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyDistanceX[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
     os << "\n" << m_option_gaussianDiscrepancyDistanceY << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyDistanceY.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyDistanceY[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
     os << "\n" << m_option_gaussianDiscrepancyDistanceZ << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyDistanceZ.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyDistanceZ[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
     os << "\n" << m_option_gaussianDiscrepancyDistanceT << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyDistanceT.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyDistanceT[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';

     os << "\n" << m_option_gaussianDiscrepancyPeriodicX << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyPeriodicX.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyPeriodicX[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
     os << "\n" << m_option_gaussianDiscrepancyPeriodicY << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyPeriodicY.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyPeriodicY[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
     os << "\n" << m_option_gaussianDiscrepancyPeriodicZ << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyPeriodicZ.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyPeriodicZ[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
     os << "\n" << m_option_gaussianDiscrepancyPeriodicT << " = {";
     for (unsigned int i = 0, size = this->m_gaussianDiscrepancyPeriodicT.size(); i != size; ++i)
       {
         os << this->m_gaussianDiscrepancyPeriodicT[i];
         if (i+1 != size)
           os << ',';
       }
     os << '}';
  os << "\n" << m_option_gaussianDiscrepancySupportThreshold << " = " << this->m_gaussianDiscrepancySupportThreshold
     << std::endl;
}

std::ostream &
operator<<(std::ostream& os, const GPMSAOptions & obj)
{
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  os << (*(obj.m_parser)) << std::endl;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
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
   const std::vector<SharedPtr<GslVector>::Type> &,
   const std::vector<typename SharedPtr<SimulationOutputMesh<GslVector> >::Type> &);


}  // End namespace QUESO
