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

#ifndef UQ_SIMULATION_MODEL_OPTIONS_H
#define UQ_SIMULATION_MODEL_OPTIONS_H

#include <queso/Environment.h>
#include <queso/BoostInputOptionsParser.h>

#define UQ_SIMULATION_MODEL_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_SIMULATION_MODEL_DATA_OUTPUT_FILE_NAME_ODV        UQ_SIMULATION_MODEL_FILENAME_FOR_NO_FILE
#define UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOW_ALL_ODV        0
#define UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOWED_SET_ODV      ""
#define UQ_SIMULATION_MODEL_P_ETA_ODV                        1
#define UQ_SIMULATION_MODEL_ZERO_RELATIVE_SINGULAR_VALUE_ODV 0.
#define UQ_SIMULATION_MODEL_CDF_THRESHOLD_FOR_P_ETA_ODV      1.
#define UQ_SIMULATION_MODEL_A_W_ODV                          0.
#define UQ_SIMULATION_MODEL_B_W_ODV                          0.
#define UQ_SIMULATION_MODEL_A_RHO_W_ODV                      0.
#define UQ_SIMULATION_MODEL_B_RHO_W_ODV                      0.
#define UQ_SIMULATION_MODEL_A_ETA_ODV                        0.
#define UQ_SIMULATION_MODEL_B_ETA_ODV                        0.
#define UQ_SIMULATION_MODEL_A_S_ODV                          0.
#define UQ_SIMULATION_MODEL_B_S_ODV                          0.

namespace boost {
  namespace program_options {
    class options_description;
  }
}

namespace QUESO {

class SmOptionsValues
{
public:
  SmOptionsValues();
  SmOptionsValues(const BaseEnvironment * env, const char * prefix);
  SmOptionsValues(const SmOptionsValues& src);
  SmOptionsValues& operator=(const SmOptionsValues& rhs);
  virtual ~SmOptionsValues();

  std::string m_prefix;

  std::string m_dataOutputFileName;
  bool m_dataOutputAllowAll;
  std::set<unsigned int> m_dataOutputAllowedSet;
  unsigned int m_p_eta;
  double m_zeroRelativeSingularValue;
  double m_cdfThresholdForPEta;
  double m_a_w;
  double m_b_w;
  double m_a_rho_w;
  double m_b_rho_w;
  double m_a_eta;
  double m_b_eta;
  double m_a_s;
  double m_b_s;

private:
  BoostInputOptionsParser * m_parser;

  std::string m_option_help;
  std::string m_option_dataOutputFileName;
  std::string m_option_dataOutputAllowAll;
  std::string m_option_dataOutputAllowedSet;
  std::string m_option_p_eta;
  std::string m_option_zeroRelativeSingularValue;
  std::string m_option_cdfThresholdForPEta;
  std::string m_option_a_w;
  std::string m_option_b_w;
  std::string m_option_a_rho_w;
  std::string m_option_b_rho_w;
  std::string m_option_a_eta;
  std::string m_option_b_eta;
  std::string m_option_a_s;
  std::string m_option_b_s;

  void copy(const SmOptionsValues& src);
};

////////////////////////////////////////////////////////////////////////////////////////////

class SimulationModelOptions
{
public:
  SimulationModelOptions(const BaseEnvironment& env, const char* prefix);
  SimulationModelOptions(const BaseEnvironment& env, const char* prefix, const SmOptionsValues& alternativeOptionsValues);
 ~SimulationModelOptions();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  SmOptionsValues        m_ov;
  std::string                   m_prefix;

private:
  void   defineMyOptions  (boost::program_options::options_description& optionsDesc) const;
  void   getMyOptionValues(boost::program_options::options_description& optionsDesc);

  const BaseEnvironment& m_env;

  boost::program_options::options_description*      m_optionsDesc;
  std::string                   m_option_help;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowAll;
  std::string                   m_option_dataOutputAllowedSet;
  std::string                   m_option_p_eta;
  std::string                   m_option_zeroRelativeSingularValue;
  std::string                   m_option_cdfThresholdForPEta;
  std::string                   m_option_a_w;
  std::string                   m_option_b_w;
  std::string                   m_option_a_rho_w;
  std::string                   m_option_b_rho_w;
  std::string                   m_option_a_eta;
  std::string                   m_option_b_eta;
  std::string                   m_option_a_s;
  std::string                   m_option_b_s;
};

std::ostream& operator<<(std::ostream& os, const SimulationModelOptions& obj);

}  // End namespace QUESO

#endif // UQ_SIMULATION_MODEL_OPTIONS_H
