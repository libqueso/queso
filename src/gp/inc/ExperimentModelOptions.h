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

#ifndef UQ_EXPERIMENT_MODEL_OPTIONS_H
#define UQ_EXPERIMENT_MODEL_OPTIONS_H

#include <queso/Environment.h>
#include <queso/BoostInputOptionsParser.h>

// _ODV = option default value
#define UQ_EXPERIMENT_MODEL_G_VALUES_ODV ""
#define UQ_EXPERIMENT_MODEL_A_V_ODV      0.
#define UQ_EXPERIMENT_MODEL_B_V_ODV      0.
#define UQ_EXPERIMENT_MODEL_A_RHO_V_ODV  0.
#define UQ_EXPERIMENT_MODEL_B_RHO_V_ODV  0.
#define UQ_EXPERIMENT_MODEL_A_Y_ODV      0.
#define UQ_EXPERIMENT_MODEL_B_Y_ODV      0.

namespace boost {
  namespace program_options {
    class options_description;
  }
}

namespace QUESO {

class EmOptionsValues
{
public:
  EmOptionsValues();
  EmOptionsValues(const BaseEnvironment * env, const char * prefix);
  EmOptionsValues(const EmOptionsValues& src);
  EmOptionsValues& operator=(const EmOptionsValues& rhs);
  virtual ~EmOptionsValues();

  std::string m_prefix;

  std::vector<unsigned int> m_Gvalues;
  double m_a_v;
  double m_b_v;
  double m_a_rho_v;
  double m_b_rho_v;
  double m_a_y;
  double m_b_y;

private:
  BoostInputOptionsParser * m_parser;

  std::string m_option_help;
  std::string m_option_Gvalues;
  std::string m_option_a_v;
  std::string m_option_b_v;
  std::string m_option_a_rho_v;
  std::string m_option_b_rho_v;
  std::string m_option_a_y;
  std::string m_option_b_y;

  void copy(const EmOptionsValues& src);
};

class ExperimentModelOptions
{
public:
  ExperimentModelOptions(const BaseEnvironment& env, const char* prefix);
  ExperimentModelOptions(const BaseEnvironment& env, const char* prefix, const EmOptionsValues& alternativeOptionsValues);
 ~ExperimentModelOptions();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  EmOptionsValues   m_ov;
  std::string              m_prefix;

private:
  void   defineMyOptions  (boost::program_options::options_description& optionsDesc) const;
  void   getMyOptionValues(boost::program_options::options_description& optionsDesc);

  const BaseEnvironment& m_env;

  boost::program_options::options_description* m_optionsDesc;
  std::string              m_option_help;
  std::string              m_option_Gvalues;
  std::string              m_option_a_v;
  std::string              m_option_b_v;
  std::string              m_option_a_rho_v;
  std::string              m_option_b_rho_v;
  std::string              m_option_a_y;
  std::string              m_option_b_y;
};

std::ostream& operator<<(std::ostream& os, const ExperimentModelOptions& obj);

}  // End namespace QUESO

#endif // UQ_EXPERIMENT_MODEL_OPTIONS_H
