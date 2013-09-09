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
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqExperimentModelOptions.h>
#include <uqMiscellaneous.h>

namespace QUESO {

uqEmOptionsValuesClass::uqEmOptionsValuesClass()
  :
  m_Gvalues(0),
  m_a_v    (UQ_EXPERIMENT_MODEL_A_V_ODV    ),
  m_b_v    (UQ_EXPERIMENT_MODEL_B_V_ODV    ),
  m_a_rho_v(UQ_EXPERIMENT_MODEL_A_RHO_V_ODV),
  m_b_rho_v(UQ_EXPERIMENT_MODEL_B_RHO_V_ODV),
  m_a_y    (UQ_EXPERIMENT_MODEL_A_Y_ODV    ),
  m_b_y    (UQ_EXPERIMENT_MODEL_B_Y_ODV    )
{
}

uqEmOptionsValuesClass::~uqEmOptionsValuesClass()
{
}

uqEmOptionsValuesClass::uqEmOptionsValuesClass(const uqEmOptionsValuesClass& src)
{
  this->copy(src);
}

uqEmOptionsValuesClass&
uqEmOptionsValuesClass::operator=(const uqEmOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}

void
uqEmOptionsValuesClass::copy(const uqEmOptionsValuesClass& src)
{
  m_Gvalues = src.m_Gvalues;
  m_a_v     = src.m_a_v;
  m_b_v     = src.m_b_v;
  m_a_rho_v = src.m_a_rho_v;
  m_b_rho_v = src.m_b_rho_v;
  m_a_y     = src.m_a_y;
  m_b_y     = src.m_b_y;

  return;
}

uqExperimentModelOptionsClass::uqExperimentModelOptionsClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_ov            (),
  m_prefix        ((std::string)(prefix) + "em_"),
  m_env           (env),
  m_optionsDesc   (new po::options_description("Experiment model options")),
  m_option_help   (m_prefix + "help"   ),
  m_option_Gvalues(m_prefix + "Gvalues"),
  m_option_a_v    (m_prefix + "a_v"    ),
  m_option_b_v    (m_prefix + "b_v"    ),
  m_option_a_rho_v(m_prefix + "a_rho_v"),
  m_option_b_rho_v(m_prefix + "b_rho_v"),
  m_option_a_y    (m_prefix + "a_y"    ),
  m_option_b_y    (m_prefix + "b_y"    )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "uqExperimentModelOptionsClass::constructor(1)",
                      "this constructor is incompatible with the abscense of an options input file");
}

uqExperimentModelOptionsClass::uqExperimentModelOptionsClass(
  const uqBaseEnvironmentClass&  env,
  const char*                    prefix,
  const uqEmOptionsValuesClass& alternativeOptionsValues)
  :
  m_ov            (alternativeOptionsValues),
  m_prefix        ((std::string)(prefix) + "em_"),
  m_env           (env),
  m_optionsDesc   (NULL),
  m_option_help   (m_prefix + "help"   ),
  m_option_Gvalues(m_prefix + "Gvalues"),
  m_option_a_v    (m_prefix + "a_v"    ),
  m_option_b_v    (m_prefix + "b_v"    ),
  m_option_a_rho_v(m_prefix + "a_rho_v"),
  m_option_b_rho_v(m_prefix + "b_rho_v"),
  m_option_a_y    (m_prefix + "a_y"    ),
  m_option_b_y    (m_prefix + "b_y"    )
{
  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "uqExperimentModelOptionsClass::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqExperimentModelOptionsClass::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

uqExperimentModelOptionsClass::~uqExperimentModelOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
} 

void
uqExperimentModelOptionsClass::scanOptionsValues()
{
  UQ_FATAL_TEST_MACRO(m_optionsDesc == NULL,
                      m_env.worldRank(),
                      "uqExperimentModelOptionsClass::scanOptionsValues()",
                      "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In uqExperimentModelOptionsClass::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
uqExperimentModelOptionsClass::defineMyOptions(po::options_description& optionsDesc) const
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                "produce help message for experiment model options")
    (m_option_Gvalues.c_str(), po::value<std::string >()->default_value(UQ_EXPERIMENT_MODEL_G_VALUES_ODV), "G values"                                         )
    (m_option_a_v.c_str(),     po::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_V_ODV     ), "a_v"                                              )
    (m_option_b_v.c_str(),     po::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_V_ODV     ), "b_v"                                              )
    (m_option_a_rho_v.c_str(), po::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_RHO_V_ODV ), "a_rho_v"                                          )
    (m_option_b_rho_v.c_str(), po::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_RHO_V_ODV ), "b_rho_v"                                          )
    (m_option_a_y.c_str(),     po::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_Y_ODV     ), "a_y"                                              )
    (m_option_b_y.c_str(),     po::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_Y_ODV     ), "b_y"                                              )
  ;

  return;
}

void
uqExperimentModelOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  std::vector<double> tmpValues(0,0.);
  if (m_env.allOptionsMap().count(m_option_Gvalues)) {
    std::string inputString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_Gvalues]).as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpValues);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqExperimentModelOptionsClass::getMyOptionValues(): tmpValues =";
    //  for (unsigned int i = 0; i < tmpValues.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpValues[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}
    unsigned int tmpSize = tmpValues.size();
    m_ov.m_Gvalues.clear();
    m_ov.m_Gvalues.resize(tmpSize,0);
    for (unsigned int i = 0; i < tmpSize; ++i) {
      m_ov.m_Gvalues[i] = (unsigned int) tmpValues[i];
    }
  }

  if (m_env.allOptionsMap().count(m_option_a_v)) {
    m_ov.m_a_v = ((const po::variable_value&) m_env.allOptionsMap()[m_option_a_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_v)) {
    m_ov.m_b_v = ((const po::variable_value&) m_env.allOptionsMap()[m_option_b_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_rho_v)) {
    m_ov.m_a_rho_v = ((const po::variable_value&) m_env.allOptionsMap()[m_option_a_rho_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_rho_v)) {
    m_ov.m_b_rho_v = ((const po::variable_value&) m_env.allOptionsMap()[m_option_b_rho_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_y)) {
    m_ov.m_a_y = ((const po::variable_value&) m_env.allOptionsMap()[m_option_a_y]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_y)) {
    m_ov.m_b_y = ((const po::variable_value&) m_env.allOptionsMap()[m_option_b_y]).as<double>();
  }

  return;
}

void
uqExperimentModelOptionsClass::print(std::ostream& os) const
{
  os << "\n" << m_option_Gvalues << " = ";
  for (unsigned int i = 0; i < m_ov.m_Gvalues.size(); ++i) {
    os << m_ov.m_Gvalues[i] << " ";
  }
  os << "\n" << m_option_a_v     << " = " << m_ov.m_a_v
     << "\n" << m_option_b_v     << " = " << m_ov.m_b_v
     << "\n" << m_option_a_rho_v << " = " << m_ov.m_a_rho_v
     << "\n" << m_option_b_rho_v << " = " << m_ov.m_b_rho_v
     << "\n" << m_option_a_y     << " = " << m_ov.m_a_y
     << "\n" << m_option_b_y     << " = " << m_ov.m_b_y
     << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqExperimentModelOptionsClass& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO
