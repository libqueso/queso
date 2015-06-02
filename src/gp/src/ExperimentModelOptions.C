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

#include <queso/ExperimentModelOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

EmOptionsValues::EmOptionsValues()
  :
    m_prefix("em_"),
    m_Gvalues(0),
    m_a_v(UQ_EXPERIMENT_MODEL_A_V_ODV),
    m_b_v(UQ_EXPERIMENT_MODEL_B_V_ODV),
    m_a_rho_v(UQ_EXPERIMENT_MODEL_A_RHO_V_ODV),
    m_b_rho_v(UQ_EXPERIMENT_MODEL_B_RHO_V_ODV),
    m_a_y(UQ_EXPERIMENT_MODEL_A_Y_ODV),
    m_b_y(UQ_EXPERIMENT_MODEL_B_Y_ODV),
    m_parser(NULL),
    m_option_help(m_prefix + "help"),
    m_option_Gvalues(m_prefix + "Gvalues"),
    m_option_a_v(m_prefix + "a_v"),
    m_option_b_v(m_prefix + "b_v"),
    m_option_a_rho_v(m_prefix + "a_rho_v"),
    m_option_b_rho_v(m_prefix + "b_rho_v"),
    m_option_a_y(m_prefix + "a_y"),
    m_option_b_y(m_prefix + "b_y")
{
}

EmOptionsValues::EmOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
    m_prefix((std::string)(prefix) + "em_"),
    m_Gvalues(0),
    m_a_v(UQ_EXPERIMENT_MODEL_A_V_ODV),
    m_b_v(UQ_EXPERIMENT_MODEL_B_V_ODV),
    m_a_rho_v(UQ_EXPERIMENT_MODEL_A_RHO_V_ODV),
    m_b_rho_v(UQ_EXPERIMENT_MODEL_B_RHO_V_ODV),
    m_a_y(UQ_EXPERIMENT_MODEL_A_Y_ODV),
    m_b_y(UQ_EXPERIMENT_MODEL_B_Y_ODV),
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
    m_option_help(m_prefix + "help"),
    m_option_Gvalues(m_prefix + "Gvalues"),
    m_option_a_v(m_prefix + "a_v"),
    m_option_b_v(m_prefix + "b_v"),
    m_option_a_rho_v(m_prefix + "a_rho_v"),
    m_option_b_rho_v(m_prefix + "b_rho_v"),
    m_option_a_y(m_prefix + "a_y"),
    m_option_b_y(m_prefix + "b_y")
{
  m_parser->registerOption(m_option_help,                                                                                "produce help message for experiment model options");
  m_parser->registerOption<std::string >(m_option_Gvalues, UQ_EXPERIMENT_MODEL_G_VALUES_ODV, "G values"                                         );
  m_parser->registerOption<double      >(m_option_a_v,     UQ_EXPERIMENT_MODEL_A_V_ODV     , "a_v"                                              );
  m_parser->registerOption<double      >(m_option_b_v,     UQ_EXPERIMENT_MODEL_B_V_ODV     , "b_v"                                              );
  m_parser->registerOption<double      >(m_option_a_rho_v, UQ_EXPERIMENT_MODEL_A_RHO_V_ODV , "a_rho_v"                                          );
  m_parser->registerOption<double      >(m_option_b_rho_v, UQ_EXPERIMENT_MODEL_B_RHO_V_ODV , "b_rho_v"                                          );
  m_parser->registerOption<double      >(m_option_a_y,     UQ_EXPERIMENT_MODEL_A_Y_ODV     , "a_y"                                              );
  m_parser->registerOption<double      >(m_option_b_y,     UQ_EXPERIMENT_MODEL_B_Y_ODV     , "b_y"                                              );

  m_parser->scanInputFile();

  m_parser->getOption<std::vector<unsigned int> >(m_option_Gvalues, m_Gvalues);
  m_parser->getOption<double      >(m_option_a_v,     m_a_v);
  m_parser->getOption<double      >(m_option_b_v,     m_b_v);
  m_parser->getOption<double      >(m_option_a_rho_v, m_a_rho_v);
  m_parser->getOption<double      >(m_option_b_rho_v, m_b_rho_v);
  m_parser->getOption<double      >(m_option_a_y,     m_a_y);
  m_parser->getOption<double      >(m_option_b_y,     m_b_y);
}

EmOptionsValues::~EmOptionsValues()
{
}

EmOptionsValues::EmOptionsValues(const EmOptionsValues& src)
{
  this->copy(src);
}

EmOptionsValues&
EmOptionsValues::operator=(const EmOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}

// void
// EmOptionsValues::defineOptions()
// {
//   (*m_optionsDescription).add_options()
//     (m_option_help.c_str(),                                                                                "produce help message for experiment model options")
//     (m_option_Gvalues.c_str(), boost::program_options::value<std::string >()->default_value(UQ_EXPERIMENT_MODEL_G_VALUES_ODV), "G values"                                         )
//     (m_option_a_v.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_V_ODV     ), "a_v"                                              )
//     (m_option_b_v.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_V_ODV     ), "b_v"                                              )
//     (m_option_a_rho_v.c_str(), boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_RHO_V_ODV ), "a_rho_v"                                          )
//     (m_option_b_rho_v.c_str(), boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_RHO_V_ODV ), "b_rho_v"                                          )
//     (m_option_a_y.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_Y_ODV     ), "a_y"                                              )
//     (m_option_b_y.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_Y_ODV     ), "b_y"                                              )
//   ;
// }
//
// void
// EmOptionsValues::getOptionValues()
// {
//   if ((*m_optionsMap).count(m_option_help)) {
//     if (m_env->subDisplayFile()) {
//       *m_env->subDisplayFile() << *m_optionsDescription
//                               << std::endl;
//     }
//   }
//
//   std::vector<double> tmpValues(0,0.);
//   if ((*m_optionsMap).count(m_option_Gvalues)) {
//     std::string inputString = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_Gvalues]).as<std::string>();
//     MiscReadDoublesFromString(inputString,tmpValues);
//     //if (m_env->subDisplayFile()) {
//     //  *m_env->subDisplayFile() << "In ExperimentModelOptions::getMyOptionValues(): tmpValues =";
//     //  for (unsigned int i = 0; i < tmpValues.size(); ++i) {
//     //    *m_env->subDisplayFile() << " " << tmpValues[i];
//     //  }
//     //  *m_env->subDisplayFile() << std::endl;
//     //}
//     unsigned int tmpSize = tmpValues.size();
//     m_Gvalues.clear();
//     m_Gvalues.resize(tmpSize,0);
//     for (unsigned int i = 0; i < tmpSize; ++i) {
//       m_Gvalues[i] = (unsigned int) tmpValues[i];
//     }
//   }
//
//   if ((*m_optionsMap).count(m_option_a_v)) {
//     m_a_v = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_v]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_v)) {
//     m_b_v = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_v]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_a_rho_v)) {
//     m_a_rho_v = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_rho_v]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_rho_v)) {
//     m_b_rho_v = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_rho_v]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_a_y)) {
//     m_a_y = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_y]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_y)) {
//     m_b_y = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_y]).as<double>();
//   }
// }

void
EmOptionsValues::copy(const EmOptionsValues& src)
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

ExperimentModelOptions::ExperimentModelOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_ov            (),
  m_prefix        ((std::string)(prefix) + "em_"),
  m_env           (env),
  m_optionsDesc   (new boost::program_options::options_description("Experiment model options")),
  m_option_help   (m_prefix + "help"   ),
  m_option_Gvalues(m_prefix + "Gvalues"),
  m_option_a_v    (m_prefix + "a_v"    ),
  m_option_b_v    (m_prefix + "b_v"    ),
  m_option_a_rho_v(m_prefix + "a_rho_v"),
  m_option_b_rho_v(m_prefix + "b_rho_v"),
  m_option_a_y    (m_prefix + "a_y"    ),
  m_option_b_y    (m_prefix + "b_y"    )
{
  queso_deprecated();

  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the abscense of an options input file");
}

ExperimentModelOptions::ExperimentModelOptions(
  const BaseEnvironment&  env,
  const char*                    prefix,
  const EmOptionsValues& alternativeOptionsValues)
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
  queso_deprecated();

  queso_require_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In ExperimentModelOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

ExperimentModelOptions::~ExperimentModelOptions()
{
  queso_deprecated();

  if (m_optionsDesc) delete m_optionsDesc;
}

void
ExperimentModelOptions::scanOptionsValues()
{
  queso_deprecated();

  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In ExperimentModelOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
ExperimentModelOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                "produce help message for experiment model options")
    (m_option_Gvalues.c_str(), boost::program_options::value<std::string >()->default_value(UQ_EXPERIMENT_MODEL_G_VALUES_ODV), "G values"                                         )
    (m_option_a_v.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_V_ODV     ), "a_v"                                              )
    (m_option_b_v.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_V_ODV     ), "b_v"                                              )
    (m_option_a_rho_v.c_str(), boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_RHO_V_ODV ), "a_rho_v"                                          )
    (m_option_b_rho_v.c_str(), boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_RHO_V_ODV ), "b_rho_v"                                          )
    (m_option_a_y.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_A_Y_ODV     ), "a_y"                                              )
    (m_option_b_y.c_str(),     boost::program_options::value<double      >()->default_value(UQ_EXPERIMENT_MODEL_B_Y_ODV     ), "b_y"                                              )
  ;

  return;
}

void
ExperimentModelOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  std::vector<double> tmpValues(0,0.);
  if (m_env.allOptionsMap().count(m_option_Gvalues)) {
    std::string inputString = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_Gvalues]).as<std::string>();
    MiscReadDoublesFromString(inputString,tmpValues);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In ExperimentModelOptions::getMyOptionValues(): tmpValues =";
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
    m_ov.m_a_v = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_v)) {
    m_ov.m_b_v = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_rho_v)) {
    m_ov.m_a_rho_v = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_rho_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_rho_v)) {
    m_ov.m_b_rho_v = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_rho_v]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_y)) {
    m_ov.m_a_y = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_y]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_y)) {
    m_ov.m_b_y = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_y]).as<double>();
  }

  return;
}

void
ExperimentModelOptions::print(std::ostream& os) const
{
  queso_deprecated();

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

std::ostream& operator<<(std::ostream& os, const ExperimentModelOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
