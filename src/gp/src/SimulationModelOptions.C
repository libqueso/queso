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

#include <queso/SimulationModelOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

SmOptionsValues::SmOptionsValues()
  :
    m_prefix("sm_"),
    m_dataOutputFileName(UQ_SIMULATION_MODEL_DATA_OUTPUT_FILE_NAME_ODV),
    m_dataOutputAllowAll(UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOW_ALL_ODV),
    //m_dataOutputAllowedSet(),
    m_p_eta(UQ_SIMULATION_MODEL_P_ETA_ODV),
    m_zeroRelativeSingularValue(UQ_SIMULATION_MODEL_ZERO_RELATIVE_SINGULAR_VALUE_ODV),
    m_cdfThresholdForPEta(UQ_SIMULATION_MODEL_CDF_THRESHOLD_FOR_P_ETA_ODV),
    m_a_w(UQ_SIMULATION_MODEL_A_W_ODV),
    m_b_w(UQ_SIMULATION_MODEL_B_W_ODV),
    m_a_rho_w(UQ_SIMULATION_MODEL_A_RHO_W_ODV),
    m_b_rho_w(UQ_SIMULATION_MODEL_B_RHO_W_ODV),
    m_a_eta(UQ_SIMULATION_MODEL_A_ETA_ODV),
    m_b_eta(UQ_SIMULATION_MODEL_B_ETA_ODV),
    m_a_s(UQ_SIMULATION_MODEL_A_S_ODV),
    m_b_s(UQ_SIMULATION_MODEL_B_S_ODV),
    m_parser(NULL),
    m_option_help(m_prefix + "help"),
    m_option_dataOutputFileName(m_prefix + "dataOutputFileName"),
    m_option_dataOutputAllowAll(m_prefix + "dataOutputAllowAll"),
    m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
    m_option_p_eta(m_prefix + "p_eta"),
    m_option_zeroRelativeSingularValue(m_prefix + "zeroRelativeSingularValue"),
    m_option_cdfThresholdForPEta(m_prefix + "cdfThresholdForPEta"),
    m_option_a_w(m_prefix + "a_w"),
    m_option_b_w(m_prefix + "b_w"),
    m_option_a_rho_w(m_prefix + "a_rho_w"),
    m_option_b_rho_w(m_prefix + "b_rho_w"),
    m_option_a_eta(m_prefix + "a_eta"),
    m_option_b_eta(m_prefix + "b_eta"),
    m_option_a_s(m_prefix + "a_s"),
    m_option_b_s(m_prefix + "b_s")
{
}

SmOptionsValues::SmOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
    m_prefix((std::string)(prefix) + "sm_"),
    m_dataOutputFileName(UQ_SIMULATION_MODEL_DATA_OUTPUT_FILE_NAME_ODV),
    m_dataOutputAllowAll(UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOW_ALL_ODV),
    //m_dataOutputAllowedSet(),
    m_p_eta(UQ_SIMULATION_MODEL_P_ETA_ODV),
    m_zeroRelativeSingularValue(UQ_SIMULATION_MODEL_ZERO_RELATIVE_SINGULAR_VALUE_ODV),
    m_cdfThresholdForPEta(UQ_SIMULATION_MODEL_CDF_THRESHOLD_FOR_P_ETA_ODV),
    m_a_w(UQ_SIMULATION_MODEL_A_W_ODV),
    m_b_w(UQ_SIMULATION_MODEL_B_W_ODV),
    m_a_rho_w(UQ_SIMULATION_MODEL_A_RHO_W_ODV),
    m_b_rho_w(UQ_SIMULATION_MODEL_B_RHO_W_ODV),
    m_a_eta(UQ_SIMULATION_MODEL_A_ETA_ODV),
    m_b_eta(UQ_SIMULATION_MODEL_B_ETA_ODV),
    m_a_s(UQ_SIMULATION_MODEL_A_S_ODV),
    m_b_s(UQ_SIMULATION_MODEL_B_S_ODV),
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
    m_option_help(m_prefix + "help"),
    m_option_dataOutputFileName(m_prefix + "dataOutputFileName"),
    m_option_dataOutputAllowAll(m_prefix + "dataOutputAllowAll"),
    m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
    m_option_p_eta(m_prefix + "p_eta"),
    m_option_zeroRelativeSingularValue(m_prefix + "zeroRelativeSingularValue"),
    m_option_cdfThresholdForPEta(m_prefix + "cdfThresholdForPEta"),
    m_option_a_w(m_prefix + "a_w"),
    m_option_b_w(m_prefix + "b_w"),
    m_option_a_rho_w(m_prefix + "a_rho_w"),
    m_option_b_rho_w(m_prefix + "b_rho_w"),
    m_option_a_eta(m_prefix + "a_eta"),
    m_option_b_eta(m_prefix + "b_eta"),
    m_option_a_s(m_prefix + "a_s"),
    m_option_b_s(m_prefix + "b_s")
{
  m_parser->registerOption(m_option_help,                                                                                                                      "produce help message for simulation model options");
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,        UQ_SIMULATION_MODEL_DATA_OUTPUT_FILE_NAME_ODV       , "name of data output file"                         );
  m_parser->registerOption<bool        >(m_option_dataOutputAllowAll,        UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOW_ALL_ODV       , "allow all or not"                                 );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,      UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOWED_SET_ODV     , "subEnvs that will write to data output file"      );
  m_parser->registerOption<unsigned int>(m_option_p_eta,                     UQ_SIMULATION_MODEL_P_ETA_ODV                       , "p_eta"                                            );
  m_parser->registerOption<double      >(m_option_zeroRelativeSingularValue, UQ_SIMULATION_MODEL_ZERO_RELATIVE_SINGULAR_VALUE_ODV, "zeroRelativeSingularValue"                        );
  m_parser->registerOption<double      >(m_option_cdfThresholdForPEta,       UQ_SIMULATION_MODEL_CDF_THRESHOLD_FOR_P_ETA_ODV     , "cdfThresholdForPEta"                              );
  m_parser->registerOption<double      >(m_option_a_w,                       UQ_SIMULATION_MODEL_A_W_ODV                         , "a_w"                                              );
  m_parser->registerOption<double      >(m_option_b_w,                       UQ_SIMULATION_MODEL_B_W_ODV                         , "b_w"                                              );
  m_parser->registerOption<double      >(m_option_a_rho_w,                   UQ_SIMULATION_MODEL_A_RHO_W_ODV                     , "a_rho_w"                                          );
  m_parser->registerOption<double      >(m_option_b_rho_w,                   UQ_SIMULATION_MODEL_B_RHO_W_ODV                     , "b_rho_w"                                          );
  m_parser->registerOption<double      >(m_option_a_eta,                     UQ_SIMULATION_MODEL_A_ETA_ODV                       , "a_eta"                                            );
  m_parser->registerOption<double      >(m_option_b_eta,                     UQ_SIMULATION_MODEL_B_ETA_ODV                       , "b_eta"                                            );
  m_parser->registerOption<double      >(m_option_a_s,                       UQ_SIMULATION_MODEL_A_S_ODV                         , "a_s"                                              );
  m_parser->registerOption<double      >(m_option_b_s,                       UQ_SIMULATION_MODEL_B_S_ODV                         , "b_s"                                              );

  m_parser->scanInputFile();

  m_parser->getOption<std::string >(m_option_dataOutputFileName,        m_dataOutputFileName);
  m_parser->getOption<bool        >(m_option_dataOutputAllowAll,        m_dataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,      m_dataOutputAllowedSet);
  m_parser->getOption<unsigned int>(m_option_p_eta,                     m_p_eta);
  m_parser->getOption<double      >(m_option_zeroRelativeSingularValue, m_zeroRelativeSingularValue);
  m_parser->getOption<double      >(m_option_cdfThresholdForPEta,       m_cdfThresholdForPEta);
  m_parser->getOption<double      >(m_option_a_w,                       m_a_w);
  m_parser->getOption<double      >(m_option_b_w,                       m_b_w);
  m_parser->getOption<double      >(m_option_a_rho_w,                   m_a_rho_w);
  m_parser->getOption<double      >(m_option_b_rho_w,                   m_b_rho_w);
  m_parser->getOption<double      >(m_option_a_eta,                     m_a_eta);
  m_parser->getOption<double      >(m_option_b_eta,                     m_b_eta);
  m_parser->getOption<double      >(m_option_a_s,                       m_a_s);
  m_parser->getOption<double      >(m_option_b_s,                       m_b_s);
}

SmOptionsValues::~SmOptionsValues()
{
}

SmOptionsValues::SmOptionsValues(const SmOptionsValues& src)
{
  this->copy(src);
}

SmOptionsValues&
SmOptionsValues::operator=(const SmOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}

// void
// SmOptionsValues::defineOptions()
// {
//   (*m_optionsDescription).add_options()
//     (m_option_help.c_str(),                                                                                                                      "produce help message for simulation model options")
//     (m_option_dataOutputFileName.c_str(),        boost::program_options::value<std::string >()->default_value(UQ_SIMULATION_MODEL_DATA_OUTPUT_FILE_NAME_ODV       ), "name of data output file"                         )
//     (m_option_dataOutputAllowAll.c_str(),        boost::program_options::value<bool        >()->default_value(UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOW_ALL_ODV       ), "allow all or not"                                 )
//     (m_option_dataOutputAllowedSet.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOWED_SET_ODV     ), "subEnvs that will write to data output file"      )
//     (m_option_p_eta.c_str(),                     boost::program_options::value<unsigned int>()->default_value(UQ_SIMULATION_MODEL_P_ETA_ODV                       ), "p_eta"                                            )
//     (m_option_zeroRelativeSingularValue.c_str(), boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_ZERO_RELATIVE_SINGULAR_VALUE_ODV), "zeroRelativeSingularValue"                        )
//     (m_option_cdfThresholdForPEta.c_str(),       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_CDF_THRESHOLD_FOR_P_ETA_ODV     ), "cdfThresholdForPEta"                              )
//     (m_option_a_w.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_W_ODV                         ), "a_w"                                              )
//     (m_option_b_w.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_W_ODV                         ), "b_w"                                              )
//     (m_option_a_rho_w.c_str(),                   boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_RHO_W_ODV                     ), "a_rho_w"                                          )
//     (m_option_b_rho_w.c_str(),                   boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_RHO_W_ODV                     ), "b_rho_w"                                          )
//     (m_option_a_eta.c_str(),                     boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_ETA_ODV                       ), "a_eta"                                            )
//     (m_option_b_eta.c_str(),                     boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_ETA_ODV                       ), "b_eta"                                            )
//     (m_option_a_s.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_S_ODV                         ), "a_s"                                              )
//     (m_option_b_s.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_S_ODV                         ), "b_s"                                              )
//   ;
// }
//
// void
// SmOptionsValues::getOptionValues()
// {
//   if ((*m_optionsMap).count(m_option_help)) {
//     if (m_env->subDisplayFile()) {
//       *m_env->subDisplayFile() << (*m_optionsDescription)
//                               << std::endl;
//     }
//   }
//
//   if ((*m_optionsMap).count(m_option_dataOutputFileName)) {
//     m_dataOutputFileName = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_dataOutputFileName]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_dataOutputAllowAll)) {
//     m_dataOutputAllowAll = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_dataOutputAllowAll]).as<bool>();
//   }
//
//   if (m_dataOutputAllowAll) {
//     m_dataOutputAllowedSet.insert(m_env->subId());
//   }
//   else if ((*m_optionsMap).count(m_option_dataOutputAllowedSet)) {
//     m_dataOutputAllowedSet.clear();
//     std::vector<double> tmpAllow(0,0.);
//     std::string inputString = (*m_optionsMap)[m_option_dataOutputAllowedSet].as<std::string>();
//     MiscReadDoublesFromString(inputString,tmpAllow);
//
//     if (tmpAllow.size() > 0) {
//       for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
//         m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
//       }
//     }
//   }
//
//   if ((*m_optionsMap).count(m_option_p_eta)) {
//     m_p_eta = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_p_eta]).as<unsigned int>();
//   }
//
//   if ((*m_optionsMap).count(m_option_zeroRelativeSingularValue)) {
//     m_zeroRelativeSingularValue = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_zeroRelativeSingularValue]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_cdfThresholdForPEta)) {
//     m_cdfThresholdForPEta = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_cdfThresholdForPEta]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_a_w)) {
//     m_a_w = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_w]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_w)) {
//     m_b_w = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_w]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_a_rho_w)) {
//     m_a_rho_w = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_rho_w]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_rho_w)) {
//     m_b_rho_w = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_rho_w]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_a_eta)) {
//     m_a_eta = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_eta]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_eta)) {
//     m_b_eta = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_eta]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_a_s)) {
//     m_a_s = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_a_s]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_b_s)) {
//     m_b_s = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_b_s]).as<double>();
//   }
// }

void
SmOptionsValues::copy(const SmOptionsValues& src)
{
  m_dataOutputFileName        = src.m_dataOutputFileName;
  m_dataOutputAllowAll        = src.m_dataOutputAllowAll;
  m_dataOutputAllowedSet      = src.m_dataOutputAllowedSet;
  m_p_eta                     = src.m_p_eta;
  m_zeroRelativeSingularValue = src.m_zeroRelativeSingularValue;
  m_cdfThresholdForPEta       = src.m_cdfThresholdForPEta;
  m_a_w                       = src.m_a_w;
  m_b_w                       = src.m_b_w;
  m_a_rho_w                   = src.m_a_rho_w;
  m_b_rho_w                   = src.m_b_rho_w;
  m_a_eta                     = src.m_a_eta;
  m_b_eta                     = src.m_b_eta;
  m_a_s                       = src.m_a_s;
  m_b_s                       = src.m_b_s;

  return;
}

SimulationModelOptions::SimulationModelOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_ov                              (),
  m_prefix                          ((std::string)(prefix) + "sm_"),
  m_env                             (env),
  m_optionsDesc                     (new boost::program_options::options_description("Simulation model options")),
  m_option_help                     (m_prefix + "help"                     ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"       ),
  m_option_dataOutputAllowAll       (m_prefix + "dataOutputAllowAll"       ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"     ),
  m_option_p_eta                    (m_prefix + "p_eta"                    ),
  m_option_zeroRelativeSingularValue(m_prefix + "zeroRelativeSingularValue"),
  m_option_cdfThresholdForPEta      (m_prefix + "cdfThresholdForPEta"      ),
  m_option_a_w                      (m_prefix + "a_w"                      ),
  m_option_b_w                      (m_prefix + "b_w"                      ),
  m_option_a_rho_w                  (m_prefix + "a_rho_w"                  ),
  m_option_b_rho_w                  (m_prefix + "b_rho_w"                  ),
  m_option_a_eta                    (m_prefix + "a_eta"                    ),
  m_option_b_eta                    (m_prefix + "b_eta"                    ),
  m_option_a_s                      (m_prefix + "a_s"                      ),
  m_option_b_s                      (m_prefix + "b_s"                      )
{
  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the abscense of an options input file");
}

SimulationModelOptions::SimulationModelOptions(
  const BaseEnvironment&  env,
  const char*                    prefix,
  const SmOptionsValues& alternativeOptionsValues)
  :
  m_ov                              (alternativeOptionsValues),
  m_prefix                          ((std::string)(prefix) + "sm_"),
  m_env                             (env),
  m_optionsDesc                     (NULL),
  m_option_help                     (m_prefix + "help"                     ),
  m_option_dataOutputFileName       (m_prefix + "dataOutputFileName"       ),
  m_option_dataOutputAllowAll       (m_prefix + "dataOutputAllowAll"       ),
  m_option_dataOutputAllowedSet     (m_prefix + "dataOutputAllowedSet"     ),
  m_option_p_eta                    (m_prefix + "p_eta"                    ),
  m_option_zeroRelativeSingularValue(m_prefix + "zeroRelativeSingularValue"),
  m_option_cdfThresholdForPEta      (m_prefix + "cdfThresholdForPEta"      ),
  m_option_a_w                      (m_prefix + "a_w"                      ),
  m_option_b_w                      (m_prefix + "b_w"                      ),
  m_option_a_rho_w                  (m_prefix + "a_rho_w"                  ),
  m_option_b_rho_w                  (m_prefix + "b_rho_w"                  ),
  m_option_a_eta                    (m_prefix + "a_eta"                    ),
  m_option_b_eta                    (m_prefix + "b_eta"                    ),
  m_option_a_s                      (m_prefix + "a_s"                      ),
  m_option_b_s                      (m_prefix + "b_s"                      )
{
  queso_deprecated();

  queso_require_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In SimulationModelOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

SimulationModelOptions::~SimulationModelOptions()
{
  queso_deprecated();

  if (m_optionsDesc) delete m_optionsDesc;
}

void
SimulationModelOptions::scanOptionsValues()
{
  queso_deprecated();

  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  //std::cout << "scan 000\n"
  //          << std::endl;
  getMyOptionValues              (*m_optionsDesc);
  //std::cout << "scan 001\n"
  //          << std::endl;

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In SimulationModelOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
SimulationModelOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                      "produce help message for simulation model options")
    (m_option_dataOutputFileName.c_str(),        boost::program_options::value<std::string >()->default_value(UQ_SIMULATION_MODEL_DATA_OUTPUT_FILE_NAME_ODV       ), "name of data output file"                         )
    (m_option_dataOutputAllowAll.c_str(),        boost::program_options::value<bool        >()->default_value(UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOW_ALL_ODV       ), "allow all or not"                                 )
    (m_option_dataOutputAllowedSet.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_SIMULATION_MODEL_DATA_OUTPUT_ALLOWED_SET_ODV     ), "subEnvs that will write to data output file"      )
    (m_option_p_eta.c_str(),                     boost::program_options::value<unsigned int>()->default_value(UQ_SIMULATION_MODEL_P_ETA_ODV                       ), "p_eta"                                            )
    (m_option_zeroRelativeSingularValue.c_str(), boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_ZERO_RELATIVE_SINGULAR_VALUE_ODV), "zeroRelativeSingularValue"                        )
    (m_option_cdfThresholdForPEta.c_str(),       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_CDF_THRESHOLD_FOR_P_ETA_ODV     ), "cdfThresholdForPEta"                              )
    (m_option_a_w.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_W_ODV                         ), "a_w"                                              )
    (m_option_b_w.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_W_ODV                         ), "b_w"                                              )
    (m_option_a_rho_w.c_str(),                   boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_RHO_W_ODV                     ), "a_rho_w"                                          )
    (m_option_b_rho_w.c_str(),                   boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_RHO_W_ODV                     ), "b_rho_w"                                          )
    (m_option_a_eta.c_str(),                     boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_ETA_ODV                       ), "a_eta"                                            )
    (m_option_b_eta.c_str(),                     boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_ETA_ODV                       ), "b_eta"                                            )
    (m_option_a_s.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_A_S_ODV                         ), "a_s"                                              )
    (m_option_b_s.c_str(),                       boost::program_options::value<double      >()->default_value(UQ_SIMULATION_MODEL_B_S_ODV                         ), "b_s"                                              )
  ;

  return;
}

void
SimulationModelOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowAll)) {
    m_ov.m_dataOutputAllowAll = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputAllowAll]).as<bool>();
  }

  if (m_ov.m_dataOutputAllowAll) {
    m_ov.m_dataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_ov.m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_p_eta)) {
    m_ov.m_p_eta = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_p_eta]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_zeroRelativeSingularValue)) {
    m_ov.m_zeroRelativeSingularValue = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_zeroRelativeSingularValue]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_cdfThresholdForPEta)) {
    m_ov.m_cdfThresholdForPEta = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_cdfThresholdForPEta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_w)) {
    m_ov.m_a_w = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_w]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_w)) {
    m_ov.m_b_w = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_w]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_rho_w)) {
    m_ov.m_a_rho_w = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_rho_w]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_rho_w)) {
    m_ov.m_b_rho_w = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_rho_w]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_eta)) {
    m_ov.m_a_eta = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_eta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_eta)) {
    m_ov.m_b_eta = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_eta]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_a_s)) {
    m_ov.m_a_s = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_a_s]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_b_s)) {
    m_ov.m_b_s = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_b_s]).as<double>();
  }

  return;
}

void
SimulationModelOptions::print(std::ostream& os) const
{
  queso_deprecated();

  os << "\n" << m_option_dataOutputFileName        << " = " << m_ov.m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowAll        << " = " << m_ov.m_dataOutputAllowAll
     << "\n" << m_option_dataOutputAllowedSet      << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_p_eta                     << " = " << m_ov.m_p_eta
     << "\n" << m_option_zeroRelativeSingularValue << " = " << m_ov.m_zeroRelativeSingularValue
     << "\n" << m_option_cdfThresholdForPEta       << " = " << m_ov.m_cdfThresholdForPEta
     << "\n" << m_option_a_w                       << " = " << m_ov.m_a_w
     << "\n" << m_option_b_w                       << " = " << m_ov.m_b_w
     << "\n" << m_option_a_rho_w                   << " = " << m_ov.m_a_rho_w
     << "\n" << m_option_b_rho_w                   << " = " << m_ov.m_b_rho_w
     << "\n" << m_option_a_eta                     << " = " << m_ov.m_a_eta
     << "\n" << m_option_b_eta                     << " = " << m_ov.m_b_eta
     << "\n" << m_option_a_s                       << " = " << m_ov.m_a_s
     << "\n" << m_option_b_s                       << " = " << m_ov.m_b_s
     << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const SimulationModelOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
