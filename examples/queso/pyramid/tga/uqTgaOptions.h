/* uq/examples/queso/pyramid/uqTgaOptions.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_OPTIONS_H__
#define __UQ_TGA_OPTIONS_H__

#include <uqEnvironment.h>

// _ODV = option default value
#define UQ_TGA_OPTION_RUN_TESTS_ODV                    0
#define UQ_TGA_OPTION_CREATE_REFERENCE_ODV             0
#define UQ_TGA_OPTION_REF_W1_ODV                       1.0
#define UQ_TGA_OPTION_REF_W2_ODV                       0.
#define UQ_TGA_OPTION_REF_W3_ODV                       0.
#define UQ_TGA_OPTION_REF_A1_ODV                       2.6000e+11
#define UQ_TGA_OPTION_REF_E1_ODV                       2.0000e+05
#define UQ_TGA_OPTION_REF_A2_ODV                       2.5900e+11
#define UQ_TGA_OPTION_REF_E2_ODV                       2.0010e+05
#define UQ_TGA_OPTION_REF_A3_ODV                       2.6000e+11
#define UQ_TGA_OPTION_REF_E3_ODV                       2.0000e+05
#define UQ_TGA_OPTION_REF_A4_ODV                       2.6000e+11
#define UQ_TGA_OPTION_REF_E4_ODV                       2.0000e+05
#define UQ_TGA_OPTION_REF_TEMPERATURE_PROFILE_ID_ODV   0
#define UQ_TGA_OPTION_REF_MAX_TIME_ODV                 0.
#define UQ_TGA_OPTION_REF_MAX_TIME_STEP_ODV            1.
#define UQ_TGA_OPTION_REF_TREAT_DATA_AS_CONTINUOUS_ODV 0
#define UQ_TGA_OPTION_REF_NUM_DISCRETE_SAMPLES_ODV     12
#define UQ_TGA_OPTION_W_MAX_TIME_STEP_ODV              .1
#define UQ_TGA_OPTION_LAMBDA_MAX_TIME_STEP_ODV         .1
#define UQ_TGA_OPTION_INTEGRALS_NUM_INTERVALS_ODV      1000

class uqTgaOptionsClass
{
public:
  uqTgaOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
 ~uqTgaOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  bool         m_runTests;
  bool         m_createReference;
  double       m_refW1;
  double       m_refW2;
  double       m_refW3;
  double       m_refA1;
  double       m_refE1;
  double       m_refA2;
  double       m_refE2;
  double       m_refA3;
  double       m_refE3;
  double       m_refA4;
  double       m_refE4;
  unsigned int m_refTemperatureProfileId;
  double       m_refMaxTime;
  double       m_refMaxTimeStep;
  bool         m_refTreatDataAsContinuous;
  unsigned int m_refNumDiscreteSamples;
  double       m_wMaxTimeStep;
  double       m_lambdaMaxTimeStep;
  unsigned int m_integralsNumIntervals;

private:
  void   defineMyOptions  (po::options_description& optionsDesc);
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  std::string              m_prefix;
  po::options_description* m_optionsDesc;

  std::string              m_option_help;
  std::string              m_option_runTests;
  std::string              m_option_createReference;
  std::string              m_option_refW1;
  std::string              m_option_refW2;
  std::string              m_option_refW3;
  std::string              m_option_refA1;
  std::string              m_option_refE1;
  std::string              m_option_refA2;
  std::string              m_option_refE2;
  std::string              m_option_refA3;
  std::string              m_option_refE3;
  std::string              m_option_refA4;
  std::string              m_option_refE4;
  std::string              m_option_refTemperatureProfileId;
  std::string              m_option_refMaxTime;
  std::string              m_option_refMaxTimeStep;
  std::string              m_option_refTreatDataAsContinuous;
  std::string              m_option_refNumDiscreteSamples;
  std::string              m_option_wMaxTimeStep;
  std::string              m_option_lambdaMaxTimeStep;
  std::string              m_option_integralsNumIntervals;
};

std::ostream& operator<<(std::ostream& os, const uqTgaOptionsClass& obj);

uqTgaOptionsClass::uqTgaOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_env                            (env),
  m_prefix                         ((std::string)(prefix) + "option_"),
  m_optionsDesc                    (new po::options_description("TGA options")),
  m_option_help                    (m_prefix + "help"                    ),
  m_option_runTests                (m_prefix + "runTests"                ),
  m_option_createReference         (m_prefix + "createReference"         ),
  m_option_refW1                   (m_prefix + "refW1"                   ),
  m_option_refW2                   (m_prefix + "refW2"                   ),
  m_option_refW3                   (m_prefix + "refW3"                   ),
  m_option_refA1                   (m_prefix + "refA1"                   ),
  m_option_refE1                   (m_prefix + "refE1"                   ),
  m_option_refA2                   (m_prefix + "refA2"                   ),
  m_option_refE2                   (m_prefix + "refE2"                   ),
  m_option_refA3                   (m_prefix + "refA3"                   ),
  m_option_refE3                   (m_prefix + "refE3"                   ),
  m_option_refA4                   (m_prefix + "refA4"                   ),
  m_option_refE4                   (m_prefix + "refE4"                   ),
  m_option_refTemperatureProfileId (m_prefix + "refTemperatureProfileId" ),
  m_option_refMaxTime              (m_prefix + "refMaxTime"              ),
  m_option_refMaxTimeStep          (m_prefix + "refMaxTimeStep"          ),
  m_option_refTreatDataAsContinuous(m_prefix + "refTreatDataAsContinuous"),
  m_option_refNumDiscreteSamples   (m_prefix + "refNumDiscreteSamples"   ),
  m_option_wMaxTimeStep            (m_prefix + "wMaxTimeStep"            ),
  m_option_lambdaMaxTimeStep       (m_prefix + "lambdaMaxTimeStep"       ),
  m_option_integralsNumIntervals   (m_prefix + "integralsNumIntervals"   )
{
}

uqTgaOptionsClass::~uqTgaOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
}

void
uqTgaOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqTgaOptionsClass::scanOptionsValues()"
                                   << ": after getting values of options with prefix '" << m_prefix
                                   << "', state of  object is:"
                                   << "\n" << *this
                                   << std::endl;

  return;
};

void
uqTgaOptionsClass::defineMyOptions(po::options_description& optionsDesc)
{
  optionsDesc.add_options()     
    (m_option_help.c_str(),                                                                                                              "produce help message for TGA options"     )
    (m_option_runTests.c_str(),                po::value<bool        >()->default_value(UQ_TGA_OPTION_RUN_TESTS_ODV),                    "run tests"                                )
    (m_option_createReference.c_str(),         po::value<bool        >()->default_value(UQ_TGA_OPTION_CREATE_REFERENCE_ODV),             "create reference"                         )
    (m_option_refW1.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_W1_ODV),                       "reference W1"                             )
    (m_option_refW2.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_W2_ODV),                       "reference W2"                             )
    (m_option_refW3.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_W3_ODV),                       "reference W3"                             )
    (m_option_refA1.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_A1_ODV),                       "reference A1"                             )
    (m_option_refE1.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_E1_ODV),                       "reference E1"                             )
    (m_option_refA2.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_A2_ODV),                       "reference A2"                             )
    (m_option_refE2.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_E2_ODV),                       "reference E2"                             )
    (m_option_refA3.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_A3_ODV),                       "reference A3"                             )
    (m_option_refE3.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_E3_ODV),                       "reference E3"                             )
    (m_option_refA4.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_A4_ODV),                       "reference A4"                             )
    (m_option_refE4.c_str(),                   po::value<double      >()->default_value(UQ_TGA_OPTION_REF_E4_ODV),                       "reference E4"                             )
    (m_option_refTemperatureProfileId.c_str(), po::value<unsigned int>()->default_value(UQ_TGA_OPTION_REF_TEMPERATURE_PROFILE_ID_ODV),   "refTemperatureProfileId"                  )
    (m_option_refMaxTime.c_str(),              po::value<double      >()->default_value(UQ_TGA_OPTION_REF_MAX_TIME_ODV),                 "refMaxTime"                               )
    (m_option_refMaxTimeStep.c_str(),          po::value<double      >()->default_value(UQ_TGA_OPTION_REF_MAX_TIME_STEP_ODV),            "refMaxTimeStep"                           )
    (m_option_refTreatDataAsContinuous.c_str(),po::value<bool        >()->default_value(UQ_TGA_OPTION_REF_TREAT_DATA_AS_CONTINUOUS_ODV), "refTreatDataAsContinuous"                 )
    (m_option_refNumDiscreteSamples.c_str(),   po::value<unsigned int>()->default_value(UQ_TGA_OPTION_REF_NUM_DISCRETE_SAMPLES_ODV),     "refNumDiscreteSamples"                    )
    (m_option_wMaxTimeStep.c_str(),            po::value<double      >()->default_value(UQ_TGA_OPTION_W_MAX_TIME_STEP_ODV),              "wMaxTimeStep"                             )
    (m_option_lambdaMaxTimeStep.c_str(),       po::value<double      >()->default_value(UQ_TGA_OPTION_LAMBDA_MAX_TIME_STEP_ODV),         "lambdaMaxTimeStep"                        )
    (m_option_integralsNumIntervals.c_str(),   po::value<unsigned int>()->default_value(UQ_TGA_OPTION_INTEGRALS_NUM_INTERVALS_ODV),      "integralsNumIntervals"                    )
  ;

  return;
}

void
uqTgaOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_runTests.c_str())) {
    m_runTests = m_env.allOptionsMap()[m_option_runTests.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_createReference.c_str())) {
    m_createReference = m_env.allOptionsMap()[m_option_createReference.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_refW1.c_str())) {
    m_refW1 = m_env.allOptionsMap()[m_option_refW1.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refW2.c_str())) {
    m_refW2 = m_env.allOptionsMap()[m_option_refW2.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refW3.c_str())) {
    m_refW3 = m_env.allOptionsMap()[m_option_refW3.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA1.c_str())) {
    m_refA1 = m_env.allOptionsMap()[m_option_refA1.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE1.c_str())) {
    m_refE1 = m_env.allOptionsMap()[m_option_refE1.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA2.c_str())) {
    m_refA2 = m_env.allOptionsMap()[m_option_refA2.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE2.c_str())) {
    m_refE2 = m_env.allOptionsMap()[m_option_refE2.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA3.c_str())) {
    m_refA3 = m_env.allOptionsMap()[m_option_refA3.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE3.c_str())) {
    m_refE3 = m_env.allOptionsMap()[m_option_refE3.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA4.c_str())) {
    m_refA4 = m_env.allOptionsMap()[m_option_refA4.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE4.c_str())) {
    m_refE4 = m_env.allOptionsMap()[m_option_refE4.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refTemperatureProfileId.c_str())) {
    m_refTemperatureProfileId = m_env.allOptionsMap()[m_option_refTemperatureProfileId.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_refMaxTime.c_str())) {
    m_refMaxTime = m_env.allOptionsMap()[m_option_refMaxTime.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refMaxTimeStep.c_str())) {
    m_refMaxTimeStep = m_env.allOptionsMap()[m_option_refMaxTimeStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refTreatDataAsContinuous.c_str())) {
    m_refTreatDataAsContinuous = m_env.allOptionsMap()[m_option_refTreatDataAsContinuous.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_refNumDiscreteSamples.c_str())) {
    m_refNumDiscreteSamples = m_env.allOptionsMap()[m_option_refNumDiscreteSamples.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_wMaxTimeStep.c_str())) {
    m_wMaxTimeStep = m_env.allOptionsMap()[m_option_wMaxTimeStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_lambdaMaxTimeStep.c_str())) {
    m_lambdaMaxTimeStep = m_env.allOptionsMap()[m_option_lambdaMaxTimeStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_integralsNumIntervals.c_str())) {
    m_integralsNumIntervals = m_env.allOptionsMap()[m_option_integralsNumIntervals.c_str()].as<unsigned int>();
  }

  return;
}

void
uqTgaOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_runTests                 << " = " << m_runTests
     << "\n" << m_option_createReference          << " = " << m_createReference
     << "\n" << m_option_refW1                    << " = " << m_refW1
     << "\n" << m_option_refW2                    << " = " << m_refW2
     << "\n" << m_option_refW3                    << " = " << m_refW3
     << "\n" << m_option_refA1                    << " = " << m_refA1
     << "\n" << m_option_refE1                    << " = " << m_refE1
     << "\n" << m_option_refA2                    << " = " << m_refA2
     << "\n" << m_option_refE2                    << " = " << m_refE2
     << "\n" << m_option_refA3                    << " = " << m_refA3
     << "\n" << m_option_refE3                    << " = " << m_refE3
     << "\n" << m_option_refA4                    << " = " << m_refA4
     << "\n" << m_option_refE4                    << " = " << m_refE4
     << "\n" << m_option_refTemperatureProfileId  << " = " << m_refTemperatureProfileId
     << "\n" << m_option_refMaxTime               << " = " << m_refMaxTime
     << "\n" << m_option_refMaxTimeStep           << " = " << m_refMaxTimeStep
     << "\n" << m_option_refTreatDataAsContinuous << " = " << m_refTreatDataAsContinuous
     << "\n" << m_option_refNumDiscreteSamples    << " = " << m_refNumDiscreteSamples
     << "\n" << m_option_wMaxTimeStep             << " = " << m_wMaxTimeStep
     << "\n" << m_option_lambdaMaxTimeStep        << " = " << m_lambdaMaxTimeStep
     << "\n" << m_option_integralsNumIntervals    << " = " << m_integralsNumIntervals
     << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqTgaOptionsClass& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_TGA_OPTIONS_H__
