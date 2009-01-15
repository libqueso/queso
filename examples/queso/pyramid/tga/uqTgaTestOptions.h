/* uq/examples/queso/pyramid/uqTgaTestOptions.h
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

#ifndef __UQ_TGA_TEST_OPTIONS_H__
#define __UQ_TGA_TEST_OPTIONS_H__

#include <uqEnvironment.h>

// _ODV = option default value
#define UQ_TGA_TEST_OUTER_PREFIX_NAME_ODV            "case0_"
#define UQ_TGA_TEST_RUN_TEMP_TIME_TEST_ODV           0
#define UQ_TGA_TEST_CREATE_REFERENCE_ODV             0
#define UQ_TGA_TEST_RUN_TIMING_TEST_ODV              0
#define UQ_TGA_TEST_RUN_GRAD_TEST_ODV                0
#define UQ_TGA_TEST_RUN_OPTIMIZATION_TEST_ODV        0
#define UQ_TGA_TEST_REF_W1_ODV                       1.0
#define UQ_TGA_TEST_REF_W2_ODV                       0.
#define UQ_TGA_TEST_REF_W3_ODV                       0.
#define UQ_TGA_TEST_REF_A1_ODV                       2.6000e+11
#define UQ_TGA_TEST_REF_E1_ODV                       2.0000e+05
#define UQ_TGA_TEST_REF_A2_ODV                       2.5900e+11
#define UQ_TGA_TEST_REF_E2_ODV                       2.0010e+05
#define UQ_TGA_TEST_REF_A3_ODV                       2.6000e+11
#define UQ_TGA_TEST_REF_E3_ODV                       2.0000e+05
#define UQ_TGA_TEST_REF_A4_ODV                       2.6000e+11
#define UQ_TGA_TEST_REF_E4_ODV                       2.0000e+05
#define UQ_TGA_TEST_GUESS_A_ODV                      2.5910e+11
#define UQ_TGA_TEST_GUESS_E_ODV                      2.0090e+05
#define UQ_TGA_TEST_REF_MAX_TIME_STEP_ODV            1.
#define UQ_TGA_TEST_REF_TREAT_DATA_AS_CONTINUOUS_ODV 0
#define UQ_TGA_TEST_REF_NUM_DISCRETE_SAMPLES_ODV     12
#define UQ_TGA_TEST_COMPUTE_HESSIAN_ODV              1
#define UQ_TGA_TEST_W_MAX_TIME_STEP_ODV              .1
#define UQ_TGA_TEST_LAMBDA_MAX_TIME_STEP_ODV         .1
#define UQ_TGA_TEST_INTEGRALS_NUM_INTERVALS_ODV      1000
#define UQ_TGA_TEST_RELATIVE_FD_STEP_ODV             1.e-8
#define UQ_TGA_TEST_WRITE_OUTPUT_ODV                 1

class uqTgaTestOptionsClass
{
public:
  uqTgaTestOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
 ~uqTgaTestOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  std::string  outerPrefixName;
  bool         runTempTimeTest;
  bool         createReference;
  bool         runTimingTest;
  bool         runGradTest;
  bool         runOptimizationTest;
  double       refW1;
  double       refW2;
  double       refW3;
  double       refA1;
  double       refE1;
  double       refA2;
  double       refE2;
  double       refA3;
  double       refE3;
  double       refA4;
  double       refE4;
  double       guessA;
  double       guessE;
  double       refMaxTimeStep;
  bool         refTreatDataAsContinuous;
  unsigned int refNumDiscreteSamples;
  bool         computeHessian;
  double       wMaxTimeStep;
  double       lambdaMaxTimeStep;
  unsigned int integralsNumIntervals;
  double       relativeFDStep;
  bool         writeOutput;

private:
  void   defineMyOptions  (po::options_description& optionsDesc);
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  std::string              m_prefix;
  po::options_description* m_optionsDesc;

  std::string              m_option_help;
  std::string              m_option_outerPrefixName;
  std::string              m_option_runTempTimeTest;
  std::string              m_option_createReference;
  std::string              m_option_runTimingTest;
  std::string              m_option_runGradTest;
  std::string              m_option_runOptimizationTest;
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
  std::string              m_option_guessA;
  std::string              m_option_guessE;
  std::string              m_option_refMaxTimeStep;
  std::string              m_option_refTreatDataAsContinuous;
  std::string              m_option_refNumDiscreteSamples;
  std::string              m_option_computeHessian;
  std::string              m_option_wMaxTimeStep;
  std::string              m_option_lambdaMaxTimeStep;
  std::string              m_option_integralsNumIntervals;
  std::string              m_option_relativeFDStep;
  std::string              m_option_writeOutput;
};

std::ostream& operator<<(std::ostream& os, const uqTgaTestOptionsClass& obj);

uqTgaTestOptionsClass::uqTgaTestOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_env                            (env),
  m_prefix                         ((std::string)(prefix) + "tests_"),
  m_optionsDesc                    (new po::options_description("TGA test options")),
  m_option_help                    (m_prefix + "help"                    ),
  m_option_outerPrefixName         (m_prefix + "outerPrefixName"         ),
  m_option_runTempTimeTest         (m_prefix + "runTempTimeTest"         ),
  m_option_createReference         (m_prefix + "createReference"         ),
  m_option_runTimingTest           (m_prefix + "runTimingTest"           ),
  m_option_runGradTest             (m_prefix + "runGradTest"             ),
  m_option_runOptimizationTest     (m_prefix + "runOptimizationTest"     ),
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
  m_option_guessA                  (m_prefix + "guessA"                  ),
  m_option_guessE                  (m_prefix + "guessE"                  ),
  m_option_refMaxTimeStep          (m_prefix + "refMaxTimeStep"          ),
  m_option_refTreatDataAsContinuous(m_prefix + "refTreatDataAsContinuous"),
  m_option_refNumDiscreteSamples   (m_prefix + "refNumDiscreteSamples"   ),
  m_option_computeHessian          (m_prefix + "computeHessian"          ),
  m_option_wMaxTimeStep            (m_prefix + "wMaxTimeStep"            ),
  m_option_lambdaMaxTimeStep       (m_prefix + "lambdaMaxTimeStep"       ),
  m_option_integralsNumIntervals   (m_prefix + "integralsNumIntervals"   ),
  m_option_relativeFDStep          (m_prefix + "relativeFDStep"          ),
  m_option_writeOutput             (m_prefix + "writeOutput"             )
{
}

uqTgaTestOptionsClass::~uqTgaTestOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
}

void
uqTgaTestOptionsClass::scanOptionsValues()
{
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqTgaTestOptionsClass::scanOptionsValues()"
                                   << ": after getting values of options with prefix '" << m_prefix
                                   << "', state of  object is:"
                                   << "\n" << *this
                                   << std::endl;

  return;
};

void
uqTgaTestOptionsClass::defineMyOptions(po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                            "produce help message for TGA test options")
    (m_option_outerPrefixName.c_str(),         po::value<std::string >()->default_value(UQ_TGA_TEST_OUTER_PREFIX_NAME_ODV),            "prefix on output variables"               )
    (m_option_runTempTimeTest.c_str(),         po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_TEMP_TIME_TEST_ODV),           "run tempTime() test"                      )
    (m_option_createReference.c_str(),         po::value<bool        >()->default_value(UQ_TGA_TEST_CREATE_REFERENCE_ODV),             "create reference"                         )
    (m_option_runTimingTest.c_str(),           po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_TIMING_TEST_ODV),              "run timing() test"                        )
    (m_option_runGradTest.c_str(),             po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_GRAD_TEST_ODV),                "run grad() test"                          )
    (m_option_runOptimizationTest.c_str(),     po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_OPTIMIZATION_TEST_ODV),        "run optimization() test"                  )
    (m_option_refW1.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_W1_ODV),                       "reference W1"                             )
    (m_option_refW2.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_W2_ODV),                       "reference W2"                             )
    (m_option_refW3.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_W3_ODV),                       "reference W3"                             )
    (m_option_refA1.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_A1_ODV),                       "reference A1"                             )
    (m_option_refE1.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_E1_ODV),                       "reference E1"                             )
    (m_option_refA2.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_A2_ODV),                       "reference A2"                             )
    (m_option_refE2.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_E2_ODV),                       "reference E2"                             )
    (m_option_refA3.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_A3_ODV),                       "reference A3"                             )
    (m_option_refE3.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_E3_ODV),                       "reference E3"                             )
    (m_option_refA4.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_A4_ODV),                       "reference A4"                             )
    (m_option_refE4.c_str(),                   po::value<double      >()->default_value(UQ_TGA_TEST_REF_E4_ODV),                       "reference E4"                             )
    (m_option_guessA.c_str(),                  po::value<double      >()->default_value(UQ_TGA_TEST_GUESS_A_ODV),                      "guess A"                                  )
    (m_option_guessE.c_str(),                  po::value<double      >()->default_value(UQ_TGA_TEST_GUESS_E_ODV),                      "guess E"                                  )
    (m_option_refMaxTimeStep.c_str(),          po::value<double      >()->default_value(UQ_TGA_TEST_REF_MAX_TIME_STEP_ODV),            "refMaxTimeStep"                           )
    (m_option_refTreatDataAsContinuous.c_str(),po::value<bool        >()->default_value(UQ_TGA_TEST_REF_TREAT_DATA_AS_CONTINUOUS_ODV), "refTreatDataAsContinuous"                 )
    (m_option_refNumDiscreteSamples.c_str(),   po::value<unsigned int>()->default_value(UQ_TGA_TEST_REF_NUM_DISCRETE_SAMPLES_ODV),     "refNumDiscreteSamples"                    )
    (m_option_computeHessian.c_str(),          po::value<bool        >()->default_value(UQ_TGA_TEST_COMPUTE_HESSIAN_ODV),              "computeHessian"                           )
    (m_option_wMaxTimeStep.c_str(),            po::value<double      >()->default_value(UQ_TGA_TEST_W_MAX_TIME_STEP_ODV),              "wMaxTimeStep"                             )
    (m_option_lambdaMaxTimeStep.c_str(),       po::value<double      >()->default_value(UQ_TGA_TEST_LAMBDA_MAX_TIME_STEP_ODV),         "lambdaMaxTimeStep"                        )
    (m_option_integralsNumIntervals.c_str(),   po::value<unsigned int>()->default_value(UQ_TGA_TEST_INTEGRALS_NUM_INTERVALS_ODV),      "integralsNumIntervals"                    )
    (m_option_relativeFDStep.c_str(),          po::value<double      >()->default_value(UQ_TGA_TEST_RELATIVE_FD_STEP_ODV),             "relativeFDStep"                           )
    (m_option_writeOutput.c_str(),             po::value<bool        >()->default_value(UQ_TGA_TEST_WRITE_OUTPUT_ODV),                 "writeOutput"                              )
  ;

  return;
}

void
uqTgaTestOptionsClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_outerPrefixName.c_str())) {
    outerPrefixName = m_env.allOptionsMap()[m_option_outerPrefixName.c_str()].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_runTempTimeTest.c_str())) {
    runTempTimeTest = m_env.allOptionsMap()[m_option_runTempTimeTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_createReference.c_str())) {
    createReference = m_env.allOptionsMap()[m_option_createReference.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_runTimingTest.c_str())) {
    runTimingTest = m_env.allOptionsMap()[m_option_runTimingTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_runGradTest.c_str())) {
    runGradTest = m_env.allOptionsMap()[m_option_runGradTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_runOptimizationTest.c_str())) {
    runOptimizationTest = m_env.allOptionsMap()[m_option_runOptimizationTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_refW1.c_str())) {
    refW1 = m_env.allOptionsMap()[m_option_refW1.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refW2.c_str())) {
    refW2 = m_env.allOptionsMap()[m_option_refW2.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refW3.c_str())) {
    refW3 = m_env.allOptionsMap()[m_option_refW3.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA1.c_str())) {
    refA1 = m_env.allOptionsMap()[m_option_refA1.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE1.c_str())) {
    refE1 = m_env.allOptionsMap()[m_option_refE1.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA2.c_str())) {
    refA2 = m_env.allOptionsMap()[m_option_refA2.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE2.c_str())) {
    refE2 = m_env.allOptionsMap()[m_option_refE2.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA3.c_str())) {
    refA3 = m_env.allOptionsMap()[m_option_refA3.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE3.c_str())) {
    refE3 = m_env.allOptionsMap()[m_option_refE3.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refA4.c_str())) {
    refA4 = m_env.allOptionsMap()[m_option_refA4.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refE4.c_str())) {
    refE4 = m_env.allOptionsMap()[m_option_refE4.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_guessA.c_str())) {
    guessA = m_env.allOptionsMap()[m_option_guessA.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_guessE.c_str())) {
    guessE = m_env.allOptionsMap()[m_option_guessE.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refMaxTimeStep.c_str())) {
    refMaxTimeStep = m_env.allOptionsMap()[m_option_refMaxTimeStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_refTreatDataAsContinuous.c_str())) {
    refTreatDataAsContinuous = m_env.allOptionsMap()[m_option_refTreatDataAsContinuous.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_refNumDiscreteSamples.c_str())) {
    refNumDiscreteSamples = m_env.allOptionsMap()[m_option_refNumDiscreteSamples.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_computeHessian.c_str())) {
    computeHessian = m_env.allOptionsMap()[m_option_computeHessian.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_wMaxTimeStep.c_str())) {
    wMaxTimeStep = m_env.allOptionsMap()[m_option_wMaxTimeStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_lambdaMaxTimeStep.c_str())) {
    lambdaMaxTimeStep = m_env.allOptionsMap()[m_option_lambdaMaxTimeStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_integralsNumIntervals.c_str())) {
    integralsNumIntervals = m_env.allOptionsMap()[m_option_integralsNumIntervals.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_relativeFDStep.c_str())) {
    relativeFDStep = m_env.allOptionsMap()[m_option_relativeFDStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_writeOutput.c_str())) {
    writeOutput = m_env.allOptionsMap()[m_option_writeOutput.c_str()].as<bool>();
  }

  return;
}

void
uqTgaTestOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_outerPrefixName          << " = " << outerPrefixName
     << "\n" << m_option_runTempTimeTest          << " = " << runTempTimeTest
     << "\n" << m_option_createReference          << " = " << createReference
     << "\n" << m_option_runTimingTest            << " = " << runTimingTest
     << "\n" << m_option_runGradTest              << " = " << runGradTest
     << "\n" << m_option_runOptimizationTest      << " = " << runOptimizationTest
     << "\n" << m_option_refW1                    << " = " << refW1
     << "\n" << m_option_refW2                    << " = " << refW2
     << "\n" << m_option_refW3                    << " = " << refW3
     << "\n" << m_option_refA1                    << " = " << refA1
     << "\n" << m_option_refE1                    << " = " << refE1
     << "\n" << m_option_refA2                    << " = " << refA2
     << "\n" << m_option_refE2                    << " = " << refE2
     << "\n" << m_option_refA3                    << " = " << refA3
     << "\n" << m_option_refE3                    << " = " << refE3
     << "\n" << m_option_refA4                    << " = " << refA4
     << "\n" << m_option_refE4                    << " = " << refE4
     << "\n" << m_option_guessA                   << " = " << guessA
     << "\n" << m_option_guessE                   << " = " << guessE
     << "\n" << m_option_refMaxTimeStep           << " = " << refMaxTimeStep
     << "\n" << m_option_refTreatDataAsContinuous << " = " << refTreatDataAsContinuous
     << "\n" << m_option_refNumDiscreteSamples    << " = " << refNumDiscreteSamples
     << "\n" << m_option_computeHessian           << " = " << computeHessian
     << "\n" << m_option_wMaxTimeStep             << " = " << wMaxTimeStep
     << "\n" << m_option_lambdaMaxTimeStep        << " = " << lambdaMaxTimeStep
     << "\n" << m_option_integralsNumIntervals    << " = " << integralsNumIntervals
     << "\n" << m_option_relativeFDStep           << " = " << relativeFDStep
     << "\n" << m_option_writeOutput              << " = " << writeOutput
     << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const uqTgaTestOptionsClass& obj)
{
  obj.print(os);

  return os;
}

////////////////////////////////////////////////////////////////////
// Supporting structure
////////////////////////////////////////////////////////////////////
struct uqTgaTestVarsStruct {
  uqTgaTestVarsStruct();
 ~uqTgaTestVarsStruct();

  uqSampled1D1DFunctionClass*  refW;
  std::ofstream*               ofs;
  uqSampled1D1DFunctionClass*  continuousWeightFunction;
  uqSampled1D1DFunctionClass*  tildeContinuousWeightFunction;
  uqDeltaSeq1D1DFunctionClass* deltaWeightFunction;
  uqDeltaSeq1D1DFunctionClass* tildeDeltaWeightFunction;
};

uqTgaTestVarsStruct::uqTgaTestVarsStruct()
  :
  refW                         (NULL),
  ofs                          (NULL),
  continuousWeightFunction     (NULL),
  tildeContinuousWeightFunction(NULL),
  deltaWeightFunction          (NULL),
  tildeDeltaWeightFunction     (NULL)
{
}

uqTgaTestVarsStruct::~uqTgaTestVarsStruct()
{
  delete refW;
  delete ofs;
  delete continuousWeightFunction;
  delete tildeContinuousWeightFunction;
  delete deltaWeightFunction;
  delete tildeDeltaWeightFunction;
}
#endif // __UQ_TGA_TEST_OPTIONS_H__
