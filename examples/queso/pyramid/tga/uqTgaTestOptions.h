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
#define UQ_TGA_TEST_RUN_TIMING_TEST_ODV              0
#define UQ_TGA_TEST_RUN_GRAD_TEST_ODV                0
#define UQ_TGA_TEST_RUN_OPTIMIZATION_TEST_ODV        0
#define UQ_TGA_TEST_GUESS_A_ODV                      2.5910e+11
#define UQ_TGA_TEST_GUESS_E_ODV                      2.0090e+05
#define UQ_TGA_TEST_COMPUTE_HESSIAN_ODV              1
#define UQ_TGA_TEST_RELATIVE_FD_STEP_ODV             1.e-8
#define UQ_TGA_TEST_NEWTON_MAX_ITERS_ODV             20
#define UQ_TGA_TEST_NEWTON_ABS_TOL_ODV               1.e-7
#define UQ_TGA_TEST_CRITICAL_TEMPERATURE_ODV         0
#define UQ_TGA_TEST_WRITE_OUTPUT_ODV                 1

class uqTgaTestOptionsClass
{
public:
  uqTgaTestOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
 ~uqTgaTestOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  std::string  m_outerPrefixName;
  bool         m_runTempTimeTest;
  bool         m_runTimingTest;
  bool         m_runGradTest;
  bool         m_runOptimizationTest;
  double       m_guessA;
  double       m_guessE;
  bool         m_computeHessian;
  double       m_relativeFDStep;
  unsigned int m_NewtonMaxIters;
  double       m_NewtonAbsTol;
  double       m_criticalTemperature;
  bool         m_writeOutput;

private:
  void   defineMyOptions  (po::options_description& optionsDesc);
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  std::string              m_prefix;
  po::options_description* m_optionsDesc;

  std::string              m_option_help;
  std::string              m_option_outerPrefixName;
  std::string              m_option_runTempTimeTest;
  std::string              m_option_runTimingTest;
  std::string              m_option_runGradTest;
  std::string              m_option_runOptimizationTest;
  std::string              m_option_guessA;
  std::string              m_option_guessE;
  std::string              m_option_computeHessian;
  std::string              m_option_relativeFDStep;
  std::string              m_option_NewtonMaxIters;
  std::string              m_option_NewtonAbsTol;
  std::string              m_option_criticalTemperature;
  std::string              m_option_writeOutput;
};

std::ostream& operator<<(std::ostream& os, const uqTgaTestOptionsClass& obj);

uqTgaTestOptionsClass::uqTgaTestOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix)
  :
  m_outerPrefixName                (UQ_TGA_TEST_OUTER_PREFIX_NAME_ODV    ),
  m_runTempTimeTest                (UQ_TGA_TEST_RUN_TEMP_TIME_TEST_ODV   ),
  m_runTimingTest                  (UQ_TGA_TEST_RUN_TIMINING_TEST_ODV    ),
  m_runGradTest                    (UQ_TGA_TEST_RUN_GRAD_TEST_ODV        ),
  m_runOptimizationTest            (UQ_TGA_TEST_RUN_OPTIMIZATION_TEST_ODV),
  m_guessA                         (UQ_TGA_TEST_GUESS_A_ODV              ),
  m_guessE                         (UQ_TGA_TEST_GUESS_E_ODV              ),
  m_computeHessian                 (UQ_TGA_TEST_COMPUTE_HESSIAN_ODV      ),
  m_relativeFDStep                 (UQ_TGA_TEST_RELATIVE_FD_STEP_ODV     ),
  m_NewtonMaxIters                 (UQ_TGA_TEST_NEWTON_MAX_ITERS_ODV     ),
  m_NewtonAbsTol                   (UQ_TGA_TEST_NEWTON_ABS_TOL_ODV       ),
  m_criticalTemperature            (UQ_TGA_TEST_CRITICAL_TEMPERATURE_ODV ),
  m_writeOutput                    (UQ_TGA_TEST_WRITE_OUTPUT_ODV         ),
  m_env                            (env),
  m_prefix                         ((std::string)(prefix) + "tests_"),
  m_optionsDesc                    (new po::options_description("TGA test options")),
  m_option_help                    (m_prefix + "help"                    ),
  m_option_outerPrefixName         (m_prefix + "outerPrefixName"         ),
  m_option_runTempTimeTest         (m_prefix + "runTempTimeTest"         ),
  m_option_runTimingTest           (m_prefix + "runTimingTest"           ),
  m_option_runGradTest             (m_prefix + "runGradTest"             ),
  m_option_runOptimizationTest     (m_prefix + "runOptimizationTest"     ),
  m_option_guessA                  (m_prefix + "guessA"                  ),
  m_option_guessE                  (m_prefix + "guessE"                  ),
  m_option_computeHessian          (m_prefix + "computeHessian"          ),
  m_option_relativeFDStep          (m_prefix + "relativeFDStep"          ),
  m_option_NewtonMaxIters          (m_prefix + "NewtonMaxIters"          ),
  m_option_NewtonAbsTol            (m_prefix + "NewtonAbsTol"            ),
  m_option_criticalTemperature     (m_prefix + "criticalTemperature"     ),
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
    (m_option_runTimingTest.c_str(),           po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_TIMING_TEST_ODV),              "run timing() test"                        )
    (m_option_runGradTest.c_str(),             po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_GRAD_TEST_ODV),                "run grad() test"                          )
    (m_option_runOptimizationTest.c_str(),     po::value<bool        >()->default_value(UQ_TGA_TEST_RUN_OPTIMIZATION_TEST_ODV),        "run optimization() test"                  )
    (m_option_guessA.c_str(),                  po::value<double      >()->default_value(UQ_TGA_TEST_GUESS_A_ODV),                      "guess A"                                  )
    (m_option_guessE.c_str(),                  po::value<double      >()->default_value(UQ_TGA_TEST_GUESS_E_ODV),                      "guess E"                                  )
    (m_option_computeHessian.c_str(),          po::value<bool        >()->default_value(UQ_TGA_TEST_COMPUTE_HESSIAN_ODV),              "computeHessian"                           )
    (m_option_relativeFDStep.c_str(),          po::value<double      >()->default_value(UQ_TGA_TEST_RELATIVE_FD_STEP_ODV),             "relativeFDStep"                           )
    (m_option_NewtonMaxIters.c_str(),          po::value<unsigned int>()->default_value(UQ_TGA_TEST_NEWTON_MAX_ITERS_ODV),             "Newton max number of iterations"          )
    (m_option_NewtonAbsTol.c_str(),            po::value<double      >()->default_value(UQ_TGA_TEST_NEWTON_ABS_TOL_ODV),               "Newton absolute tolerance"                )
    (m_option_criticalTemperature.c_str(),     po::value<double      >()->default_value(UQ_TGA_TEST_CRITICAL_TEMPERATURE_ODV),         "critical temperature"                     )
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
    m_outerPrefixName = m_env.allOptionsMap()[m_option_outerPrefixName.c_str()].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_runTempTimeTest.c_str())) {
    m_runTempTimeTest = m_env.allOptionsMap()[m_option_runTempTimeTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_runTimingTest.c_str())) {
    m_runTimingTest = m_env.allOptionsMap()[m_option_runTimingTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_runGradTest.c_str())) {
    m_runGradTest = m_env.allOptionsMap()[m_option_runGradTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_runOptimizationTest.c_str())) {
    m_runOptimizationTest = m_env.allOptionsMap()[m_option_runOptimizationTest.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_guessA.c_str())) {
    m_guessA = m_env.allOptionsMap()[m_option_guessA.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_guessE.c_str())) {
    m_guessE = m_env.allOptionsMap()[m_option_guessE.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_computeHessian.c_str())) {
    m_computeHessian = m_env.allOptionsMap()[m_option_computeHessian.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_relativeFDStep.c_str())) {
    m_relativeFDStep = m_env.allOptionsMap()[m_option_relativeFDStep.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_NewtonMaxIters.c_str())) {
    m_NewtonMaxIters = m_env.allOptionsMap()[m_option_NewtonMaxIters.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_NewtonAbsTol.c_str())) {
    m_NewtonAbsTol = m_env.allOptionsMap()[m_option_NewtonAbsTol.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_criticalTemperature.c_str())) {
    m_criticalTemperature = m_env.allOptionsMap()[m_option_criticalTemperature.c_str()].as<double>();
    globalTgaCriticalTemperature = m_criticalTemperature;
  }

  if (m_env.allOptionsMap().count(m_option_writeOutput.c_str())) {
    m_writeOutput = m_env.allOptionsMap()[m_option_writeOutput.c_str()].as<bool>();
  }

  return;
}

void
uqTgaTestOptionsClass::print(std::ostream& os) const
{
  os <<         m_option_outerPrefixName          << " = " << m_outerPrefixName
     << "\n" << m_option_runTempTimeTest          << " = " << m_runTempTimeTest
     << "\n" << m_option_runTimingTest            << " = " << m_runTimingTest
     << "\n" << m_option_runGradTest              << " = " << m_runGradTest
     << "\n" << m_option_runOptimizationTest      << " = " << m_runOptimizationTest
     << "\n" << m_option_guessA                   << " = " << m_guessA
     << "\n" << m_option_guessE                   << " = " << m_guessE
     << "\n" << m_option_computeHessian           << " = " << m_computeHessian
     << "\n" << m_option_relativeFDStep           << " = " << m_relativeFDStep
     << "\n" << m_option_NewtonMaxIters           << " = " << m_NewtonMaxIters
     << "\n" << m_option_NewtonAbsTol             << " = " << m_NewtonAbsTol   
     << "\n" << m_option_criticalTemperature      << " = " << m_criticalTemperature   
     << "\n" << m_option_writeOutput              << " = " << m_writeOutput
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

  const uqBase1D1DFunctionClass* tempFunction;
  uqSampled1D1DFunctionClass*    refW;
  std::ofstream*                 ofs;
  uqSampled1D1DFunctionClass*    continuousWeightFunction;
  uqSampled1D1DFunctionClass*    tildeContinuousWeightFunction;
  uqDeltaSeq1D1DFunctionClass*   deltaWeightFunction;
  uqDeltaSeq1D1DFunctionClass*   tildeDeltaWeightFunction;
};

uqTgaTestVarsStruct::uqTgaTestVarsStruct()
  :
  tempFunction                 (NULL),
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
  delete tempFunction;
  delete refW;
  delete ofs;
  delete continuousWeightFunction;
  delete tildeContinuousWeightFunction;
  delete deltaWeightFunction;
  delete tildeDeltaWeightFunction;
}
#endif // __UQ_TGA_TEST_OPTIONS_H__
