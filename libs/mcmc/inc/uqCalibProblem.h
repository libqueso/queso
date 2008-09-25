/* uq/libs/mcmc/inc/uqCalibProblem.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_CALIB_PROBLEM_H__
#define __UQ_CALIB_PROBLEM_H__

#include <uqParamSpace.h>

#include <uqProbDensity.h>       // For substep 1 in appls with a calibration problem
#include <uqProposalDensity.h>   // For substep 3
#include <uqProposalGenerator.h> // For substep 3

#include <uqDefaultPrior.h>
#include <uqMarkovChainSG1.h>
#include <uqRandomVariable.h>

// _ODV = option default value
#define UQ_CALIB_PROBLEM_SOLVER_ODV "bayes_mc" // Bayesian formula + Markov Chain

template <class P_V,class P_M>
class uqCalibProblemClass
{
public:
  uqCalibProblemClass(const uqEnvironmentClass&                     env,
                      const char*                                   prefix,
                      const uqParamSpaceClass            <P_V,P_M>& paramSpace,
                      const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensityObj,  // Set in substep 1 in appls with a calibration prob.
                      const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunctionObj, // Set in substep 2
                      P_M*                                          proposalCovMatrix,     // Set in substep 3
                      const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensityObj,    // Set in substep 3
                      const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGeneratorObj); // Set in substep 3 // FIX ME: use such object
  uqCalibProblemClass(const uqEnvironmentClass&                     env,
                      const char*                                   prefix,
                      const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensityObj,  // Set in substep 1 in appls with a calibration prob.
                      const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunctionObj, // Set in substep 2
                      P_M*                                          proposalCovMatrix,     // Set in substep 3
                      const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensityObj,    // Set in substep 3
                      const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGeneratorObj); // Set in substep 3 // FIX ME: use such object
 ~uqCalibProblemClass();

  const uqParamSpaceClass           <P_V,P_M>& paramSpace                () const;

        void                                   prepareSolver             ();
        void                                   solve                     ();
        void                                   solve                     (const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj);

  const uqRandomVariableClass       <P_V,P_M>& solution                  () const;    
  const uqProbDensity_BaseClass     <P_V,P_M>& solutionProbDensityObj    () const;
  const uqSampleGenerator_BaseClass <P_V,P_M>& solutionSampleGeneratorObj() const;

        void                                   print                     (std::ostream& os) const;

private:
        void commonConstructor();
        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                          m_env;
        std::string                                  m_prefix;
  const uqParamSpaceClass                 <P_V,P_M>* m_paramSpace;
        bool                                         m_userSpacesAreNull;

        po::options_description*                     m_optionsDesc;
        std::string                                  m_option_help;
        std::string                                  m_option_solver;

	std::string                                  m_solverString;

  const uqProbDensity_BaseClass           <P_V,P_M>* m_priorParamDensityObj;
        bool                                         m_userPriorDensityIsNull;
        uqDefault_M2lPriorRoutine_DataType<P_V,P_M>  m_m2lPriorRoutine_Data;
        P_V*                                         m_paramPriorMus;
        P_V*                                         m_paramPriorSigmas;
  const uqProbDensity_BaseClass           <P_V,P_M>& m_likelihoodFunctionObj;
        P_M*                                         m_proposalCovMatrix;
  const uqProposalDensity_BaseClass       <P_V,P_M>* m_proposalDensityObj;
  const uqProposalGenerator_BaseClass     <P_V,P_M>* m_proposalGeneratorObj;

        uqMarkovChainSGClass              <P_V,P_M>* m_mcSampler;
        uqProbDensity_BaseClass           <P_V,P_M>* m_solutionProbDensityObj;
        uqSampleGenerator_BaseClass       <P_V,P_M>* m_solutionSampleGeneratorObj;
        uqRandomVariableClass             <P_V,P_M>* m_solutionObj;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqCalibProblemClass<P_V,P_M>& obj);

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::uqCalibProblemClass(
  const uqEnvironmentClass&                     env,
  const char*                                   prefix,
  const uqParamSpaceClass            <P_V,P_M>& paramSpace,
  const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensityObj,
  const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunctionObj,
  P_M*                                          proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensityObj,
  const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGeneratorObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "cal_"),
  m_paramSpace                (&paramSpace),
  m_userSpacesAreNull         (false),
  m_optionsDesc               (new po::options_description("UQ Calibration Problem")),
  m_option_help               (m_prefix + "help"  ),
  m_option_solver             (m_prefix + "solver"),
  m_solverString              (UQ_CALIB_PROBLEM_SOLVER_ODV),
  m_priorParamDensityObj      (priorParamDensityObj),
  m_userPriorDensityIsNull    (priorParamDensityObj == NULL),
  m_paramPriorMus             (NULL),
  m_paramPriorSigmas          (NULL),
  m_likelihoodFunctionObj     (likelihoodFunctionObj),
  m_proposalCovMatrix         (proposalCovMatrix),
  m_proposalDensityObj        (proposalDensityObj),
  m_proposalGeneratorObj      (proposalGeneratorObj),
  m_mcSampler                 (NULL),
  m_solutionProbDensityObj    (NULL),
  m_solutionSampleGeneratorObj(NULL),
  m_solutionObj               (NULL)
{
  commonConstructor();
}

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::uqCalibProblemClass(
  const uqEnvironmentClass&                     env,
  const char*                                   prefix,
  const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensityObj,
  const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunctionObj,
  P_M*                                          proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensityObj,
  const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGeneratorObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "cal_"),
  m_paramSpace                (new uqParamSpaceClass<P_V,P_M>(m_env,m_prefix.c_str())),
  m_userSpacesAreNull         (true),
  m_optionsDesc               (new po::options_description("UQ Calibration Problem")),
  m_option_help               (m_prefix + "help"  ),
  m_option_solver             (m_prefix + "solver"),
  m_solverString              (UQ_CALIB_PROBLEM_SOLVER_ODV),
  m_priorParamDensityObj      (priorParamDensityObj),
  m_userPriorDensityIsNull    (priorParamDensityObj == NULL),
  m_paramPriorMus             (NULL),
  m_paramPriorSigmas          (NULL),
  m_likelihoodFunctionObj     (likelihoodFunctionObj),
  m_proposalCovMatrix         (proposalCovMatrix),
  m_proposalDensityObj        (proposalDensityObj),
  m_proposalGeneratorObj      (proposalGeneratorObj),
  m_mcSampler                 (NULL),
  m_solutionProbDensityObj    (NULL),
  m_solutionSampleGeneratorObj(NULL),
  m_solutionObj               (NULL)
{
  commonConstructor();
}

template <class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::commonConstructor()
{
  if (m_env.rank() == 0) std::cout << "Entering uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_userPriorDensityIsNull) {
    m_paramPriorMus    = new P_V(m_paramSpace->priorMuValues   ());
    m_paramPriorSigmas = new P_V(m_paramSpace->priorSigmaValues());
    m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
    m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;

    m_priorParamDensityObj = new uqRoutineProbDensity_Class<P_V,P_M>(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                     (void *) &m_m2lPriorRoutine_Data,
                                                                     true); // the routine computes [-2.*ln(Likelihood)]
  }

  m_solutionProbDensityObj = new uqBayesianProbDensity_Class<P_V,P_M>(m_priorParamDensityObj,
                                                                     &m_likelihoodFunctionObj);

  // Instantiate the distribution calculator.
  m_mcSampler = new uqMarkovChainSGClass<P_V,P_M>(m_env,
                                                  m_prefix.c_str(),
                                                 *m_paramSpace,
                                                 *m_solutionProbDensityObj,
                                                  m_proposalCovMatrix,
                                                  m_proposalDensityObj,
                                                  m_proposalGeneratorObj);

  if (m_env.rank() == 0) std::cout << "Leaving uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;

  return;
}

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::~uqCalibProblemClass()
{
  if (m_solutionObj               ) delete m_solutionObj;
  if (m_solutionSampleGeneratorObj) delete m_solutionSampleGeneratorObj;
  if (m_solutionProbDensityObj    ) delete m_solutionProbDensityObj;
  if (m_mcSampler                 ) delete m_mcSampler;

  if (m_userPriorDensityIsNull) { 
    delete m_priorParamDensityObj;
    delete m_paramPriorSigmas;
    delete m_paramPriorMus;
  }

  if (m_optionsDesc) delete m_optionsDesc;

  if (m_userSpacesAreNull) {
    if (m_paramSpace) delete m_paramSpace;
  }
}

template<class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                         "produce help message for calibration problem")
    (m_option_solver.c_str(), po::value<std::string>()->default_value(UQ_CALIB_PROBLEM_SOLVER_ODV), "algorithm for calibration"                   )
  ;

  return;
}

template<class P_V,class P_M>
void
  uqCalibProblemClass<P_V,P_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_solver.c_str())) {
    m_solverString = m_env.allOptionsMap()[m_option_solver.c_str()].as<std::string>();
  }

  return;
}

template <class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::prepareSolver()
{
  return;
}

template <class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::solve()
{
  m_mcSampler->calculateDistributions();

  if (m_solutionSampleGeneratorObj) delete m_solutionSampleGeneratorObj;
  m_solutionSampleGeneratorObj = new uqSampleGenerator_BaseClass<P_V,P_M>(&(m_mcSampler->chain()));

  return;
}

template <class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::solve(const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj)
{
  m_mcSampler->calculateDistributions(priorParamDensityObj);

  if (m_solutionSampleGeneratorObj) delete m_solutionSampleGeneratorObj;
  m_solutionSampleGeneratorObj = new uqSampleGenerator_BaseClass<P_V,P_M>(&(m_mcSampler->chain()));

  return;
}

template <class P_V,class P_M>
const uqParamSpaceClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M>
const uqProbDensity_BaseClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::solutionProbDensityObj() const
{
  UQ_FATAL_TEST_MACRO(m_solutionProbDensityObj == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M>::solutionProbDensityObj()",
                      "posterior param density object is being requested but it has not been created yet");

  return *m_solutionProbDensityObj;
}

template <class P_V,class P_M>
const uqSampleGenerator_BaseClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::solutionSampleGeneratorObj() const
{
  UQ_FATAL_TEST_MACRO(m_solutionSampleGeneratorObj == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M>::solutionSampleGeneratorObj()",
                      "posterior param generator object is being requested but it has not been created yet");

  return *m_solutionSampleGeneratorObj;
}

template <class P_V,class P_M>
const uqRandomVariableClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::solution() const
{
  UQ_FATAL_TEST_MACRO(m_solutionObj == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M>::solutionSampleGeneratorObj()",
                      "solution is being requested but it has not been created yet");

  return *m_solutionObj;
}

template <class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_solver  << " = " << m_solverString;
  os << std::endl;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqCalibProblemClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_CALIB_PROBLEM_H__
