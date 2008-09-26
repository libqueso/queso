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
                      const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensity,  // Set in substep 1 in appls with a calibration prob.
                      const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunction, // Set in substep 2
                      P_M*                                          proposalCovMatrix,     // Set in substep 3
                      const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensity,    // Set in substep 3
                      const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGenerator); // Set in substep 3 // FIX ME: use such object
  uqCalibProblemClass(const uqEnvironmentClass&                     env,
                      const char*                                   prefix,
                      const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensity,  // Set in substep 1 in appls with a calibration prob.
                      const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunction, // Set in substep 2
                      P_M*                                          proposalCovMatrix,     // Set in substep 3
                      const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensity,    // Set in substep 3
                      const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGenerator); // Set in substep 3 // FIX ME: use such object
 ~uqCalibProblemClass();

        void                                   prepareSolver             ();
        void                                   solve                     ();
        void                                   solve                     (const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensity);

  const uqParamSpaceClass           <P_V,P_M>& paramSpace                () const;
  const uqRandomVariableClass       <P_V,P_M>& solution                  () const;    

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

  const uqProbDensity_BaseClass           <P_V,P_M>* m_priorParamDensity;
        bool                                         m_userPriorDensityIsNull;
        uqDefault_M2lPriorRoutine_DataType<P_V,P_M>  m_m2lPriorRoutine_Data;
        P_V*                                         m_paramPriorMus;
        P_V*                                         m_paramPriorSigmas;
  const uqProbDensity_BaseClass           <P_V,P_M>& m_likelihoodFunction;
        P_M*                                         m_proposalCovMatrix;
  const uqProposalDensity_BaseClass       <P_V,P_M>* m_proposalDensity;
  const uqProposalGenerator_BaseClass     <P_V,P_M>* m_proposalGenerator;

        uqMarkovChainSGClass              <P_V,P_M>* m_mcSeqGenerator;
        uqProbDensity_BaseClass           <P_V,P_M>* m_solutionProbDensity;
        uqRealizer_BaseClass              <P_V,P_M>* m_solutionRealizer;
        uqRandomVariableClass             <P_V,P_M>* m_solution;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqCalibProblemClass<P_V,P_M>& obj);

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::uqCalibProblemClass(
  const uqEnvironmentClass&                     env,
  const char*                                   prefix,
  const uqParamSpaceClass            <P_V,P_M>& paramSpace,
  const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensity,
  const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunction,
  P_M*                                          proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensity,
  const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGenerator)
  :
  m_env                   (env),
  m_prefix                ((std::string)(prefix) + "cal_"),
  m_paramSpace            (&paramSpace),
  m_userSpacesAreNull     (false),
  m_optionsDesc           (new po::options_description("UQ Calibration Problem")),
  m_option_help           (m_prefix + "help"  ),
  m_option_solver         (m_prefix + "solver"),
  m_solverString          (UQ_CALIB_PROBLEM_SOLVER_ODV),
  m_priorParamDensity     (priorParamDensity),
  m_userPriorDensityIsNull(priorParamDensity == NULL),
  m_paramPriorMus         (NULL),
  m_paramPriorSigmas      (NULL),
  m_likelihoodFunction    (likelihoodFunction),
  m_proposalCovMatrix     (proposalCovMatrix),
  m_proposalDensity       (proposalDensity),
  m_proposalGenerator     (proposalGenerator),
  m_mcSeqGenerator        (NULL),
  m_solutionProbDensity   (NULL),
  m_solutionRealizer      (NULL),
  m_solution              (NULL)
{
  commonConstructor();
}

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::uqCalibProblemClass(
  const uqEnvironmentClass&                     env,
  const char*                                   prefix,
  const uqProbDensity_BaseClass      <P_V,P_M>* priorParamDensity,
  const uqProbDensity_BaseClass      <P_V,P_M>& likelihoodFunction,
  P_M*                                          proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>* proposalDensity,
  const uqProposalGenerator_BaseClass<P_V,P_M>* proposalGenerator)
  :
  m_env                   (env),
  m_prefix                ((std::string)(prefix) + "cal_"),
  m_paramSpace            (new uqParamSpaceClass<P_V,P_M>(m_env,m_prefix.c_str())),
  m_userSpacesAreNull     (true),
  m_optionsDesc           (new po::options_description("UQ Calibration Problem")),
  m_option_help           (m_prefix + "help"  ),
  m_option_solver         (m_prefix + "solver"),
  m_solverString          (UQ_CALIB_PROBLEM_SOLVER_ODV),
  m_priorParamDensity     (priorParamDensity),
  m_userPriorDensityIsNull(priorParamDensity    == NULL),
  m_paramPriorMus         (NULL),
  m_paramPriorSigmas      (NULL),
  m_likelihoodFunction    (likelihoodFunction),
  m_proposalCovMatrix     (proposalCovMatrix),
  m_proposalDensity       (proposalDensity),
  m_proposalGenerator     (proposalGenerator),
  m_mcSeqGenerator        (NULL),
  m_solutionProbDensity   (NULL),
  m_solutionRealizer      (NULL),
  m_solution              (NULL)
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

    m_priorParamDensity = new uqRoutineProbDensity_Class<P_V,P_M>(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                     (void *) &m_m2lPriorRoutine_Data,
                                                                     true); // the routine computes [-2.*ln(Likelihood)]
  }

  m_solutionProbDensity = new uqBayesianProbDensity_Class<P_V,P_M>(m_priorParamDensity,
                                                                     &m_likelihoodFunction);

  // Instantiate the distribution calculator.
  m_mcSeqGenerator = new uqMarkovChainSGClass<P_V,P_M>(m_env,
                                                       m_prefix.c_str(),
                                                      *m_paramSpace,
                                                      *m_solutionProbDensity,
                                                       m_proposalCovMatrix,
                                                       m_proposalDensity,
                                                       m_proposalGenerator);

  if (m_env.rank() == 0) std::cout << "Leaving uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;

  return;
}

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::~uqCalibProblemClass()
{
  if (m_solution           ) delete m_solution;
  if (m_solutionRealizer   ) delete m_solutionRealizer;
  if (m_solutionProbDensity) delete m_solutionProbDensity;
  if (m_mcSeqGenerator     ) delete m_mcSeqGenerator;

  if (m_userPriorDensityIsNull) { 
    delete m_priorParamDensity;
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
  m_mcSeqGenerator->calculateDistributions();

  if (m_solutionRealizer) delete m_solutionRealizer;
  m_solutionRealizer = new uqRealizer_BaseClass<P_V,P_M>(&(m_mcSeqGenerator->chain()));

  m_solution = new uqRandomVariableClass<P_V,P_M>(m_solutionProbDensity,
                                                  m_solutionRealizer);

  return;
}

template <class P_V,class P_M>
void
uqCalibProblemClass<P_V,P_M>::solve(const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensity)
{
  m_mcSeqGenerator->calculateDistributions(priorParamDensity);

  if (m_solutionRealizer) delete m_solutionRealizer;
  m_solutionRealizer = new uqRealizer_BaseClass<P_V,P_M>(&(m_mcSeqGenerator->chain()));

  m_solution = new uqRandomVariableClass<P_V,P_M>(m_solutionProbDensity,
                                                  m_solutionRealizer);
  return;
}

template <class P_V,class P_M>
const uqParamSpaceClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::paramSpace() const
{
  return *m_paramSpace;
}

#if 0
template <class P_V,class P_M>
const uqProbDensity_BaseClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::solutionProbDensity() const
{
  UQ_FATAL_TEST_MACRO(m_solutionProbDensity == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M>::solutionProbDensity()",
                      "solution density is being requested but it has not been created yet");

  return *m_solutionProbDensity;
}

template <class P_V,class P_M>
const uqRealizer_BaseClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::solutionSampleGenerator() const
{
  UQ_FATAL_TEST_MACRO(m_solutionRealizer == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M>::solutionSampleGenerator()",
                      "solution realizer is being requested but it has not been created yet");

  return *m_solutionRealizer;
}
#endif

template <class P_V,class P_M>
const uqRandomVariableClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M>::solution() const
{
  UQ_FATAL_TEST_MACRO(m_solution == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M>::solution()",
                      "solution is being requested but it has not been created yet");

  return *m_solution;
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
