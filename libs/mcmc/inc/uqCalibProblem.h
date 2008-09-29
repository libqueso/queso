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

#include <uqProposalDensity.h>   // For substep 3
#include <uqProposalGenerator.h> // For substep 3

#include <uqDefaultPrior.h>
#include <uqMarkovChainSG1.h>
#include <uqVectorRV.h>

// _ODV = option default value
#define UQ_CALIB_PROBLEM_SOLVER_ODV "bayes_mc" // Bayesian formula + Markov Chain

template <class P_V,class P_M>
class uqCalibProblemClass
{
public:
  uqCalibProblemClass(const char*                                  prefix,
                      const uqBaseVectorRVClass         <P_V,P_M>& priorRv,
                      const uqBaseVectorProbDensityClass<P_V,P_M>& likelihoodFunction,
                            uqBaseVectorRVClass         <P_V,P_M>& postRv);
 ~uqCalibProblemClass();

        void solveWithBayesMarkovChain(void* transitionKernel);

        void print                    (std::ostream& os) const;

private:
        void defineMyOptions          (po::options_description& optionsDesc);
        void getMyOptionValues        (po::options_description& optionsDesc);

  const uqEnvironmentClass&                          m_env;
        std::string                                  m_prefix;

        po::options_description*                     m_optionsDesc;
        std::string                                  m_option_help;
        std::string                                  m_option_solver;

	std::string                                  m_solverString;

  const uqBaseVectorRVClass               <P_V,P_M>& m_priorRv;
  const uqBaseVectorProbDensityClass      <P_V,P_M>& m_likelihoodFunction;
        uqBaseVectorRVClass               <P_V,P_M>& m_postRv;

  const uqBaseVectorProbDensityClass      <P_V,P_M>* m_priorParamDensity;
        bool                                         m_userPriorDensityIsNull;
        uqDefault_M2lPriorRoutine_DataType<P_V,P_M>  m_m2lPriorRoutine_Data;
        P_V*                                         m_paramPriorMus;
        P_V*                                         m_paramPriorSigmas;
        P_M*                                         m_proposalCovMatrix;
  const uqProposalDensity_BaseClass       <P_V,P_M>* m_proposalDensity;
  const uqProposalGenerator_BaseClass     <P_V,P_M>* m_proposalGenerator;

        uqMarkovChainSGClass              <P_V,P_M>* m_mcSeqGenerator;
        uqBaseVectorProbDensityClass      <P_V,P_M>* m_solutionProbDensity;
        uqBaseVectorRealizerClass         <P_V,P_M>* m_solutionRealizer;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqCalibProblemClass<P_V,P_M>& obj);

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::uqCalibProblemClass(
  const char*                                  prefix,
  const uqBaseVectorRVClass         <P_V,P_M>& priorRv,
  const uqBaseVectorProbDensityClass<P_V,P_M>& likelihoodFunction,
        uqBaseVectorRVClass         <P_V,P_M>& postRv)
  :
  m_env                   (priorRv.env()),
  m_prefix                ((std::string)(prefix) + "cal_"),
  m_optionsDesc           (new po::options_description("UQ Calibration Problem")),
  m_option_help           (m_prefix + "help"  ),
  m_option_solver         (m_prefix + "solver"),
  m_solverString          (UQ_CALIB_PROBLEM_SOLVER_ODV),
  m_priorRv               (priorRv),
  m_likelihoodFunction    (likelihoodFunction),
  m_postRv                (postRv),
  m_priorParamDensity     (NULL),
  m_userPriorDensityIsNull(true),
  m_paramPriorMus         (NULL),
  m_paramPriorSigmas      (NULL),
  m_proposalCovMatrix     (NULL),
  m_proposalDensity       (NULL),
  m_proposalGenerator     (NULL),
  m_mcSeqGenerator        (NULL),
  m_solutionProbDensity   (NULL),
  m_solutionRealizer      (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_env.rank() == 0) std::cout << "Leaving uqCalibProblemClass<P_V,P_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;

  return;
}

template <class P_V,class P_M>
uqCalibProblemClass<P_V,P_M>::~uqCalibProblemClass()
{
  if (m_solutionRealizer   ) delete m_solutionRealizer;
  if (m_solutionProbDensity) delete m_solutionProbDensity;
  if (m_mcSeqGenerator     ) delete m_mcSeqGenerator;

  if (m_userPriorDensityIsNull) { 
    delete m_priorParamDensity;
    delete m_paramPriorSigmas;
    delete m_paramPriorMus;
  }

  if (m_optionsDesc) delete m_optionsDesc;
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
uqCalibProblemClass<P_V,P_M>::solveWithBayesMarkovChain(void* transitionKernel)
{
  if (m_solutionRealizer   ) delete m_solutionRealizer;
  if (m_solutionProbDensity) delete m_solutionProbDensity;
  if (m_mcSeqGenerator     ) delete m_mcSeqGenerator;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqCalibProblemClass<P_V,P_M>::solveWithBayesMarkovChain()"
              << ", m_priorRv contents: " << m_priorRv
              << std::endl;
  }

  // Bayesian step
  m_solutionProbDensity = new uqBayesianVectorProbDensityClass<P_V,P_M>(&m_priorRv.probDensity(),
                                                                        &m_likelihoodFunction);
  m_postRv.setProbDensity(*m_solutionProbDensity);

  m_postRv.setComponent(0, // GAMBIARRA
                        m_priorRv.component(0).minValue(),
                        m_priorRv.component(0).maxValue(),
                        0.,
                        INFINITY);
  m_postRv.setComponent(1, // GAMBIARRA
                        m_priorRv.component(1).minValue(),
                        m_priorRv.component(1).maxValue(),
                        0.,
                        INFINITY);

  // Markov Chain step
  m_mcSeqGenerator = new uqMarkovChainSGClass<P_V,P_M>(m_prefix.c_str(),
                                                       m_postRv,
                                                       m_proposalCovMatrix,
                                                       m_proposalDensity,
                                                       m_proposalGenerator);
  m_mcSeqGenerator->generateSequence(m_postRv.chain());
  m_solutionRealizer = new uqBaseVectorRealizerClass<P_V,P_M>(&(m_postRv.chain()));
  m_postRv.setRealizer(*m_solutionRealizer);

  return;
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
