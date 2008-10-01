/* uq/libs/mcmc/inc/uqPropagProblem.h
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

#ifndef __UQ_PROPAG_PROBLEM_H__
#define __UQ_PROPAG_PROBLEM_H__

#include <uqVectorFunction.h> // For substep 4 (or 1) in appls. with propagation

#include <uqMonteCarloSG.h>
#include <uqVectorProbDensity.h>
#include <uqVectorRV.h>

#undef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION

// _ODV = option default value
#define UQ_PROPAG_PROBLEM_SKIP_SOLUTION_ODV 0
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
#define UQ_PROPAG_PROBLEM_SOLVER_ODV        "mc_kde" // Monte Carlo + Kernel Density Estimator
#endif

template <class P_V,class P_M,class Q_V,class Q_M>
class uqPropagProblemClass
{
public:
  uqPropagProblemClass(const char*                                   prefix,
                       const uqBaseVectorRVClass  <P_V,P_M>&         paramRv,
                       const uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction,
                             uqBaseVectorRVClass  <Q_V,Q_M>&         qoiRV);
 ~uqPropagProblemClass();

        void                                  solveWithMonteCarloKde();

        void                                  print          (std::ostream& os) const;

private:
        void commonConstructor();
        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                            m_env;
        std::string                                    m_prefix;

        po::options_description*                       m_optionsDesc;
        std::string                                    m_option_help;
	std::string                                    m_option_skipSolution;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
        std::string                                    m_option_solver;
#endif

        bool                                           m_skipSolution;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
	std::string                                    m_solverString;
#endif

  const uqBaseVectorRVClass         <P_V,P_M>&         m_paramRv;
  const uqVectorFunctionClass       <P_V,P_M,Q_V,Q_M>& m_qoiFunction;
        uqBaseVectorRVClass         <Q_V,Q_M>&         m_qoiRv;

        uqBaseVectorProbDensityClass<Q_V,Q_M>*         m_solutionProbDensity;
        uqBaseVectorRealizerClass   <Q_V,Q_M>*         m_solutionRealizer;

        uqMonteCarloSGClass         <P_V,P_M,Q_V,Q_M>* m_mcSeqGenerator;
        uqSequenceOfVectorsClass    <Q_V,Q_M>*         m_chain1;
        uqArrayOfSequencesClass     <Q_V,Q_M>*         m_chain2;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::uqPropagProblemClass(
  const char*                                   prefix,
  const uqBaseVectorRVClass  <P_V,P_M>&         paramRv,
  const uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction,
        uqBaseVectorRVClass  <Q_V,Q_M>&         qoiRv)
  :
  m_env                (paramRv.env()),
  m_prefix             ((std::string)(prefix) + "pro_"),
  m_optionsDesc        (new po::options_description("UQ Propagation Problem")),
  m_option_help        (m_prefix + "help"        ),
  m_option_skipSolution(m_prefix + "skipSolution"),
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  m_option_solver      (m_prefix + "solver"      ),
#endif
  m_skipSolution       (UQ_PROPAG_PROBLEM_SKIP_SOLUTION_ODV),
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  m_solverString       (UQ_PROPAG_PROBLEM_SOLVER_ODV),
#endif
  m_paramRv            (paramRv),
  m_qoiFunction        (qoiFunction),
  m_qoiRv              (qoiRv),
  m_solutionProbDensity(NULL),
  m_solutionRealizer   (NULL),
  m_mcSeqGenerator     (NULL),
  m_chain1             (NULL),
  m_chain2             (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_env.rank() == 0) std::cout << "Leaving uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::~uqPropagProblemClass()
{
  if (m_chain1) {
    m_chain1->clear();
    delete m_chain1;
  }
  if (m_chain2) {
    m_chain2->clear();
    delete m_chain2;
  }
  if (m_mcSeqGenerator     ) delete m_mcSeqGenerator;
  if (m_solutionRealizer   ) delete m_solutionRealizer;
  if (m_solutionProbDensity) delete m_solutionProbDensity;
  if (m_optionsDesc        ) delete m_optionsDesc;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                       "produce help message for propagation problem")
    (m_option_skipSolution.c_str(), po::value<bool       >()->default_value(UQ_PROPAG_PROBLEM_SKIP_SOLUTION_ODV), "skip solution process"                       )
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
    (m_option_solver.c_str(),       po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_SOLVER_ODV       ), "algorithm for propagation"                   )
#endif
  ;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
  uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_skipSolution.c_str())) {
    m_skipSolution = m_env.allOptionsMap()[m_option_skipSolution.c_str()].as<bool>();
  }

#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver.c_str())) {
    m_solverString = m_env.allOptionsMap()[m_option_solver.c_str()].as<std::string>();
  }
#endif

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarloKde()
{
  // Instantiate the distribution calculator.
  m_mcSeqGenerator = new uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>(m_prefix.c_str(),
                                                              m_paramRv,
                                                              m_qoiFunction,
                                                              m_qoiRv);

  m_chain1 = new uqSequenceOfVectorsClass<Q_V,Q_M>(m_qoiRv.imageSpace(),0);
  //m_chain2 = new uqArrayOfSequencesClass <Q_V,Q_M>(m_qoiRv.imageSpace(),0);

  m_mcSeqGenerator->generateSequence(*m_chain1);
  m_solutionRealizer = new uqSequentialVectorRealizerClass<Q_V,Q_M>(m_prefix.c_str(),
                                                                   *m_chain1);
  m_qoiRv.setRealizer(*m_solutionRealizer);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_skipSolution << " = " << m_skipSolution;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  os << "\n" << m_option_solver       << " = " << m_solverString;
#endif
  os << std::endl;
}

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROPAG_PROBLEM_H__
