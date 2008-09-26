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

#include <uqParamSpace.h>
#include <uqQoISpace.h>

#include <uqVectorFunction.h> // For substep 4 (or 1) in appls. with propagation

#include <uqMonteCarloSG.h>
#include <uqProbDensity.h>
#include <uqVectorRV.h>

// _ODV = option default value
#define UQ_PROPAG_PROBLEM_SOLVER_ODV "mc_kde" // Monte Carlo + Kernel Density Estimator

template <class P_V,class P_M,class Q_V,class Q_M>
class uqPropagProblemClass
{
public:
  uqPropagProblemClass(const uqEnvironmentClass&                     env,
                       const char*                                   prefix,
                       const uqParamSpaceClass    <P_V,P_M>&         paramSpace,
                       const uqQoISpaceClass      <Q_V,Q_M>&         qoiSpace,
                       const uqVectorRVClass      <P_V,P_M>*         paramRV,
                       const uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction); // Set in substep 4
  uqPropagProblemClass(const uqEnvironmentClass&                     env,
                       const char*                                   prefix,
                       const uqVectorRVClass      <P_V,P_M>*         paramRV,
                             unsigned int                            qoiSpaceDim,
                       const uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction); // Set in substep 4
 ~uqPropagProblemClass();

  const uqParamSpaceClass          <P_V,P_M>& paramSpace     () const;
  const uqQoISpaceClass            <Q_V,Q_M>& qoiSpace       () const;

        void                                  solve          ();
        void                                  solve          (const uqVectorRVClass<P_V,P_M>& paramRV);

//const uqProbDensity_BaseClass<Q_V,Q_M>& solutionProbDensity    () const;
//const uqVectorRVClass  <Q_V,Q_M>& solutionSampleGenerator() const;
  const uqVectorRVClass  <Q_V,Q_M>& solution                  () const;

        void                                  print          (std::ostream& os) const;

private:
        void commonConstructor();
        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                       m_env;
        std::string                               m_prefix;
  const uqParamSpaceClass      <P_V,P_M>*         m_paramSpace;
  const uqQoISpaceClass        <Q_V,Q_M>*         m_qoiSpace;
        bool                                      m_userSpacesAreNull;

        po::options_description*                  m_optionsDesc;
        std::string                               m_option_help;
        std::string                               m_option_solver;

	std::string                               m_solverString;

  const uqVectorRVClass        <P_V,P_M>*         m_paramRV;
  const uqVectorFunctionClass  <P_V,P_M,Q_V,Q_M>& m_qoiFunction;

        uqMonteCarloSGClass    <P_V,P_M,Q_V,Q_M>* m_mcSeqGenerator;
        uqProbDensity_BaseClass<Q_V,Q_M>*         m_solutionProbDensity;
        uqRealizer_BaseClass   <Q_V,Q_M>*         m_solutionRealizer;
        uqVectorRVClass        <Q_V,Q_M>*         m_solution;

};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::uqPropagProblemClass(
  const uqEnvironmentClass&                     env,
  const char*                                   prefix,
  const uqParamSpaceClass    <P_V,P_M>&         paramSpace,
  const uqQoISpaceClass      <Q_V,Q_M>&         qoiSpace,
  const uqVectorRVClass      <P_V,P_M>*         paramRV,
  const uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction)
  :
  m_env                (env),
  m_prefix             ((std::string)(prefix) + "pro_"),
  m_paramSpace         (&paramSpace),
  m_qoiSpace           (&qoiSpace),
  m_userSpacesAreNull  (false),
  m_optionsDesc        (new po::options_description("UQ Propagation Problem")),
  m_option_help        (m_prefix + "help"  ),
  m_option_solver      (m_prefix + "solver"),
  m_solverString       (UQ_PROPAG_PROBLEM_SOLVER_ODV),
  m_paramRV            (paramRV),
  m_qoiFunction        (qoiFunction),
  m_mcSeqGenerator     (NULL),
  m_solutionProbDensity(NULL),
  m_solutionRealizer   (NULL),
  m_solution           (NULL)
{
  commonConstructor();
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::uqPropagProblemClass(
  const uqEnvironmentClass&                     env,
  const char*                                   prefix,
  const uqVectorRVClass      <P_V,P_M>*         paramRV,
        unsigned int                            qoiSpaceDim,
  const uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction)
  :
  m_env                (env),
  m_prefix             ((std::string)(prefix) + "pro_"),
  m_paramSpace         (new uqParamSpaceClass<P_V,P_M>(m_env,m_prefix.c_str())),
  m_qoiSpace           (new uqQoISpaceClass  <Q_V,Q_M>(m_env,m_prefix.c_str())),
  m_userSpacesAreNull  (true),
  m_optionsDesc        (new po::options_description("UQ Propagation Problem")),
  m_option_help        (m_prefix + "help"  ),
  m_option_solver      (m_prefix + "solver"),
  m_solverString       (UQ_PROPAG_PROBLEM_SOLVER_ODV),
  m_paramRV            (paramRV),
  m_qoiFunction        (qoiFunction),
  m_mcSeqGenerator     (NULL),
  m_solutionProbDensity(NULL),
  m_solutionRealizer   (NULL),
  m_solution           (NULL)
{
  commonConstructor();
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::commonConstructor()
{
  if (m_env.rank() == 0) std::cout << "Entering uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  // Instantiate the distribution calculator.
  m_mcSeqGenerator = new uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                              m_prefix.c_str(),
                                                             *m_paramSpace,
                                                             *m_qoiSpace,
                                                              m_paramRV,
                                                              m_qoiFunction);

  if (m_env.rank() == 0) std::cout << "Leaving uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::~uqPropagProblemClass()
{
  if (m_solution           ) delete m_solution;
  if (m_solutionRealizer   ) delete m_solutionRealizer;
  if (m_solutionProbDensity) delete m_solutionProbDensity;
  if (m_mcSeqGenerator     ) delete m_mcSeqGenerator;

  if (m_optionsDesc        ) delete m_optionsDesc;
  if (m_userSpacesAreNull) {
    if (m_qoiSpace  ) delete m_qoiSpace;
    if (m_paramSpace) delete m_paramSpace;
  }
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                          "produce help message for propagation problem")
    (m_option_solver.c_str(), po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_SOLVER_ODV), "algorithm for propagation"                   )
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

  if (m_env.allOptionsMap().count(m_option_solver.c_str())) {
    m_solverString = m_env.allOptionsMap()[m_option_solver.c_str()].as<std::string>();
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solve()
{
  m_mcSeqGenerator->generateSequence();

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solve(const uqVectorRVClass<P_V,P_M>& paramRV)
{
  m_mcSeqGenerator->generateSequence(paramRV);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqParamSpaceClass<P_V,P_M>&
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqQoISpaceClass<Q_V,Q_M>&
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::qoiSpace() const
{
  return *m_qoiSpace;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqVectorRVClass<Q_V,Q_M>&
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solution() const
{
  UQ_FATAL_TEST_MACRO(m_solution == NULL,
                      m_env.rank(),
                      "uqPropagProblemClass<P_V,P_M>::solution()",
                      "solution is being requested but it has not been created yet");

  return *m_solution;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_solver << " = " << m_solverString;
  os << std::endl;
}

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROPAG_PROBLEM_H__
