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

#include <uqQoIFunction.h> // For substep 4 (or 1) in appls. with propagation

#include <uqMonteCarloSG.h>
#include <uqProbDensity.h>
#include <uqSampleGenerator.h>

// _ODV = option default value
#define UQ_PROPAG_PROBLEM_SOLVER_ODV "mc_kde" // Monte Carlo + Kernel Density Estimator

template <class P_V,class P_M,class Q_V,class Q_M>
class uqPropagProblemClass
{
public:
  uqPropagProblemClass(const uqEnvironmentClass&                           env,
                       const char*                                         prefix,
                       const uqParamSpaceClass          <P_V,P_M>&         paramSpace,
                       const uqQoISpaceClass            <Q_V,Q_M>&         qoiSpace,
                       const uqProbDensity_BaseClass    <P_V,P_M>*         paramDensityObj,
                       const uqSampleGenerator_BaseClass<P_V,P_M>*         paramGeneratorObj,
                       const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj); // Set in substep 4
  uqPropagProblemClass(const uqEnvironmentClass&                           env,
                       const char*                                         prefix,
                       const uqProbDensity_BaseClass    <P_V,P_M>*         paramDensityObj,
                       const uqSampleGenerator_BaseClass<P_V,P_M>*         paramGeneratorObj,
                       const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj); // Set in substep 4
 ~uqPropagProblemClass();

  const uqParamSpaceClass          <P_V,P_M>& paramSpace     () const;
  const uqQoISpaceClass            <Q_V,Q_M>& qoiSpace       () const;

        void                                  solve          ();
        void                                  solve          (const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj);

  const uqProbDensity_BaseClass    <Q_V,Q_M>& solutionProbDensityObj  () const;
  const uqSampleGenerator_BaseClass<Q_V,Q_M>& solutionSampleGeneratorObj() const;

        void                                  print          (std::ostream& os) const;

private:
        void commonConstructor();
        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                           m_env;
        std::string                                   m_prefix;
  const uqParamSpaceClass          <P_V,P_M>*         m_paramSpace;
  const uqQoISpaceClass            <Q_V,Q_M>*         m_qoiSpace;
        bool                                          m_userSpacesAreNull;

        po::options_description*                      m_optionsDesc;
        std::string                                   m_option_help;
        std::string                                   m_option_solver;

	std::string                                   m_solverString;

  const uqProbDensity_BaseClass    <P_V,P_M>*         m_paramDensityObj;
  const uqSampleGenerator_BaseClass<P_V,P_M>*         m_paramGeneratorObj;
  const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& m_qoiFunctionObj;

        uqMonteCarloSGClass        <P_V,P_M,Q_V,Q_M>* m_mcSampler;
        uqProbDensity_BaseClass    <Q_V,Q_M>*         m_solutionProbDensityObj;
        uqSampleGenerator_BaseClass<Q_V,Q_M>*         m_solutionSampleGeneratorObj;

};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::uqPropagProblemClass(
  const uqEnvironmentClass&                           env,
  const char*                                         prefix,
  const uqParamSpaceClass          <P_V,P_M>&         paramSpace,
  const uqQoISpaceClass            <Q_V,Q_M>&         qoiSpace,
  const uqProbDensity_BaseClass    <P_V,P_M>*         paramDensityObj,
  const uqSampleGenerator_BaseClass<P_V,P_M>*         paramGeneratorObj,
  const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "pro_"),
  m_paramSpace                (&paramSpace),
  m_qoiSpace                  (&qoiSpace),
  m_userSpacesAreNull         (false),
  m_optionsDesc               (new po::options_description("UQ Propagation Problem")),
  m_option_help               (m_prefix + "help"  ),
  m_option_solver             (m_prefix + "solver"),
  m_solverString              (UQ_PROPAG_PROBLEM_SOLVER_ODV),
  m_paramDensityObj           (paramDensityObj),
  m_paramGeneratorObj         (paramGeneratorObj),
  m_qoiFunctionObj            (qoiFunctionObj),
  m_mcSampler                 (NULL),
  m_solutionProbDensityObj    (NULL),
  m_solutionSampleGeneratorObj(NULL)
{
  commonConstructor();
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::uqPropagProblemClass(
  const uqEnvironmentClass&                           env,
  const char*                                         prefix,
  const uqProbDensity_BaseClass    <P_V,P_M>*         paramDensityObj,
  const uqSampleGenerator_BaseClass<P_V,P_M>*         paramGeneratorObj,
  const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "pro_"),
  m_paramSpace                (new uqParamSpaceClass<P_V,P_M>(m_env,m_prefix.c_str())),
  m_qoiSpace                  (new uqQoISpaceClass  <Q_V,Q_M>(m_env,m_prefix.c_str())),
  m_userSpacesAreNull         (true),
  m_optionsDesc               (new po::options_description("UQ Propagation Problem")),
  m_option_help               (m_prefix + "help"  ),
  m_option_solver             (m_prefix + "solver"),
  m_solverString              (UQ_PROPAG_PROBLEM_SOLVER_ODV),
  m_paramDensityObj           (paramDensityObj),
  m_paramGeneratorObj         (paramGeneratorObj),
  m_qoiFunctionObj            (qoiFunctionObj),
  m_mcSampler                 (NULL),
  m_solutionProbDensityObj    (NULL),
  m_solutionSampleGeneratorObj(NULL)
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
  m_mcSampler = new uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                         m_prefix.c_str(),
                                                        *m_paramSpace,
                                                        *m_qoiSpace,
                                                         m_paramDensityObj,
                                                         m_paramGeneratorObj,
                                                         m_qoiFunctionObj);

  if (m_env.rank() == 0) std::cout << "Leaving uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::~uqPropagProblemClass()
{
  if (m_solutionSampleGeneratorObj) delete m_solutionSampleGeneratorObj;
  if (m_solutionProbDensityObj    ) delete m_solutionProbDensityObj;
  if (m_mcSampler                 ) delete m_mcSampler;

  if (m_optionsDesc               ) delete m_optionsDesc;
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
  m_mcSampler->calculateDistributions();

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solve(const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj)
{
  m_mcSampler->calculateDistributions(paramGeneratorObj);

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
