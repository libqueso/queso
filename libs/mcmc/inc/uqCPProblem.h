/* uq/libs/mcmc/inc/uqCPProblem.h
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

#ifndef __UQ_CP_PROBLEM_H__
#define __UQ_CP_PROBLEM_H__

#include <uqParamSpace.h>
#include <uqObservableSpace.h>
#include <uqQoISpace.h>

#include <uqCalibProblem.h>
#include <uqPropagProblem.h>

#include <uqDefaultPrior.h>
#include <uqBayesianMarkovChainDC1.h>
#include <uqMonteCarloDC.h>

// _ODV = option default value

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
class uqCPProblemClass
{
public:
  uqCPProblemClass(const uqEnvironmentClass&                              env,
                   const char*                                            prefix,
                   const uqProbDensity_BaseClass       <P_V,P_M>*         m2lPriorParamDensityObj,  // Set in substep 1 in appls with a CP problem
                   const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunctionObj, // Set in substep 2
                   P_M*                                                   proposalCovMatrix,        // Set in substep 3
                   const uqProposalDensity_BaseClass   <P_V,P_M>*         proposalDensityObj,       // Set in substep 3
                   const uqProposalGenerator_BaseClass <P_V,P_M>*         proposalGeneratorObj,     // Set in substep 3 // FIX ME: use such object
                   const uqQoIFunction_BaseClass       <P_V,P_M,Q_V,Q_M>* qoiFunctionObj);          // Set in substep 4
 ~uqCPProblemClass();

  const uqParamSpaceClass         <P_V,P_M>&  paramSpace                () const;
  const uqObservableSpaceClass    <L_V,L_M>&  observableSpace           () const;
  const uqQoISpaceClass           <Q_V,Q_M>&  qoiSpace                  () const;

  const uqProbDensity_BaseClass    <P_V,P_M>& posteriorParamDensityObj  () const;
  const uqSampleGenerator_BaseClass<P_V,P_M>& posteriorParamGeneratorObj() const;

        void                                  solve                     ();

        void                                  print                     (std::ostream& os) const;

private:
  const uqEnvironmentClass&                                  m_env;
        std::string                                          m_prefix;
        uqParamSpaceClass     <P_V,P_M>*                     m_paramSpace;
        uqObservableSpaceClass<L_V,L_M>*                     m_observableSpace;
        uqQoISpaceClass       <Q_V,Q_M>*                     m_qoiSpace;

        po::options_description*                             m_optionsDesc;
        std::string                                          m_option_help;

  const uqLikelihoodFunction_BaseClass    <P_V,P_M,L_V,L_M>* m_m2lLikelihoodFunctionObj;
  const uqQoIFunction_BaseClass           <P_V,P_M,Q_V,Q_M>* m_qoiFunctionObj;

        uqBayesianMarkovChainDCClass      <P_V,P_M,L_V,L_M>* m_bmcDc;
        uqProbDensity_BaseClass           <P_V,P_M>*         m_posteriorParamDensityObj;
        uqSampleGenerator_BaseClass       <P_V,P_M>*         m_posteriorParamGeneratorObj;

        uqMonteCarloDCClass               <P_V,P_M,Q_V,Q_M>* m_mcDc;
        uqProbDensity_BaseClass           <Q_V,Q_M>*         m_qoiDensityObj;
        uqSampleGenerator_BaseClass       <Q_V,Q_M>*         m_qoiGeneratorObj;

        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);
};

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::uqCPProblemClass(
  const uqEnvironmentClass&                              env,
  const char*                                            prefix,
  const uqProbDensity_BaseClass       <P_V,P_M>*         m2lPriorParamDensityObj,
  const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunctionObj,
  P_M*                                                   proposalCovMatrix,
  const uqProposalDensity_BaseClass   <P_V,P_M>*         proposalDensityObj,
  const uqProposalGenerator_BaseClass <P_V,P_M>*         proposalGeneratorObj,
  const uqQoIFunction_BaseClass       <P_V,P_M,Q_V,Q_M>* qoiFunctionObj)
  :
  m_env                          (env),
  m_prefix                       ((std::string)(prefix) + "cp_"),
  m_paramSpace                   (new uqParamSpaceClass<P_V,P_M>(env,m_prefix.c_str())),
  m_observableSpace              (NULL),
  m_qoiSpace                     (NULL),
  m_optionsDesc                  (new po::options_description("UQ CP Problem")),
  m_option_help                  (m_prefix + "help"                    ),
  m_userPriorDensityIsNull       (m2lPriorParamDensityObj == NULL),
  m_m2lPriorParamDensityObj      (m2lPriorParamDensityObj),
  m_paramPriorMus                (NULL),
  m_paramPriorSigmas             (NULL),
  m_m2lLikelihoodFunctionObj     (m2lLikelihoodFunctionObj),
  m_proposalCovMatrix            (proposalCovMatrix),
  m_proposalDensityObj           (proposalDensityObj),
  m_proposalGeneratorObj         (proposalGeneratorObj),
  m_qoiFunctionObj               (qoiFunctionObj),
  m_bmcDc                        (NULL),
  m_posteriorParamDensityObj     (NULL),
  m_posteriorParamGeneratorObj   (NULL),
  m_mcDc                         (NULL),
  m_qoiDensityObj                (NULL),
  m_qoiGeneratorObj              (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = " << m_prefix
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (1) {
    m_observableSpace = new uqObservableSpaceClass<L_V,L_M>(m_env,
                                                            m_prefix.c_str());

    if (m_userPriorDensityIsNull) {
      m_paramPriorMus    = new P_V(m_paramSpace->priorMuValues   ());
      m_paramPriorSigmas = new P_V(m_paramSpace->priorSigmaValues());
      m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
      m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;

      m_m2lPriorParamDensityObj = new uqM2lProbDensity_Class<P_V,P_M>(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                      (void *) &m_m2lPriorRoutine_Data);
    }

    // Instantiate the distribution calculator.
    m_bmcDc = new uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>(m_env,
                                                                (m_prefix + "cal_").c_str(),
                                                               *m_paramSpace,
                                                               *m_observableSpace,
                                                               *m_m2lPriorParamDensityObj,
                                                               *m_m2lLikelihoodFunctionObj,
                                                                m_proposalCovMatrix,
                                                                m_proposalDensityObj,
                                                                m_proposalGeneratorObj);
  }

  if (1) {
    m_qoiSpace = new uqQoISpaceClass<Q_V,Q_M>(m_env,
                                              m_prefix.c_str());

    UQ_FATAL_TEST_MACRO(m_qoiFunctionObj == NULL,
                        m_env.rank(),
                        "uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()",
                        "propagation is being requested but 'm_qoiFunctionObj' is null");

    // Instantiate the distribution calculator.
    m_mcDc = new uqMonteCarloDCClass<P_V,P_M,L_V,L_M>(m_env,
                                                      (m_prefix + "pro_").c_str(),
                                                     *m_paramSpace,
                                                     *m_qoiSpace,
                                                     *m_qoiFunctionObj);
  }

  if (m_env.rank() == 0) std::cout << "Leaving uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = " << m_prefix
                                   << std::endl;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::~uqCPProblemClass()
{
  if (m_optionsDesc    ) delete m_optionsDesc;
  if (m_qoiSpace       ) delete m_qoiSpace;
  if (m_observableSpace) delete m_observableSpace;
  if (m_paramSpace     ) delete m_paramSpace;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(), "produce help message for CP problem")
  ;

  return;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
  uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::solve()
{
  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqParamSpaceClass<P_V,P_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqQoISpaceClass<Q_V,Q_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::qoiSpace() const
{
  return *m_qoiSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqProbDensity_BaseClass<P_V,P_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamDensityObj() const
{
  UQ_FATAL_TEST_MACRO(m_posteriorParamDensityObj == NULL,
                      m_env.rank(),
                      "uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamDensityObj()",
                      "posterior param density object is being requested but it has not been created yet");

  return *m_posteriorParamDensityObj;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqSampleGenerator_BaseClass<P_V,P_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamGeneratorObj() const
{
  UQ_FATAL_TEST_MACRO(m_posteriorParamGeneratorObj == NULL,
                      m_env.rank(),
                      "uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamDensityObj()",
                      "posterior param generator object is being requested but it has not been created yet");

  return *m_posteriorParamGeneratorObj;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqObservableSpaceClass<L_V,L_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::observableSpace() const
{
  return *m_observableSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_CP_PROBLEM_H__
