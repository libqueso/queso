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
#include <uqObservableSpace.h>

#include <uqBayesProbDensity.h>
#include <uqProbDensity.h>       // For substep 1 in appls with a calibration problem
#include <uqScalarLhFunction.h>  // For substep 2
#include <uqVectorLhFunction.h>
#include <uqProposalDensity.h>   // For substep 3
#include <uqProposalGenerator.h> // For substep 3

#include <uqDefaultPrior.h>
#include <uqBayesianMarkovChainDC1.h>
#include <uqSampleGenerator.h>

// _ODV = option default value
#define UQ_CALIB_PROBLEM_DISTR_CALCULATOR_ODV "BMC_DC" // Bayesian Markov Chain Distribution Calculator

template <class P_V,class P_M,class L_V,class L_M>
class uqCalibProblemClass
{
public:
  uqCalibProblemClass(const uqEnvironmentClass&                             env,
                      const char*                                           prefix,
                      const uqParamSpaceClass            <P_V,P_M>&         paramSpace,
                      const uqObservableSpaceClass       <L_V,L_M>&         observableSpace,
                      const uqProbDensity_BaseClass      <P_V,P_M>*         m2lPriorParamDensityObj, // Set in substep 1 in appls with a calibration prob.
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
                      const uqScalarLhFunction_BaseClass <P_V,P_M>&         m2lScalarLhFunctionObj,  // Set in substep 2
#else
                      const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,  // Set in substep 2
#endif
                      P_M*                                                  proposalCovMatrix,       // Set in substep 3
                      const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,      // Set in substep 3
                      const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj);   // Set in substep 3 // FIX ME: use such object
  uqCalibProblemClass(const uqEnvironmentClass&                             env,
                      const char*                                           prefix,
                      const uqProbDensity_BaseClass      <P_V,P_M>*         m2lPriorParamDensityObj, // Set in substep 1 in appls with a calibration prob.
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
                      const uqScalarLhFunction_BaseClass <P_V,P_M>&         m2lScalarLhFunctionObj,  // Set in substep 2
#else
                      const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,  // Set in substep 2
#endif
                      P_M*                                                  proposalCovMatrix,       // Set in substep 3
                      const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,      // Set in substep 3
                      const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj);   // Set in substep 3 // FIX ME: use such object
 ~uqCalibProblemClass();

  const uqParamSpaceClass           <P_V,P_M>& paramSpace                () const;
  const uqObservableSpaceClass      <L_V,L_M>& observableSpace           () const;

        void                                   solve                     ();
        void                                   solve                     (const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj);

  const uqBayesProbDensity_BaseClass<P_V,P_M>& posteriorParamDensityObj  () const;
  const uqSampleGenerator_BaseClass <P_V,P_M>& posteriorParamGeneratorObj() const;

        void                                   print                     (std::ostream& os) const;

private:
        void commonConstructor();
        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                                  m_env;
        std::string                                          m_prefix;
  const uqParamSpaceClass     <P_V,P_M>*                     m_paramSpace;
  const uqObservableSpaceClass<L_V,L_M>*                     m_observableSpace;
        bool                                                 m_userSpacesAreNull;

        po::options_description*                             m_optionsDesc;
        std::string                                          m_option_help;
        std::string                                          m_option_distrCalculator;

	std::string                                          m_distrCalculator;

  const uqProbDensity_BaseClass           <P_V,P_M>*         m_m2lPriorParamDensityObj;
        bool                                                 m_userPriorDensityIsNull;
        uqDefault_M2lPriorRoutine_DataType<P_V,P_M>          m_m2lPriorRoutine_Data;
        P_V*                                                 m_paramPriorMus;
        P_V*                                                 m_paramPriorSigmas;
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  const uqScalarLhFunction_BaseClass      <P_V,P_M>&         m_m2lScalarLhFunctionObj;
#else
  const uqVectorLhFunction_BaseClass      <P_V,P_M,L_V,L_M>& m_m2lVectorLhFunctionObj;
#endif
        P_M*                                                 m_proposalCovMatrix;
  const uqProposalDensity_BaseClass       <P_V,P_M>*         m_proposalDensityObj;
  const uqProposalGenerator_BaseClass     <P_V,P_M>*         m_proposalGeneratorObj;

        uqBayesianMarkovChainDCClass      <P_V,P_M,L_V,L_M>* m_bmcDc;
        uqBayesProbDensity_BaseClass      <P_V,P_M>*         m_posteriorParamDensityObj;
        uqSampleGenerator_BaseClass       <P_V,P_M>*         m_posteriorParamGeneratorObj;
};

template<class P_V,class P_M,class L_V,class L_M>
std::ostream& operator<<(std::ostream& os, const uqCalibProblemClass<P_V,P_M,L_V,L_M>& obj);

template <class P_V,class P_M,class L_V,class L_M>
uqCalibProblemClass<P_V,P_M,L_V,L_M>::uqCalibProblemClass(
  const uqEnvironmentClass&                             env,
  const char*                                           prefix,
  const uqParamSpaceClass            <P_V,P_M>&         paramSpace,
  const uqObservableSpaceClass       <L_V,L_M>&         observableSpace,
  const uqProbDensity_BaseClass      <P_V,P_M>*         m2lPriorParamDensityObj,
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  const uqScalarLhFunction_BaseClass <P_V,P_M>&         m2lScalarLhFunctionObj,
#else
  const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,
#endif
  P_M*                                                  proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,
  const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "cal_"),
  m_paramSpace                (&paramSpace),
  m_observableSpace           (&observableSpace),
  m_userSpacesAreNull         (false),
  m_optionsDesc               (new po::options_description("UQ Calibration Problem")),
  m_option_help               (m_prefix + "help"           ),
  m_option_distrCalculator    (m_prefix + "distrCalculator"),
  m_distrCalculator           (UQ_CALIB_PROBLEM_DISTR_CALCULATOR_ODV ),
  m_m2lPriorParamDensityObj   (m2lPriorParamDensityObj),
  m_userPriorDensityIsNull    (m2lPriorParamDensityObj == NULL),
  m_paramPriorMus             (NULL),
  m_paramPriorSigmas          (NULL),
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  m_m2lScalarLhFunctionObj    (m2lScalarLhFunctionObj),
#else
  m_m2lVectorLhFunctionObj    (m2lVectorLhFunctionObj),
#endif
  m_proposalCovMatrix         (proposalCovMatrix),
  m_proposalDensityObj        (proposalDensityObj),
  m_proposalGeneratorObj      (proposalGeneratorObj),
  m_bmcDc                     (NULL),
  m_posteriorParamDensityObj  (NULL),
  m_posteriorParamGeneratorObj(NULL)
{
  commonConstructor();
}

template <class P_V,class P_M,class L_V,class L_M>
uqCalibProblemClass<P_V,P_M,L_V,L_M>::uqCalibProblemClass(
  const uqEnvironmentClass&                            env,
  const char*                                          prefix,
  const uqProbDensity_BaseClass       <P_V,P_M>*       m2lPriorParamDensityObj,
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  const uqScalarLhFunction_BaseClass<P_V,P_M>&         m2lScalarLhFunctionObj,
#else
  const uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,
#endif
  P_M*                                                 proposalCovMatrix,
  const uqProposalDensity_BaseClass   <P_V,P_M>*       proposalDensityObj,
  const uqProposalGenerator_BaseClass <P_V,P_M>*       proposalGeneratorObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "cal_"),
  m_paramSpace                (new uqParamSpaceClass     <P_V,P_M>(m_env,m_prefix.c_str())),
  m_observableSpace           (new uqObservableSpaceClass<L_V,L_M>(m_env,m_prefix.c_str())),
  m_userSpacesAreNull         (true),
  m_optionsDesc               (new po::options_description("UQ Calibration Problem")),
  m_option_help               (m_prefix + "help"           ),
  m_option_distrCalculator    (m_prefix + "distrCalculator"),
  m_distrCalculator           (UQ_CALIB_PROBLEM_DISTR_CALCULATOR_ODV),
  m_m2lPriorParamDensityObj   (m2lPriorParamDensityObj),
  m_userPriorDensityIsNull    (m2lPriorParamDensityObj == NULL),
  m_paramPriorMus             (NULL),
  m_paramPriorSigmas          (NULL),
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  m_m2lScalarLhFunctionObj    (m2lScalarLhFunctionObj),
#else
  m_m2lVectorLhFunctionObj    (m2lVectorLhFunctionObj),
#endif
  m_proposalCovMatrix         (proposalCovMatrix),
  m_proposalDensityObj        (proposalDensityObj),
  m_proposalGeneratorObj      (proposalGeneratorObj),
  m_bmcDc                     (NULL),
  m_posteriorParamDensityObj  (NULL),
  m_posteriorParamGeneratorObj(NULL)
{
  commonConstructor();
}

template <class P_V,class P_M,class L_V,class L_M>
void
uqCalibProblemClass<P_V,P_M,L_V,L_M>::commonConstructor()
{
  if (m_env.rank() == 0) std::cout << "Entering uqCalibProblemClass<P_V,P_M,L_V,L_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqCalibProblemClass<P_V,P_M,L_V,L_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_userPriorDensityIsNull) {
    m_paramPriorMus    = new P_V(m_paramSpace->priorMuValues   ());
    m_paramPriorSigmas = new P_V(m_paramSpace->priorSigmaValues());
    m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
    m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;

    m_m2lPriorParamDensityObj = new uqM2lProbDensity_Class<P_V,P_M>(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                    (void *) &m_m2lPriorRoutine_Data);
  }

#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  m_posteriorParamDensityObj = new uqM2lBayesProbDensity_Class<P_V,P_M>(m_m2lPriorParamDensityObj,
                                                                        &m_m2lScalarLhFunctionObj);
#endif

  // Instantiate the distribution calculator.
  m_bmcDc = new uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>(m_env,
                                                              m_prefix.c_str(),
                                                             *m_paramSpace,
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
                                                             *m_posteriorParamDensityObj,
#else
                                                             *m_observableSpace,
                                                             *m_m2lPriorParamDensityObj,
                                                              m_m2lVectorLhFunctionObj,
#endif
                                                              m_proposalCovMatrix,
                                                              m_proposalDensityObj,
                                                              m_proposalGeneratorObj);

  if (m_env.rank() == 0) std::cout << "Leaving uqCalibProblemClass<P_V,P_M,L_V,L_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << ", m_userSpacesAreNull = " << m_userSpacesAreNull
                                   << std::endl;

  return;
}

template <class P_V,class P_M,class L_V,class L_M>
uqCalibProblemClass<P_V,P_M,L_V,L_M>::~uqCalibProblemClass()
{
  if (m_posteriorParamGeneratorObj) delete m_posteriorParamGeneratorObj;
  if (m_posteriorParamDensityObj  ) delete m_posteriorParamDensityObj;
  if (m_bmcDc                     ) delete m_bmcDc;

  if (m_userPriorDensityIsNull) { 
    delete m_m2lPriorParamDensityObj;
    delete m_paramPriorSigmas;
    delete m_paramPriorMus;
  }

  if (m_optionsDesc) delete m_optionsDesc;

  if (m_userSpacesAreNull) {
    if (m_observableSpace) delete m_observableSpace;
    if (m_paramSpace     ) delete m_paramSpace;
  }
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqCalibProblemClass<P_V,P_M,L_V,L_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                             "produce help message for calibration problem")
    (m_option_distrCalculator.c_str(),  po::value<std::string>()->default_value(UQ_CALIB_PROBLEM_DISTR_CALCULATOR_ODV), "algorithm for calibration"                   )
  ;

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
  uqCalibProblemClass<P_V,P_M,L_V,L_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_distrCalculator.c_str())) {
    m_distrCalculator = m_env.allOptionsMap()[m_option_distrCalculator.c_str()].as<std::string>();
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M>
void
uqCalibProblemClass<P_V,P_M,L_V,L_M>::solve()
{
  m_bmcDc->calculateDistributions();

  if (m_posteriorParamGeneratorObj) delete m_posteriorParamGeneratorObj;
  m_posteriorParamGeneratorObj = new uqSampleGenerator_BaseClass<P_V,P_M>(&(m_bmcDc->chain()));

  return;
}

template <class P_V,class P_M,class L_V,class L_M>
void
uqCalibProblemClass<P_V,P_M,L_V,L_M>::solve(const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj)
{
  m_bmcDc->calculateDistributions(priorParamDensityObj);

  if (m_posteriorParamGeneratorObj) delete m_posteriorParamGeneratorObj;
  m_posteriorParamGeneratorObj = new uqSampleGenerator_BaseClass<P_V,P_M>(&(m_bmcDc->chain()));

  return;
}

template <class P_V,class P_M,class L_V,class L_M>
const uqParamSpaceClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M,L_V,L_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M,class L_V,class L_M>
const uqBayesProbDensity_BaseClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M,L_V,L_M>::posteriorParamDensityObj() const
{
  UQ_FATAL_TEST_MACRO(m_posteriorParamDensityObj == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M,L_V,L_M>::posteriorParamDensityObj()",
                      "posterior param density object is being requested but it has not been created yet");

  return *m_posteriorParamDensityObj;
}

template <class P_V,class P_M,class L_V,class L_M>
const uqSampleGenerator_BaseClass<P_V,P_M>&
uqCalibProblemClass<P_V,P_M,L_V,L_M>::posteriorParamGeneratorObj() const
{
  UQ_FATAL_TEST_MACRO(m_posteriorParamGeneratorObj == NULL,
                      m_env.rank(),
                      "uqCalibProblemClass<P_V,P_M,L_V,L_M>::posteriorParamGeneratorObj()",
                      "posterior param generator object is being requested but it has not been created yet");

  return *m_posteriorParamGeneratorObj;
}

template <class P_V,class P_M,class L_V,class L_M>
const uqObservableSpaceClass<L_V,L_M>&
uqCalibProblemClass<P_V,P_M,L_V,L_M>::observableSpace() const
{
  return *m_observableSpace;
}

template <class P_V,class P_M,class L_V,class L_M>
void
uqCalibProblemClass<P_V,P_M,L_V,L_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_distrCalculator  << " = " << m_distrCalculator;
  os << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M>
std::ostream& operator<<(std::ostream& os, const uqCalibProblemClass<P_V,P_M,L_V,L_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_CALIB_PROBLEM_H__
