/* uq/libs/mcmc/inc/uqValidProblemStage.h
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

#ifndef __UQ_VALID_PROBLEM_STAGE_H__
#define __UQ_VALID_PROBLEM_STAGE_H__

#include <uqParamSpace.h>
#include <uqObservableSpace.h>
#include <uqQoISpace.h>

#include <uqProbDensity.h>        // For substep x.1 and substep x.4 in applications setting a validation problem stage
#include <uqLikelihoodFunction.h> // For substep x.2
#include <uqProposalDensity.h>    // For substep x.3
#include <uqProposalGenerator.h>  // For substep x.3
#include <uqSampleGenerator.h>    // For substep x.4
#include <uqQoIFunction.h>        // For substep x.5

#include <uqDefaultPrior.h>
#include <uqBayesianMarkovChainDC1.h>
#include <uqMonteCarloDC.h>

// _ODV = option default value
#define UQ_VALID_PROBLEM_STAGE_CALIB_PERFORM_ODV           0
#define UQ_VALID_PROBLEM_STAGE_CALIB_INPUT_STAGE_ID_ODV    -1
#define UQ_VALID_PROBLEM_STAGE_CALIB_DISTR_CALCULATOR_ODV  "BMC_DC" // Bayesian Markov Chain Distribution Calculator
#define UQ_VALID_PROBLEM_STAGE_PROPAG_PERFORM_ODV          0
#define UQ_VALID_PROBLEM_STAGE_PROPAG_INPUT_STAGE_ID_ODV   -1
#define UQ_VALID_PROBLEM_STAGE_PROPAG_DISTR_CALCULATOR_ODV "MC_DC" // Monte Carlo Distribution Calculator

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
class uqValidProblemStageClass
{
public:
  uqValidProblemStageClass(const uqEnvironmentClass&                              env,
                           const char*                                            prefix,
                           const uqProbDensity_BaseClass       <P_V,P_M>*         m2lPriorParamDensityObj,  // Set in substep x.1 in applications setting a validation problem stage
                           const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunctionObj, // Set in substep x.2
                           P_M*                                                   proposalCovMatrix,        // Set in substep x.3
                           const uqProposalDensity_BaseClass   <P_V,P_M>*         proposalDensityObj,       // Set in substep x.3
                           const uqProposalGenerator_BaseClass <P_V,P_M>*         proposalGeneratorObj,     // Set in substep x.3 // FIX ME: accomodate code in order to use such object
                           const uqProbDensity_BaseClass       <P_V,P_M>*         propagParamDensityObj,    // Set in substep x.4 // FIX ME: accomodate code in order to use such object
                           const uqSampleGenerator_BaseClass   <P_V,P_M>*         propagParamGeneratorObj,  // Set in substep x.4
                           const uqQoIFunction_BaseClass       <P_V,P_M,Q_V,Q_M>* qoiFunctionObj);          // Set in substep x.5
 ~uqValidProblemStageClass();

        bool                                  isCalibRequested          () const;
        int                                   calibInputStageId         () const;
        void                                  calibrateParamDistribs    ();
        void                                  calibrateParamDistribs    (const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj);
        bool                                  isPropagRequested         () const;
        int                                   propagInputStageId        () const;
        void                                  propagateParamDistribs    ();
        void                                  propagateParamDistribs    (const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj);

  const uqParamSpaceClass         <P_V,P_M>&  paramSpace                () const;
  const uqObservableSpaceClass    <L_V,L_M>&  observableSpace           () const;
  const uqQoISpaceClass           <Q_V,Q_M>&  qoiSpace                  () const;

  const uqProbDensity_BaseClass    <P_V,P_M>& posteriorParamDensityObj  () const;
  const uqSampleGenerator_BaseClass<P_V,P_M>& posteriorParamGeneratorObj() const;

        void                                  print                     (std::ostream& os) const;

private:
  const uqEnvironmentClass&                                  m_env;
        std::string                                          m_prefix;
        uqParamSpaceClass     <P_V,P_M>*                     m_paramSpace;
        uqObservableSpaceClass<L_V,L_M>*                     m_observableSpace;
        uqQoISpaceClass       <Q_V,Q_M>*                     m_qoiSpace;

        po::options_description*                             m_optionsDesc;
        std::string                                          m_option_help;
        std::string                                          m_option_calib_perform;
        std::string                                          m_option_calib_inputStageId;
        std::string                                          m_option_calib_distrCalculator;
        std::string                                          m_option_propag_perform;
        std::string                                          m_option_propag_inputStageId;
        std::string                                          m_option_propag_distrCalculator;

        bool                                                 m_calibPerform;
        int                                                  m_calibInputStageId;
	std::string                                          m_calibDistrCalculator;
        bool                                                 m_propagPerform;
        int                                                  m_propagInputStageId;
	std::string                                          m_propagDistrCalculator;

        bool                                                 m_userPriorDensityIsNull;
  const uqProbDensity_BaseClass           <P_V,P_M>*         m_m2lPriorParamDensityObj;
        uqDefault_M2lPriorRoutine_DataType<P_V,P_M>          m_m2lPriorRoutine_Data;
        P_V*                                                 m_paramPriorMus;
        P_V*                                                 m_paramPriorSigmas;
  const uqLikelihoodFunction_BaseClass    <P_V,P_M,L_V,L_M>* m_m2lLikelihoodFunctionObj;
        P_M*                                                 m_proposalCovMatrix;
  const uqProposalDensity_BaseClass       <P_V,P_M>*         m_proposalDensityObj;
  const uqProposalGenerator_BaseClass     <P_V,P_M>*         m_proposalGeneratorObj;
  const uqProbDensity_BaseClass           <P_V,P_M>*         m_propagParamDensityObj;
  const uqSampleGenerator_BaseClass       <P_V,P_M>*         m_propagParamGeneratorObj;
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
std::ostream& operator<<(std::ostream& os, const uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::uqValidProblemStageClass(
  const uqEnvironmentClass&                              env,
  const char*                                            prefix,
  const uqProbDensity_BaseClass       <P_V,P_M>*         m2lPriorParamDensityObj,
  const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunctionObj,
  P_M*                                                   proposalCovMatrix,
  const uqProposalDensity_BaseClass   <P_V,P_M>*         proposalDensityObj,
  const uqProposalGenerator_BaseClass <P_V,P_M>*         proposalGeneratorObj,
  const uqProbDensity_BaseClass       <P_V,P_M>*         propagParamDensityObj,
  const uqSampleGenerator_BaseClass   <P_V,P_M>*         propagParamGeneratorObj,
  const uqQoIFunction_BaseClass       <P_V,P_M,Q_V,Q_M>* qoiFunctionObj)
  :
  m_env                          (env),
  m_prefix                       (prefix),
  m_paramSpace                   (new uqParamSpaceClass<P_V,P_M>(env,m_prefix.c_str())),
  m_observableSpace              (NULL),
  m_qoiSpace                     (NULL),
  m_optionsDesc                  (new po::options_description("UQ Validation Problem Stage")),
  m_option_help                  (m_prefix + "help"                       ),
  m_option_calib_perform         (m_prefix + "calib_"  + "perform"        ),
  m_option_calib_inputStageId    (m_prefix + "calib_"  + "inputStageId"   ),
  m_option_calib_distrCalculator (m_prefix + "calib_"  + "distrCalculator"),
  m_option_propag_perform        (m_prefix + "propag_" + "perform"        ),
  m_option_propag_inputStageId   (m_prefix + "propag_" + "inputStageId"   ),
  m_option_propag_distrCalculator(m_prefix + "propag_" + "distrCalculator"),
  m_calibPerform                 (UQ_VALID_PROBLEM_STAGE_CALIB_PERFORM_ODV          ),
  m_calibInputStageId            (UQ_VALID_PROBLEM_STAGE_CALIB_INPUT_STAGE_ID_ODV   ),
  m_calibDistrCalculator         (UQ_VALID_PROBLEM_STAGE_CALIB_DISTR_CALCULATOR_ODV ),
  m_propagPerform                (UQ_VALID_PROBLEM_STAGE_PROPAG_PERFORM_ODV         ),
  m_propagInputStageId           (UQ_VALID_PROBLEM_STAGE_PROPAG_INPUT_STAGE_ID_ODV  ),
  m_propagDistrCalculator        (UQ_VALID_PROBLEM_STAGE_PROPAG_DISTR_CALCULATOR_ODV),
  m_userPriorDensityIsNull       (m2lPriorParamDensityObj == NULL),
  m_m2lPriorParamDensityObj      (m2lPriorParamDensityObj),
  m_paramPriorMus                (NULL),
  m_paramPriorSigmas             (NULL),
  m_m2lLikelihoodFunctionObj     (m2lLikelihoodFunctionObj),
  m_proposalCovMatrix            (proposalCovMatrix),
  m_proposalDensityObj           (proposalDensityObj),
  m_proposalGeneratorObj         (proposalGeneratorObj),
  m_propagParamDensityObj        (propagParamDensityObj),
  m_propagParamGeneratorObj      (propagParamGeneratorObj),
  m_qoiFunctionObj               (qoiFunctionObj),
  m_bmcDc                        (NULL),
  m_posteriorParamDensityObj     (NULL),
  m_posteriorParamGeneratorObj   (NULL),
  m_mcDc                         (NULL),
  m_qoiDensityObj                (NULL),
  m_qoiGeneratorObj              (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_calibPerform) {
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
                                                                (m_prefix + "calib_").c_str(),
                                                               *m_paramSpace,
                                                               *m_observableSpace,
                                                               *m_m2lPriorParamDensityObj,
                                                               *m_m2lLikelihoodFunctionObj,
                                                                m_proposalCovMatrix,
                                                                m_proposalDensityObj,
                                                                m_proposalGeneratorObj);
  }

  if (m_propagPerform) {
    m_qoiSpace = new uqQoISpaceClass<Q_V,Q_M>(m_env,
                                              m_prefix.c_str());

    UQ_FATAL_TEST_MACRO(m_qoiFunctionObj == NULL,
                        m_env.rank(),
                        "uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()",
                        "propagation is being requested but 'm_qoiFunctionObj' is null");

    // Instantiate the distribution calculator.
    m_mcDc = new uqMonteCarloDCClass<P_V,P_M,L_V,L_M>(m_env,
                                                      (m_prefix + "propag_").c_str(),
                                                     *m_paramSpace,
                                                     *m_qoiSpace,
                                                      m_propagParamDensityObj,
                                                      m_propagParamGeneratorObj,
                                                     *m_qoiFunctionObj);
  }

  if (m_env.rank() == 0) std::cout << "Leaving uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << std::endl;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::~uqValidProblemStageClass()
{
  if (m_qoiGeneratorObj           ) delete m_qoiGeneratorObj;
  if (m_qoiDensityObj             ) delete m_qoiDensityObj;
  if (m_mcDc                      ) delete m_mcDc;
  if (m_posteriorParamGeneratorObj) delete m_posteriorParamGeneratorObj;
  if (m_posteriorParamDensityObj  ) delete m_posteriorParamDensityObj;
  if (m_bmcDc                     ) delete m_bmcDc;

  if (m_userPriorDensityIsNull) { 
    delete m_m2lPriorParamDensityObj;
    delete m_paramPriorSigmas;
    delete m_paramPriorMus;
  }

  if (m_optionsDesc    ) delete m_optionsDesc;
  if (m_qoiSpace       ) delete m_qoiSpace;
  if (m_observableSpace) delete m_observableSpace;
  if (m_paramSpace     ) delete m_paramSpace;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                "produce help message for validation problem stage")
    (m_option_calib_perform.c_str(),          po::value<bool       >()->default_value(UQ_VALID_PROBLEM_STAGE_CALIB_PERFORM_ODV          ), "calibrate parameters"                             )
    (m_option_calib_inputStageId.c_str(),     po::value<int        >()->default_value(UQ_VALID_PROBLEM_STAGE_CALIB_INPUT_STAGE_ID_ODV   ), "eventual input stage for calibration"             )
    (m_option_calib_distrCalculator.c_str(),  po::value<std::string>()->default_value(UQ_VALID_PROBLEM_STAGE_CALIB_DISTR_CALCULATOR_ODV ), "algorithm for calibration"                        )
    (m_option_propag_perform.c_str(),         po::value<bool       >()->default_value(UQ_VALID_PROBLEM_STAGE_PROPAG_PERFORM_ODV         ), "propagate distributions"                          )
    (m_option_propag_inputStageId.c_str(),    po::value<int        >()->default_value(UQ_VALID_PROBLEM_STAGE_PROPAG_INPUT_STAGE_ID_ODV  ), "eventual input stage for propagation"             )
    (m_option_propag_distrCalculator.c_str(), po::value<std::string>()->default_value(UQ_VALID_PROBLEM_STAGE_PROPAG_DISTR_CALCULATOR_ODV), "algorithm for propagation"                        )
  ;

  return;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
  uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_calib_perform.c_str())) {
    m_calibPerform = m_env.allOptionsMap()[m_option_calib_perform.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_calib_inputStageId.c_str())) {
    m_calibInputStageId = m_env.allOptionsMap()[m_option_calib_inputStageId.c_str()].as<int>();
  }

  if (m_env.allOptionsMap().count(m_option_calib_distrCalculator.c_str())) {
    m_calibDistrCalculator = m_env.allOptionsMap()[m_option_calib_distrCalculator.c_str()].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_propag_perform.c_str())) {
    m_propagPerform = m_env.allOptionsMap()[m_option_propag_perform.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_propag_inputStageId.c_str())) {
    m_propagInputStageId = m_env.allOptionsMap()[m_option_propag_inputStageId.c_str()].as<int>();
  }

  if (m_env.allOptionsMap().count(m_option_propag_distrCalculator.c_str())) {
    m_propagDistrCalculator = m_env.allOptionsMap()[m_option_propag_distrCalculator.c_str()].as<std::string>();
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
bool
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::isCalibRequested() const
{
  return m_calibPerform;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
int
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::calibInputStageId() const
{
  return m_calibInputStageId;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::calibrateParamDistribs()
{
  m_bmcDc->calculateDistributions();

  if (m_posteriorParamGeneratorObj) delete m_posteriorParamGeneratorObj;
  m_posteriorParamGeneratorObj = new uqSampleGenerator_BaseClass<P_V,P_M>(&(m_bmcDc->chain()));

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::calibrateParamDistribs(const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj)
{
  m_bmcDc->calculateDistributions(priorParamDensityObj);

  if (m_posteriorParamGeneratorObj) delete m_posteriorParamGeneratorObj;
  m_posteriorParamGeneratorObj = new uqSampleGenerator_BaseClass<P_V,P_M>(&(m_bmcDc->chain()));

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
bool
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::isPropagRequested() const
{
  return m_propagPerform;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
int
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::propagInputStageId() const
{
  return m_propagInputStageId;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::propagateParamDistribs()
{
  m_mcDc->calculateDistributions();

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::propagateParamDistribs(const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj)
{
  m_mcDc->calculateDistributions(paramGeneratorObj);

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqParamSpaceClass<P_V,P_M>&
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqQoISpaceClass<Q_V,Q_M>&
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::qoiSpace() const
{
  return *m_qoiSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqProbDensity_BaseClass<P_V,P_M>&
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamDensityObj() const
{
  UQ_FATAL_TEST_MACRO(m_posteriorParamDensityObj == NULL,
                      m_env.rank(),
                      "uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamDensityObj()",
                      "posterior param density object is being requested but it has not been created yet");

  return *m_posteriorParamDensityObj;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqSampleGenerator_BaseClass<P_V,P_M>&
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamGeneratorObj() const
{
  UQ_FATAL_TEST_MACRO(m_posteriorParamGeneratorObj == NULL,
                      m_env.rank(),
                      "uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::posteriorParamDensityObj()",
                      "posterior param generator object is being requested but it has not been created yet");

  return *m_posteriorParamGeneratorObj;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqObservableSpaceClass<L_V,L_M>&
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::observableSpace() const
{
  return *m_observableSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_calib_perform          << " = " << m_calibPerform
     << "\n" << m_option_calib_inputStageId     << " = " << m_calibInputStageId
     << "\n" << m_option_calib_distrCalculator  << " = " << m_calibDistrCalculator
     << "\n" << m_option_propag_perform         << " = " << m_propagPerform
     << "\n" << m_option_propag_inputStageId    << " = " << m_propagInputStageId
     << "\n" << m_option_propag_distrCalculator << " = " << m_propagDistrCalculator;
  os << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqValidProblemStageClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_VALID_PROBLEM_STAGE_H__
