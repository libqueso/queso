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

#include <uqCalibProblem.h>
#include <uqPropagProblem.h>

// _ODV = option default value
#define UQ_VALID_PROBLEM_STAGE_CALIB_PERFORM_ODV         0
#define UQ_VALID_PROBLEM_STAGE_CALIB_INPUT_STAGE_ID_ODV  0
#define UQ_VALID_PROBLEM_STAGE_PROPAG_PERFORM_ODV        0
#define UQ_VALID_PROBLEM_STAGE_PROPAG_INPUT_STAGE_ID_ODV 0

template <class P_V,class P_M,class Q_V,class Q_M>
class uqValidProblemStageClass
{
public:
  uqValidProblemStageClass(const uqEnvironmentClass&                             env,
                           const char*                                           prefix,
                           const char*                                           sufix,
                           const uqProbDensity_BaseClass      <P_V,P_M>*         priorParamDensityObj,    // Set in substep x.1 in appls with a validation problem stage
                           const uqProbDensity_BaseClass      <P_V,P_M>*         likelihoodFunctionObj,   // Set in substep x.2
                                 P_M*                                            proposalCovMatrix,       // Set in substep x.3
                           const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,      // Set in substep x.3
                           const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj,    // Set in substep x.3 // FIX ME: accomodate code in order to use such object
                           const uqProbDensity_BaseClass      <P_V,P_M>*         propagParamDensityObj,   // Set in substep x.4 // FIX ME: accomodate code in order to use such object
                           const uqSampleGenerator_BaseClass  <P_V,P_M>*         propagParamGeneratorObj, // Set in substep x.4
                           const uqQoIFunction_BaseClass      <P_V,P_M,Q_V,Q_M>* qoiFunctionObj);         // Set in substep x.5
 ~uqValidProblemStageClass();

  const uqParamSpaceClass   <P_V,P_M>&         paramSpace        () const;
  const uqQoISpaceClass     <Q_V,Q_M>&         qoiSpace          () const;

        bool                                   isCalibRequested  () const;
        unsigned int                           calibInputStageId () const;
        void                                   solveCalibration  ();
        void                                   solveCalibration  (const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj);

        bool                                   isPropagRequested () const;
        unsigned int                           propagInputStageId() const;
        void                                   solvePropagation  ();
        void                                   solvePropagation  (const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj);

  const uqCalibProblemClass <P_V,P_M>&         calibProblem      () const;
  const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& propagProblem     () const;

        void                                   print             (std::ostream& os) const;

private:
        void                                   defineMyOptions   (po::options_description& optionsDesc);
        void                                   getMyOptionValues (po::options_description& optionsDesc);

  const uqEnvironmentClass&                    m_env;
        std::string                            m_prefix;
        uqParamSpaceClass   <P_V,P_M>*         m_paramSpace;
        uqQoISpaceClass     <Q_V,Q_M>*         m_qoiSpace;

        po::options_description*               m_optionsDesc;
        std::string                            m_option_help;
        std::string                            m_option_calib_perform;
        std::string                            m_option_calib_inputStageId;
        std::string                            m_option_propag_perform;
        std::string                            m_option_propag_inputStageId;

        bool                                   m_calibPerform;
        int                                    m_calibInputStageId;
        bool                                   m_propagPerform;
        int                                    m_propagInputStageId;

        uqCalibProblemClass <P_V,P_M>*         m_calibProblem;
        uqPropagProblemClass<P_V,P_M,Q_V,Q_M>* m_propagProblem;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::uqValidProblemStageClass(
  const uqEnvironmentClass&                             env,
  const char*                                           prefix,
  const char*                                           sufix,
  const uqProbDensity_BaseClass      <P_V,P_M>*         priorParamDensityObj,
  const uqProbDensity_BaseClass      <P_V,P_M>*         likelihoodFunctionObj,
        P_M*                                            proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,
  const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj,
  const uqProbDensity_BaseClass      <P_V,P_M>*         propagParamDensityObj,
  const uqSampleGenerator_BaseClass  <P_V,P_M>*         propagParamGeneratorObj,
  const uqQoIFunction_BaseClass      <P_V,P_M,Q_V,Q_M>* qoiFunctionObj)
  :
  m_env                       (env),
  m_prefix                    ((std::string)(prefix) + "s" + (std::string)(sufix)),
  m_paramSpace                (new uqParamSpaceClass<P_V,P_M>(env,m_prefix.c_str())),
  m_qoiSpace                  (new uqQoISpaceClass  <Q_V,Q_M>(env,m_prefix.c_str())),
  m_optionsDesc               (new po::options_description("UQ Validation Problem Stage")),
  m_option_help               (m_prefix + "help"              ),
  m_option_calib_perform      (m_prefix + "performCalibration"),
  m_option_calib_inputStageId (m_prefix + "inputStageIdForCal"),
  m_option_propag_perform     (m_prefix + "performPropagation"),
  m_option_propag_inputStageId(m_prefix + "inputStageIdForPro"),
  m_calibPerform              (UQ_VALID_PROBLEM_STAGE_CALIB_PERFORM_ODV        ),
  m_calibInputStageId         (UQ_VALID_PROBLEM_STAGE_CALIB_INPUT_STAGE_ID_ODV ),
  m_propagPerform             (UQ_VALID_PROBLEM_STAGE_PROPAG_PERFORM_ODV       ),
  m_propagInputStageId        (UQ_VALID_PROBLEM_STAGE_PROPAG_INPUT_STAGE_ID_ODV),
  m_calibProblem              (NULL),
  m_propagProblem             (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = " << m_prefix
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_calibPerform) {
    m_calibProblem = new uqCalibProblemClass<P_V,P_M>(m_env,
                                                      m_prefix.c_str(),
                                                     *m_paramSpace,
                                                      priorParamDensityObj,
                                                     *likelihoodFunctionObj,
                                                      proposalCovMatrix,
                                                      proposalDensityObj,
                                                      proposalGeneratorObj);
  }

  if (m_propagPerform) {
    m_propagProblem = new uqPropagProblemClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                                m_prefix.c_str(),
                                                               *m_paramSpace,
                                                               *m_qoiSpace,
                                                                NULL,
                                                                NULL,
                                                               *qoiFunctionObj);
  }

  if (m_env.rank() == 0) std::cout << "Leaving uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = " << m_prefix
                                   << std::endl;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::~uqValidProblemStageClass()
{
  if (m_propagProblem) delete m_propagProblem;
  if (m_calibProblem ) delete m_calibProblem;
  if (m_optionsDesc  ) delete m_optionsDesc;
  if (m_qoiSpace     ) delete m_qoiSpace;
  if (m_paramSpace   ) delete m_paramSpace;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                    "produce help message for validation problem stage")
    (m_option_calib_perform.c_str(),       po::value<bool>()->default_value(UQ_VALID_PROBLEM_STAGE_CALIB_PERFORM_ODV        ), "calibrate parameters"                             )
    (m_option_calib_inputStageId.c_str(),  po::value<int >()->default_value(UQ_VALID_PROBLEM_STAGE_CALIB_INPUT_STAGE_ID_ODV ), "eventual input stage for calibration"             )
    (m_option_propag_perform.c_str(),      po::value<bool>()->default_value(UQ_VALID_PROBLEM_STAGE_PROPAG_PERFORM_ODV       ), "propagate distributions"                          )
    (m_option_propag_inputStageId.c_str(), po::value<int >()->default_value(UQ_VALID_PROBLEM_STAGE_PROPAG_INPUT_STAGE_ID_ODV), "eventual input stage for propagation"             )
  ;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
  uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::getMyOptionValues(
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

  if (m_env.allOptionsMap().count(m_option_propag_perform.c_str())) {
    m_propagPerform = m_env.allOptionsMap()[m_option_propag_perform.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_propag_inputStageId.c_str())) {
    m_propagInputStageId = m_env.allOptionsMap()[m_option_propag_inputStageId.c_str()].as<int>();
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
bool
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::isCalibRequested() const
{
  return m_calibPerform;
}

template <class P_V,class P_M,class Q_V,class Q_M>
unsigned int
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::calibInputStageId() const
{
  return m_calibInputStageId;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::solveCalibration()
{
  if (m_calibPerform) {
    m_calibProblem->solve();
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::solveCalibration(const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj)
{
  if (m_calibPerform) {
    m_calibProblem->solve(priorParamDensityObj);
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
bool
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::isPropagRequested() const
{
  return m_propagPerform;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqCalibProblemClass<P_V,P_M>&
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::calibProblem() const
{
  return *m_calibProblem;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>&
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::propagProblem() const
{
  return *m_propagProblem;
}

template <class P_V,class P_M,class Q_V,class Q_M>
unsigned int
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::propagInputStageId() const
{
  return m_propagInputStageId;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::solvePropagation()
{
  if (m_propagPerform) {
    m_propagProblem->solve();
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::solvePropagation(const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj)
{
  if (m_propagPerform) {
    m_propagProblem->solve(paramGeneratorObj);
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqParamSpaceClass<P_V,P_M>&
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqQoISpaceClass<Q_V,Q_M>&
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::qoiSpace() const
{
  return *m_qoiSpace;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_calib_perform       << " = " << m_calibPerform
     << "\n" << m_option_calib_inputStageId  << " = " << m_calibInputStageId
     << "\n" << m_option_propag_perform      << " = " << m_propagPerform
     << "\n" << m_option_propag_inputStageId << " = " << m_propagInputStageId;
  os << std::endl;
}

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqValidProblemStageClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_VALID_PROBLEM_STAGE_H__
