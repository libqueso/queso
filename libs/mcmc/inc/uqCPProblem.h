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

// _ODV = option default value
#define UQ_CP_PROBLEM_CALIB_PERFORM_ODV  0
#define UQ_CP_PROBLEM_PROPAG_PERFORM_ODV 0

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
class uqCPProblemClass
{
public:
  uqCPProblemClass(const uqEnvironmentClass&                             env,
                   const char*                                           prefix,
                   const uqProbDensity_BaseClass      <P_V,P_M>*         m2lPriorParamDensityObj, // Set in substep 1 in appls with a CP problem
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
                   const uqScalarLhFunction_BaseClass <P_V,P_M>&         m2lScalarLhFunctionObj,  // Set in substep 2
#else
                   const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,  // Set in substep 2
#endif
                   P_M*                                                  proposalCovMatrix,       // Set in substep 3
                   const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,      // Set in substep 3
                   const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj,    // Set in substep 3 // FIX ME: use such object
                   const uqQoIFunction_BaseClass      <P_V,P_M,Q_V,Q_M>& qoiFunctionObj);         // Set in substep 4
 ~uqCPProblemClass();

  const uqParamSpaceClass     <P_V,P_M>&         paramSpace       () const;
  const uqObservableSpaceClass<L_V,L_M>&         observableSpace  () const;
  const uqQoISpaceClass       <Q_V,Q_M>&         qoiSpace         () const;

        void                                     solve            ();

  const uqCalibProblemClass   <P_V,P_M,L_V,L_M>& calibProblem     () const;
  const uqPropagProblemClass  <P_V,P_M,Q_V,Q_M>& propagProblem    () const;

        void                                     print            (std::ostream& os) const;

private:
        void                                     defineMyOptions  (po::options_description& optionsDesc);
        void                                     getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                      m_env;
        std::string                              m_prefix;
        uqParamSpaceClass     <P_V,P_M>*         m_paramSpace;
        uqObservableSpaceClass<L_V,L_M>*         m_observableSpace;
        uqQoISpaceClass       <Q_V,Q_M>*         m_qoiSpace;

        po::options_description*                 m_optionsDesc;
        std::string                              m_option_help;
        std::string                              m_option_calib_perform;
        std::string                              m_option_propag_perform;

        bool                                     m_calibPerform;
        bool                                     m_propagPerform;

        uqCalibProblemClass   <P_V,P_M,L_V,L_M>* m_calibProblem;
        uqPropagProblemClass  <P_V,P_M,Q_V,Q_M>* m_propagProblem;
};

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::uqCPProblemClass(
  const uqEnvironmentClass&                             env,
  const char*                                           prefix,
  const uqProbDensity_BaseClass      <P_V,P_M>*         m2lPriorParamDensityObj,
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  const uqScalarLhFunction_BaseClass <P_V,P_M>&         m2lScalarLhFunctionObj,
#else
  const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,
#endif
  P_M*                                                  proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,
  const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj,
  const uqQoIFunction_BaseClass      <P_V,P_M,Q_V,Q_M>& qoiFunctionObj)
  :
  m_env                  (env),
  m_prefix               ((std::string)(prefix) + "cp_"),
  m_paramSpace           (new uqParamSpaceClass     <P_V,P_M>(env,m_prefix.c_str())),
  m_observableSpace      (new uqObservableSpaceClass<P_V,P_M>(env,m_prefix.c_str())),
  m_qoiSpace             (new uqQoISpaceClass       <P_V,P_M>(env,m_prefix.c_str())),
  m_optionsDesc          (new po::options_description("UQ CP Problem")),
  m_option_help          (m_prefix + "help"            ),
  m_option_calib_perform (m_prefix + "cal_" + "perform"),
  m_option_propag_perform(m_prefix + "pro_" + "perform"),
  m_calibPerform         (UQ_CP_PROBLEM_CALIB_PERFORM_ODV ),
  m_propagPerform        (UQ_CP_PROBLEM_PROPAG_PERFORM_ODV),
  m_calibProblem         (NULL),
  m_propagProblem        (NULL)
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

  if (m_calibPerform) {
    m_calibProblem = new uqCalibProblemClass<P_V,P_M,L_V,L_M>(m_env,
                                                              m_prefix.c_str(),
                                                             *m_paramSpace,
                                                             *m_observableSpace,
                                                              m2lPriorParamDensityObj,
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
                                                              m2lScalarLhFunctionObj,
#else
                                                              m2lVectorLhFunctionObj,
#endif
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
                                                                qoiFunctionObj);
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
    (m_option_help.c_str(),                                                                               "produce help message for CP problem")
    (m_option_calib_perform.c_str(),  po::value<bool>()->default_value(UQ_CP_PROBLEM_CALIB_PERFORM_ODV ), "calibrate parameters"               )
    (m_option_propag_perform.c_str(), po::value<bool>()->default_value(UQ_CP_PROBLEM_PROPAG_PERFORM_ODV), "propagate distributions"            )
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

  if (m_env.allOptionsMap().count(m_option_calib_perform.c_str())) {
    m_calibPerform = m_env.allOptionsMap()[m_option_calib_perform.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_propag_perform.c_str())) {
    m_propagPerform = m_env.allOptionsMap()[m_option_propag_perform.c_str()].as<bool>();
  }


  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqParamSpaceClass<P_V,P_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::paramSpace() const
{
  return *m_paramSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqObservableSpaceClass<L_V,L_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::observableSpace() const
{
  return *m_observableSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqQoISpaceClass<Q_V,Q_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::qoiSpace() const
{
  return *m_qoiSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::solve()
{
  if (m_calibPerform) {
    m_calibProblem->solve();

    if (m_propagPerform) {
      m_propagProblem->solve(m_calibProblem->posteriorParamGeneratorObj());
    }
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqCalibProblemClass<P_V,P_M,L_V,L_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::calibProblem() const
{
  return *m_calibProblem;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>&
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::propagProblem() const
{
  return *m_propagProblem;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_calib_perform  << " = " << m_calibPerform
     << "\n" << m_option_propag_perform << " = " << m_propagPerform;
  os << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqCPProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_CP_PROBLEM_H__
