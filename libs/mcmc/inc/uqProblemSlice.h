/* uq/libs/mcmc/inc/uqProblemSlice.h
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

#ifndef __UQ_PROBLEM_SLICE_H__
#define __UQ_PROBLEM_SLICE_H__

#include <uqParamSpace.h>
#include <uqObservableSpace.h>
//#include <uqQoISpace.h>
#include <uqProbDensity.h>
#include <uqLikelihoodFunction.h>
//#include <uqProposalDensity.h>
//#include <uqProposalGenerator.h>
//#include <uqQoIFunction.h>
#include <uqDefaultPrior.h>
#include <uqDRAM_MarkovChainGenerator.h>

// _ODV = option default value
#define UQ_PROBLEM_SLICE_CALIBRATE_ODV          0
#define UQ_PROBLEM_SLICE_CALIB_INPUT_SLICE_ODV  -1
#define UQ_PROBLEM_SLICE_PROPAGATE_ODV          0
#define UQ_PROBLEM_SLICE_PROPAG_INPUT_SLICE_ODV -1

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
class uqProblemSliceClass
{
public:
  uqProblemSliceClass(const uqEnvironmentClass&                               env,
                      const char*                                             prefix,
                      const uq_ProbDensity_BaseClass       <P_V,P_M>*         m2lPriorProbDensity_Obj,
                      const uq_LikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunction_Obj,
                      P_M*                                                    proposalCovMatrix,
                      void*                                                   userProposalDensity_Obj,
                      void*                                                   userProposalGenerator_Obj,
                      const uq_ProbDensity_BaseClass       <P_V,P_M>*         qoiInputProbDensity_Obj,
                      void*                                                   qoiFunction_Obj);
 ~uqProblemSliceClass();

        bool                             isCalibRequested ();
        void                             calibrate        ();
        bool                             isPropagRequested();
        void                             propagate        ();

  const uqParamSpaceClass     <P_V,P_M>& paramSpace       () const;
  const uqObservableSpaceClass<L_V,L_M>& observableSpace  () const;

        void                             print            (std::ostream& os) const;

private:
  const uqEnvironmentClass&              m_env;
        std::string                      m_prefix;
        uqParamSpaceClass     <P_V,P_M>* m_paramSpace;
        uqObservableSpaceClass<L_V,L_M>* m_observableSpace;
//      uqQoISpaceClass       <Q_V,Q_M>* m_qoiSpace;

        po::options_description*         m_optionsDesc;
        std::string                      m_option_help;
        std::string                      m_option_calibrate;
        std::string                      m_option_calibInputSlice;
        std::string                      m_option_propagate;
        std::string                      m_option_propagInputSlice;

        bool                             m_calibrate;
        int                              m_calibInputSlice;
        bool                             m_propagate;
        int                              m_propagInputSlice;

        bool                                                 m_userPriorDensityIsNull;
  const uq_ProbDensity_BaseClass          <P_V,P_M>*         m_m2lPriorProbDensity_Obj;
        uqDefault_M2lPriorRoutine_DataType<P_V,P_M>          m_m2lPriorRoutine_Data;
        P_V*                                                 m_paramPriorMus;
        P_V*                                                 m_paramPriorSigmas;
  const uq_LikelihoodFunction_BaseClass   <P_V,P_M,L_V,L_M>* m_m2lLikelihoodFunction_Obj;
  const P_M*                                                 m_proposalCovMatrix;
//const uq_ProposalDensity_BaseClass      <P_V,P_M>*         m_proposalDensity_Obj;
//const uq_ProposalGenerator_BaseClass    <P_V,P_M>*         m_proposalGenerator_Obj;
//const uq_QoIFunction_BaseClass          <P_V,P_M,Q_V,Q_M>* m_qoiFunction_Obj;

  uqDRAM_MarkovChainGeneratorClass<P_V,P_M,L_V,L_M>*         m_mcg;

        void defineMyOptions  (po::options_description& optionsDesc);
        void getMyOptionValues(po::options_description& optionsDesc);
};

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::uqProblemSliceClass(
  const uqEnvironmentClass&                               env,
  const char*                                             prefix,
  const uq_ProbDensity_BaseClass       <P_V,P_M>*         m2lPriorProbDensity_Obj,
  const uq_LikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunction_Obj,
  P_M*                                                    proposalCovMatrix,
  void*                                                   userProposalDensity_Obj,
  void*                                                   userProposalGenerator_Obj,
  const uq_ProbDensity_BaseClass       <P_V,P_M>*         qoiInputProbDensity_Obj,
  void*                                                   qoiFunction_Obj)
  :
  m_env                      (env),
  m_prefix                   (prefix),
  m_paramSpace               (new uqParamSpaceClass     <P_V,P_M>(env,m_prefix.c_str())),
  m_observableSpace          (new uqObservableSpaceClass<L_V,L_M>(env,m_prefix.c_str())),
  m_optionsDesc              (new po::options_description("UQ Problem Slice")),
  m_option_help              (m_prefix + "help"            ),
  m_option_calibrate         (m_prefix + "calibrate"       ),
  m_option_calibInputSlice   (m_prefix + "calibInputSlice" ),
  m_option_propagate         (m_prefix + "propagate"       ),
  m_option_propagInputSlice  (m_prefix + "propagInputSlice"),
  m_calibrate                (UQ_PROBLEM_SLICE_CALIBRATE_ODV         ),
  m_calibInputSlice          (UQ_PROBLEM_SLICE_CALIB_INPUT_SLICE_ODV ),
  m_propagate                (UQ_PROBLEM_SLICE_PROPAGATE_ODV         ),
  m_propagInputSlice         (UQ_PROBLEM_SLICE_PROPAG_INPUT_SLICE_ODV),
  m_userPriorDensityIsNull   (false),
  m_m2lPriorProbDensity_Obj  (NULL),
  m_paramPriorMus            (NULL),
  m_paramPriorSigmas         (NULL),
  m_m2lLikelihoodFunction_Obj(NULL),
  m_proposalCovMatrix        (NULL),
  m_mcg                      (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  m_userPriorDensityIsNull  = (m2lPriorProbDensity_Obj == NULL);
  if (m_userPriorDensityIsNull) {
    m_paramPriorMus    = new P_V(m_paramSpace->priorMuValues   ());
    m_paramPriorSigmas = new P_V(m_paramSpace->priorSigmaValues());
    m_m2lPriorRoutine_Data.paramPriorMus    = m_paramPriorMus;
    m_m2lPriorRoutine_Data.paramPriorSigmas = m_paramPriorSigmas;

    m_m2lPriorProbDensity_Obj = new uq_M2lProbDensity_Class<P_V,P_M>(uqDefault_M2lPriorRoutine<P_V,P_M>, // use default prior() routine
                                                                     (void *) &m_m2lPriorRoutine_Data);
  }
  else {
    m_m2lPriorProbDensity_Obj = m2lPriorProbDensity_Obj;
  }

  m_m2lLikelihoodFunction_Obj = m2lLikelihoodFunction_Obj;
  m_proposalCovMatrix         = proposalCovMatrix;

  // Define the Markov chain generator.
  m_mcg = new uqDRAM_MarkovChainGeneratorClass<P_V,P_M,L_V,L_M>(m_env,
                                                                m_prefix.c_str(),
                                                                *m_paramSpace,
                                                                *m_observableSpace,
                                                                *m_m2lPriorProbDensity_Obj,
                                                                *m_m2lLikelihoodFunction_Obj);

  if (m_env.rank() == 0) std::cout << "Leaving uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << std::endl;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::~uqProblemSliceClass()
{
  if (m_mcg                   ) delete m_mcg;
  if (m_userPriorDensityIsNull) { 
    delete m_m2lPriorProbDensity_Obj;
    delete m_paramPriorSigmas;
    delete m_paramPriorMus;
  }
  if (m_optionsDesc           ) delete m_optionsDesc;
  if (m_observableSpace       ) delete m_observableSpace;
  if (m_paramSpace            ) delete m_paramSpace;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                        "produce help message for UQ slice"    )
    (m_option_calibrate.c_str(),        po::value<bool>()->default_value(UQ_PROBLEM_SLICE_CALIBRATE_ODV),          "calibrate"                            )
    (m_option_calibInputSlice.c_str(),  po::value<int >()->default_value(UQ_PROBLEM_SLICE_CALIB_INPUT_SLICE_ODV),  "eventual input slice for calibration)")
    (m_option_propagate.c_str(),        po::value<bool>()->default_value(UQ_PROBLEM_SLICE_PROPAGATE_ODV),          "propagate"                            )
    (m_option_propagInputSlice.c_str(), po::value<int >()->default_value(UQ_PROBLEM_SLICE_PROPAG_INPUT_SLICE_ODV), "eventual input slice for propagation)")
  ;

  return;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
  uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_calibrate.c_str())) {
    m_calibrate = m_env.allOptionsMap()[m_option_calibrate.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_calibInputSlice.c_str())) {
    m_calibInputSlice = m_env.allOptionsMap()[m_option_calibInputSlice.c_str()].as<int>();
  }

  if (m_env.allOptionsMap().count(m_option_propagate.c_str())) {
    m_propagate = m_env.allOptionsMap()[m_option_propagate.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_propagInputSlice.c_str())) {
    m_propagInputSlice = m_env.allOptionsMap()[m_option_propagInputSlice.c_str()].as<int>();
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
bool
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::isCalibRequested()
{
  return m_calibrate;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::calibrate()
{
  m_mcg->generateChains(m_proposalCovMatrix);

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
bool
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::isPropagRequested()
{
  return m_propagate;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::propagate()
{
  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqParamSpaceClass<P_V,P_M>&
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::paramSpace() const
{
  return *m_paramSpace;
}


template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqObservableSpaceClass<L_V,L_M>&
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::observableSpace() const
{
  return *m_observableSpace;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_calibrate        << " = " << m_calibrate
     << "\n" << m_option_calibInputSlice  << " = " << m_calibInputSlice
     << "\n" << m_option_propagate        << " = " << m_propagate
     << "\n" << m_option_propagInputSlice << " = " << m_propagInputSlice;
  os << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROBLEM_SLICE_H__
