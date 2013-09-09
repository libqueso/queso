//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_VALIDATION_CYCLE_H__
#define __UQ_VALIDATION_CYCLE_H__

#include <uqStatisticalInverseProblem.h>
#include <uqStatisticalForwardProblem.h>

namespace QUESO {

/*! \file uqValidationCycle.h
 * \brief A templated class for validation cycle of the examples validationCycle and validationCycle2. 
 *
 * \class ValidationCycle
 * \brief A templated class for validation cycle of the examples validationCycle and validationCycle2. 
 * 
 * It has two stages: calibration and validation. First, in the calibration stage, the inverse problem 
 * solution (posterior RV) is the input parameter vector RV for the forward problem. Then, in the 
 * Validation stage, the posterior vector RV from calibration stage (solution of the forward problem of 
 * the calibration stage) is the prior vector RV for an inverse problem. Then, the solution of this
 * inverse problem is once more the input parameter vector RV for the (validation) forward problem.\n
 * 
 * The examples validationCycle and validationCycle2 use the present class to solve the same TGA problem, 
 * and they only differ in implementation styles. */

template <class P_V,class P_M,class Q_V,class Q_M>
class ValidationCycle
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  ValidationCycle(const BaseEnvironment&      env,
                         const char*                        prefix,
                         const VectorSpace<P_V,P_M>& paramSpace,
                         const VectorSpace<P_V,P_M>& qoiSpace);
 
  //! Destructor.
  ~ValidationCycle();
  //@}
  
  //! @name Misc methods
  //@{
  //! Access to the environment variable (m_env).
  const BaseEnvironment& env() const;
  //@}
  
  //! @name Statistical methods
  //@{
  //! Instantiate an inverse problem for the calibration stage.
  void instantiateCalIP(SipOptionsValues*                  optionsValues,
                        const BaseVectorRV      <P_V,P_M>& priorRv,
                        const BaseScalarFunction<P_V,P_M>& likelihoodFunctionObj);
  //double (*likelihoodRoutinePtr)(const P_V& paramValues, const void* routineDataPtr),
  //const void* likelihoodRoutineDataPtr,
  //bool routineComputesMinus2LogOfDensity);
  
  //! Instantiate a forward problem for the calibration stage.  
  void instantiateCalFP(SfpOptionsValues* optionsValues,
                        void (*qoiRoutinePtr)(const P_V&                    domainVector,
                                              const P_V*                    domainDirection,
                                              const void*                   functionDataPtr,
                                                    Q_V&                    imageVector,
                                                    DistArray<P_V*>* gradVectors,
                                                    DistArray<P_M*>* hessianMatrices,
                                                    DistArray<P_V*>* hessianEffects),
                        const void* qoiRoutineDataPtr);

  //! Inverse problem of the calibration stage (const) . 
  /*! It is an instance of class StatisticalInverseProblem<>.*/
  const StatisticalInverseProblem<P_V,P_M>& calIP() const;
  
  //! Inverse problem of the calibration stage (non-const) . 
  /*! It is an instance of class StatisticalInverseProblem<>.*/
  StatisticalInverseProblem<P_V,P_M>& calIP();

  //! Forward problem of the calibration stage (const) . 
  /*! It is an instance of class StatisticalForwardProblem<>.*/
  const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& calFP() const;
  
  //! Forward problem of the calibration stage (non-const) . 
  /*! It is an instance of class StatisticalForwardProblem<>.*/
  StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& calFP();

  //! Instantiate an inverse problem for the validation stage.
  void instantiateValIP(SipOptionsValues*                  optionsValues,
                        const BaseScalarFunction<P_V,P_M>& likelihoodFunctionObj);
  //double (*likelihoodRoutinePtr)(const P_V& paramValues, const void* routineDataPtr),
  //const void* likelihoodRoutineDataPtr,
  //bool routineComputesMinus2LogOfDensity);

  //! Instantiate a forward problem for the validation stage.
  void instantiateValFP(SfpOptionsValues* optionsValues,
                        void (*qoiRoutinePtr)(const P_V&                    domainVector,
                                              const P_V*                    domainDirection,
                                              const void*                   functionDataPtr,
                                                    Q_V&                    imageVector,
                                                    DistArray<P_V*>* gradVectors,
                                                    DistArray<P_M*>* hessianMatrices,
                                                    DistArray<P_V*>* hessianEffects),
                        const void* qoiRoutineDataPtr);

  //! Inverse problem of the validation stage (const) . 
  /*! It is an instance of class StatisticalInverseProblem<>.*/
  const StatisticalInverseProblem<P_V,P_M>& valIP() const;
  
  //! Inverse problem of the validation stage (non-const) . 
  /*! It is an instance of class StatisticalInverseProblem<>.*/
  StatisticalInverseProblem<P_V,P_M>& valIP();

   //! Forward problem of the validation stage (const) . 
  /*! It is an instance of class StatisticalForwardProblem<>.*/
  const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& valFP() const;
  
    //! Forward problem of the validation stage (non-const) . 
  /*! It is an instance of class StatisticalForwardProblem<>.*/
  StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& valFP();
  //@}
private:
  const BaseEnvironment&                            m_env;
        std::string                                        m_prefix;
  const VectorSpace<P_V,P_M>&                       m_paramSpace;
  const VectorSpace<Q_V,Q_M>&                       m_qoiSpace;

  const BaseVectorRV             <P_V,P_M>*         m_calPriorRv;               // instantiated outside this class!!
  const	BaseScalarFunction       <P_V,P_M>*         m_calLikelihoodFunctionObj; // instantiated outside this class!!
        GenericVectorRV          <P_V,P_M>*         m_calPostRv;
        StatisticalInverseProblem<P_V,P_M>*         m_calIP;

        GenericVectorFunction    <P_V,P_M,Q_V,Q_M>* m_calQoiFunctionObj;
        GenericVectorRV          <Q_V,Q_M>*         m_calQoiRv;
        StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>* m_calFP;

  const BaseScalarFunction       <P_V,P_M>*         m_valLikelihoodFunctionObj; // instantiated outside this class!!
        GenericVectorRV          <P_V,P_M>*         m_valPostRv;
        StatisticalInverseProblem<P_V,P_M>*         m_valIP;

        GenericVectorFunction    <P_V,P_M,Q_V,Q_M>* m_valQoiFunctionObj;
        GenericVectorRV          <Q_V,Q_M>*         m_valQoiRv;
        StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>* m_valFP;
};

// Constructor -------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
ValidationCycle<P_V,P_M,Q_V,Q_M>::ValidationCycle(
  const BaseEnvironment&      env,
  const char*                        prefix,
  const VectorSpace<P_V,P_M>& paramSpace,
  const VectorSpace<P_V,P_M>& qoiSpace)
  :
  m_env                     (env),
  m_prefix                  ((std::string)(prefix) + "cycle_"),
  m_paramSpace              (paramSpace),
  m_qoiSpace                (qoiSpace),
  m_calLikelihoodFunctionObj(NULL),
  m_calPostRv               (NULL),
  m_calIP                   (NULL),
  m_calQoiFunctionObj       (NULL),
  m_calQoiRv                (NULL),
  m_calFP                   (NULL),
  m_valLikelihoodFunctionObj(NULL),
  m_valPostRv               (NULL),
  m_valIP                   (NULL),
  m_valQoiFunctionObj       (NULL),
  m_valQoiRv                (NULL),
  m_valFP                   (NULL)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering ValidationCycle<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving ValidationCycle<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  return;
}
// Destructor---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
ValidationCycle<P_V,P_M,Q_V,Q_M>::~ValidationCycle()
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering ValidationCycle::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_valFP)             delete m_valFP;
  if (m_valQoiRv)          delete m_valQoiRv;
  if (m_valQoiFunctionObj) delete m_valQoiFunctionObj;
  if (m_valIP)             delete m_valIP;
  if (m_valPostRv)         delete m_valPostRv;
  if (m_calFP)             delete m_calFP;
  if (m_calQoiRv)          delete m_calQoiRv;
  if (m_calQoiFunctionObj) delete m_calQoiFunctionObj;
  if (m_calIP)             delete m_calIP;
  if (m_calPostRv)         delete m_calPostRv;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving ValidationCycle::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Misc methods-------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const BaseEnvironment&
ValidationCycle<P_V,P_M,Q_V,Q_M>::env() const
{
  return m_env;
}
// Statistical methods------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
ValidationCycle<P_V,P_M,Q_V,Q_M>::instantiateCalIP(
  SipOptionsValues*                  optionsValues,
  const BaseVectorRV      <P_V,P_M>& priorRv,
  const BaseScalarFunction<P_V,P_M>& likelihoodFunctionObj)
  //double (*likelihoodRoutinePtr)(const P_V& paramValues, const void* routineDataPtr),
  //const void* likelihoodRoutineDataPtr,
  //bool routineComputesMinus2LogOfDensity)
{
  // Calibration stage: Prior vector RV
  m_calPriorRv = &priorRv;

  // Calibration stage: Likelihood function object (e.g., ln[likelihood])
  m_calLikelihoodFunctionObj = &likelihoodFunctionObj;

  // Calibration stage: Posterior vector RV
  m_calPostRv = new GenericVectorRV<P_V,P_M> ("cal_post_", // Extra prefix before the default "RV_" prefix
                                                     m_paramSpace);

  // Calibration stage: Inverse problem
  m_calIP = new StatisticalInverseProblem<P_V,P_M> ((m_prefix+"cal_").c_str(), // Extra prefix before the default "ip_" prefix
                                                           optionsValues,
                                                           *m_calPriorRv,
                                                           *m_calLikelihoodFunctionObj,
                                                           *m_calPostRv);

  return;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const StatisticalInverseProblem<P_V,P_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::calIP() const
{
  return *m_calIP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
StatisticalInverseProblem<P_V,P_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::calIP()
{
  return *m_calIP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
ValidationCycle<P_V,P_M,Q_V,Q_M>::instantiateCalFP(
  SfpOptionsValues* optionsValues,
  void (*qoiRoutinePtr)(const P_V&                    domainVector,
                        const P_V*                    domainDirection,
                        const void*                   functionDataPtr,
                              Q_V&                    imageVector,
                              DistArray<P_V*>* gradVectors,
                              DistArray<P_M*>* hessianMatrices,
                              DistArray<P_V*>* hessianEffects),
  const void* qoiRoutineDataPtr)
{
  // Calibration stage: Input parameter vector RV for forward = output posterior vector RV of inverse

  // Calibration stage: QoI function object
  m_calQoiFunctionObj = new GenericVectorFunction<P_V,P_M,Q_V,Q_M> ("cal_qoi_", // Extra prefix before the default "func_" prefix
                                                                           m_paramSpace,
                                                                           m_qoiSpace,
                                                                           qoiRoutinePtr,
                                                                           qoiRoutineDataPtr);

  // Calibration stage: QoI vector RV
  m_calQoiRv = new GenericVectorRV<Q_V,Q_M> ("cal_qoi_", // Extra prefix before the default "RV_" prefix
                                                    m_qoiSpace);

  // Calibration stage: Forward problem
  m_calFP = new StatisticalForwardProblem<P_V,P_M,Q_V,Q_M> ((m_prefix+"cal_").c_str(), // Extra prefix before the default "fp_" prefix
                                                                   optionsValues,
                                                                   *m_calPostRv, // forward input = inverse output
                                                                   *m_calQoiFunctionObj,
                                                                   *m_calQoiRv);

  return;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::calFP() const
{
  return *m_calFP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::calFP()
{
  return *m_calFP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
ValidationCycle<P_V,P_M,Q_V,Q_M>::instantiateValIP(
  SipOptionsValues*                  optionsValues,
  const BaseScalarFunction<P_V,P_M>& likelihoodFunctionObj)
  //double (*likelihoodRoutinePtr)(const P_V& paramValues, const void* routineDataPtr),
  //const void* likelihoodRoutineDataPtr,
  //bool routineComputesMinus2LogOfDensity)
{
  // Validation stage: Prior vector RV = posterior vector RV from calibration stage

  // Validation stage: Likelihood function object (e.g., ln[likelihood])
  m_valLikelihoodFunctionObj = &likelihoodFunctionObj;

  // Validation stage: Posterior vector RV
  m_valPostRv = new GenericVectorRV<P_V,P_M> ("val_post_", // Extra prefix before the default "RV_" prefix
                                                     m_paramSpace);

  // Validation stage: Inverse problem
  m_valIP = new StatisticalInverseProblem<P_V,P_M> ((m_prefix+"val_").c_str(), // Extra prefix before the default "ip_" prefix
                                                           optionsValues,
                                                           *m_calPostRv, // 'validation stage' inverse input = 'calibration stage' inverse output
                                                           *m_valLikelihoodFunctionObj,
                                                           *m_valPostRv);

  return;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const StatisticalInverseProblem<P_V,P_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::valIP() const
{
  return *m_valIP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
StatisticalInverseProblem<P_V,P_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::valIP()
{
  return *m_valIP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
ValidationCycle<P_V,P_M,Q_V,Q_M>::instantiateValFP(
  SfpOptionsValues* optionsValues,
  void (*qoiRoutinePtr)(const P_V&                    domainVector,
                        const P_V*                    domainDirection,
                        const void*                   functionDataPtr,
                              Q_V&                    imageVector,
                              DistArray<P_V*>* gradVectors,
                              DistArray<P_M*>* hessianMatrices,
                              DistArray<P_V*>* hessianEffects),
  const void* qoiRoutineDataPtr)
{
  // Validation stage: Input parameter vector RV for forward = output posterior vector RV of inverse

  // Validation stage: QoI function object
  m_valQoiFunctionObj = new GenericVectorFunction<P_V,P_M,Q_V,Q_M> ("val_qoi_", // Extra prefix before the default "func_" prefix
                                                                           m_paramSpace,
                                                                           m_qoiSpace,
                                                                           qoiRoutinePtr,
                                                                           qoiRoutineDataPtr);

  // Validation stage: QoI vector RV
  m_valQoiRv = new GenericVectorRV<Q_V,Q_M> ("val_qoi_", // Extra prefix before the default "RV_" prefix
                                                    m_qoiSpace);

  // Validation stage: Forward problem
  m_valFP = new StatisticalForwardProblem<P_V,P_M,Q_V,Q_M> ((m_prefix+"val_").c_str(),       // Extra prefix before the default "fp_" prefix
                                                                   optionsValues,
                                                                   *m_valPostRv, // forward input = inverse output
                                                                   *m_valQoiFunctionObj,
                                                                   *m_valQoiRv);

  return;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::valFP() const
{
  return *m_valFP;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>&
ValidationCycle<P_V,P_M,Q_V,Q_M>::valFP()
{
  return *m_valFP;
}

}  // End namespace QUESO

#endif // __UQ_VALIDATION_CYCLE_H__
