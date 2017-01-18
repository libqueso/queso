//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef UQ_VALIDATION_CYCLE_H
#define UQ_VALIDATION_CYCLE_H

#include <queso/VectorFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file ValidationCycle.h
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

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
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

}  // End namespace QUESO

#endif // UQ_VALIDATION_CYCLE_H
