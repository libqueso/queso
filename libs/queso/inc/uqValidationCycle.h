/* uq/libs/queso/inc/uqValidationCycle.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_VALIDATION_CYCLE_H__
#define __UQ_VALIDATION_CYCLE_H__

#include <uqStatisticalInverseProblem.h>
#include <uqStatisticalForwardProblem.h>

template <class P_V,class P_M,class Q_V,class Q_M>
class uqValidationCycleClass
{
public:
  uqValidationCycleClass(const uqEnvironmentClass&          env,
                         const char*                        prefix,
                         const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                         const uqVectorSpaceClass<P_V,P_M>& qoiSpace);
 ~uqValidationCycleClass();

  void setCalIP(const uqBaseVectorRVClass <P_V,P_M>& priorRv,
                double (*likelihoodRoutinePtr)(const P_V& paramValues, const void* routineDataPtr),
                const void* likelihoodRoutineDataPtr,
                bool routineComputesMinus2LogOfDensity);

  void setCalFP(void (*qoiRoutinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector),
                const void* qoiRoutineDataPtr);

  const uqStatisticalInverseProblemClass<P_V,P_M>& calIP() const;
        uqStatisticalInverseProblemClass<P_V,P_M>& calIP();

  const uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>& calFP() const;
        uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>& calFP();

private:
  const uqEnvironmentClass&          m_env;
        std::string                  m_prefix;
  const uqVectorSpaceClass<P_V,P_M>& m_paramSpace;
  const uqVectorSpaceClass<Q_V,Q_M>& m_qoiSpace;

  const uqBaseVectorRVClass             <P_V,P_M>* m_calPriorRv;
	uqGenericVectorPdfClass         <P_V,P_M>* m_calLikelihoodFunctionObj;
        uqGenericVectorRVClass          <P_V,P_M>* m_calPostRv;
        uqStatisticalInverseProblemClass<P_V,P_M>* m_calIP;

        uqGenericVectorFunctionClass    <P_V,P_M,Q_V,Q_M>* m_calQoiFunctionObj;
        uqGenericVectorRVClass          <Q_V,Q_M>*         m_calQoiRv;
        uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>* m_calFP;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::uqValidationCycleClass(
  const uqEnvironmentClass&          env,
  const char*                        prefix,
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const uqVectorSpaceClass<P_V,P_M>& qoiSpace)
  :
  m_env       (env),
  m_prefix    ((std::string)(prefix) + ""),
  m_paramSpace(paramSpace),
  m_qoiSpace  (qoiSpace)
{
  if (m_env.rank() == 0) std::cout << "Entering uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;

  if (m_env.rank() == 0) std::cout << "Leaving uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::~uqValidationCycleClass()
{
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::setCalIP(
  const uqBaseVectorRVClass <P_V,P_M>& priorRv,
  double (*likelihoodRoutinePtr)(const P_V& paramValues, const void* routineDataPtr),
  const void* likelihoodRoutineDataPtr,
  bool routineComputesMinus2LogOfDensity)
{
  // Calibration stage: Prior vector rv
  m_calPriorRv = &priorRv;

  // Calibration stage: Likelihood function object (e.g., -2*ln[likelihood])
  m_calLikelihoodFunctionObj = new uqGenericVectorPdfClass<P_V,P_M> ("cal_like_", // Extra prefix before the default "genpd_" prefix
                                                                     m_paramSpace,
                                                                     likelihoodRoutinePtr,
                                                                     likelihoodRoutineDataPtr,
                                                                     routineComputesMinus2LogOfDensity);

  // Calibration stage: Posterior vector rv
  m_calPostRv = new uqGenericVectorRVClass<P_V,P_M> ("cal_post_", // Extra prefix before the default "rv_" prefix
                                                     m_paramSpace);

  // Calibration stage: Inverse problem
  m_calIP = new uqStatisticalInverseProblemClass<P_V,P_M> ("cal_", // Extra prefix before the default "ip_" prefix
                                                           *m_calPriorRv,
                                                           *m_calLikelihoodFunctionObj,
                                                           *m_calPostRv);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqStatisticalInverseProblemClass<P_V,P_M>&
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::calIP() const
{
  return *m_calIP;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqStatisticalInverseProblemClass<P_V,P_M>&
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::calIP()
{
  return *m_calIP;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::setCalFP(
  void (*qoiRoutinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector),
  const void* qoiRoutineDataPtr)
{
  // Calibration stage: Input param vector rv for forward = output posterior vector rv of inverse

  // Calibration stage: Qoi function object
  m_calQoiFunctionObj = new uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M> ("cal_qoi_", // Extra prefix before the default "func_" prefix
                                                                           m_paramSpace,
                                                                           m_qoiSpace,
                                                                           qoiRoutinePtr,
                                                                           qoiRoutineDataPtr);

  // Calibration stage: Qoi vector rv
  m_calQoiRv = new uqGenericVectorRVClass<Q_V,Q_M> ("cal_qoi_", // Extra prefix before the default "rv_" prefix
                                                    m_qoiSpace);

  // Calibration stage: Forward problem
  m_calFP = new uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M> ("cal_",       // Extra prefix before the default "fp_" prefix
                                                                   *m_calPostRv, // forward input = inverse output
                                                                   *m_calQoiFunctionObj,
                                                                   *m_calQoiRv);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>&
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::calFP() const
{
  return *m_calFP;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>&
uqValidationCycleClass<P_V,P_M,Q_V,Q_M>::calFP()
{
  return *m_calFP;
}

#endif // __UQ_VALIDATION_CYCLE_H__
