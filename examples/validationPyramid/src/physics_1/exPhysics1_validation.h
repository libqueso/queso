/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __EX_PHYSICS_1_VALIDATION_H__
#define __EX_PHYSICS_1_VALIDATION_H__

#include <uq1D1DFunction.h>
#include <exPhysics1_likelihood.h>
#include <exPhysics1_qoi.h>
#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>

//********************************************************
// Class "exPhysics1Validation", instantiated by main()
//********************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
class exPhysics1ValidationClass : public uqModelValidationClass<P_V,P_M,Q_V,Q_M>
{
public:
  exPhysics1ValidationClass(const uqBaseEnvironmentClass& env,
                            const char*                   prefix);
 ~exPhysics1ValidationClass();

  void  run                   ();

private:
  void  runCalibrationStage   ();
  void  runValidationStage    ();
  void  runComparisonStage    ();

  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_cycle;

  uqAsciiTableClass<P_V,P_M>*                            m_paramsTable;
  const EpetraExt::DistArray<std::string>*               m_paramNames;         // instantiated outside this class!!
  P_V*                                                   m_paramMinValues;     // instantiated outside this class!!
  P_V*                                                   m_paramMaxValues;     // instantiated outside this class!!
  P_V*                                                   m_paramInitialValues; // instantiated outside this class!!
  uqVectorSpaceClass<P_V,P_M>*                           m_paramSpace;
  uqVectorSetClass<P_V,P_M>*                             m_paramDomain;

  uqAsciiTableClass<P_V,P_M>*                            m_qoisTable;
  const EpetraExt::DistArray<std::string>*               m_qoiNames; // instantiated outside this class!!
  uqVectorSpaceClass<Q_V,Q_M>*                           m_qoiSpace;

  double                                                 m_predBeta;
  double                                                 m_predInitialTemp;
  double                                                 m_predCriticalW;
  double                                                 m_predCriticalTime;

  uqBaseVectorRVClass<P_V,P_M>*                          m_calPriorRv;
  std::vector<exPhysics1LikelihoodInfoStruct<P_V,P_M>* > m_calLikelihoodInfoVector;
  uqBaseScalarFunctionClass<P_V,P_M>*                    m_calLikelihoodFunctionObj;
  exPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>*              m_calQoiRoutineInfo;

  std::vector<exPhysics1LikelihoodInfoStruct<P_V,P_M>* > m_valLikelihoodInfoVector;
  uqBaseScalarFunctionClass<P_V,P_M>*                    m_valLikelihoodFunctionObj;
  exPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>*              m_valQoiRoutineInfo;
};

template <class P_V,class P_M,class Q_V,class Q_M>
exPhysics1ValidationClass<P_V,P_M,Q_V,Q_M>::exPhysics1ValidationClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  uqModelValidationClass<P_V,P_M,Q_V,Q_M>(env,prefix),
  m_paramsTable             (NULL),
  m_paramNames              (NULL),
  m_paramMinValues          (NULL),
  m_paramMaxValues          (NULL),
  m_paramInitialValues      (NULL),
  m_paramSpace              (NULL),
  m_paramDomain             (NULL),
  m_qoisTable               (NULL),
  m_qoiNames                (NULL),
  m_qoiSpace                (NULL),
  m_calPriorRv              (NULL),
  m_calLikelihoodInfoVector (0),
  m_calLikelihoodFunctionObj(NULL),
  m_calQoiRoutineInfo       (NULL),
  m_valLikelihoodInfoVector (0),
  m_valLikelihoodFunctionObj(NULL),
  m_valQoiRoutineInfo       (NULL)
{
  if (m_env.fullRank() == 0) {
    std::cout << "Entering exPhysics1Validation::constructor()\n"
              << std::endl;
  }

  // Read Ascii file with information on parameters
  m_paramsTable = new uqAsciiTableClass<P_V,P_M> (m_env,
                                                  2,    // # of rows
                                                  3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                                  NULL, // All extra columns are of 'double' type
                                                  "physics1InputData/params.tab");

  m_paramNames = &(m_paramsTable->stringColumn(0));
  m_paramMinValues     = new P_V(m_paramsTable->doubleColumn(1));
  m_paramMaxValues     = new P_V(m_paramsTable->doubleColumn(2));
  m_paramInitialValues = new P_V(m_paramsTable->doubleColumn(3));

  m_paramSpace = new uqVectorSpaceClass<P_V,P_M>(m_env,
                                                 "param_", // Extra prefix before the default "space_" prefix
                                                 m_paramsTable->numRows(),
                                                 NULL);//m_paramNames);

  m_paramDomain = new uqBoxSubsetClass<P_V,P_M>("param_",
                                                *m_paramSpace,
                                                *m_paramMinValues,
                                                *m_paramMaxValues);

  // Read Ascii file with information on qois
  m_qoisTable = new uqAsciiTableClass<P_V,P_M>(m_env,
                                               1,    // # of rows
                                               0,    // # of cols after 'parameter name': none
                                               NULL, // All extra columns are of 'double' type
                                               "physics1InputData/qois.tab");

  m_qoiNames = &(m_qoisTable->stringColumn(0));

  m_qoiSpace = new uqVectorSpaceClass<Q_V,Q_M>(m_env,
                                               "qoi_", // Extra prefix before the default "space_" prefix
                                               m_qoisTable->numRows(),
                                               NULL);//m_qoiNames);

  // Instantiate the validation cycle
  m_cycle = new uqValidationCycleClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                        m_prefix.c_str(), // Use the prefix passed above
                                                        *m_paramSpace,
                                                        *m_qoiSpace);

  m_predBeta         = 250.; // Should /60. --> COMPATIBILITY WITH OLD VERSION
  m_predInitialTemp  = 0.1;
  m_predCriticalW    = 0.;
  m_predCriticalTime = 3.9;

  if (m_env.fullRank() == 0) {
    std::cout << "Leaving exPhysics1Validation::constructor()\n"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
exPhysics1ValidationClass<P_V,P_M,Q_V,Q_M>::~exPhysics1ValidationClass()
{
  if (m_env.fullRank() == 0) {
    std::cout << "Entering exPhysics1Validation::destructor()"
              << std::endl;
  }

  if (m_valQoiRoutineInfo)         delete m_valQoiRoutineInfo;
  if (m_valLikelihoodFunctionObj)  delete m_valLikelihoodFunctionObj;
  for (unsigned int i = 0; i < m_valLikelihoodInfoVector.size(); ++i) {
    if (m_valLikelihoodInfoVector[i]) delete m_valLikelihoodInfoVector[i];
  }
  if (m_calQoiRoutineInfo)         delete m_calQoiRoutineInfo;
  if (m_calLikelihoodFunctionObj)  delete m_calLikelihoodFunctionObj;
  for (unsigned int i = 0; i < m_calLikelihoodInfoVector.size(); ++i) {
    if (m_calLikelihoodInfoVector[i]) delete m_calLikelihoodInfoVector[i];
  }
  if (m_calPriorRv)                delete m_calPriorRv;
  if (m_qoiSpace)                  delete m_qoiSpace;
//if (m_qoiNames)                  delete m_qoiNames; // instantiated outside this class!!
  if (m_qoisTable)                 delete m_qoisTable;
  if (m_paramDomain)               delete m_paramDomain;
  if (m_paramSpace)                delete m_paramSpace;
//if (m_paramInitialValues)        delete m_paramInitialValues; // instantiated outside this class!!
//if (m_paramMaxValues)            delete m_paramMaxValues;     // instantiated outside this class!!
//if (m_paramMinValues)            delete m_paramMinValues;     // instantiated outside this class!!
//if (m_paramNames)                delete m_paramNames;         // instantiated outside this class!!
  if (m_paramsTable)               delete m_paramsTable;

  if (m_env.fullRank() == 0) {
    std::cout << "Leaving exPhysics1Validation::destructor()"
              << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
exPhysics1ValidationClass<P_V,P_M,Q_V,Q_M>::run()
{
  if (m_env.fullRank() == 0) {
    std::cout << "Entering exPhysics1Validation::run()"
              << std::endl;
  }

  runCalibrationStage();
  runValidationStage();
  runComparisonStage();

  if (m_env.fullRank() == 0) {
    std::cout << "Leaving exPhysics1Validation::run()"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
exPhysics1ValidationClass<P_V,P_M,Q_V,Q_M>::runCalibrationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.fullRank() == 0) {
    std::cout << "Entering exPhysics1Validation::runCalibrationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_calPriorRv = new uqUniformVectorRVClass<P_V,P_M> ("cal_prior_", // Extra prefix before the default "rv_" prefix
                                                      *m_paramDomain);

  m_calLikelihoodInfoVector.resize(3,NULL);
  m_calLikelihoodInfoVector[0] = new exPhysics1LikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"physics1InputData/scenario_5_K_min.dat" );
  m_calLikelihoodInfoVector[1] = new exPhysics1LikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"physics1InputData/scenario_25_K_min.dat");
  m_calLikelihoodInfoVector[2] = new exPhysics1LikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"physics1InputData/scenario_50_K_min.dat");

  m_calLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("cal_like_",
                                                                         *m_paramDomain,
                                                                         exPhysics1LikelihoodRoutine<P_V,P_M>,
                                                                         (void *) &m_calLikelihoodInfoVector,
                                                                         true); // the routine computes [ln(function)]

  m_cycle->instantiateCalIP(NULL,
                            *m_calPriorRv,
                            *m_calLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* calProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newProposalMatrix(NULL,
                                                                                                   m_paramInitialValues);
  m_cycle->calIP().solveWithBayesMetropolisHastings(NULL,
						    *m_paramInitialValues,
                                                    calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Deal with forward problem
  m_calQoiRoutineInfo = new exPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>(*m_paramSpace,
                                                                     m_predInitialTemp,
                                                                     m_predBeta,
                                                                     m_predCriticalW,
                                                                     m_predCriticalTime);

  m_cycle->instantiateCalFP(NULL,
                            exPhysics1QoiRoutine<P_V,P_M,Q_V,Q_M>,
                            (void *) m_calQoiRoutineInfo);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->calFP().solveWithMonteCarlo(NULL); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.fullRank() == 0) {
    std::cout << "Leaving exPhysics1Validation::runCalibrationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total exPhysics1Validation::runCalibrationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
exPhysics1ValidationClass<P_V,P_M,Q_V,Q_M>::runValidationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.fullRank() == 0) {
    std::cout << "Entering exPhysics1Validation::runValidationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_valLikelihoodInfoVector.resize(1,NULL);
  m_valLikelihoodInfoVector[0] = new exPhysics1LikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"physics1InputData/scenario_100_K_min.dat");

  m_valLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("val_like_",
                                                                         *m_paramDomain,
                                                                         exPhysics1LikelihoodRoutine<P_V,P_M>,
                                                                         (void *) &m_valLikelihoodInfoVector,
                                                                         true); // the routine computes [ln(function)]

  m_cycle->instantiateValIP(NULL,
                            *m_valLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  const uqSequentialVectorRealizerClass<P_V,P_M>* tmpRealizer = dynamic_cast< const uqSequentialVectorRealizerClass<P_V,P_M>* >(&(m_cycle->calIP().postRv().realizer()));
  P_M* valProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newProposalMatrix(&tmpRealizer->unifiedSampleVarVector(),  // Use 'realizer()' because post. rv was computed with Markov Chain
                                                                                                   &tmpRealizer->unifiedSampleExpVector()); // Use these values as the initial values
  m_cycle->valIP().solveWithBayesMetropolisHastings(NULL,
                                                    tmpRealizer->unifiedSampleExpVector(),
                                                    valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Deal with forward problem
  m_valQoiRoutineInfo = new exPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>(*m_paramSpace,
                                                                     m_predInitialTemp,
                                                                     m_predBeta,
                                                                     m_predCriticalW,
                                                                     m_predCriticalTime);

  m_cycle->instantiateValFP(NULL,
                            exPhysics1QoiRoutine<P_V,P_M,Q_V,Q_M>,
                            (void *) m_valQoiRoutineInfo);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->valFP().solveWithMonteCarlo(NULL); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.fullRank() == 0) {
    std::cout << "Leaving exPhysics1Validation::runValidationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total exPhysics1Validation::runValidationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
exPhysics1ValidationClass<P_V,P_M,Q_V,Q_M>::runComparisonStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.fullRank() == 0) {
    std::cout << "Entering exPhysics1Validation::runComparisonStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  if (m_cycle->calFP().computeSolutionFlag() &&
      m_cycle->valFP().computeSolutionFlag()) {
    Q_V* epsilonVec = m_cycle->calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    Q_V cdfDistancesVec(m_cycle->calFP().qoiRv().imageSet().vectorSpace().zeroVector());
    horizontalDistances(m_cycle->calFP().qoiRv().unifiedCdf(),
                        m_cycle->valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().fullRank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(m_cycle->valFP().qoiRv().unifiedCdf(),
                        m_cycle->calFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().fullRank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(m_cycle->calFP().qoiRv().unifiedCdf(),
                        m_cycle->valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().fullRank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(m_cycle->calFP().qoiRv().unifiedCdf(),
                        m_cycle->valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().fullRank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(m_cycle->calFP().qoiRv().unifiedCdf(),
                        m_cycle->valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().fullRank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(m_cycle->calFP().qoiRv().unifiedCdf(),
                        m_cycle->valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().fullRank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.fullRank() == 0) {
    std::cout << "Leaving exPhysics1Validation::runComparisonStage() at " << ctime(&timevalNow.tv_sec)
              << "Total exPhysics1Validation::runComparisonStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

#endif // __EX_PHYSICS_1_VALIDATION_H__
