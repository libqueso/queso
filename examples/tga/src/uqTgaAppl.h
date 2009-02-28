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

#ifndef __UQ_TGA_APPL_H__
#define __UQ_TGA_APPL_H__

#include <uqTgaLikelihood.h>
#include <uqTgaQoi.h>
#include <uqValidationCycle.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

//Just declaration: actual code is below
template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAppl_ComparisonStage(uqValidationCycleClass<P_V,P_M,Q_V,Q_M>& cycle);

//********************************************************
// The driving routine "uqAppl()": called by main()
// There are 5 main stages:
// 1) initilization
// 2) the 'calibration stage'
// 3) the 'validation stage'
// 4) the 'comparison stage'
// 5) release memory
// Stages 2, 3 and 4 constitute the actual validation cycle.
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAppl(const uqBaseEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqTgaExample' example\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  //******************************************************
  // Task 1 of 5: instantiation of basic classes
  //******************************************************

  // Read Ascii file with information on parameters
  uqAsciiTableClass<P_V,P_M> paramsTable(env,
                                         2,    // # of rows
                                         3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                         NULL, // All extra columns are of 'double' type
                                         "inputData/params.tab");

  const EpetraExt::DistArray<std::string>& paramNames = paramsTable.stringColumn(0);
  P_V                                      paramMinValues    (paramsTable.doubleColumn(1));
  P_V                                      paramMaxValues    (paramsTable.doubleColumn(2));
  P_V                                      paramInitialValues(paramsTable.doubleColumn(3));

  uqVectorSpaceClass<P_V,P_M> paramSpace(env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         paramsTable.numRows(),
                                         &paramNames);

  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",
                                        paramSpace,
                                        paramMinValues,
                                        paramMaxValues);

  // Read Ascii file with information on qois
  uqAsciiTableClass<P_V,P_M> qoisTable(env,
                                       1,    // # of rows
                                       0,    // # of cols after 'qoi name': none
                                       NULL, // All extra columns are of 'double' type
                                       "inputData/qois.tab");

  const EpetraExt::DistArray<std::string>& qoiNames = qoisTable.stringColumn(0);

  double beta_prediction         = 250.;
  double criticalMass_prediction = 0.;
  double criticalTime_prediction = 3.9;

  uqVectorSpaceClass<Q_V,Q_M> qoiSpace(env,
                                       "qoi_", // Extra prefix before the default "space_" prefix
                                       qoisTable.numRows(),
                                       &qoiNames);

  // Instantiate the validation cycle
  uqValidationCycleClass<P_V,P_M,Q_V,Q_M> cycle(env,
                                                "", // No extra prefix
                                                paramSpace,
                                                qoiSpace);

  //********************************************************
  // Task 2 of 5: calibration stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning 'calibration stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Inverse problem: instantiate the prior rv
  uqUniformVectorRVClass<P_V,P_M> calPriorRv("cal_prior_", // Extra prefix before the default "rv_" prefix
                                             paramDomain);

  // Inverse problem: instantiate the likelihood function object (data + routine)
  likelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(env,
                                                                 "inputData/scenario_5_K_min.dat",
                                                                 "inputData/scenario_25_K_min.dat",
                                                                 "inputData/scenario_50_K_min.dat");

  uqGenericScalarFunctionClass<P_V,P_M> calLikelihoodFunctionObj("cal_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &calLikelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(function)]

  // Inverse problem: instantiate it (posterior rv is instantiated internally)
  cycle.instantiateCalIP(calPriorRv,
                         calLikelihoodFunctionObj);

  // Inverse problem: solve it, that is, set 'pdf' and 'realizer' of the posterior rv
  P_M* calProposalCovMatrix = cycle.calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(cycle.calIP().priorRv().pdf().domainVarVector(),
                                                                                                paramInitialValues);
  cycle.calIP().solveWithBayesMarkovChain(paramInitialValues,
                                          calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Forward problem: instantiate it (parameter rv = posterior rv of inverse problem; qoi rv is instantiated internally)
  qoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> calQoiRoutine_Data;
  calQoiRoutine_Data.m_beta         = beta_prediction;
  calQoiRoutine_Data.m_criticalMass = criticalMass_prediction;
  calQoiRoutine_Data.m_criticalTime = criticalTime_prediction;

  cycle.instantiateCalFP(qoiRoutine<P_V,P_M,Q_V,Q_M>,
                         (void *) &calQoiRoutine_Data);

  // Forward problem: solve it, that is, set 'realizer' and 'cdf' of the qoi rv
  cycle.calFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'calibration stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'calibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds\n"
              << std::endl;
  }

  //********************************************************
  // Task 3 of 5: validation stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning 'validation stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Inverse problem: no need to instantiate the prior rv (= posterior rv of calibration inverse problem)

  // Inverse problem: instantiate the likelihood function object (data + routine)
  likelihoodRoutine_DataClass<P_V,P_M> valLikelihoodRoutine_Data(env,
                                                                 "inputData/scenario_100_K_min.dat",
                                                                 NULL,
                                                                 NULL);

  uqGenericScalarFunctionClass<P_V,P_M> valLikelihoodFunctionObj("val_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &valLikelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(function)]

  // Inverse problem: instantiate it (posterior rv is instantiated internally)
  cycle.instantiateValIP(valLikelihoodFunctionObj);

  // Inverse problem: solve it, that is, set 'pdf' and 'realizer' of the posterior rv
  P_M* valProposalCovMatrix = cycle.calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(cycle.calIP().postRv().realizer().imageVarVector(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                                cycle.calIP().postRv().realizer().imageExpVector()); // Use these values as the initial values
  cycle.valIP().solveWithBayesMarkovChain(cycle.calIP().postRv().realizer().imageExpVector(),
                                          valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Forward problem: instantiate it (parameter rv = posterior rv of inverse problem; qoi rv is instantiated internally)
  qoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> valQoiRoutine_Data;
  valQoiRoutine_Data.m_beta         = beta_prediction;
  valQoiRoutine_Data.m_criticalMass = criticalMass_prediction;
  valQoiRoutine_Data.m_criticalTime = criticalTime_prediction;

  cycle.instantiateValFP(qoiRoutine<P_V,P_M,Q_V,Q_M>,
                         (void *) &valQoiRoutine_Data);

  // Forward problem: solve it, that is, set 'realizer' and 'cdf' of the qoi rv
  cycle.valFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'validation stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'validation stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds\n"
              << std::endl;
  }

  //********************************************************
  // Task 4 of 5: comparison stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning 'comparison stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  uqAppl_ComparisonStage(cycle);

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'comparison stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'comparison stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds\n"
              << std::endl;
  }

  //******************************************************
  // Task 5 of 5: release memory before leaving routine.
  //******************************************************

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqTgaExample' example"
              << std::endl;
  }

  return;
}

//********************************************************
// The 'comparison stage' of the driving routine "uqAppl()"
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAppl_ComparisonStage(uqValidationCycleClass<P_V,P_M,Q_V,Q_M>& cycle)
{
  if (cycle.calFP().computeSolutionFlag() &&
      cycle.valFP().computeSolutionFlag()) {
    Q_V cdfDistancesVec(cycle.calFP().qoiRv().imageSet().vectorSpace().zeroVector());

    // Epsilon = 0.02
    Q_V* epsilonVec = cycle.calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(cycle.valFP().qoiRv().cdf(),
                        cycle.calFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  return;
}
#endif // __UQ_TGA_APPL_H__
