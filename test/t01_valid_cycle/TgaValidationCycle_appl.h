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

#ifndef EX_TGA_VALIDATION_CYCLE_APPL_H
#define EX_TGA_VALIDATION_CYCLE_APPL_H

#include <TgaValidationCycle_likelihood.h>
#include <TgaValidationCycle_qoi.h>
#include <queso/ValidationCycle.h>
#include <queso/VectorSubset.h>
#include <queso/AsciiTable.h>
#include <queso/GenericScalarFunction.h>
#include <queso/UniformVectorRV.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <cstdlib>

//Just declaration: actual code is below
template<class P_V,class P_M,class Q_V,class Q_M>
void
uqAppl_LocalComparisonStage(QUESO::ValidationCycle<P_V,P_M,Q_V,Q_M>& cycle);

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqAppl_UnifiedComparisonStage(QUESO::ValidationCycle<P_V,P_M,Q_V,Q_M>& cycle);

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
uqAppl(const QUESO::BaseEnvironment& env)
{
  if (env.fullRank() == 0) {
    std::cout << "Beginning run of 'uqTgaExample' example\n"
              << std::endl;
  }

  const std::string inputDataDir = std::string(std::getenv("srcdir"))
    + "/t01_valid_cycle/inputData";

  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  //******************************************************
  // Task 1 of 5: instantiation of basic classes
  //******************************************************

  // Instantiate the parameter space
  std::vector<std::string> paramNames(2,"");
  paramNames[0] = "A_param";
  paramNames[1] = "E_param";
  QUESO::VectorSpace<P_V,P_M> paramSpace(env,"param_",paramNames.size(),&paramNames);

  // Instantiate the parameter domain
  P_V paramMinValues    (paramSpace.zeroVector());
  paramMinValues[0] = 2.40e+11;
  paramMinValues[1] = 1.80e+05;
  P_V paramMaxValues    (paramSpace.zeroVector());
  paramMaxValues[0] = 2.80e+11;
  paramMaxValues[1] = 2.20e+05;
  QUESO::BoxSubset<P_V,P_M> paramDomain("param_",
                                        paramSpace,
                                        paramMinValues,
                                        paramMaxValues);

  // Instantiate the qoi space
  std::vector<std::string> qoiNames(1,"");
  qoiNames[0] = "TimeFor25PercentOfMass";
  QUESO::VectorSpace<Q_V,Q_M> qoiSpace(env,"qoi_",qoiNames.size(),&qoiNames);

  // Instantiate the validation cycle
  QUESO::ValidationCycle<P_V,P_M,Q_V,Q_M> cycle(env,
                                                "", // No extra prefix
                                                paramSpace,
                                                qoiSpace);

  //********************************************************
  // Task 2 of 5: calibration stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (iRC) {}; // just to remove compiler warning
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'calibration stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Inverse problem: instantiate the prior rv
  QUESO::UniformVectorRV<P_V,P_M> calPriorRv("cal_prior_", // Extra prefix before the default "rv_" prefix
                                             paramDomain);

  // Inverse problem: instantiate the likelihood function object (data + routine)
  likelihoodRoutine_Data<P_V,P_M> calLikelihoodRoutine_Data(env,
                                                                 (inputDataDir+"/scenario_5_K_min.dat").c_str(),
                                                                 (inputDataDir+"/scenario_25_K_min.dat").c_str(),
                                                                 (inputDataDir+"/scenario_50_K_min.dat").c_str());

  QUESO::GenericScalarFunction<P_V,P_M> calLikelihoodFunctionObj("cal_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &calLikelihoodRoutine_Data,
                                                                 true); // the routine computes [ln(function)]

  // Inverse problem: instantiate it (posterior rv is instantiated internally)
  cycle.instantiateCalIP(NULL,
                         calPriorRv,
                         calLikelihoodFunctionObj);

  // Inverse problem: solve it, that is, set 'pdf' and 'realizer' of the posterior rv
  P_V paramInitialValues(paramSpace.zeroVector());
  if (env.numSubEnvironments() == 1) {
    // For regression test purposes
    paramInitialValues[0] = 2.41e+11;
    paramInitialValues[1] = 2.19e+05;
  }
  else {
    calPriorRv.realizer().realization(paramInitialValues);
  }

  P_M* calProposalCovMatrix = cycle.calIP().postRv().imageSet().vectorSpace().newProposalMatrix(NULL,&paramInitialValues);
  cycle.calIP().solveWithBayesMetropolisHastings(NULL,
                                                 paramInitialValues,
                                                 calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Forward problem: instantiate it (parameter rv = posterior rv of inverse problem; qoi rv is instantiated internally)
  double beta_prediction         = 250.;
  double criticalMass_prediction = 0.;
  double criticalTime_prediction = 3.9;

  qoiRoutine_Data<P_V,P_M,Q_V,Q_M> calQoiRoutine_Data;
  calQoiRoutine_Data.m_beta         = beta_prediction;
  calQoiRoutine_Data.m_criticalMass = criticalMass_prediction;
  calQoiRoutine_Data.m_criticalTime = criticalTime_prediction;

  cycle.instantiateCalFP(NULL,
                         qoiRoutine<P_V,P_M,Q_V,Q_M>,
                         (void *) &calQoiRoutine_Data);

  // Forward problem: solve it, that is, set 'realizer' and 'cdf' of the qoi rv
  cycle.calFP().solveWithMonteCarlo(NULL); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Ending 'calibration stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'calibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds\n"
              << std::endl;
  }

  //********************************************************
  // Task 3 of 5: validation stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'validation stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Inverse problem: no need to instantiate the prior rv (= posterior rv of calibration inverse problem)

  // Inverse problem: instantiate the likelihood function object (data + routine)
  likelihoodRoutine_Data<P_V,P_M> valLikelihoodRoutine_Data(env,
                                                                 (inputDataDir+"/scenario_100_K_min.dat").c_str(),
                                                                 NULL,
                                                                 NULL);

  QUESO::GenericScalarFunction<P_V,P_M> valLikelihoodFunctionObj("val_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &valLikelihoodRoutine_Data,
                                                                 true); // the routine computes [ln(function)]

  // Inverse problem: instantiate it (posterior rv is instantiated internally)
  cycle.instantiateValIP(NULL,
                         valLikelihoodFunctionObj);

  // Inverse problem: solve it, that is, set 'pdf' and 'realizer' of the posterior rv
  const QUESO::SequentialVectorRealizer<P_V,P_M>* tmpRealizer = dynamic_cast< const QUESO::SequentialVectorRealizer<P_V,P_M>* >(&(cycle.calIP().postRv().realizer()));
  P_M* valProposalCovMatrix = cycle.calIP().postRv().imageSet().vectorSpace().newProposalMatrix(&tmpRealizer->unifiedSampleVarVector(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                                &tmpRealizer->unifiedSampleExpVector()); // Use these values as the initial values
  cycle.valIP().solveWithBayesMetropolisHastings(NULL,
                                                 tmpRealizer->unifiedSampleExpVector(),
                                                 valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Forward problem: instantiate it (parameter rv = posterior rv of inverse problem; qoi rv is instantiated internally)
  qoiRoutine_Data<P_V,P_M,Q_V,Q_M> valQoiRoutine_Data;
  valQoiRoutine_Data.m_beta         = beta_prediction;
  valQoiRoutine_Data.m_criticalMass = criticalMass_prediction;
  valQoiRoutine_Data.m_criticalTime = criticalTime_prediction;

  cycle.instantiateValFP(NULL,
                         qoiRoutine<P_V,P_M,Q_V,Q_M>,
                         (void *) &valQoiRoutine_Data);

  // Forward problem: solve it, that is, set 'realizer' and 'cdf' of the qoi rv
  cycle.valFP().solveWithMonteCarlo(NULL); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Ending 'validation stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'validation stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds\n"
              << std::endl;
  }

  //********************************************************
  // Task 4 of 5: comparison stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'comparison stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  uqAppl_LocalComparisonStage(cycle);
  if (env.numSubEnvironments() > 1) {
    uqAppl_UnifiedComparisonStage(cycle);
  }

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Ending 'comparison stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'comparison stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds\n"
              << std::endl;
  }

  //******************************************************
  // Task 5 of 5: release memory before leaving routine.
  //******************************************************

  if (env.fullRank() == 0) {
    std::cout << "Finishing run of 'uqTgaExample' example"
              << std::endl;
  }

  return;
}

//********************************************************
// The 'local comparison stage' of the driving routine "uqAppl()"
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void
uqAppl_LocalComparisonStage(QUESO::ValidationCycle<P_V,P_M,Q_V,Q_M>& cycle)
{
  if (cycle.calFP().computeSolutionFlag() &&
      cycle.valFP().computeSolutionFlag()) {
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    Q_V cdfDistancesVec(cycle.calFP().qoiRv().imageSet().vectorSpace().zeroVector());

    // Epsilon = 0.02
    Q_V* epsilonVec = cycle.calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    horizontalDistances(cycle.calFP().qoiRv().unifiedCdf(),
                        cycle.valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().subDisplayFile()) {
      *cycle.env().subDisplayFile() << "For epsilonVec = "    << *epsilonVec
                                   << ", cdfDistancesVec = " << cdfDistancesVec
                                   << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(cycle.valFP().qoiRv().unifiedCdf(),
                        cycle.calFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().subDisplayFile()) {
      *cycle.env().subDisplayFile() << "For epsilonVec = "                             << *epsilonVec
                                   << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                                   << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(cycle.calFP().qoiRv().unifiedCdf(),
                        cycle.valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().subDisplayFile()) {
      *cycle.env().subDisplayFile() << "For epsilonVec = "    << *epsilonVec
                                   << ", cdfDistancesVec = " << cdfDistancesVec
                                   << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(cycle.calFP().qoiRv().unifiedCdf(),
                        cycle.valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().subDisplayFile()) {
      *cycle.env().subDisplayFile() << "For epsilonVec = "    << *epsilonVec
                                   << ", cdfDistancesVec = " << cdfDistancesVec
                                   << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(cycle.calFP().qoiRv().unifiedCdf(),
                        cycle.valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().subDisplayFile()) {
      *cycle.env().subDisplayFile() << "For epsilonVec = "    << *epsilonVec
                                   << ", cdfDistancesVec = " << cdfDistancesVec
                                   << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(cycle.calFP().qoiRv().unifiedCdf(),
                        cycle.valFP().qoiRv().unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().subDisplayFile()) {
      *cycle.env().subDisplayFile() << "For epsilonVec = "    << *epsilonVec
                                   << ", cdfDistancesVec = " << cdfDistancesVec
                                   << std::endl;
    }

    delete epsilonVec;
#endif
  }

  return;
}

//********************************************************
// The 'unified comparison stage' of the driving routine "uqAppl()"
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void
uqAppl_UnifiedComparisonStage(QUESO::ValidationCycle<P_V,P_M,Q_V,Q_M>& cycle)
{
  if (cycle.calFP().computeSolutionFlag() &&
      cycle.valFP().computeSolutionFlag()) {
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    Q_V cdfDistancesVec(cycle.calFP().qoiRv().imageSet().vectorSpace().zeroVector());

    // Epsilon = 0.02
    Q_V* epsilonVec = cycle.calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    horizontalDistances(cycle.calFP().qoiRv_unifiedCdf(),
                        cycle.valFP().qoiRv_unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().fullRank() == 0) {
      std::cout << "For epsilonVec = "           << *epsilonVec
                << ", unifiedCdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(cycle.valFP().qoiRv_unifiedCdf(),
                        cycle.calFP().qoiRv_unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().fullRank() == 0) {
      std::cout << "For epsilonVec = "                                    << *epsilonVec
                << ", unifiedCdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(cycle.calFP().qoiRv_unifiedCdf(),
                        cycle.valFP().qoiRv_unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().fullRank() == 0) {
      std::cout << "For epsilonVec = "           << *epsilonVec
                << ", unifiedCdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(cycle.calFP().qoiRv_unifiedCdf(),
                        cycle.valFP().qoiRv_unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().fullRank() == 0) {
      std::cout << "For epsilonVec = "           << *epsilonVec
                << ", unifiedCdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(cycle.calFP().qoiRv_unifiedCdf(),
                        cycle.valFP().qoiRv_unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().fullRank() == 0) {
      std::cout << "For epsilonVec = "           << *epsilonVec
                << ", unifiedCdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(cycle.calFP().qoiRv_unifiedCdf(),
                        cycle.valFP().qoiRv_unifiedCdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().fullRank() == 0) {
      std::cout << "For epsilonVec = "           << *epsilonVec
                << ", unifiedCdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
#endif
  }

  return;
}
#endif // EX_TGA_VALIDATION_CYCLE_APPL_H
