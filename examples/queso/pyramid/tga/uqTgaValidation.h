/* uq/examples/queso/pyramid/uqTgaValidation.h
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

#ifndef __UQ_TGA_VALIDATION_H__
#define __UQ_TGA_VALIDATION_H__

#include <uq1D1DFunction.h>
#include <uqTgaLikelihood.h>
#include <uqTgaQoi.h>
#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>

//********************************************************
// Class "uqTgaValidation", instantiated by main()
//********************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
class uqTgaValidationClass : public uqModelValidationClass<P_V,P_M,Q_V,Q_M>
{
public:
  uqTgaValidationClass(const uqBaseEnvironmentClass& env,
                       const char*                   prefix);
 ~uqTgaValidationClass();

  void  run           ();
  void  runGradTest   (const std::string& outerPrefixName,
                       double refA,
                       double refE,
                       bool   treatReferenceDataAsContinuous,
                       double guessA,
                       double guessE);

private:
  void  runCalibrationStage();
  void  runValidationStage ();
  void  runComparisonStage ();

  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_cycle;

  uqAsciiTableClass<P_V,P_M>*                       m_paramsTable;
  const EpetraExt::DistArray<std::string>*          m_paramNames;         // instantiated outside this class!!
  P_V*                                              m_paramMinValues;     // instantiated outside this class!!
  P_V*                                              m_paramMaxValues;     // instantiated outside this class!!
  P_V*                                              m_paramInitialValues; // instantiated outside this class!!
  uqVectorSpaceClass<P_V,P_M>*                      m_paramSpace;
  uqVectorSetClass<P_V,P_M>*                        m_paramDomain;

  uqAsciiTableClass<P_V,P_M>*                       m_qoisTable;
  const EpetraExt::DistArray<std::string>*          m_qoiNames; // instantiated outside this class!!
  uqVectorSpaceClass<Q_V,Q_M>*                      m_qoiSpace;

  double                                            m_predBeta;
  double                                            m_predInitialTemp;
  double                                            m_predCriticalW;
  double                                            m_predCriticalTime;

  uqBaseVectorRVClass<P_V,P_M>*                     m_calPriorRv;
  std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M>* > m_calLikelihoodInfoVector;
  uqBaseScalarFunctionClass<P_V,P_M>*               m_calLikelihoodFunctionObj;
  uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>*              m_calQoiRoutineInfo;

  std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M>* > m_valLikelihoodInfoVector;
  uqBaseScalarFunctionClass<P_V,P_M>*               m_valLikelihoodFunctionObj;
  uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>*              m_valQoiRoutineInfo;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::uqTgaValidationClass(
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
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::constructor()\n"
              << std::endl;
  }

  // Read Ascii file with information on parameters
  m_paramsTable = new uqAsciiTableClass<P_V,P_M> (m_env,
                                                  2,    // # of rows
                                                  3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                                  NULL, // All extra columns are of 'double' type
                                                  "tga/params.tab");

  m_paramNames = &(m_paramsTable->stringColumn(0));
  m_paramMinValues     = new P_V(m_paramsTable->doubleColumn(1));
  m_paramMaxValues     = new P_V(m_paramsTable->doubleColumn(2));
  m_paramInitialValues = new P_V(m_paramsTable->doubleColumn(3));

  m_paramSpace = new uqVectorSpaceClass<P_V,P_M>(m_env,
                                                 "param_", // Extra prefix before the default "space_" prefix
                                                 m_paramsTable->numRows(),
                                                 m_paramNames);

  m_paramDomain = new uqBoxSubsetClass<P_V,P_M>("param_",
                                                *m_paramSpace,
                                                *m_paramMinValues,
                                                *m_paramMaxValues);

  // Read Ascii file with information on qois
  m_qoisTable = new uqAsciiTableClass<P_V,P_M>(m_env,
                                               1,    // # of rows
                                               0,    // # of cols after 'parameter name': none
                                               NULL, // All extra columns are of 'double' type
                                               "tga/qois.tab");

  m_qoiNames = &(m_qoisTable->stringColumn(0));

  m_qoiSpace = new uqVectorSpaceClass<Q_V,Q_M>(m_env,
                                               "qoi_", // Extra prefix before the default "space_" prefix
                                               m_qoisTable->numRows(),
                                               m_qoiNames);

  // Instantiate the validation cycle
  m_cycle = new uqValidationCycleClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                        m_prefix.c_str(), // Use the prefix passed above
                                                        *m_paramSpace,
                                                        *m_qoiSpace);

  m_predBeta         = 250.; // Should /60. --> COMPATIBILITY WITH OLD VERSION
  m_predInitialTemp  = 0.1;
  m_predCriticalW    = 0.;
  m_predCriticalTime = 3.9;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::constructor()\n"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::~uqTgaValidationClass()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::destructor()"
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

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::destructor()"
              << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::run()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::run()"
              << std::endl;
  }

  runCalibrationStage();
  runValidationStage();
  runComparisonStage();

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::run()"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runCalibrationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runCalibrationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_calPriorRv = new uqUniformVectorRVClass<P_V,P_M> ("cal_prior_", // Extra prefix before the default "rv_" prefix
                                                      *m_paramDomain);

  m_calLikelihoodInfoVector.resize(3,NULL);
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_5_K_min.dat");
  m_calLikelihoodInfoVector[1] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_25_K_min.dat");
  m_calLikelihoodInfoVector[2] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_50_K_min.dat");

  m_calLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("cal_like_",
                                                                         *m_paramDomain,
                                                                         uqTgaLikelihoodRoutine<P_V,P_M>,
                                                                         (void *) &m_calLikelihoodInfoVector,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setCalIP(*m_calPriorRv,
                    *m_calLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* calProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().priorRv().pdf().domainVarVector(),
                                                                                                  *m_paramInitialValues);
  m_cycle->calIP().solveWithBayesMarkovChain(*m_paramInitialValues,
                                             calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Deal with forward problem
  m_calQoiRoutineInfo = new uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>(*m_paramSpace,
                                                                m_predInitialTemp,
                                                                m_predBeta,
                                                                m_predCriticalW,
                                                                m_predCriticalTime);

  m_cycle->setCalFP(uqTgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_calQoiRoutineInfo);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->calFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runCalibrationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqTgaValidation::runCalibrationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runValidationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runValidationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_valLikelihoodInfoVector.resize(1,NULL);
  m_valLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_100_K_min.dat");

  m_valLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("val_like_",
                                                                         *m_paramDomain,
                                                                         uqTgaLikelihoodRoutine<P_V,P_M>,
                                                                         (void *) &m_valLikelihoodInfoVector,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setValIP(*m_valLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* valProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().postRv().realizer().imageVarVector(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                                   m_cycle->calIP().postRv().realizer().imageExpVector()); // Use these values as the initial values
  m_cycle->valIP().solveWithBayesMarkovChain(m_cycle->calIP().postRv().realizer().imageExpVector(),
                                             valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Deal with forward problem
  m_valQoiRoutineInfo = new uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>(*m_paramSpace,
                                                                m_predInitialTemp,
                                                                m_predBeta,
                                                                m_predCriticalW,
                                                                m_predCriticalTime);

  m_cycle->setValFP(uqTgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_valQoiRoutineInfo);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->valFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runValidationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqTgaValidation::runValidationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runComparisonStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runComparisonStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  if (m_cycle->calFP().computeSolutionFlag() &&
      m_cycle->valFP().computeSolutionFlag()) {
    Q_V* epsilonVec = m_cycle->calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    Q_V cdfDistancesVec(m_cycle->calFP().qoiRv().imageSet().vectorSpace().zeroVector());
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(m_cycle->valFP().qoiRv().cdf(),
                        m_cycle->calFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runComparisonStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqTgaValidation::runComparisonStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runGradTest(
  const std::string& outerPrefixName,
  double refA,
  double refE,
  bool   treatReferenceDataAsContinuous,
  double guessA,
  double guessE)
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runGradTest()"
              << std::endl;
  }

  // Read parameters for temperature function (which is linear wrt time)
  // Read discrete measurements (which might be replaced by continuous measurements)
  m_calLikelihoodInfoVector.resize(1,NULL);
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_5_K_min.dat");

  // Miscellaneous
  std::ofstream ofs( (outerPrefixName+".m").c_str(), std::ofstream::out | std::ofstream::trunc); // Output file
  P_V paramValues(m_paramSpace->zeroVector());

  // Compute w(refA,refE), using derivative wrt temperature

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(refA,refE) using derivative wrt temperature..."
              << std::endl;
  }

  paramValues[0] = refA;
  paramValues[1] = refE;

  uqTgaComputableWClass<P_V,P_M> tmpW(*m_paramSpace,
                                      *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));
  tmpW.computeUsingTemp(paramValues,
                        822., // maximumTemp
                        NULL, // referenceW
                        NULL);
  std::string tmpPrefixName = outerPrefixName + "withTempW_";
  tmpW.printForMatlab(ofs,tmpPrefixName);


  //////////////////////////////////////////////////////////////////
  // Step 1 of 4: Compute referenceW
  //////////////////////////////////////////////////////////////////

  // Compute w(refA,refE), using derivative wrt time

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(refA,refE) using derivative wrt time..."
              << std::endl;
  }

  tmpW.computeUsingTime(paramValues,
                        false, // computeWAndLambdaGradsAlso
                        NULL,  // referenceW
                        NULL,
                        NULL);

  uqTgaStorageClass<P_V,P_M> refW(tmpW.times(),
                                  tmpW.temps(),
                                  tmpW.ws(),
                                  NULL,  // use all variances = 1.
                                  treatReferenceDataAsContinuous);

  //////////////////////////////////////////////////////////////////
  // Step 2 of 4: Compute gradient(misfit) and Hessian(misfit) wrt A and E, using Lagrangian multipliers
  //////////////////////////////////////////////////////////////////

  // Change the reference data to be used inside the likelihood routine,
  // from discrete (read from file) to continuous (computed and saved in 'refW')
  m_calLikelihoodInfoVector[0]->setReferenceW(refW);

  // Create the objects where the likelihood routine will store its calculations.
  // These objects will then be plotted in Matlab
  uqTgaStorageClass<P_V,P_M> refWForPlot   (refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> wForPlot      (refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> wAForPlot     (refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> wEForPlot     (refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> misfitForPlot (refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> lambdaForPlot (refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> lambdaAForPlot(refW.dataIsContinuousWithTime());
  uqTgaStorageClass<P_V,P_M> lambdaEForPlot(refW.dataIsContinuousWithTime());
  m_calLikelihoodInfoVector[0]->setPlottingVariables(&refWForPlot,
                                                     &wForPlot,
                                                     &wAForPlot,
                                                     &wEForPlot,
                                                     &misfitForPlot,
                                                     &lambdaForPlot,
                                                     &lambdaAForPlot,
                                                     &lambdaEForPlot);

  paramValues[0] = guessA;
  paramValues[1] = guessE;

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute gradient(misfit) and Hessian(misfit) w.r.t. parameters A and E, using adjoint..."
              << std::endl;
  }

  P_V gradWithLM(m_paramSpace->zeroVector());
  P_M* hessianWithLM = m_paramSpace->newMatrix();
  double guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                       (const void *)&m_calLikelihoodInfoVector,
                                                       &gradWithLM,
                                                       hessianWithLM,
                                                       NULL);
  std::cout << "guessMisfit = "     << guessMisfit
            << "\ngradWithLM = "    << gradWithLM
            << "\nhessianWithLM = " << *hessianWithLM
            << std::endl;
  delete hessianWithLM;

  tmpPrefixName = outerPrefixName + "refW_";
  refWForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "w_";
  wForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "wA_";
  wAForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "wE_";
  wEForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "misfit_";
  misfitForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "lambda_";
  lambdaForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "lambdaA_";
  lambdaAForPlot.printForMatlab(ofs,tmpPrefixName);

  tmpPrefixName = outerPrefixName + "lambdaE_";
  lambdaEForPlot.printForMatlab(ofs,tmpPrefixName);

  //////////////////////////////////////////////////////////////////
  // Step 3 of 4: Compute gradient(misfit) wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute gradient(misfit) w.r.t. parameters A and E, using finite differences..."
              << std::endl;
  }

  double deltaA = guessA * 1.e-6;

  paramValues[0] = guessA-deltaA;
  paramValues[1] = guessE;
  double valueAm = 0.;
  tmpW.computeUsingTime(paramValues,
                        false, // computeWAndLambdaGradsAlso
                        &refW,
                        &valueAm,
                        NULL);
  std::cout << "valueAm = " << valueAm << std::endl;

  paramValues[0] = guessA+deltaA;
  paramValues[1] = guessE;
  double valueAp = 0.;
  tmpW.computeUsingTime(paramValues,
                        false, // computeWAndLambdaGradsAlso
                        &refW,
                        &valueAp,
                        NULL);
  std::cout << "valueAp = " << valueAp << std::endl;

  double deltaE = guessE * 1.e-6;

  paramValues[0] = guessA;
  paramValues[1] = guessE-deltaE;
  double valueEm = 0.;
  tmpW.computeUsingTime(paramValues,
                        false, // computeWAndLambdaGradsAlso
                        &refW,
                        &valueEm,
                        NULL);
  std::cout << "valueEm = " << valueEm << std::endl;

  paramValues[0] = guessA;
  paramValues[1] = guessE+deltaE;
  double valueEp = 0.;
  tmpW.computeUsingTime(paramValues,
                        false, // computeWAndLambdaGradsAlso
                        &refW,
                        &valueEp,
                        NULL);
  std::cout << "valueEp = " << valueEp << std::endl;

  P_V gradWithFD(m_paramSpace->zeroVector());
  gradWithFD[0] = (valueAp-valueAm)/2./deltaA;
  gradWithFD[1] = (valueEp-valueEm)/2./deltaE;
  std::cout << "gradWithFD = " << gradWithFD
            << "\ngrad relative error = " << (gradWithFD - gradWithLM).norm2()/gradWithFD.norm2()
            << std::endl;

  //////////////////////////////////////////////////////////////////
  // Step 4 of 4: Compute Hessian(misfit), wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////

  // Compute \grad(misfit), using finite differences
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute Hessian(misfit) w.r.t. parameters A and E, using finite differences..."
              << std::endl;
  }

  paramValues[0] = guessA-deltaA;
  paramValues[1] = guessE;
  P_V grad_Am(m_paramSpace->zeroVector());
  double tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                     (const void *)&m_calLikelihoodInfoVector,
                                                     &grad_Am,
                                                     NULL, // Hessian
                                                     NULL);
  std::cout << "grad_Am = " << grad_Am << std::endl;

  paramValues[0] = guessA+deltaA;
  paramValues[1] = guessE;
  P_V grad_Ap(m_paramSpace->zeroVector());
  tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                              (const void *)&m_calLikelihoodInfoVector,
                                              &grad_Ap,
                                              NULL, // Hessian
                                              NULL);
  std::cout << "grad_Ap = " << grad_Ap << std::endl;

  paramValues[0] = guessA;
  paramValues[1] = guessE-deltaE;
  P_V grad_Em(m_paramSpace->zeroVector());
  tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                              (const void *)&m_calLikelihoodInfoVector,
                                              &grad_Em,
                                              NULL, // Hessian
                                              NULL);
  std::cout << "grad_Em = " << grad_Em << std::endl;

  paramValues[0] = guessA;
  paramValues[1] = guessE+deltaE;
  P_V grad_Ep(m_paramSpace->zeroVector());
  tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                              (const void *)&m_calLikelihoodInfoVector,
                                              &grad_Ep,
                                              NULL, // Hessian
                                              NULL);
  std::cout << "grad_Ep = " << grad_Ep << std::endl;

  P_V column1((.5/deltaA)*(grad_Ap-grad_Am));
  P_V column2((.5/deltaE)*(grad_Ep-grad_Em));
  P_M* hessianWithFD = m_paramSpace->newMatrix();
  (*hessianWithFD)(0,0)=column1[0];
  (*hessianWithFD)(1,0)=column1[1];
  (*hessianWithFD)(0,1)=column2[0];
  (*hessianWithFD)(1,1)=column2[1];
  std::cout << "\nhessianWithFD = " << *hessianWithFD
            << std::endl;
  delete hessianWithFD;

  // Finish
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_VALIDATION_H__

#if 0
  // Compute w(guessA,guessE)

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(guessA,guessE)..."
              << std::endl;
  }

  uqTgaStorageClass<P_V,P_M> weigthedMisfitData(refW.dataIsContinuousWithTime());
  tmpW.computeUsingTime(paramValues,
                        false, // computeWAndLambdaGradsAlso
                        &refW,
                        NULL,
                        &weigthedMisfitData);
  tmpPrefixName = outerPrefixName + "guessW_";
  tmpW.printForMatlab(ofs,tmpPrefixName);

  // AQUI: print misfitData as well 

  // Compute lambda(guessA,guessE)

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute lambda(guessA,guessE)..."
              << std::endl;
  }

  uqTgaLambdaClass<P_V,P_M> guessLambda(*m_paramSpace,
                                        *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));
  guessLambda.compute(paramValues,
                      false, // computeWAndLambdaGradsAlso
                      weigthedMisfitData,
                      tmpW);
  tmpPrefixName = outerPrefixName + "guessLambda_";
  guessLambda.printForMatlab(ofs,tmpPrefixName);
#endif
