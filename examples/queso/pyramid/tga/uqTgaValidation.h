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
  void  runGradTest   ();

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
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runGradTest()
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
  std::ofstream ofs("nada1.m", std::ofstream::out | std::ofstream::trunc); // Output file
  P_V params(m_paramSpace->zeroVector());

  // Compute w(A_reference,E_reference), using derivative wrt temperature

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(A_reference,E_reference) using derivative wrt temperature..."
              << std::endl;
  }

  params[0] = 2.6090e+11; // Reference _A
  params[1] = 1.9910e+05; // Reference _E

  uqTgaComputableWClass<P_V,P_M> tmpW(*m_paramSpace,
                                      *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));
  tmpW.computeUsingTemp(params,
                        822., // maximumTemp
                        NULL, // referenceW
                        NULL);

  unsigned int tmpSize = tmpW.times().size();
  ofs << "tempTemps = zeros(" << tmpSize
      << ",1);"
      << "\ntempWs = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\ntempTemps(" << i+1 << ",1) = " << tmpW.temps()[i] << ";"
        << "\ntempWs("    << i+1 << ",1) = " << tmpW.ws()   [i] << ";";
  }

  // Compute w(A_reference,E_reference), using derivative wrt time

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(A_reference,E_reference) using derivative wrt time..."
              << std::endl;
  }

  tmpW.computeUsingTime(params,
                        false, // computeWAndLambdaGradsAlso
                        NULL,  // referenceW
                        NULL,
                        NULL);

  uqTgaStorageWClass<P_V,P_M> referenceW(tmpW.times(),
                                         tmpW.temps(),
                                         tmpW.ws(),
                                         NULL,  // use all variances = 1.
                                         true); // treat data as continuous with time

  tmpSize = referenceW.times().size();
  ofs << "referenceTemps = zeros(" << tmpSize
      << ",1);"
      << "\nreferenceWs = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nreferenceTemps(" << i+1 << ",1) = " << referenceW.temps() [i] << ";"
        << "\nreferenceWs("    << i+1 << ",1) = " << referenceW.values()[i] << ";";
  }

  // Change measured data, from discrete (read from file) to continuous (computed with 'uqTgaComputableWClass' object)
  m_calLikelihoodInfoVector[0]->setReferenceW(referenceW);

  // Compute w(A_guess,E_guess)

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(A_guess,E_guess)..."
              << std::endl;
  }

  params[0] = 2.6000e+11; // Guess _A
  params[1] = 2.0000e+05; // Guess _E
  uqTgaStorageWClass<P_V,P_M> weigthedMisfitData(referenceW.continuous());
  tmpW.computeUsingTime(params,
                        false, // computeWAndLambdaGradsAlso
                        &referenceW,
                        NULL,
                        &weigthedMisfitData);
 
  tmpSize = tmpW.times().size();
  ofs << "\nguessWTemps = zeros(" << tmpSize
      << ",1);"
      << "\nguessWs = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nguessWTemps(" << i+1 << ",1) = " << tmpW.temps()[i] << ";"
        << "\nguessWs("     << i+1 << ",1) = " << tmpW.ws()   [i] << ";";
  }

  // Compute lambda(A_guess,E_guess)

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute lambda(A_guess,E_guess)..."
              << std::endl;
  }

  uqTgaLambdaClass<P_V,P_M> guessLambda(*m_paramSpace,
                                        *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));
  guessLambda.compute(params,
                      false, // computeWAndLambdaGradsAlso
                      weigthedMisfitData,
                      tmpW);

  tmpSize = guessLambda.times().size();
  ofs << "\nguessLambdaTemps = zeros(" << tmpSize
      << ",1);"
      << "\nguessLambdas = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nguessLambdaTemps(" << i+1 << ",1) = " << guessLambda.temps()[i]   << ";"
        << "\nguessLambdas("     << i+1 << ",1) = " << guessLambda.lambdas()[i] << ";";
  }

  //////////////////////////////////////////////////////////////////
  // Compute Gradient(misfit) wrt A and E, using Lagrangian multipliers
  //////////////////////////////////////////////////////////////////

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute misfit gradient w.r.t. parameters A and E, using adjoint..."
              << std::endl;
  }

  P_V misfitGrad(m_paramSpace->zeroVector());
  double valueCenter = uqTgaLikelihoodRoutine<P_V,P_M>(params,
                                                       (const void *)&m_calLikelihoodInfoVector,
                                                       &misfitGrad,
                                                       NULL, // Hessian
                                                       NULL);
  std::cout << "valueCenter = "  << valueCenter
            << "\nmisfitGrad = " << misfitGrad
            << std::endl;

  //////////////////////////////////////////////////////////////////
  // Compute Gradient(misfit) wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////
  double deltaA = params[0] * 1.e-6;

  params[0] = 2.6000e+11-deltaA; // A
  params[1] = 2.0000e+05;        // E
  double valueAm = 0.;
  tmpW.computeUsingTime(params,
                        false, // computeWAndLambdaGradsAlso
                        &referenceW,
                        &valueAm,
                        NULL);
  std::cout << "valueAm = " << valueAm << std::endl;

  params[0] = 2.6000e+11+deltaA; // A
  params[1] = 2.0000e+05;        // E
  double valueAp = 0.;
  tmpW.computeUsingTime(params,
                        false, // computeWAndLambdaGradsAlso
                        &referenceW,
                        &valueAp,
                        NULL);
  std::cout << "valueAp = " << valueAp << std::endl;

  double deltaE = params[1] * 1.e-6;

  params[0] = 2.6000e+11;        // A
  params[1] = 2.0000e+05-deltaE; // E
  double valueEm = 0.;
  tmpW.computeUsingTime(params,
                        false, // computeWAndLambdaGradsAlso
                        &referenceW,
                        &valueEm,
                        NULL);
  std::cout << "valueEm = " << valueEm << std::endl;

  params[0] = 2.6000e+11;        // A
  params[1] = 2.0000e+05+deltaE; // E
  double valueEp = 0.;
  tmpW.computeUsingTime(params,
                        false, // computeWAndLambdaGradsAlso
                        &referenceW,
                        &valueEp,
                        NULL);
  std::cout << "valueEp = " << valueEp << std::endl;

  std::cout << "So, grad (finite diff) = " << (valueAp-valueAm)/2./deltaA << ", " << (valueEp-valueEm)/2./deltaE << std::endl;

  //////////////////////////////////////////////////////////////////
  // Compute Hessian(misfit) wrt A and E, using Lagrangian multipliers
  //////////////////////////////////////////////////////////////////
#if 0
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute misfit gradient w.r.t. parameters A and E, using adjoint..."
              << std::endl;
  }

  P_M* misfitHessian = m_paramSpace->newMatrix();
  valueCenter = uqTgaLikelihoodRoutine<P_V,P_M>(params,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                NULL,
                                                misfitHessian,
                                                NULL);
  std::cout << "valueCenter = "     << valueCenter
            << "\nmisfitHessian = " << *misfitHessian
            << std::endl;
  delete misfitHessian;
#endif
  //////////////////////////////////////////////////////////////////
  // Compute Hessian(misfit), wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////

#if 0
  std::vector<double> misfitVarRatios(0);
  std::vector<double> allWTemps(0);
  std::vector<double> allWs(0);
  double valueCenter = 0.;
  valueCenter = tgaConstraintEquation<P_V,P_M>(params,
                                               *(m_calLikelihoodInfoVector[0]),
                                               false,
                                               &misfitVarRatios,
                                               &allWTemps,
                                               &allWs);
  std::cout << "In runGradTest()"
            << ": valueCenter = " << valueCenter
            << std::endl;

  tgaAdjointEquation<P_V,P_M>(params,
                              *(m_calLikelihoodInfoVector[0]),
                              allWTemps[allWTemps.size()-1],
                              misfitVarRatios,
                              &allLambdaTemps,
                              &allLambdas);

  // Compute \grad(misfit), using adjoint
  P_V LagrangianGradWrtParams(m_paramSpace->zeroVector());
  tgaDesignEquation<P_V,P_M>(params,
                             *(m_calLikelihoodInfoVector[0]),
                             allWTemps,
                             allWs,
                             allLambdaTemps,
                             allLambdas,
			     LagrangianGradWrtParams,
                             NULL);

  std::cout << "lagGrad = " << LagrangianGradWrtParams << std::endl;
#endif
  // Compute \grad(misfit), using finite differences
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": computing misfit gradient w.r.t. parameters A and E, using finite differences..."
              << std::endl;
  }

  // Finish
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_VALIDATION_H__
