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
  void  runGradTest   (char* argv[]);

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
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_5_K_min.dat", NULL,NULL,NULL);
  m_calLikelihoodInfoVector[1] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_25_K_min.dat",NULL,NULL,NULL);
  m_calLikelihoodInfoVector[2] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_50_K_min.dat",NULL,NULL,NULL);

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
  m_valLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_100_K_min.dat",NULL,NULL,NULL);

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

// Example for argv
// ./uqPyramid_gsl -i uqPyramid.inp Prefix RefA       RefE       GuessA     GuessE     RefDt RefContinuous RefNumDiscreteSamples Hessian Wdt Ldt IntegralNumIntervals FD    writeOut
// ./uqPyramid_gsl -i uqPyramid.inp case1_ 2.6000e+11 2.0000e+05 2.5910e+11 2.0090e+05 10.   0             1000                  1       1.  1.  1000                 1.e-8 0
template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runGradTest(char* argv[])
{
  std::string  outerPrefixName                  (argv[ 3]     ); // case1_
  double       refA                     = strtod(argv[ 4],NULL); // refA
  double       refE                     = strtod(argv[ 5],NULL); // refE
  double       guessA                   = strtod(argv[ 6],NULL); // guessA
  double       guessE                   = strtod(argv[ 7],NULL); // guessE
  double       refMaxTimeStep           = strtod(argv[ 8],NULL); // 10.
  bool         treatRefDataAsContinuous = atoi  (argv[ 9]     ); // 0
  unsigned int refNumDiscreteSamples    = strtod(argv[10],NULL); // 1000
  bool         computeHessian           = atoi  (argv[11]     ); // 1
  double       wMaxTimeStep             = strtod(argv[12],NULL); // .11
  double       lambdaMaxTimeStep        = strtod(argv[13],NULL); // .55
  unsigned int integralsNumIntervals    = atoi  (argv[14]     ); // 1001
  double       relativeFDStep           = strtod(argv[15],NULL); // 1.1e-7
  bool         writeOutput              = atoi  (argv[16]     ); // 0

  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runGradTest() with:"
              << "\n outerPrefixName          = " << outerPrefixName
              << "\n refA                     = " << refA
              << "\n refE                     = " << refE
              << "\n guessA                   = " << guessA
              << "\n guessE                   = " << guessE
              << "\n refMaxTimeStep           = " << refMaxTimeStep
              << "\n treatRefDataAsContinuous = " << treatRefDataAsContinuous
              << "\n refNumDiscreteSamples    = " << refNumDiscreteSamples
              << "\n computeHessian           = " << computeHessian
              << "\n wMaxTimeStep             = " << wMaxTimeStep
              << "\n lambdaMaxTimeStep        = " << lambdaMaxTimeStep
              << "\n integralsNumIntervals    = " << integralsNumIntervals
              << "\n relativeFDStep           = " << relativeFDStep
              << "\n writeOutput              = " << writeOutput
              << std::endl;
  }

  // Read parameters for temperature function (which is linear wrt time)
  // Read discrete measurements (which might be replaced by continuous measurements)
  m_calLikelihoodInfoVector.resize(1,NULL);
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,
                                                                        "tga/scenario_5_K_min.dat",
                                                                        &wMaxTimeStep,
                                                                        &lambdaMaxTimeStep,
                                                                        &integralsNumIntervals);
  m_calLikelihoodInfoVector[0]->setCheckingVariables(true,&relativeFDStep); // IMPORTANT

  // Miscellaneous
  std::ofstream ofs( (outerPrefixName+".m").c_str(), std::ofstream::out | std::ofstream::trunc); // Output file
  P_V paramValues(m_paramSpace->zeroVector());
  std::string tmpPrefixName;

  //////////////////////////////////////////////////////////////////
  // Step 1 of 4: Compute data for referenceW
  //////////////////////////////////////////////////////////////////

  paramValues[0] = refA;
  paramValues[1] = refE;

  uqTgaWClass<P_V,P_M> tmpW(*m_paramSpace,
                            *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  // Compute w(refA,refE), using derivative wrt temperature

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(refA,refE) using derivative wrt temperature..."
              << std::endl;
  }

  tmpW.computeUsingTemp(paramValues,
                        822., // maximumTemp
                        NULL, // referenceW
                        NULL);
  tmpPrefixName = outerPrefixName + "withTempW_";
  tmpW.printForMatlab(ofs,tmpPrefixName);
#endif

  // Compute w(refA,refE), using derivative wrt time

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute w(refA,refE) using derivative wrt time..."
              << std::endl;
  }

  tmpW.compute(paramValues,
               0.,             // maximumTime
               refMaxTimeStep, // maximum delta time
               false,          // computeWAndLambdaGradsAlso
               NULL,           // referenceW
               NULL,           // weightFunction
               NULL,           // misfitValue
               NULL);          // diffFunction

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": tmpW.times().size() = " << tmpW.times().size()
              << ", tmpW.times()[0] = "     << tmpW.times()[0]
              << ", tmpW.values()[0] = "    << tmpW.ws   ()[0]
              << ", tmpW.times()[max] = "   << tmpW.times()[tmpW.times().size()-1]
              << ", tmpW.values()[max] = "  << tmpW.ws   ()[tmpW.times().size()-1]
              << std::endl;
  }

  //////////////////////////////////////////////////////////////////
  // Step 2 of 4: Compute gradient(misfit) and Hessian(misfit) wrt A and E, using Lagrangian multipliers
  //////////////////////////////////////////////////////////////////

  //std::vector<double> constantTmpWs(tmpW.times().size(),.5);

  // Change the reference data to be used inside the likelihood routine,
  // from discrete (read from file) to continuous (computed and saved in 'refW')
  uqSampled1D1DFunctionClass* refW = NULL;
  if (treatRefDataAsContinuous) {
    refW = new uqSampled1D1DFunctionClass(tmpW.times(),
                                          tmpW.ws()); // tmpW.ws() or constantTmpWs
  }
  else {
    unsigned int tmpSize = tmpW.times().size();
    double minTime = tmpW.times()[0];
    double maxTime = tmpW.times()[tmpSize-1];
    double refDeltaTime = (maxTime-minTime)/((double) (refNumDiscreteSamples-1));

    std::vector<double> vecOfTimes (refNumDiscreteSamples,0.);
    std::vector<double> vecOfValues(refNumDiscreteSamples,0.);
    unsigned int startingId = 0;
    for (unsigned int j = 0; j < refNumDiscreteSamples; ++j) {
      if (j == 0) {
        vecOfTimes [j] = minTime;
        vecOfValues[j] = tmpW.ws()[0];
      }
      else if (j == (refNumDiscreteSamples - 1)) {
        vecOfTimes [j] = maxTime;
        vecOfValues[j] = tmpW.ws()[tmpSize-1];
      }
      else {
        vecOfTimes [j] = minTime + ((double) j)*refDeltaTime;
        tmpW.interpolate(vecOfTimes[j],
                         startingId,
                         &vecOfValues[j],
                         NULL,
                         NULL);
      }
    }

    refW = new uqSampled1D1DFunctionClass(vecOfTimes,
                                          vecOfValues);
  }
  m_calLikelihoodInfoVector[0]->changeReferenceW(*refW);

  // Change the weight functions
  unsigned int refSize = refW->domainValues().size();
  double       maxTime = refW->maxDomainValue();
  std::vector<double> vecOfDeltas (refSize,1.);
  std::vector<double> vecOfWeights(refSize,1.);
  std::vector<double> aux1        (refSize,0.);
  std::vector<double> aux2        (refSize,0.);
  std::vector<double> aux3        (refSize,0.);
  for (unsigned int j = 0; j < refSize; ++j) {
    aux1[j] = maxTime-refW->domainValues()[refSize-1-j];
    aux2[j] = vecOfDeltas                 [refSize-1-j];
    aux3[j] = vecOfWeights                [refSize-1-j];
  }

  uqSampled1D1DFunctionClass*  continuousWeightFunction      = NULL;
  uqSampled1D1DFunctionClass*  tildeContinuousWeightFunction = NULL;
  uqDeltaSeq1D1DFunctionClass* deltaWeightFunction           = NULL;
  uqDeltaSeq1D1DFunctionClass* tildeDeltaWeightFunction      = NULL;
  if (treatRefDataAsContinuous) {
    continuousWeightFunction = new uqSampled1D1DFunctionClass(refW->domainValues(),
                                                              vecOfWeights);
    tildeContinuousWeightFunction = new uqSampled1D1DFunctionClass(aux1,
                                                                   aux3);

    m_calLikelihoodInfoVector[0]->changeWeightFunctions(continuousWeightFunction,tildeContinuousWeightFunction); // or NULL,NULL
  }
  else {
    deltaWeightFunction = new uqDeltaSeq1D1DFunctionClass(refW->domainValues(),
                                                          vecOfDeltas,
                                                          vecOfWeights);
    tildeDeltaWeightFunction = new uqDeltaSeq1D1DFunctionClass(aux1,
                                                               aux2,
                                                               aux3);

    m_calLikelihoodInfoVector[0]->changeWeightFunctions(deltaWeightFunction,tildeDeltaWeightFunction); // or NULL,NULL
  }

  // Perform calculations
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute gradient(misfit) and Hessian(misfit) w.r.t. parameters A and E, using adjoint..."
              << std::endl;
  }

  paramValues[0] = guessA;
  paramValues[1] = guessE;

  P_V gradWithLM(m_paramSpace->zeroVector());
  P_M* hessianWithLM = NULL;
  if (computeHessian) hessianWithLM = m_paramSpace->newMatrix();
  double guessMisfit = 0.;
  guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                &gradWithLM,
                                                hessianWithLM,
                                                NULL);
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": guessMisfit = " << guessMisfit
              << std::endl;
  }
  delete hessianWithLM;

  if (writeOutput) {
    ofs << "\n" << outerPrefixName << "hessianIsAvailable = " << computeHessian << ";" << std::endl;
    ofs << "\n" << outerPrefixName << "refA = "               << refA           << ";" << std::endl;
    ofs << "\n" << outerPrefixName << "refE = "               << refE           << ";" << std::endl;
    ofs << "\n" << outerPrefixName << "guessA = "             << guessA         << ";" << std::endl;
    ofs << "\n" << outerPrefixName << "guessE = "             << guessE         << ";" << std::endl;

    tmpPrefixName = outerPrefixName + "refW_";
    m_calLikelihoodInfoVector[0]->m_refWPtrForPlot->printForMatlab(ofs,tmpPrefixName);

    tmpPrefixName = outerPrefixName + "w_";
    m_calLikelihoodInfoVector[0]->m_wPtrForPlot->printForMatlab(ofs,tmpPrefixName);

    if (computeHessian) {
      tmpPrefixName = outerPrefixName + "wA_";
      m_calLikelihoodInfoVector[0]->m_wAPtrForPlot->printForMatlab(ofs,tmpPrefixName);

      tmpPrefixName = outerPrefixName + "wE_";
      m_calLikelihoodInfoVector[0]->m_wEPtrForPlot->printForMatlab(ofs,tmpPrefixName);
    }

    tmpPrefixName = outerPrefixName + "diff_";
    m_calLikelihoodInfoVector[0]->m_diffPtrForPlot->printForMatlab(ofs,tmpPrefixName);

    tmpPrefixName = outerPrefixName + "lambda_";
    m_calLikelihoodInfoVector[0]->m_lambdaPtrForPlot->printForMatlab(ofs,tmpPrefixName);

    if (computeHessian) {
      tmpPrefixName = outerPrefixName + "lambdaA_";
      m_calLikelihoodInfoVector[0]->m_lambdaAPtrForPlot->printForMatlab(ofs,tmpPrefixName);

      tmpPrefixName = outerPrefixName + "lambdaE_";
      m_calLikelihoodInfoVector[0]->m_lambdaEPtrForPlot->printForMatlab(ofs,tmpPrefixName);
    }
  }

  // Finish
  delete tildeDeltaWeightFunction;
  delete deltaWeightFunction;
  delete tildeContinuousWeightFunction;
  delete continuousWeightFunction;
  delete refW;
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_VALIDATION_H__
