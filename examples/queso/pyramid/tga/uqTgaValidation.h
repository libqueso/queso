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

struct uqTgaTestInputStruct {
  std::string  outerPrefixName;
  double       refA;
  double       refE;
  double       guessA;
  double       guessE;
  double       refMaxTimeStep;
  bool         treatRefDataAsContinuous;
  unsigned int refNumDiscreteSamples;
  bool         computeHessian;
  double       wMaxTimeStep;
  double       lambdaMaxTimeStep;
  unsigned int integralsNumIntervals;
  double       relativeFDStep;
  bool         writeOutput;
};

struct uqTgaTestVarsStruct {
  uqTgaTestVarsStruct();
 ~uqTgaTestVarsStruct();

  uqSampled1D1DFunctionClass*  refW;
  std::ofstream*               ofs;
  uqSampled1D1DFunctionClass*  continuousWeightFunction;
  uqSampled1D1DFunctionClass*  tildeContinuousWeightFunction;
  uqDeltaSeq1D1DFunctionClass* deltaWeightFunction;
  uqDeltaSeq1D1DFunctionClass* tildeDeltaWeightFunction;
};

uqTgaTestVarsStruct::uqTgaTestVarsStruct()
  :
  refW                         (NULL),
  ofs                          (NULL),
  continuousWeightFunction     (NULL),
  tildeContinuousWeightFunction(NULL),
  deltaWeightFunction          (NULL),
  tildeDeltaWeightFunction     (NULL)
{
}

uqTgaTestVarsStruct::~uqTgaTestVarsStruct()
{
  delete refW;
  delete ofs;
  delete continuousWeightFunction;
  delete tildeContinuousWeightFunction;
  delete deltaWeightFunction;
  delete tildeDeltaWeightFunction;
}

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
  void  runTests      (char* argv[]);

private:
  void  runCalibrationStage();
  void  runValidationStage ();
  void  runComparisonStage ();

  void  runTempTimeTest    ();
  void  prepareRefForTests ();
  void  runGradTest        ();
  void  runOptimizationTest();

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

  uqTgaTestInputStruct                              m_testInput;
  uqTgaTestVarsStruct                               m_testVars;
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
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runTests(char* argv[])
{
  m_testInput.outerPrefixName          =       (argv[ 3]     ); // case1_
  m_testInput.refA                     = strtod(argv[ 4],NULL); // refA
  m_testInput.refE                     = strtod(argv[ 5],NULL); // refE
  m_testInput.guessA                   = strtod(argv[ 6],NULL); // guessA
  m_testInput.guessE                   = strtod(argv[ 7],NULL); // guessE
  m_testInput.refMaxTimeStep           = strtod(argv[ 8],NULL); // 10.
  m_testInput.treatRefDataAsContinuous = atoi  (argv[ 9]     ); // 0
  m_testInput.refNumDiscreteSamples    = strtod(argv[10],NULL); // 1000
  m_testInput.computeHessian           = atoi  (argv[11]     ); // 1
  m_testInput.wMaxTimeStep             = strtod(argv[12],NULL); // .11
  m_testInput.lambdaMaxTimeStep        = strtod(argv[13],NULL); // .55
  m_testInput.integralsNumIntervals    = atoi  (argv[14]     ); // 1001
  m_testInput.relativeFDStep           = strtod(argv[15],NULL); // 1.1e-7
  m_testInput.writeOutput              = atoi  (argv[16]     ); // 0

  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runTests() with:"
              << "\n outerPrefixName          = " << m_testInput.outerPrefixName
              << "\n refA                     = " << m_testInput.refA
              << "\n refE                     = " << m_testInput.refE
              << "\n guessA                   = " << m_testInput.guessA
              << "\n guessE                   = " << m_testInput.guessE
              << "\n refMaxTimeStep           = " << m_testInput.refMaxTimeStep
              << "\n treatRefDataAsContinuous = " << m_testInput.treatRefDataAsContinuous
              << "\n refNumDiscreteSamples    = " << m_testInput.refNumDiscreteSamples
              << "\n computeHessian           = " << m_testInput.computeHessian
              << "\n wMaxTimeStep             = " << m_testInput.wMaxTimeStep
              << "\n lambdaMaxTimeStep        = " << m_testInput.lambdaMaxTimeStep
              << "\n integralsNumIntervals    = " << m_testInput.integralsNumIntervals
              << "\n relativeFDStep           = " << m_testInput.relativeFDStep
              << "\n writeOutput              = " << m_testInput.writeOutput
              << std::endl;
  }

  // Read parameters for temperature function (which is linear wrt time)
  // Read discrete measurements (which might be replaced by continuous measurements)
  m_calLikelihoodInfoVector.resize(1,NULL);
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,
                                                                        "tga/scenario_5_K_min.dat",
                                                                        &m_testInput.wMaxTimeStep,
                                                                        &m_testInput.lambdaMaxTimeStep,
                                                                        &m_testInput.integralsNumIntervals);

  // Output file
  m_testVars.ofs = new std::ofstream( (m_testInput.outerPrefixName+".m").c_str(), std::ofstream::out | std::ofstream::trunc);

  //runTempTimeTest    ();
  prepareRefForTests ();
  runGradTest        ();
  //runOptimizationTest();

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runTempTimeTest()
{
  P_V paramValues(m_paramSpace->zeroVector());
  std::string tmpPrefixName;

  paramValues[0] = 2.6000e+11;
  paramValues[1] = 2.0000e+05;

  uqTgaWClass<P_V,P_M> tmpW(*m_paramSpace,
                            *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  // Compute w using derivative wrt temperature
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runTempTimeTest()"
              << ": computing w using derivative wrt temperature..."
              << std::endl;
  }

  tmpW.computeUsingTemp(paramValues,
                        822., // maximumTemp
                        NULL, // referenceW
                        NULL);

  tmpPrefixName = m_testInput.outerPrefixName + "withTempW_";
  tmpW.printForMatlab(*m_testVars.ofs,tmpPrefixName);
#endif

  // Compute w using derivative wrt time
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runTempTimeTest()"
              << ": computing w using derivative wrt time..."
              << std::endl;
  }

  tmpW.compute(paramValues,
               0.,                         // maximumTime
               m_testInput.refMaxTimeStep, // maximum delta time
               false,                      // computeWAndLambdaGradsAlso
               NULL,                       // referenceW
               NULL,                       // weightFunction
               NULL,                       // misfitValue
               NULL);                      // diffFunction

  tmpPrefixName = m_testInput.outerPrefixName + "withTimeW_";
  tmpW.printForMatlab(*m_testVars.ofs,tmpPrefixName);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::prepareRefForTests()
{
  //////////////////////////////////////////////////////////////////
  // Compute data for refW
  //////////////////////////////////////////////////////////////////
  P_V paramValues(m_paramSpace->zeroVector());

  paramValues[0] = m_testInput.refA;
  paramValues[1] = m_testInput.refE;

  uqTgaWClass<P_V,P_M> tmpW(*m_paramSpace,
                            *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

  tmpW.compute(paramValues,
               0.,                         // maximumTime
               m_testInput.refMaxTimeStep, // maximum delta time
               false,                      // computeWAndLambdaGradsAlso
               NULL,                       // referenceW
               NULL,                       // weightFunction
               NULL,                       // misfitValue
               NULL);                      // diffFunction

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::prepareRefForTests()"
              << ": tmpW.times().size() = " << tmpW.times().size()
              << ", tmpW.times()[0] = "     << tmpW.times()[0]
              << ", tmpW.values()[0] = "    << tmpW.ws   ()[0]
              << ", tmpW.times()[max] = "   << tmpW.times()[tmpW.times().size()-1]
              << ", tmpW.values()[max] = "  << tmpW.ws   ()[tmpW.times().size()-1]
              << std::endl;
  }

  //////////////////////////////////////////////////////////////////
  // Set refW
  //////////////////////////////////////////////////////////////////

  //std::vector<double> constantTmpWs(tmpW.times().size(),.5);
  if (m_testInput.treatRefDataAsContinuous) {
    m_testVars.refW = new uqSampled1D1DFunctionClass(tmpW.times(),
                                                     tmpW.ws()); // tmpW.ws() or constantTmpWs
  }
  else {
    unsigned int tmpSize = tmpW.times().size();
    double minTime = tmpW.times()[0];
    double maxTime = tmpW.times()[tmpSize-1];
    double refDeltaTime = (maxTime-minTime)/((double) (m_testInput.refNumDiscreteSamples-1));

    std::vector<double> vecOfTimes (m_testInput.refNumDiscreteSamples,0.);
    std::vector<double> vecOfValues(m_testInput.refNumDiscreteSamples,0.);
    unsigned int startingId = 0;
    for (unsigned int j = 0; j < m_testInput.refNumDiscreteSamples; ++j) {
      if (j == 0) {
        vecOfTimes [j] = minTime;
        vecOfValues[j] = tmpW.ws()[0];
      }
      else if (j == (m_testInput.refNumDiscreteSamples - 1)) {
        vecOfTimes [j] = maxTime;
        vecOfValues[j] = tmpW.ws()[tmpSize-1];
      }
      else {
        vecOfTimes [j] = minTime + ((double) j)*refDeltaTime;
        tmpW.interpolate(vecOfTimes[j],
                         startingId,
                         1.,
                         &vecOfValues[j],
                         NULL,
                         NULL);
      }
    }

    m_testVars.refW = new uqSampled1D1DFunctionClass(vecOfTimes,
                                                     vecOfValues);
  }

  //////////////////////////////////////////////////////////////////
  // Prepare weight functions
  //////////////////////////////////////////////////////////////////
  unsigned int refSize = m_testVars.refW->domainValues().size();
  double       maxTime = m_testVars.refW->maxDomainValue();
  std::vector<double> vecOfDeltas (refSize,INFINITY);
  std::vector<double> vecOfWeights(refSize,1.);
  std::vector<double> aux1        (refSize,0.);
  std::vector<double> aux2        (refSize,0.);
  std::vector<double> aux3        (refSize,0.);
  for (unsigned int j = 0; j < refSize; ++j) {
    aux1[j] = maxTime-m_testVars.refW->domainValues()[refSize-1-j];
    aux2[j] = vecOfDeltas                 [refSize-1-j];
    aux3[j] = vecOfWeights                [refSize-1-j];
  }

  if (m_testInput.treatRefDataAsContinuous) {
    m_testVars.continuousWeightFunction = new uqSampled1D1DFunctionClass(m_testVars.refW->domainValues(),
                                                                         vecOfWeights);
    m_testVars.tildeContinuousWeightFunction = new uqSampled1D1DFunctionClass(aux1,
                                                                              aux3);
  }
  else {
    m_testVars.deltaWeightFunction = new uqDeltaSeq1D1DFunctionClass(m_testVars.refW->domainValues(),
                                                                     vecOfDeltas,
                                                                     vecOfWeights);
    m_testVars.tildeDeltaWeightFunction = new uqDeltaSeq1D1DFunctionClass(aux1,
                                                                          aux2,
                                                                          aux3);
  }

  //////////////////////////////////////////////////////////////////
  // Set likelihood information appropriately
  //////////////////////////////////////////////////////////////////

  // Change the reference data to be used inside the likelihood routine,
  // from the (discrete) one read from file
  // to the (discrete or continuous) one computed and saved in 'refW'.
  m_calLikelihoodInfoVector[0]->changeReferenceW(*m_testVars.refW);

  if (m_testInput.treatRefDataAsContinuous) {
    m_calLikelihoodInfoVector[0]->changeWeightFunctions(m_testVars.continuousWeightFunction,m_testVars.tildeContinuousWeightFunction); // or NULL,NULL
  }
  else {
    m_calLikelihoodInfoVector[0]->changeWeightFunctions(m_testVars.deltaWeightFunction,m_testVars.tildeDeltaWeightFunction); // or NULL,NULL
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runGradTest()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runGradTest()..."
              << std::endl;
  }

  int iRC;
  struct timeval timeval0;
  struct timeval timeval1;

  P_V paramValues(m_paramSpace->zeroVector());
  paramValues[0] = m_testInput.guessA;
  paramValues[1] = m_testInput.guessE;

  P_V gradWithLM(m_paramSpace->zeroVector());
  P_M* hessianWithLM = NULL;
  if (m_testInput.computeHessian) hessianWithLM = m_paramSpace->newMatrix();
  double guessMisfit = 0.;

#if 1
  // Run without checking, just in order to measure time
  iRC = gettimeofday(&timeval0, NULL);
  guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                &gradWithLM,
                                                hessianWithLM,
                                                NULL);
  iRC = gettimeofday(&timeval1, NULL);
  double total_usecs = (timeval1.tv_sec * 1.e+6 + timeval1.tv_usec) - (timeval0.tv_sec * 1.e+6 + timeval0.tv_usec);
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": total_usecs = " << total_usecs
              << std::endl;
  }
#else
  // Run with checking against finite differences
  m_calLikelihoodInfoVector[0]->setCheckingVariables(true,&m_testInput.relativeFDStep); // IMPORTANT
  guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                &gradWithLM,
                                                hessianWithLM,
                                                NULL);
  m_calLikelihoodInfoVector[0]->setCheckingVariables(false,&m_testInput.relativeFDStep); // IMPORTANT
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": guessMisfit = " << guessMisfit
              << std::endl;
  }

  if (m_testInput.writeOutput) {
    std::string tmpPrefixName;

    *m_testVars.ofs << "\n" << m_testInput.outerPrefixName << "hessianIsAvailable = " << m_testInput.computeHessian << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testInput.outerPrefixName << "refA = "               << m_testInput.refA           << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testInput.outerPrefixName << "refE = "               << m_testInput.refE           << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testInput.outerPrefixName << "guessA = "             << m_testInput.guessA         << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testInput.outerPrefixName << "guessE = "             << m_testInput.guessE         << ";" << std::endl;

    tmpPrefixName = m_testInput.outerPrefixName + "refW_";
    m_calLikelihoodInfoVector[0]->m_refWPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    tmpPrefixName = m_testInput.outerPrefixName + "w_";
    m_calLikelihoodInfoVector[0]->m_wPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    if (m_testInput.computeHessian) {
      tmpPrefixName = m_testInput.outerPrefixName + "wA_";
      m_calLikelihoodInfoVector[0]->m_wAPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

      tmpPrefixName = m_testInput.outerPrefixName + "wE_";
      m_calLikelihoodInfoVector[0]->m_wEPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);
    }

    tmpPrefixName = m_testInput.outerPrefixName + "diff_";
    m_calLikelihoodInfoVector[0]->m_diffPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    tmpPrefixName = m_testInput.outerPrefixName + "lambda_";
    m_calLikelihoodInfoVector[0]->m_lambdaPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    if (m_testInput.computeHessian) {
      tmpPrefixName = m_testInput.outerPrefixName + "lambdaA_";
      m_calLikelihoodInfoVector[0]->m_lambdaAPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

      tmpPrefixName = m_testInput.outerPrefixName + "lambdaE_";
      m_calLikelihoodInfoVector[0]->m_lambdaEPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);
    }
  }
#endif

  delete hessianWithLM;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runOptimizationTest()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runOptimizationTest()"
              << std::endl;
  }

  P_V paramValues(m_paramSpace->zeroVector());
  paramValues[0] = m_testInput.guessA;
  paramValues[1] = m_testInput.guessE;

  double objFunctionValue = 0.;
  P_V    objFunctionGradient(m_paramSpace->zeroVector());
  P_M*   objFunctionHessian = m_paramSpace->newMatrix();
  P_V    paramsStep(m_paramSpace->zeroVector());

  // Apply Newton method
  bool newtonSucceeded = false;
  int newtonFailureReason = 1;
  unsigned int newtonLoopId = 0;
  while (newtonLoopId < 10) {
    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "Beggining newtonLoopId = " << newtonLoopId
                << " with params = "           << stringA << " " << stringE
                << std::endl;
    }

    // Compute objective function: value, gradient and Hessian
    objFunctionValue = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                       (const void *)&m_calLikelihoodInfoVector,
                                                       &objFunctionGradient,
                                                       objFunctionHessian,
                                                       NULL);
    double gradNorm = objFunctionGradient.norm2();

    double a = (*objFunctionHessian)(0,0);
    double b = (*objFunctionHessian)(0,1);
    double c = (*objFunctionHessian)(1,0);
    double d = (*objFunctionHessian)(1,1);
    double determinant = a*d - b*c;
    double frobNorm = sqrt(a*a + b*b + c*c + d*d);
    double e2 = .5*(a + d + sqrt( (a-d)*(a-d) + 4*b*c ));
    double e1 = .5*(a + d - sqrt( (a-d)*(a-d) + 4*b*c ));

    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In newtonLoopId = "    << newtonLoopId
                << ", with params = "      << stringA << " " << stringE
                << ": objFunctionValue = " << objFunctionValue
                << ", objFunctionGrad = "  << objFunctionGradient
                << ", gradNorm = "         << gradNorm
                << ", Hessian = \n"        << *objFunctionHessian
                << ", determinant = "      << determinant
                << ", frobNorm = "         << frobNorm
                << ", e2 = "               << e2
                << ", e1 = "               << e1
                << ", e2*e1 = "            << e2*e1
                << std::endl;
    }

    unsigned int tauId = 0;
    double tau = 0.;
    while (determinant <= 0.) {
      tauId++;
      tau += frobNorm;
      double newA = a + tau;
      double newD = d + tau;
      determinant = newA*newD - b*c;
    }

    if (tauId > 0) {
      (*objFunctionHessian)(0,0) += tau;
      (*objFunctionHessian)(1,1) += tau;

      a = (*objFunctionHessian)(0,0);
      b = (*objFunctionHessian)(0,1);
      c = (*objFunctionHessian)(1,0);
      d = (*objFunctionHessian)(1,1);
      determinant = a*d - b*c;
      frobNorm = sqrt(a*a + b*b + c*c + d*d);
      e2 = .5*(a + d + sqrt( (a-d)*(a-d) + 4*b*c ));
      e1 = .5*(a + d - sqrt( (a-d)*(a-d) + 4*b*c ));

      if (m_env.rank() == 0) {
        char stringA[64];
        char stringE[64];
        sprintf(stringA,"%12.6e",paramValues[0]);
        sprintf(stringE,"%12.6e",paramValues[1]);
        std::cout << "In newtonLoopId = "   << newtonLoopId
                  << ", with params = "     << stringA << " " << stringE
                  << ": tauId = "           << tauId
                  << ", tau = "             << tau
                  << ", new Hessian = \n"   << *objFunctionHessian
                  << ", new determinant = " << determinant
                  << ", new frobNorm = "    << frobNorm
                  << ", new e2 = "          << e2
                  << ", new e1 = "          << e1
                  << ", e2*e1 = "           << e2*e1
                  << std::endl;
      }
    }

    // Check sufficiently small gradient norm
    if (gradNorm < 1.e-7) {
      newtonSucceeded = true;
      break;
    }

    // Compute Newton step
    objFunctionHessian->invertMultiply(-1.*objFunctionGradient,paramsStep);
    double directionalDerivative = scalarProduct(objFunctionGradient,paramsStep);

    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In newtonLoopId = "         << newtonLoopId
                << " with params = "            << stringA << " " << stringE
                << ": paramsStep = "            << paramsStep
                << ", directionalDerivative = " << directionalDerivative
                << std::endl;
    }

    // Perform line search
    bool lineSearchSucceeded = false;
    unsigned int lineSearchLoopId = 0;
    double alphaFactor = .5;
    double alpha = 1./alphaFactor;
    while (lineSearchLoopId < 30) {
      double c1 = 1.e-4;
      //double c2 = .9;
      P_V attemptedParamValues(m_paramSpace->zeroVector());

      alpha *= alphaFactor;
      double lineSearchThreshold = objFunctionValue + alpha * c1 * directionalDerivative;
      attemptedParamValues = paramValues + alpha*paramsStep;

      double attemptedObjFunctionValue = uqTgaLikelihoodRoutine<P_V,P_M>(attemptedParamValues,
                                                                         (const void *)&m_calLikelihoodInfoVector,
                                                                         NULL,
                                                                         NULL,
                                                                         NULL);
      if (m_env.rank() == 0) {
        char stringA[64];
        char stringE[64];
        sprintf(stringA,"%12.6e",paramValues[0]);
        sprintf(stringE,"%12.6e",paramValues[1]);
        std::cout << "In newtonLoopId = "             << newtonLoopId
                  << " with params = "                << stringA << " " << stringE
                  << ", alpha = "                     << alpha
                  << ": attemptedObjFunctionValue = " << attemptedObjFunctionValue
                  << ", lineSearchThreshold = "       << lineSearchThreshold
                  << std::endl;
      }

      // Check sufficient decrease
      if (attemptedObjFunctionValue <= lineSearchThreshold) {
        lineSearchSucceeded = true;
        // Update position
        paramValues = attemptedParamValues;
        break;
      }

      lineSearchLoopId++;
    }

    if (!lineSearchSucceeded) {
      newtonFailureReason = 2;
      break;
    }

    newtonLoopId++;
  }
  delete objFunctionHessian;

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runOptimizationTest():";
    if (newtonSucceeded) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << " solution = " << stringA << " " << stringE;
    }
    else {
      std::cout << " Newton failed: reason = " << newtonFailureReason;
    }
    std::cout << std::endl;
  }

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runOptimizationTest()"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_VALIDATION_H__
