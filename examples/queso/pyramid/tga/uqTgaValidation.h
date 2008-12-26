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
#include <uqTgaComputableW.h>
#include <uqTgaOdes.h>
#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
tgaLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr, P_V* gradVector, P_M* hessianMatrix, P_V* hessianEffect)
{
  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering tgaLikelihoodRoutine()..." << std::endl;
  }

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "In tgaLikelihoodRoutine()"
              << ", A = " << paramValues[0]
              << ", E = " << paramValues[1]
              << std::endl;
  }

  bool computeLambda = false;
  bool computeWAndLambdaGradsAlso = false;

  if (gradVector) computeLambda = true;
  if (hessianMatrix || hessianEffect) {
    computeLambda              = true;
    computeWAndLambdaGradsAlso = true;
  }

  double resultValue = 0.;
  const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *>& info = *((const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *> *) functionDataPtr);

  /////////////////////////////////////////////////////////////////////////////
  // Loop on scenarios
  /////////////////////////////////////////////////////////////////////////////
  for (unsigned int i = 0; i < info.size(); ++i) {
#ifdef QUESO_USE_NEW_W_CLASS
    /////////////////////////////////////////////////////////////////////////////
    // Compute W
    /////////////////////////////////////////////////////////////////////////////
    uqTgaComputableWClass<P_V,P_M> wObj(info[i]->m_paramSpace,
                                        *(info[i]->m_temperatureFunctionObj),
                                        false); // useOdeWithDerivativeWrtTime // COMPATIBILITY WITH OLD VERSION
    wObj.compute(paramValues,
                 computeWAndLambdaGradsAlso,
                 info[i]->m_referenceW,
                 1900.); // COMPATIBILITY WITH OLD VERSION

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to the likelihood ("tmpResultValue")
    /////////////////////////////////////////////////////////////////////////////
    unsigned int numMisfits = info[i]->m_referenceW->times().size();
    std::vector<double> twiceDiffVec(numMisfits,0.);
    double tmpResultValue = 0.; // COMPATIBILITY WITH OLD VERSION
    for (unsigned int j = 0; j < numMisfits; ++j) {
      double diff = wObj.diffsForMisfit()[j];
      double ratio = diff/info[i]->m_referenceW->variances()[j];
      twiceDiffVec[j] = 2.*ratio;
      tmpResultValue += diff*diff/info[i]->m_referenceW->variances()[j]; //COMPATIBILITY WITH OLD VERSION

      if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
        std::cout << "In tgaLikelihoodRoutine()"
                  << ", misfitId = "     << j
                  << ", measuredTemp = " << info[i]->m_referenceW->temps()[j]
                  << ", measuredW = "    << info[i]->m_referenceW->ws()[j]
                  << ": computedW = "    << info[i]->m_referenceW->ws()[j] + wObj.diffsForMisfit()[j]
                  << ", misfitValue = "  << diff
                  << std::endl;
      }
    }
    resultValue += tmpResultValue;

    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In tgaLikelihoodRoutine()"
                << ", A = "                              << stringA
                << ", E = "                              << stringE
                << ", beta = "                           << info[i]->m_temperatureFunctionObj->deriv(0.)
                << ": finished ode loop after "          << wObj.ws().size()
                << " iterations, with weigthedMisfit = " << tmpResultValue
                << std::endl;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute Lambda, if necessary
    /////////////////////////////////////////////////////////////////////////////
#if 0
    uqTgaLambdaClass<P_V,P_M> lambdaObj(info[i]->m_paramSpace,
                                        *(info[i]->m_temperatureFunctionObj),
                                        false); // useOdeWithDerivativeWrtTime
    if (computeLambda) {
      lambdaObj.compute(paramValues,
                        computeWAndLambdaGradsAlso,
                        twiceDiffVec,
                        info[i]->m_referenceW->continuous(),
                        wObj);
    }
#endif
    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to gradVector
    /////////////////////////////////////////////////////////////////////////////
    if (gradVector) {
      P_V tmpGradVector(info[i]->m_paramSpace.zeroVector());
      *gradVector += tmpGradVector;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to hessianMatrix and/or hessianEffect
    /////////////////////////////////////////////////////////////////////////////
    if (hessianMatrix || hessianEffect) {
      P_M* tmpMatrix = info[i]->m_paramSpace.newMatrix();

      if (hessianMatrix) *hessianMatrix += *tmpMatrix;
      if (hessianEffect) *hessianEffect += ((*tmpMatrix) * paramValues);

      delete tmpMatrix;
    }
#else
    // Compute likelihood for scenario
    resultValue += tgaConstraintEquation<P_V,P_M>(paramValues,*(info[i]),true,NULL,NULL,NULL);
#endif
  }

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving tgaLikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}

//********************************************************
// QoI function object for both forward problems of the validation cycle.
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data class that carries the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
uqTgaQoiInfoStruct
{
  uqTgaQoiInfoStruct(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                     double                             initialTemp,
                     double                             beta,
                     double                             criticalW,
                     double                             criticalTime);
 ~uqTgaQoiInfoStruct();

  const uqVectorSpaceClass<P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*           m_temperatureFunctionObj;
  double                             m_criticalW;
  double                             m_criticalTime;
#ifdef QUESO_USE_NEW_W_CLASS
#else
  bool                               m_useTimeAsDomainVariable;
#endif
};

template<class P_V,class P_M,class Q_V,class Q_M>
  uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>::uqTgaQoiInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  double                             initialTemp,
  double                             beta,
  double                             criticalW,
  double                             criticalTime)
  :
  m_paramSpace            (paramSpace),
  m_temperatureFunctionObj(new uqLinear1D1DFunctionClass(-INFINITY,INFINITY,0.,initialTemp,beta)),
  m_criticalW             (criticalW),
  m_criticalTime          (criticalTime)
#ifdef QUESO_USE_NEW_W_CLASS
#else
  ,
  m_useTimeAsDomainVariable(false) // COMPATIBILITY WITH OLD VERSION
#endif
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>::~uqTgaQoiInfoStruct()
{
  delete m_temperatureFunctionObj;
}

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void tgaQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  const uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>& info = *((const uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M> *) functionDataPtr);

#ifdef QUESO_USE_NEW_W_CLASS
  uqTgaComputableWClass<P_V,P_M> wObj(info.m_paramSpace,
                                      *(info.m_temperatureFunctionObj),
                                      false); // useOdeWithDerivativeWrtTime = false // COMPATIBILITY WITH OLD VERSION
  wObj.compute(paramValues,
               false, // computeGradAlso
               NULL,  // referenceW
               info.m_criticalTime*info.m_temperatureFunctionObj->deriv(0.)); // Should add initialTemp --> COMPATIBILITY WITH OLD VERSION

  unsigned int tmpSize = wObj.ws().size();
  qoiValues[0] = wObj.ws()[tmpSize-1]; // QoI = mass fraction remaining at critical time
#else
  double A = paramValues[0];
  double E = paramValues[1];

  bool   useTimeAsDomainVariable = info.m_useTimeAsDomainVariable;
  double beta                    = info.m_temperatureFunctionObj->deriv(0.);
  double initialTemp             = info.m_temperatureFunctionObj->value(0.);
  double criticalW               = info.m_criticalW;
  double criticalTime            = info.m_criticalTime;

  double stateTimeDotOdeParameters[]={A,E,beta,initialTemp};
  double stateTempDotOdeParameters[]={A,E,beta};
      	
  // integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,1);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sysTime = {tgaStateTimeDotOdeFunction, NULL, 1, (void *)stateTimeDotOdeParameters};
        gsl_odeiv_system     sysTemp = {tgaStateTempDotOdeFunction, NULL, 1, (void *)stateTempDotOdeParameters}; 
	
  double previousTemp = 0.;
  double currentTemp  = initialTemp;
  double deltaTemp    = 1e-3;
  double maximumTemp  = criticalTime*beta;
  double crossingTemp = 0.;

  double currentTime  = 0.;
  double deltaTime    = 1e-3;
  double maximumTime  = (maximumTemp-initialTemp)/beta; // CONVERSION TO TIME

  double currentW [1];
  double previousW[1];
  currentW [0]=1.;
  previousW[0]=1.;
	
  unsigned int loopId = 0;
  while ((currentTemp < maximumTemp) &&
         (currentW[0] > criticalW  )) {
    loopId++;
    int status = 0;
    if (useTimeAsDomainVariable) {
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, maximumTime, &deltaTime, currentW);
      currentTemp = initialTemp + beta*currentTime;
    }
    else {
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, maximumTemp, &deltaTemp, currentW);
    }
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        paramValues.env().rank(),
                        "tgaQoiRoutine()",
                        "gsl_odeiv_evolve_apply() failed");

    if (currentW[0] <= criticalW) {
      crossingTemp = previousTemp + (currentTemp - previousTemp) * (previousW[0]-criticalW)/(previousW[0]-currentW[0]);
    }
		
    previousTemp = currentTemp;
    previousW[0] = currentW[0];
  }

  if (criticalW    > 0.) qoiValues[0] = (crossingTemp-initialTemp)/beta; // QoI = time to achieve critical mass
  if (criticalTime > 0.) qoiValues[0] = currentW[0];                     // QoI = mass fraction remaining at critical time
	
  if ((paramValues.env().verbosity() >= 3) && (paramValues.env().rank() == 0)) {
    std::cout << "In tgaQoiRoutine()"
              << ", A = "            << A
              << ", E = "            << E
              << ", beta = "         << beta
              << ", loopSize = "     << loopId
              << ", criticalTime = " << criticalTime
              << ", criticalW = "    << criticalW
              << ": qoi = "          << qoiValues[0]
              << std::endl;
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);
#endif

  return;
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

  void  run            ();
  void  runGradTest    ();
  //void  runMinimization();

private:
  void  runCalibrationStage();
  void  runValidationStage ();
  void  runComparisonStage ();

  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_cycle;

  uqAsciiTableClass<P_V,P_M>*               m_paramsTable;
  const EpetraExt::DistArray<std::string>*  m_paramNames;         // instantiated outside this class!!
  P_V*                                      m_paramMinValues;     // instantiated outside this class!!
  P_V*                                      m_paramMaxValues;     // instantiated outside this class!!
  P_V*                                      m_paramInitialValues; // instantiated outside this class!!
  uqVectorSpaceClass<P_V,P_M>*              m_paramSpace;
  uqVectorSetClass<P_V,P_M>*                m_paramDomain;

  uqAsciiTableClass<P_V,P_M>*               m_qoisTable;
  const EpetraExt::DistArray<std::string>*  m_qoiNames; // instantiated outside this class!!
  uqVectorSpaceClass<Q_V,Q_M>*              m_qoiSpace;

  double                                    m_predBeta;
  double                                    m_predInitialTemp;
  double                                    m_predCriticalW;
  double                                    m_predCriticalTime;

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
                                                                         tgaLikelihoodRoutine<P_V,P_M>,
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

  m_cycle->setCalFP(tgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
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
                                                                         tgaLikelihoodRoutine<P_V,P_M>,
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

  m_cycle->setValFP(tgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
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

  // Read discrete measurements
  m_calLikelihoodInfoVector.resize(1,NULL);
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_5_K_min.dat");

  // Miscellaneous
  std::ofstream ofs("nada1.m", std::ofstream::out | std::ofstream::trunc); // Output file
  P_V params(m_paramSpace->zeroVector());

  // Compute w(A_reference,E_reference)
  params[0] = 2.6090e+11; // Reference _A
  params[1] = 1.9910e+05; // Reference _E
  uqTgaComputableWClass<P_V,P_M> referenceW(m_env, *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj), false); // useOdeWithDerivativeWrtTime = false
  referenceW.compute(params,
                     true, // computeGradAlso
                     NULL, // referenceW
                     0.);  // maximumTemp

  unsigned int tmpSize = referenceW.times().size();
  ofs << "referenceTemps = zeros(" << tmpSize
      << ",1);"
      << "\nreferenceWs = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nreferenceTemps(" << i+1 << ",1) = " << referenceW.temps()[i] << ";"
        << "\nreferenceWs("    << i+1 << ",1) = " << referenceW.ws()   [i] << ";";
  }

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute misfit gradient w.r.t. parameters A and E, using adjoint..."
              << std::endl;
  }

  // Compute w(A_guess,E_guess)
  params[0] = 2.6000e+11; // Guess _A
  params[1] = 2.0000e+05; // Guess _E
  uqTgaComputableWClass<P_V,P_M> guessW(m_env, *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj), false); // useOdeWithDerivativeWrtTime = false
  guessW.compute(params,
                 true, // computeGradAlso
                 &referenceW,
                 0.);  // maximumTemp

  tmpSize = guessW.times().size();
  ofs << "\nguessWTemps = zeros(" << tmpSize
      << ",1);"
      << "\nguessWs = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nguessWTemps(" << i+1 << ",1) = " << guessW.temps()[i] << ";"
        << "\nguessWs("     << i+1 << ",1) = " << guessW.ws()   [i] << ";";
  }

  // Compute lambda(A_guess,E_guess)
#if 0
  uqTgaLambdaClass<P_V,P_M> guessLambda(m_env, *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj), false); // useOdeWithDerivativeWrtTime = false
  guessLambda.compute(params,
                      true, // computeGradAlso
                      guessW);

  tmpSize = guessLambda.times().size();
  ofs << "\nguessLambdaTemps = zeros(" << tmpSize
      << ",1);"
      << "\nguessLambdas = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nguessLambdaTemps(" << i+1 << ",1) = " << guessLambda.temps[i]  << ";"
        << "\nguessLambdas("     << i+1 << ",1) = " << guessLambda.lamdas[i] << ";";
  }
#endif
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

#ifdef QUESO_RUN_SIMPLIFIED_W_CASE
  return;
#endif

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

#ifdef QUESO_USE_NEW_W_CLASS
#else
  double deltaA = params[0] * 1.e-6;

  params[0] = 2.6000e+11-deltaA; // A
  params[1] = 2.0000e+05;        // E
  double valueAm = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodInfoVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueAm = " << valueAm << std::endl;

  params[0] = 2.6000e+11+deltaA; // A
  params[1] = 2.0000e+05;        // E
  double valueAp = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodInfoVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueAp = " << valueAp << std::endl;

  double deltaE = params[1] * 1.e-6;

  params[0] = 2.6000e+11;        // A
  params[1] = 2.0000e+05-deltaE; // E
  double valueEm = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodInfoVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueEm = " << valueEm << std::endl;

  params[0] = 2.6000e+11;        // A
  params[1] = 2.0000e+05+deltaE; // E
  double valueEp = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodInfoVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueEp = " << valueEp << std::endl;

  std::cout << "So, grad (finite diff) = " << (valueAp-valueAm)/2./deltaA << ", " << (valueEp-valueEm)/2./deltaE << std::endl;
#endif
  // Finish
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}
#if 0
template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runMinimization()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runMinimization()"
              << std::endl;
  }

  P_V params(m_paramSpace->zeroVector());
  params[0] = 2.6000e+11; // A
  params[1] = 2.0000e+05; // E
  P_V externalLagGrad(m_paramSpace->zeroVector());
  double externalLagGradNorm = externalLagGrad.norm2();

  unsigned int externalLoopId = 0;
  do {
    std::cout << "Beggining externalLoopId = " << externalLoopId
              << " with params = "             << params
              << std::endl;

    m_calLikelihoodInfoVector.resize(1,NULL);
    m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,"tga/scenario_5_K_min.dat");

    std::vector<double> misfitVarRatios(0);
    std::vector<double> allWTemps(0);
    std::vector<double> allWs(0);
    double tmpValue = 0.;
    tmpValue = tgaConstraintEquation<P_V,P_M>(params,
                                              *(m_calLikelihoodInfoVector[0]),
                                              false,
                                              &misfitVarRatios,
                                              &allWTemps,
                                              &allWs);
    //for (unsigned int i = 0; i < allWTemps.size(); ++i) {
    //  std::cout << allWTemps[i] << " " << allWs[i] << std::endl;
    //}

    std::vector<double> allLambdaTemps(0);
    std::vector<double> allLambdas(0);
    tgaAdjointEquation<P_V,P_M>(params,
                                *(m_calLikelihoodInfoVector[0]),
                                allWTemps[allWTemps.size()-1],
                                misfitVarRatios,
                                &allLambdaTemps,
                                &allLambdas);

    //for (unsigned int i = 0; i < allLambdaTemps.size(); ++i) {
    //  std::cout << allLambdaTemps[i] << " " << allLambdas[i] << std::endl;
    //}

    double extraTerm;
    tgaDesignEquation<P_V,P_M>(params,
                               *(m_calLikelihoodInfoVector[0]),
                               allWTemps,
                               allWs,
                               allLambdaTemps,
                               allLambdas,
                               externalLagGrad,
                               &extraTerm);
    std::cout << "externalLagGrad = " << externalLagGrad << std::endl;
    externalLagGradNorm = externalLagGrad.norm2();
    std::cout << "At externalLoopId = "     << externalLoopId
              << " with params = "          << params
              << ": misfit = "              << tmpValue
              << ", externalLagGradNorm = " << externalLagGradNorm
              << std::endl;

    if (1.e-7 < externalLagGradNorm) {
      P_V newtonLagGrad(externalLagGrad);
      double newtonLagGradNorm = newtonLagGrad.norm2();

      // Apply Newton method
      unsigned int newtonLoopId = 0;
      do {
        // Compute Newton step
        P_V paramsStep(m_paramSpace->zeroVector());
        double A = params[0];
        P_M mat(params,0.);
        mat(0,0)=0.;
        mat(0,1)=newtonLagGrad[1]/A;
        mat(1,0)=newtonLagGrad[1]/A;
        mat(1,1)=extraTerm;
        paramsStep *= 0.;
        mat.invertMultiply(-1.*newtonLagGrad,paramsStep);
        std::cout << "In newtonLoopId = " << newtonLoopId
                  << " with params = "    << params
                  << ": paramsStep = "    << paramsStep
                  << std::endl;

        // Perform line search
        unsigned int lineSearchLoopId = 0;
        double directionalDerivative = -newtonLagGrad.norm2Sq();
        double c = 0.1;
        double alpha = 2.0;
        double currentMeritFunction = .5*newtonLagGrad.norm2Sq();
        double newMeritFunction = 0.;
        double lineSearchThreshold = 0.;
        do {
          alpha *= .5;
          lineSearchThreshold = currentMeritFunction + alpha * c * directionalDerivative;
          params += alpha*paramsStep;
          tgaDesignEquation<P_V,P_M>(params,
                                     *(m_calLikelihoodInfoVector[0]),
                                     allWTemps,
                                     allWs,
                                     allLambdaTemps,
                                     allLambdas,
                                     newtonLagGrad,
                                     NULL);
          newtonLagGradNorm = newtonLagGrad.norm2();
          newMeritFunction = .5*newtonLagGrad.norm2Sq();
          std::cout << "In lineSearchLoopId = "   << lineSearchLoopId
                    << ", with params = "         << params
                    << ": newMeritFunction = "    << newMeritFunction
                    << ", lineSearchThreshold = " << lineSearchThreshold
                    << std::endl;
          lineSearchLoopId++;
        } while ((newMeritFunction > lineSearchThreshold) && (lineSearchLoopId < 5));
        std::cout << "In newtonLoopId = "                << newtonLoopId
                  << ", with params = "                  << params
                  << ": exited line search after "       <<  lineSearchLoopId
                  << " iterations; newtonLagGradNorm = " << newtonLagGradNorm
                  << std::endl;

        newtonLoopId++;
      } while ((1.e-7 < newtonLagGradNorm) && (newtonLoopId < 10));
    }

    std::cout << "Ending externalLoopId = " << externalLoopId
              << "\n"
              << std::endl;
    externalLoopId++;
  } while ((1.e-8 < externalLagGradNorm) && (externalLoopId < 10));

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runMinimization()"
              << std::endl;
  }

  return;
}
#endif
#endif // __UQ_TGA_VALIDATION_H__
