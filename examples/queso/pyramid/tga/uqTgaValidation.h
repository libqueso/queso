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

#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472
#define GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO 1.0e-6
#define QUESO_USE_FABRICATED_MEASUREMENT
#undef  QUESO_USE_OBSERVED_MINUS_COMPUTED
#undef  QUESO_RUN_SIMPLIFIED_W_CASE

// The state dot ODE function
int tgaStateTimeDotOdeFunction(double time, const double w[], double f[], void *info)
{
  double* odeParameters = (double *)info;
  double A           = odeParameters[0];
  double E           = odeParameters[1];
  double beta        = odeParameters[2];
  double initialTemp = odeParameters[3];
  double temp        = initialTemp + beta*time;

  f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp)); // No division by 'beta' (CONVERSION TO TIME)

  return GSL_SUCCESS;
}

int tgaStateTempDotOdeFunction(double temp, const double w[], double f[], void *info)
{
  double* odeParameters = (double *)info;
  double A    = odeParameters[0];
  double E    = odeParameters[1];
  double beta = odeParameters[2];

  f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp))/beta;

  return GSL_SUCCESS;
}

// The Lagrange multiplier dot ODE function
typedef struct
{
  double A;
  double E;
  double beta;
  double initialTemp;
  double nextContributionTime;
  double nextContributionTemp;
  double nextContributionValue;
  const std::vector<double>* fabricatedTemps;
  const std::vector<double>* fabricatedWs;
  const std::vector<double>* simulatedWs;
} tgaLambdaDotOdeInfo_st;

int tgaLambdaTimeDotOdeFunction(double time, const double L[], double f[], void *info)
{
  tgaLambdaDotOdeInfo_st* odeInfo = (tgaLambdaDotOdeInfo_st *)info;
  double A           = odeInfo->A;
  double E           = odeInfo->E;
  double beta        = odeInfo->beta;
  double initialTemp = odeInfo->initialTemp;
  double temp        = initialTemp + beta*time;

  if (odeInfo->fabricatedTemps) {
    const std::vector<double>& fabricatedTemps = *(odeInfo->fabricatedTemps);
    const std::vector<double>& fabricatedWs    = *(odeInfo->fabricatedWs   );
    const std::vector<double>& simulatedWs     = *(odeInfo->simulatedWs    );

    double tempTilda = fabricatedTemps[0] + fabricatedTemps[fabricatedTemps.size()-1] - temp;
    f[0] = -A*L[0]*exp(-E/(R_CONSTANT*tempTilda));

    unsigned int tmpId = 0;
    double contribution = 0.;
    if (tempTilda > fabricatedTemps[0]) while(tmpId < fabricatedTemps.size()) {
      if (tempTilda < fabricatedTemps[tmpId]) {
        double ratio = (tempTilda - fabricatedTemps[tmpId-1])/(fabricatedTemps[tmpId] - fabricatedTemps[tmpId-1]);
        double fabricatedW = fabricatedWs[tmpId-1]+ratio*(fabricatedWs[tmpId]-fabricatedWs[tmpId-1]);
        double simulatedW  = simulatedWs [tmpId-1]+ratio*(simulatedWs [tmpId]-simulatedWs [tmpId-1]);
        contribution = 2.*(simulatedW - fabricatedW);
        break;
      }
      else if (tempTilda == fabricatedTemps[tmpId]) {
        contribution = 2.*(simulatedWs[tmpId] - fabricatedWs[tmpId]);
        break;
      }
      tmpId++;
    }
    UQ_FATAL_TEST_MACRO((tmpId == fabricatedTemps.size()),
                        0,
                        "tgaLambdaTempDotOdeFunction, in uqTgaValidation.h",
                        "tmpId got too large");
    //std::cout << "In tgaLambdaTempDotOdeFunction()"
    //          << ": temp = "                << temp
    //          << ", adding contribution = " << contribution
    //         << std::endl;
    f[0] -= contribution;
  }
  else if (temp == odeInfo->nextContributionTemp) {
    f[0] = A*L[0]*exp(-E/(R_CONSTANT*temp));
    std::cout << "In tgaLambdaTimeDotOdeFunction()"
              << ": time = "                         << time
              << ", nextContributionTime = "         << odeInfo->nextContributionTime
              << ", adding nextContributionValue = " << odeInfo->nextContributionValue
              << std::endl;
    f[0] += odeInfo->nextContributionValue;
  }
  //else {
  //  std::cout << "In tgaLambdaTimeDotOdeFunction()"
  //            << ": temp = "                 << temp
  //            << ", nextContributionTime = " << odeInfo->nextContributionTime
  //            << ", adding 0."
  //            << std::endl;
  //}
  // No division by 'beta' (CONVERSION TO TIME)

  return GSL_SUCCESS;
}

int tgaLambdaTempDotOdeFunction(double temp, const double L[], double f[], void *info)
{
  //std::cout << "Entering tgaLambdaTempDotOdeFunction()"
  //          << ": temp = " << temp
  //          << std::endl;
  tgaLambdaDotOdeInfo_st* odeInfo = (tgaLambdaDotOdeInfo_st *)info;
  double A    = odeInfo->A;
  double E    = odeInfo->E;
  double beta = odeInfo->beta;

  if (odeInfo->fabricatedTemps) {
    const std::vector<double>& fabricatedTemps = *(odeInfo->fabricatedTemps);
    const std::vector<double>& fabricatedWs    = *(odeInfo->fabricatedWs   );
    const std::vector<double>& simulatedWs     = *(odeInfo->simulatedWs    );

    double tempTilda = fabricatedTemps[0] + fabricatedTemps[fabricatedTemps.size()-1] - temp;
    f[0] = -A*L[0]*exp(-E/(R_CONSTANT*tempTilda));

    unsigned int tmpId = 0;
    double contribution = 0.;
    if (tempTilda > fabricatedTemps[0]) while(tmpId < fabricatedTemps.size()) {
      if (tempTilda < fabricatedTemps[tmpId]) {
        double ratio = (tempTilda - fabricatedTemps[tmpId-1])/(fabricatedTemps[tmpId] - fabricatedTemps[tmpId-1]);
        double fabricatedW = fabricatedWs[tmpId-1]+ratio*(fabricatedWs[tmpId]-fabricatedWs[tmpId-1]);
        double simulatedW  = simulatedWs [tmpId-1]+ratio*(simulatedWs [tmpId]-simulatedWs [tmpId-1]);
        contribution = 2.*(simulatedW - fabricatedW);
        break;
      }
      else if (tempTilda == fabricatedTemps[tmpId]) {
        contribution = 2.*(simulatedWs[tmpId] - fabricatedWs[tmpId]);
        break;
      }
      tmpId++;
    }
    UQ_FATAL_TEST_MACRO((tmpId == fabricatedTemps.size()),
                        0,
                        "tgaLambdaTempDotOdeFunction, in uqTgaValidation.h",
                        "tmpId got too large");
    //std::cout << "In tgaLambdaTempDotOdeFunction()"
    //          << ": temp = "                << temp
    //          << ", adding contribution = " << contribution
    //         << std::endl;
    f[0] -= contribution;
  }
  else if (temp == odeInfo->nextContributionTemp) {
    f[0] = -A*L[0]*exp(-E/(R_CONSTANT*temp));
    std::cout << "In tgaLambdaTempDotOdeFunction()"
              << ": temp = "                         << temp
              << ", nextContributionTemp = "         << odeInfo->nextContributionTemp
              << ", adding nextContributionValue = " << odeInfo->nextContributionValue
              << std::endl;
    f[0] += odeInfo->nextContributionValue;
  }
  //else {
  //  std::cout << "In tgaLambdaTempDotOdeFunction()"
  //            << ": temp = "                         << temp
  //            << ", nextContributionTemp = "         << odeInfo->nextContributionTemp
  //            << ", adding 0."
  //            << std::endl;
  //}
  f[0] /= beta;

  //std::cout << "Leaving tgaLambdaTempDotOdeFunction()"
  //          << ": temp = " << temp
  //          << std::endl;

  return GSL_SUCCESS;
}

//********************************************************
// Likelihood function object for both inverse problems of the validation cycle.
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************

// The (user defined) data class that carries the data needed by the (user defined) likelihood routine
template<class P_V, class P_M>
struct
tgaLikelihoodRoutine_DataClass
{
  tgaLikelihoodRoutine_DataClass(const uqBaseEnvironmentClass& env,
                                 const std::string&            inpName);
 ~tgaLikelihoodRoutine_DataClass();

  bool                m_useTimeAsDomainVariable;
  double              m_beta;
  double              m_initialTemp;
  std::vector<double> m_measuredTemps; // temperatures
  std::vector<double> m_measuredWs;    // relative masses
  std::vector<double> m_measurementVs; // variances
  std::vector<double> m_fabricatedTemps;
  std::vector<double> m_fabricatedWs;
  std::vector<double> m_simulatedWs;
};

template<class P_V, class P_M>
tgaLikelihoodRoutine_DataClass<P_V,P_M>::tgaLikelihoodRoutine_DataClass(
  const uqBaseEnvironmentClass& env,
  const std::string&            inpName)
  :
  m_useTimeAsDomainVariable(true), // IMPORTANT
  m_beta           (0.),
  m_initialTemp    (0.),
  m_measuredTemps  (0),
  m_measuredWs     (0),
  m_measurementVs  (0),
  m_fabricatedTemps(0),
  m_fabricatedWs   (0),
  m_simulatedWs    (0)
{
  // Read experimental data
  if (env.rank() == 0) {
    std::cout << "In tgaLikelihoodRoutine_DataClass(), reading file '"
              << inpName << "'\n"
              << std::endl;
  }

  // Open input file on experimental data
  FILE *inp;
  inp = fopen(inpName.c_str(),"r");

  // Read kinetic parameters and convert heating rate to K/s
  unsigned int numMeasurements;
  fscanf(inp,"%lf %lf %d",&m_beta,&m_initialTemp,&numMeasurements);
  m_beta /= 60.;
  m_measuredTemps.resize(numMeasurements,0.);
  m_measuredWs.resize   (numMeasurements,0.);
  m_measurementVs.resize(numMeasurements,0.);
#ifdef QUESO_USE_FABRICATED_MEASUREMENT
  m_useTimeAsDomainVariable = true; // IMPORTANT
  m_fabricatedTemps.resize(100000,0.);
  m_fabricatedWs.resize   (100000,0.);
  m_simulatedWs.resize    (100000,0.);
  double A = 2.6090e+11; // Reference _A
  double E = 1.9910e+05; // Reference _E
  double stateTimeDotOdeParameters[]={A,E,m_beta,m_initialTemp};
  double stateTempDotOdeParameters[]={A,E,m_beta};
    	
  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,1);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sysTime = {tgaStateTimeDotOdeFunction, NULL, 1, (void *)stateTimeDotOdeParameters};
        gsl_odeiv_system     sysTemp = {tgaStateTempDotOdeFunction, NULL, 1, (void *)stateTempDotOdeParameters}; 

  double currentTemp = m_initialTemp;
  double deltaTemp   = 1e-1;

  double currentTime = 0.;
  double deltaTime   = 1e-1;

  double currentW[1];
  currentW[0]=1.;

  unsigned int loopId = 0;
  m_fabricatedTemps[loopId] = currentTemp;
  m_fabricatedWs   [loopId] = currentW[0];
  while (currentW[0] > 0.) {
    int status = 0;
    if (m_useTimeAsDomainVariable) {
      deltaTime = .1;
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, currentTime+deltaTime, &deltaTime, currentW);
      currentTemp = m_initialTemp + m_beta*currentTime;
    }
    else {
      deltaTemp = .1;
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, currentTemp+deltaTemp, &deltaTemp, currentW);
    }
    if (currentW[0] < 1.e-6) currentW[0] = 0.; // IMPORTANT

    loopId++;
    m_fabricatedTemps[loopId] = currentTemp;
    m_fabricatedWs   [loopId] = currentW[0];
  }
  m_fabricatedTemps.resize(loopId+1);
  m_fabricatedWs.resize   (loopId+1);
  m_simulatedWs.resize    (loopId+1);

  std::cout << "Fabricated signal has " << m_fabricatedTemps.size() << " samples" << std::endl;

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);
#endif
  
  unsigned int whileSize = 0;
  double tmpTemp;
  double tmpW;
  double tmpV;
  while (fscanf(inp,"%lf %lf %lf",&tmpTemp,&tmpW,&tmpV) != EOF) {
    UQ_FATAL_TEST_MACRO((whileSize >= numMeasurements),
                        env.rank(),
                        "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
                        "input file 1 has too many measurements");
    m_measuredTemps[whileSize] = tmpTemp;
    m_measuredWs   [whileSize] = tmpW;
    m_measurementVs[whileSize] = tmpV;
    whileSize++;
  }
  UQ_FATAL_TEST_MACRO((whileSize != numMeasurements),
                      env.rank(),
                      "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
                      "input file 1 has a smaller number of measurements than expected");

  // Close input file on experimental data
  fclose(inp);
}

template<class P_V, class P_M>
tgaLikelihoodRoutine_DataClass<P_V,P_M>::~tgaLikelihoodRoutine_DataClass()
{
}

// The TGA constraint equation
template<class P_V,class P_M>
double
tgaConstraintEquation(
  const P_V&                                     paramValues,
  tgaLikelihoodRoutine_DataClass<P_V,P_M>& info, // NO CONST
  bool                                           justComputeMisfit,
  std::vector<double>*                           misfitVarRatios,
  std::vector<double>*                           allWTemps,
  std::vector<double>*                           allWs)
{
  double resultValue = 0.;

    double A = paramValues[0];
    double E = paramValues[1];
    bool   useTimeAsDomainVariable             = info.m_useTimeAsDomainVariable;
    double beta                                = info.m_beta;
    double initialTemp                         = info.m_initialTemp;
    const std::vector<double>& measuredTemps   = info.m_measuredTemps;
    const std::vector<double>& measuredWs      = info.m_measuredWs;
    const std::vector<double>& measurementVs   = info.m_measurementVs;
    const std::vector<double>& fabricatedTemps = info.m_fabricatedTemps;
    const std::vector<double>& fabricatedWs    = info.m_fabricatedWs;
          std::vector<double>& simulatedWs     = info.m_simulatedWs;

    std::cout << "In tgaConstraintEquation()"
              << ": fabricatedTemps.size() = " << fabricatedTemps.size()
              << ", A = " << A
              << ", E = " << E
              << std::endl;

    double stateTimeDotOdeParameters[]={A,E,beta,initialTemp};
    double stateTempDotOdeParameters[]={A,E,beta};
      	
    // Integration
    const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,1);
          gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
          gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sysTime = {tgaStateTimeDotOdeFunction, NULL, 1, (void *)stateTimeDotOdeParameters};
          gsl_odeiv_system     sysTemp = {tgaStateTempDotOdeFunction, NULL, 1, (void *)stateTempDotOdeParameters}; 

    double previousTemp = 0.;
    double currentTemp  = initialTemp;
    double deltaTemp    = 1e-3;
    double maximumTemp  = 1900.;//2.*measuredTemps[measuredTemps.size()-1];

    double currentTime  = 0.;
    double deltaTime    = 1e-3;
    double maximumTime  = (maximumTemp-initialTemp)/beta; // CONVERSION TO TIME

    double previousW[1];
    double currentW [1];
    previousW[0]=1.;
    currentW [0]=1.;

    if (misfitVarRatios) misfitVarRatios->resize(measuredWs.size(),0.);
    std::vector<double> computedWs  (measuredWs.size(),0.);
    std::vector<double> misfitValues(measuredWs.size(),0.);

    unsigned int loopId = 0;
    if (allWTemps) {
      allWTemps->resize(allWTemps->size()+1000,0.);
      allWs->resize(allWs->size()+1000,0.);
      (*allWTemps)[loopId] = currentTemp;
      (*allWs)[loopId] = currentW[0];
    }

    unsigned int fabricatedId = 0;
    simulatedWs[fabricatedId] = currentW[0];
    unsigned int misfitId = 0;
    bool continueOnWhile = true;
    if (fabricatedTemps.size() > 0) {
      continueOnWhile = (fabricatedId <= (fabricatedTemps.size()-2));
    }
    else if (justComputeMisfit) {
      continueOnWhile = (currentTemp < maximumTemp) && (misfitId < measuredWs.size());
    }
    else {
      continueOnWhile = (1.e-16 < currentW[0]);
    }
    while (continueOnWhile) {
      loopId++;
      int status = 0;
      if (fabricatedTemps.size() > 0) {
	//std::cout << "In tgaConstraintEquation(), fabricated case"
        //          << ": loopId = "       << loopId
        //          << ", fabricatedId = " << fabricatedId
        //          << ", currentTemp = "  << currentTemp
        //          << ", limitTemp = "    << fabricatedTemps[fabricatedId+1]
        //          << ", currentTime = "  << (currentTemp-initialTemp)/beta
        //          << ", currentW[0] = "  << currentW[0]
        //          << std::endl;
        status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, fabricatedTemps[fabricatedId+1], &deltaTemp, currentW);
      }
      else if (useTimeAsDomainVariable) {
	//std::cout << "currentTime = " << currentTime << std::endl;
        status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, maximumTime, &deltaTime, currentW);
        currentTemp = initialTemp + beta*currentTime;
      }
      else {
        status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, maximumTemp, &deltaTemp, currentW);
      }
      //std::cout << "In tgaConstraintEquation()"
      //          << ": currentTemp = " << currentTemp
      //          << std::endl;
      if (currentW[0] < 1.e-6) { // IMPORTANT
        //std::cout << "In tgaConstraintEquation()"
        //          << ": currentTemp = " << currentTemp
        //          << ", currentW[0] = " << currentW[0]
        //          << std::endl;
        currentW[0] = 0.;
      }
      UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                          paramValues.env().rank(),
                          "tgaConstraintEquation()",
                          "gsl_odeiv_evolve_apply() failed");
      if (allWTemps) {
        if (loopId >= allWTemps->size()) {
          allWTemps->resize(allWTemps->size()+1000,0.);
          allWs->resize(allWs->size()+1000,0.);
        }
        (*allWTemps)[loopId] = currentTemp;
	(*allWs)[loopId] = currentW[0];
      }
		
      if (fabricatedTemps.size() > 0) {
        if (currentTemp == fabricatedTemps[fabricatedId+1]) {
          fabricatedId++;
          simulatedWs[fabricatedId] = currentW[0];
#ifdef QUESO_USE_OBSERVED_MINUS_COMPUTED
          double tmpMisfitValue = fabricatedWs[fabricatedId] - simulatedWs[fabricatedId];
#else
          double tmpMisfitValue = simulatedWs[fabricatedId] - fabricatedWs[fabricatedId];
#endif
          double tmpToAdd = tmpMisfitValue*tmpMisfitValue*(currentTemp-previousTemp)/beta;
          //std::cout << "In tgaConstraintEquation()"
          //          << ": previousTemp = "   << previousTemp
          //          << ", currentTemp = "    << currentTemp
          //          << ", tmpMisfitValue = " << tmpMisfitValue
          //          << ", tmpToAdd = "       << tmpToAdd
          //          << std::endl;
          resultValue += tmpToAdd;
        }
      }
      else while ( (misfitId < measuredWs.size()) && (previousTemp <= measuredTemps[misfitId]) && (measuredTemps[misfitId] <= currentTemp) ) {
        computedWs[misfitId] = (measuredTemps[misfitId]-previousTemp)*(currentW[0]-previousW[0])/(currentTemp-previousTemp) + previousW[0];
#ifdef QUESO_USE_OBSERVED_MINUS_COMPUTED
        misfitValues[misfitId] = measuredWs[misfitId]-computedWs[misfitId];
#else
        misfitValues[misfitId] = computedWs[misfitId]-measuredWs[misfitId];
#endif
        if (misfitVarRatios) (*misfitVarRatios)[misfitId] = misfitValues[misfitId]/measurementVs[misfitId];
        resultValue += (misfitValues[misfitId]*misfitValues[misfitId])/measurementVs[misfitId];
        if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
          std::cout << "In tgaConstraintEquation()"
                    << ", misfitId = "     << misfitId
                    << ", measuredTemp = " << measuredTemps[misfitId]
                    << ", measuredW = "    << measuredWs[misfitId]
                    << ": computedW = "    << computedWs[misfitId]
                    << ", misfitValue = "  << misfitValues[misfitId]
                    << std::endl;
        }
        misfitId++;
      }
		
      previousTemp = currentTemp;
      previousW[0] = currentW[0];

      if (fabricatedTemps.size() > 0) {
        continueOnWhile = (fabricatedId <= (fabricatedTemps.size()-2));
      }
      else if (justComputeMisfit) {
        continueOnWhile = (currentTemp < maximumTemp) && (misfitId < measuredWs.size());
      }
      else {
        continueOnWhile = (1.e-16 < currentW[0]);
      }
    }
    if (allWTemps) {
      allWTemps->resize(loopId+1);
      allWs->resize(loopId+1);
    }
	
    if ((paramValues.env().verbosity() >= 0) && (paramValues.env().rank() == 0)) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",A);
      sprintf(stringE,"%12.6e",E);
      std::cout << "In tgaConstraintEquation()"
                << ", A = "                              << stringA
                << ", E = "                              << stringE
                << ", beta = "                           << beta
                << ": finished ode loop after "          << loopId
                << " iterations, with weigthedMisfit = " << resultValue
                << std::endl;
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);

  return resultValue;
}

// The TGA adjoint equation
template<class P_V,class P_M>
void
tgaAdjointEquation(
  const P_V&                                     paramValues,
  const tgaLikelihoodRoutine_DataClass<P_V,P_M>& info,
  double                                         maximumTemp,
  const std::vector<double>&                     misfitVarRatios,
  std::vector<double>*                           allLambdaTemps,
  std::vector<double>*                           allLambdas)
{
  bool useTimeAsDomainVariable = info.m_useTimeAsDomainVariable;
  const std::vector<double>& measuredTemps = info.m_measuredTemps;

  UQ_FATAL_TEST_MACRO((misfitVarRatios.size() != measuredTemps.size()),
                      paramValues.env().rank(),
                      "tgaAdjointEquation()",
                      "vectors 'misfitRatios' and 'measuredTemps' have different sizes");

  UQ_FATAL_TEST_MACRO((misfitVarRatios.size() == 0),
                      paramValues.env().rank(),
                      "tgaAdjointEquation()",
                      "vector 'misfitRatios' has size zero");

  tgaLambdaDotOdeInfo_st odeInfo;
  odeInfo.A           = paramValues[0];
  odeInfo.E           = paramValues[1];
  odeInfo.beta        = info.m_beta;
  odeInfo.initialTemp = info.m_initialTemp;
  odeInfo.nextContributionTime  = 0.;
  odeInfo.nextContributionTemp  = 0.;
  odeInfo.nextContributionValue = 0.;
  odeInfo.fabricatedTemps       = &(info.m_fabricatedTemps);
  odeInfo.fabricatedWs          = &(info.m_fabricatedWs);
  odeInfo.simulatedWs           = &(info.m_simulatedWs);

  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,1);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sysTime = {tgaLambdaTimeDotOdeFunction, NULL, 1, (void *)&odeInfo};
        gsl_odeiv_system     sysTemp = {tgaLambdaTempDotOdeFunction, NULL, 1, (void *)&odeInfo}; 

  double currentTemp = odeInfo.initialTemp;
  double deltaTemp   = 1e-3;

  double currentTime = 0.;
  double deltaTime   = 1e-3;
  double maximumTime = (maximumTemp-odeInfo.initialTemp)/odeInfo.beta; // CONVERSION

  double currentLambda[1];
  currentLambda[0]=0.;

  unsigned int loopId = 0;
  unsigned int nextContributionId = 0;
  if (allLambdaTemps) {
    allLambdaTemps->resize(allLambdaTemps->size()+1000,0.);
    allLambdas->resize(allLambdas->size()+1000,0.);
    (*allLambdaTemps)[loopId] = currentTemp;
    (*allLambdas)[loopId] = currentLambda[0];
  }

  std::cout << "In tgaAdjointEquation()"
            << ", maximumTime = " << maximumTime
            << ", currentTime = " << currentTime
            << ", maximumTemp = " << maximumTemp
            << ", currentTemp = " << currentTemp
            << std::endl;

  double evolveTimeLimit = maximumTime;
  double evolveTempLimit = maximumTemp;

  while (currentTemp < maximumTemp) {
    loopId++;
    int status = 0;
    if (useTimeAsDomainVariable) {
      if (nextContributionId < measuredTemps.size()) {
        evolveTimeLimit = (measuredTemps[nextContributionId]-odeInfo.initialTemp)/odeInfo.beta;
        if (currentTime == evolveTimeLimit) {
          odeInfo.nextContributionTime  = (measuredTemps[nextContributionId]-odeInfo.initialTemp)/odeInfo.beta; // CONVERSION TO TIME
          odeInfo.nextContributionValue = 2.*misfitVarRatios[nextContributionId];
          nextContributionId++;
          if (nextContributionId >= measuredTemps.size()) {
            evolveTimeLimit = maximumTime;
          }
          else {
            evolveTimeLimit = (measuredTemps[nextContributionId]-odeInfo.initialTemp)/odeInfo.beta; // CONVERSION TO TIME
          }
        }
      }
      else {
        odeInfo.nextContributionTime  = maximumTime;
        odeInfo.nextContributionValue = 0.;
      }
    }
    else {
      if (nextContributionId < measuredTemps.size()) {
        evolveTempLimit = measuredTemps[nextContributionId];
        if (currentTemp == evolveTempLimit) {
          odeInfo.nextContributionTemp  = measuredTemps[nextContributionId];
          odeInfo.nextContributionValue = 2.*misfitVarRatios[nextContributionId];
          nextContributionId++;
          if (nextContributionId >= measuredTemps.size()) {
            evolveTempLimit = maximumTemp;
          }
          else {
            evolveTempLimit = measuredTemps[nextContributionId];
          }
        }
      }
      else {
        odeInfo.nextContributionTemp  = maximumTemp;
        odeInfo.nextContributionValue = 0.;
      }
    }
    //std::cout << "In tgaAdjointEquation(), before evolve()"
    //          << ": loopId = "          << loopId
    //          << ", evolveTempLimit = " << evolveTempLimit
    //          << ", currentTemp = "     << currentTemp
    //          << ", currentLambda = "   << currentLambda[0]
    //          << std::endl;
    if (useTimeAsDomainVariable) {
      //std::cout << "currentTime = " << currentTime << std::endl;
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, evolveTimeLimit, &deltaTime, currentLambda);
      currentTemp = odeInfo.initialTemp + odeInfo.beta*currentTime;
    }
    else {
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, evolveTempLimit, &deltaTemp, currentLambda);
    }
    //std::cout << "In tgaAdjointEquation(), after evolve()"
    //          << ": loopId = "          << loopId
    //          << ", evolveTempLimit = " << evolveTempLimit
    //          << ", currentTemp = "     << currentTemp
    //          << ", currentLambda = "   << currentLambda[0]
    //          << std::endl;
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        paramValues.env().rank(),
                        "tgaAdjointEquation()",
                        "gsl_odeiv_evolve_apply() failed");
    if (allLambdaTemps) {
      if (loopId >= allLambdaTemps->size()) {
        allLambdaTemps->resize(allLambdaTemps->size()+1000,0.);
        allLambdas->resize(allLambdas->size()+1000,0.);
      }
      (*allLambdaTemps)[loopId] = currentTemp;
      (*allLambdas)[loopId] = currentLambda[0];
    }
  }
  if (allLambdaTemps) {
    allLambdaTemps->resize(loopId+1);
    allLambdas->resize(loopId+1);
  }

  unsigned int pos1 = 0;
  unsigned int pos2 = allLambdas->size()-1;
  double minTemp = (*allLambdaTemps)[0];
  double maxTemp = (*allLambdaTemps)[allLambdaTemps->size()-1];
  while (pos2 >= pos1) {
    double tmpTemp = (*allLambdaTemps)[pos1];
    (*allLambdaTemps)[pos1] = minTemp + maxTemp - (*allLambdaTemps)[pos2];
    (*allLambdaTemps)[pos2] = minTemp + maxTemp - tmpTemp;

    double tmpLambda = (*allLambdas)[pos1];
    (*allLambdas)[pos1] = (*allLambdas)[pos2];
    (*allLambdas)[pos2] = tmpLambda;

    pos1++;
    pos2--;
  }
	
  if ((paramValues.env().verbosity() >= 0) && (paramValues.env().rank() == 0)) {
    char stringA[64];
    char stringE[64];
    sprintf(stringA,"%12.6e",odeInfo.A);
    sprintf(stringE,"%12.6e",odeInfo.E);
    std::cout << "In tgaAdjointEquation()"
              << ", A = "                     << stringA
              << ", E = "                     << stringE
              << ", beta = "                  << odeInfo.beta
              << ": finished ode loop after " << loopId
              << " iterations"
              << std::endl;
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

// The TGA design equation
template<class P_V,class P_M>
void
tgaDesignEquation(
  const P_V&                                     paramValues,
  const tgaLikelihoodRoutine_DataClass<P_V,P_M>& info,
  std::vector<double>&                           allWTemps,
  std::vector<double>&                           allWs,
  std::vector<double>&                           allLambdaTemps,
  std::vector<double>&                           allLambdas,
  P_V&                                           LagrangianGradWrtParams,
  double*                                        extraTerm)
{
  std::cout << "Entering tgaDesignEquation()"
            << ", allWTemps[0] = " << allWTemps[0]
            << ", allLambdaTemps[0] = " << allLambdaTemps[0]
            << std::endl;
#if 0
  UQ_FATAL_TEST_MACRO((allWTemps[0] != allLambdaTemps[0]),
                      paramValues.env().rank(),
                      "tgaDesignEquation()",
                      "allWTemps[0] and allLambdaTemps[0] are different");

  UQ_FATAL_TEST_MACRO((allWTemps[allWTemps.size()-1] != allLambdaTemps[allLambdaTemps.size()-1]),
                      paramValues.env().rank(),
                      "tgaDesignEquation()",
                      "allWTemps[last] and allLambdaTemps[last] are different");
#endif
  double A = paramValues[0];
  double E = paramValues[1];

  unsigned int currentWIntervalId      = 0;
  unsigned int currentLambdaIntervalId = 0;

  unsigned int numIntervals = 1000;
  double tempIntervalSize = (allWTemps[allWTemps.size()-1]-allWTemps[0])/((double)numIntervals);
  double firstTemp = allWTemps[0]+.5*tempIntervalSize;

  std::cout << "In tgaAdjointEquation()"
            << ": beginning integration loop on temperature interval "
            << "["          << allWTemps[0]
            << ", "         << allWTemps[allWTemps.size()-1]
            << "], with = " << numIntervals
            << " subintervals"
            << "; allWTemps.size() = "      << allWTemps.size()
            << ", allLambdaTemps.size() = " << allLambdaTemps.size()
            << std::endl;

  LagrangianGradWrtParams *= 0.;
  if (extraTerm) *extraTerm = 0.;

  for (unsigned int i = 0; i < numIntervals; ++i) {
    double temp = firstTemp + ((double) i)*tempIntervalSize;

    while (allWTemps[currentWIntervalId+1] <= temp) {
      currentWIntervalId++;
    }
    if ((allWTemps[currentWIntervalId] <= temp ) &&
        (temp < allWTemps[currentWIntervalId+1])) {
      // Ok
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          paramValues.env().rank(),
                          "tgaDesignEquation()",
                          "invalid situation with allWTemps");
    }
    double wTempDelta = allWTemps[currentWIntervalId+1] - allWTemps[currentWIntervalId];
    double wTempRatio = (temp - allWTemps[currentWIntervalId])/wTempDelta;
    double wDelta = allWs[currentWIntervalId+1] - allWs[currentWIntervalId];
    double w      = allWs[currentWIntervalId] + wTempRatio * wDelta;

    while (allLambdaTemps[currentLambdaIntervalId+1] < temp) {
      currentLambdaIntervalId++;
    }
    if ((allLambdaTemps[currentLambdaIntervalId] <= temp ) &&
        (temp < allLambdaTemps[currentLambdaIntervalId+1])) {
      // Ok
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          paramValues.env().rank(),
                          "tgaDesignEquation()",
                          "invalid situation with allLambdaTemps");
    }
    double lambdaTempDelta = allLambdaTemps[currentLambdaIntervalId+1] - allLambdaTemps[currentLambdaIntervalId];
    double lambdaTempRatio = (temp - allLambdaTemps[currentLambdaIntervalId])/lambdaTempDelta;
    double lambdaDelta = allLambdas[currentLambdaIntervalId+1] - allLambdas[currentLambdaIntervalId];
    double lambda      = allLambdas[currentLambdaIntervalId] + lambdaTempRatio * lambdaDelta;

    LagrangianGradWrtParams[0] +=                          lambda * w * exp(-E/(R_CONSTANT*temp));
    LagrangianGradWrtParams[1] += -(A/(R_CONSTANT*temp)) * lambda * w * exp(-E/(R_CONSTANT*temp));
    if (extraTerm) *extraTerm  +=  (A/(R_CONSTANT*temp*R_CONSTANT*temp)) * lambda * w * exp(-E/(R_CONSTANT*temp));
  }

  std::cout << "In tgaAdjointEquation()"
            << ": finsihed integration loop with"
            << " currentWIntervalId = " << currentWIntervalId
            << ", currentLambdaIntervalId = " << currentLambdaIntervalId
            << std::endl;

#if 0
  UQ_FATAL_TEST_MACRO((currentWIntervalId != (allWTemps.size()-2)),
                      paramValues.env().rank(),
                      "tgaDesignEquation()",
                      "currentWIntervalId finished with an invalid value");

  UQ_FATAL_TEST_MACRO((currentLambdaIntervalId != (allLambdaTemps.size()-2)),
                      paramValues.env().rank(),
                      "tgaDesignEquation()",
                      "currentLambdaIntervalId finished with an invalid value");
#endif

  LagrangianGradWrtParams *= (tempIntervalSize/info.m_beta); // Division by beta due to the integration w.r.t. temperature

  return;
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
tgaLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering tgaLikelihoodRoutine()..." << std::endl;
  }

  double resultValue = 0.;
  const std::vector<tgaLikelihoodRoutine_DataClass<P_V,P_M> *>& info = *((const std::vector<tgaLikelihoodRoutine_DataClass<P_V,P_M> *> *) functionDataPtr);

  for (unsigned int i = 0; i < info.size(); ++i) {
    // Compute likelihood for scenario
    resultValue += tgaConstraintEquation(paramValues,*(info[i]),true,NULL,NULL,NULL);
  }

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving tgaLikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}

// The (user defined) grad likelihood routine
template<class P_V,class P_M>
void
tgaLikelihoodGradRoutine(const P_V& paramValues, const void* functionDataPtr, P_V& gradVector)
{
  gradVector *= 0.;

  return;
}

// The (user defined) grad likelihood routine
template<class P_V,class P_M>
void
tgaLikelihoodHessianRoutine(const P_V& paramValues, const void* functionDataPtr, P_M& hessianMatrix)
{
  hessianMatrix *= 0.;

  return;
}

//********************************************************
// QoI function object for both forward problems of the validation cycle.
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data class that carries the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
tgaQoiRoutine_DataClass
{
  bool   m_useTimeAsDomainVariable;
  double m_beta;
  double m_initialTemp;
  double m_criticalW;
  double m_criticalTime;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void tgaQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  double A            = paramValues[0];
  double E            = paramValues[1];

  const tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>& info = *((const tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr);
  bool   useTimeAsDomainVariable = info.m_useTimeAsDomainVariable;
  double beta                    = info.m_beta;
  double initialTemp             = info.m_initialTemp;
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

  void run            ();
  void runGradTest    ();
  void runMinimization();

private:
  void  runCalibrationStage();
  void  runValidationStage ();
  void  runComparisonStage ();

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

  double                                                 m_predUseTimeAsDomainVariable;
  double                                                 m_predBeta;
  double                                                 m_predInitialTemp;
  double                                                 m_predCriticalW;
  double                                                 m_predCriticalTime;

  uqBaseVectorRVClass<P_V,P_M>*                          m_calPriorRv;
  std::vector<tgaLikelihoodRoutine_DataClass<P_V,P_M>* > m_calLikelihoodRoutine_DataVector;
  uqBaseScalarFunctionClass<P_V,P_M>*                    m_calLikelihoodFunctionObj;
  tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>*              m_calQoiRoutine_Data;

  std::vector<tgaLikelihoodRoutine_DataClass<P_V,P_M>* > m_valLikelihoodRoutine_DataVector;
  uqBaseScalarFunctionClass<P_V,P_M>*                    m_valLikelihoodFunctionObj;
  tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>*              m_valQoiRoutine_Data;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::uqTgaValidationClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  uqModelValidationClass<P_V,P_M,Q_V,Q_M>(env,prefix),
  m_paramsTable                    (NULL),
  m_paramNames                     (NULL),
  m_paramMinValues                 (NULL),
  m_paramMaxValues                 (NULL),
  m_paramInitialValues             (NULL),
  m_paramSpace                     (NULL),
  m_paramDomain                    (NULL),
  m_qoisTable                      (NULL),
  m_qoiNames                       (NULL),
  m_qoiSpace                       (NULL),
  m_calPriorRv                     (NULL),
  m_calLikelihoodRoutine_DataVector(0),
  m_calLikelihoodFunctionObj       (NULL),
  m_calQoiRoutine_Data             (NULL),
  m_valLikelihoodRoutine_DataVector(0),
  m_valLikelihoodFunctionObj       (NULL),
  m_valQoiRoutine_Data             (NULL)
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

  m_predUseTimeAsDomainVariable = false; // IMPORTANT
  m_predBeta                    = 250.;
  m_predInitialTemp             = 0.1;
  m_predCriticalW               = 0.;
  m_predCriticalTime            = 3.9;

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

  if (m_valQoiRoutine_Data)        delete m_valQoiRoutine_Data;
  if (m_valLikelihoodFunctionObj)  delete m_valLikelihoodFunctionObj;
  for (unsigned int i = 0; i < m_valLikelihoodRoutine_DataVector.size(); ++i) {
    if (m_valLikelihoodRoutine_DataVector[i]) delete m_valLikelihoodRoutine_DataVector[i];
  }
  if (m_calQoiRoutine_Data)        delete m_calQoiRoutine_Data;
  if (m_calLikelihoodFunctionObj)  delete m_calLikelihoodFunctionObj;
  for (unsigned int i = 0; i < m_calLikelihoodRoutine_DataVector.size(); ++i) {
    if (m_calLikelihoodRoutine_DataVector[i]) delete m_calLikelihoodRoutine_DataVector[i];
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

  m_calLikelihoodRoutine_DataVector.resize(3,NULL);
  m_calLikelihoodRoutine_DataVector[0] = new tgaLikelihoodRoutine_DataClass<P_V,P_M>(m_env,"tga/scenario_5_K_min.dat");
  m_calLikelihoodRoutine_DataVector[1] = new tgaLikelihoodRoutine_DataClass<P_V,P_M>(m_env,"tga/scenario_25_K_min.dat");
  m_calLikelihoodRoutine_DataVector[2] = new tgaLikelihoodRoutine_DataClass<P_V,P_M>(m_env,"tga/scenario_50_K_min.dat");

  m_calLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("cal_like_",
                                                                         *m_paramDomain,
                                                                         tgaLikelihoodRoutine<P_V,P_M>,
                                                                         tgaLikelihoodGradRoutine<P_V,P_M>,
                                                                         tgaLikelihoodHessianRoutine<P_V,P_M>,
                                                                         (void *) &m_calLikelihoodRoutine_DataVector,
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
  m_calQoiRoutine_Data = new tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>();
  m_calQoiRoutine_Data->m_useTimeAsDomainVariable = m_predUseTimeAsDomainVariable;
  m_calQoiRoutine_Data->m_beta                    = m_predBeta;
  m_calQoiRoutine_Data->m_initialTemp             = m_predInitialTemp;
  m_calQoiRoutine_Data->m_criticalW               = m_predCriticalW;
  m_calQoiRoutine_Data->m_criticalTime            = m_predCriticalTime;

  m_cycle->setCalFP(tgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_calQoiRoutine_Data);

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
  m_valLikelihoodRoutine_DataVector.resize(1,NULL);
  m_valLikelihoodRoutine_DataVector[0] = new tgaLikelihoodRoutine_DataClass<P_V,P_M> (m_env,"tga/scenario_100_K_min.dat");

  m_valLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("val_like_",
                                                                         *m_paramDomain,
                                                                         tgaLikelihoodRoutine<P_V,P_M>,
                                                                         tgaLikelihoodGradRoutine<P_V,P_M>,
                                                                         tgaLikelihoodHessianRoutine<P_V,P_M>,
                                                                         (void *) &m_valLikelihoodRoutine_DataVector,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setValIP(*m_valLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* valProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().postRv().realizer().imageVarVector(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                                   m_cycle->calIP().postRv().realizer().imageExpVector()); // Use these values as the initial values
  m_cycle->valIP().solveWithBayesMarkovChain(m_cycle->calIP().postRv().realizer().imageExpVector(),
                                             valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Deal with forward problem
  m_valQoiRoutine_Data = new tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>();
  m_valQoiRoutine_Data->m_useTimeAsDomainVariable = m_predUseTimeAsDomainVariable;
  m_valQoiRoutine_Data->m_beta                    = m_predBeta;
  m_valQoiRoutine_Data->m_initialTemp             = m_predInitialTemp;
  m_valQoiRoutine_Data->m_criticalW               = m_predCriticalW;
  m_valQoiRoutine_Data->m_criticalTime            = m_predCriticalTime;

  m_cycle->setValFP(tgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_valQoiRoutine_Data);

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

  m_calLikelihoodRoutine_DataVector.resize(1,NULL);
  m_calLikelihoodRoutine_DataVector[0] = new tgaLikelihoodRoutine_DataClass<P_V,P_M>(m_env,"tga/scenario_5_K_min.dat");
  P_V params(m_paramSpace->zeroVector());

  // Open file
  std::ofstream ofs("nada1.m", std::ofstream::out | std::ofstream::trunc);
  unsigned int tmpSize = m_calLikelihoodRoutine_DataVector[0]->m_fabricatedTemps.size();
  ofs << "fabricatedTemps = zeros(" << tmpSize
      << ",1);"
      << "\nfabricatedWs = zeros(" << tmpSize
      << ",1);";
  for (unsigned int i = 0; i < tmpSize; ++i) {
    ofs << "\nfabricatedTemps(" << i+1 << ",1) = " << m_calLikelihoodRoutine_DataVector[0]->m_fabricatedTemps[i] << ";"
        << "\nfabricatedWs("    << i+1 << ",1) = " << m_calLikelihoodRoutine_DataVector[0]->m_fabricatedWs[i]    << ";";
  }

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": compute misfit gradient w.r.t. parameters A and E, using adjoint..."
              << std::endl;
  }

  params[0] = 2.6000e+11; // Initial guess _A
  params[1] = 2.0000e+05; // Initial guess _E

  std::vector<double> misfitVarRatios(0);
  std::vector<double> allWTemps(0);
  std::vector<double> allWs(0);
  double valueCenter = 0.;
  valueCenter = tgaConstraintEquation<P_V,P_M>(params,
                                               *(m_calLikelihoodRoutine_DataVector[0]),
                                               false,
                                               &misfitVarRatios,
                                               &allWTemps,
                                               &allWs);
  std::cout << "In runGradTest()"
            << ": valueCenter = " << valueCenter
            << std::endl;

  ofs << "\nallWTemps = zeros(" << allWTemps.size()
      << ",1);"
      << "\nallWs = zeros(" << allWTemps.size()
      << ",1);";
  for (unsigned int i = 0; i < allWTemps.size(); ++i) {
    ofs << "\nallWTemps(" << i+1 << ",1) = " << allWTemps[i] << ";"
        << "\nallWs("     << i+1 << ",1) = " << allWs[i]     << ";";
  }

  std::vector<double> allLambdaTemps(0);
  std::vector<double> allLambdas(0);
  tgaAdjointEquation<P_V,P_M>(params,
                              *(m_calLikelihoodRoutine_DataVector[0]),
                              allWTemps[allWTemps.size()-1],
                              misfitVarRatios,
                              &allLambdaTemps,
                              &allLambdas);

  ofs << "\nallLambdaTemps = zeros(" << allLambdaTemps.size()
      << ",1);"
      << "\nallLambdas = zeros(" << allLambdaTemps.size()
      << ",1);";
  for (unsigned int i = 0; i < allLambdaTemps.size(); ++i) {
    ofs << "\nallLambdaTemps(" << i+1 << ",1) = " << allLambdaTemps[i] << ";"
        << "\nallLambdas("     << i+1 << ",1) = " << allLambdas[i]     << ";";
  }

#ifdef QUESO_RUN_SIMPLIFIED_W_CASE
  return;
#endif

  P_V LagrangianGradWrtParams(m_paramSpace->zeroVector());
  tgaDesignEquation<P_V,P_M>(params,
                             *(m_calLikelihoodRoutine_DataVector[0]),
                             allWTemps,
                             allWs,
                             allLambdaTemps,
                             allLambdas,
			     LagrangianGradWrtParams,
                             NULL);

  std::cout << "lagGrad = " << LagrangianGradWrtParams << std::endl;

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": computing misfit gradient w.r.t. parameters A and E, using finite differences..."
              << std::endl;
  }

  double deltaA = params[0] * 1.e-6;

  params[0] = 2.6000e+11-deltaA; // A
  params[1] = 2.0000e+05;        // E
  double valueAm = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodRoutine_DataVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueAm = " << valueAm << std::endl;

  params[0] = 2.6000e+11+deltaA; // A
  params[1] = 2.0000e+05;        // E
  double valueAp = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodRoutine_DataVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueAp = " << valueAp << std::endl;

  double deltaE = params[1] * 1.e-6;

  params[0] = 2.6000e+11;        // A
  params[1] = 2.0000e+05-deltaE; // E
  double valueEm = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodRoutine_DataVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueEm = " << valueEm << std::endl;

  params[0] = 2.6000e+11;        // A
  params[1] = 2.0000e+05+deltaE; // E
  double valueEp = tgaConstraintEquation<P_V,P_M>(params,
                                                  *(m_calLikelihoodRoutine_DataVector[0]),
                                                  true,
                                                  NULL,
                                                  NULL,
                                                  NULL);
  std::cout << "valueEp = " << valueEp << std::endl;

  std::cout << "So, grad (finite diff) = " << (valueAp-valueAm)/2./deltaA << ", " << (valueEp-valueEm)/2./deltaE << std::endl;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}

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

    m_calLikelihoodRoutine_DataVector.resize(1,NULL);
    m_calLikelihoodRoutine_DataVector[0] = new tgaLikelihoodRoutine_DataClass<P_V,P_M>(m_env,"tga/scenario_5_K_min.dat");

    std::vector<double> misfitVarRatios(0);
    std::vector<double> allWTemps(0);
    std::vector<double> allWs(0);
    double tmpValue = 0.;
    tmpValue = tgaConstraintEquation<P_V,P_M>(params,
                                              *(m_calLikelihoodRoutine_DataVector[0]),
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
                                *(m_calLikelihoodRoutine_DataVector[0]),
                                allWTemps[allWTemps.size()-1],
                                misfitVarRatios,
                                &allLambdaTemps,
                                &allLambdas);

    //for (unsigned int i = 0; i < allLambdaTemps.size(); ++i) {
    //  std::cout << allLambdaTemps[i] << " " << allLambdas[i] << std::endl;
    //}

    double extraTerm;
    tgaDesignEquation<P_V,P_M>(params,
                               *(m_calLikelihoodRoutine_DataVector[0]),
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
                                     *(m_calLikelihoodRoutine_DataVector[0]),
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
#endif // __UQ_TGA_VALIDATION_H__
