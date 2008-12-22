/* uq/examples/queso/pyramid/uqTgaOdes.h
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

#ifndef __UQ_TGA_ODES_H__
#define __UQ_TGA_ODES_H__

#include <uqTgaDataClass.h>

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

#endif // __UQ_TGA_ODES_H__
