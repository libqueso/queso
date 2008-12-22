/* uq/examples/queso/pyramid/uqTgaDataClass.h
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

#ifndef __UQ_TGA_DATA_CLASS_H__
#define __UQ_TGA_DATA_CLASS_H__

#include <uqTgaDefines.h>
#include <uqEnvironment.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

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

//********************************************************
// Likelihood data object for both inverse problems of the validation cycle.
// A likelihood function object is provided by user and is called by the UQ library.
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

#endif // __UQ_TGA_VALIDATION_H__
