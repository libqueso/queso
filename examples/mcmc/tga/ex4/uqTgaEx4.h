/* uq/examples/mcmc/tga/uqTgaEx4.h
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_EX4_H__
#define __UQ_TGA_EX4_H__

//#include <uqStateSpace.h>
#include <uqValidProblem.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472

// The ODE (state dot) function
int func(double t, const double Mass[], double f[], void *info)
{
  double* params = (double *)info;
  double A    = params[0];
  double E    = params[1];
  double beta = params[2];

  f[0] = -A*Mass[0]*exp(-E/(R_CONSTANT*t))/beta;

  return GSL_SUCCESS;
}

//********************************************************
// Likelihood function object for the first validation problem stage (with prefix "s1_").
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************

// The (user defined) data type for the data needed by the (user defined) likelihood routine
template<class P_V, class P_M>
struct
s1_likelihoodRoutine_DataType
{
  double               beta1;
  double               variance1;
  std::vector<double>* Te1; // temperatures
  std::vector<double>* Me1; // relative masses
};

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
s1_likelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  double A                       = paramValues[0];
  double E                       = paramValues[1];
  double beta1                   =  ((s1_likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta1;
  double variance1               =  ((s1_likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->variance1;
  const std::vector<double>& Te1 = *((s1_likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Te1;
  const std::vector<double>& Me1 = *((s1_likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Me1;
  std::vector<double> Mt1(Me1.size(),0.);

  double params[]={A,E,beta1};
      	
  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
        gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
        gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sys = {func, NULL, 1, (void *)params}; 
	
  double t = 0.1, t1 = 800.;
  double h = 1e-3;
  double Mass[1];
  Mass[0]=1.;
  
  unsigned int i = 0;
  double t_old = 0.;
  double M_old[1];
  M_old[0]=1.;
	
  double misfit1=0.;
  //unsigned int loopSize = 0;
  while (t < t1) {
    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, Mass);
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        paramValues.env().rank(),
                        "s1_likelihoodRoutine()",
                        "gsl_odeiv_evolve_apply() failed");
    //printf("t = %6.1lf, mass = %10.4lf\n",t,Mass[0]);
    //loopSize++;
		
    while ( (t_old <= Te1[i]) && (Te1[i] <= t) ) {
      Mt1[i] = (Te1[i]-t_old)*(Mass[0]-M_old[0])/(t-t_old) + M_old[0];
      misfit1 += (Me1[i]-Mt1[i])*(Me1[i]-Mt1[i]);
      //printf("%i %lf %lf %lf %lf\n",i,Te1[i],Me1[i],Mt1[i],misfit1);
      i++;
    }
		
    t_old=t;
    M_old[0]=Mass[0];
  }
  double resultValue = misfit1/variance1;
	
  //printf("loopSize = %d\n",loopSize);
  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    printf("In s1_likelihoodRoutine(), A = %g, E = %g, beta1 = %.3lf: misfit1 = %lf, likelihood1 = %lf.\n",A,E,beta1,misfit1,resultValue);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return resultValue;
}

//********************************************************
// QoI function object for the first validation problem stage (with prefix "s1_").
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data type for the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
s1_qoiRoutine_DataType
{
  double beta1;
  double criticalMass1;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void s1_qoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  double A             = paramValues[0];
  double E             = paramValues[1];
  double beta1         = ((s1_qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->beta1;
  double criticalMass1 = ((s1_qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->criticalMass1;

  double params[]={A,E,beta1};
      	
  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
        gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
        gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sys = {func, NULL, 1, (void *)params}; 
	
  double t = 0.1, t1 = 800.;
  double h = 1e-3;
  double Mass[1];
  Mass[0]=1.;
  
  double t_old = 0.;
  double M_old[1];
  M_old[0]=1.;
	
  double crossingTemperature = 0.;
  //unsigned int loopSize = 0;
  while ((t       < t1           ) &&
         (Mass[0] > criticalMass1)) {
    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, Mass);
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        paramValues.env().rank(),
                        "s1_qoiRoutine()",
                        "gsl_odeiv_evolve_apply() failed");
    //printf("t = %6.1lf, mass = %10.4lf\n",t,Mass[0]);
    //loopSize++;

    if (Mass[0] <= criticalMass1) {
      crossingTemperature = t_old + (t - t_old) * (M_old[0]-criticalMass1)/(M_old[0]-Mass[0]);
    }
		
    t_old=t;
    M_old[0]=Mass[0];
  }

  qoiValues[0] = crossingTemperature/beta1;
	
  //printf("loopSize = %d\n",loopSize);
  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    printf("In s1_qoiRoutine(), A = %g, E = %g, beta1 = %.3lf, criticalMass1 = %.3lf: qoi = %lf.\n",A,E,beta1,criticalMass1,qoiValues[0]);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

//********************************************************
// The validation problem driving routine "uqAppl()": called by main()
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqTgaEx4' example\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //*****************************************************
  // Step 1 of 4: Instantiate the validation problem
  //*****************************************************
  uqValidProblemClass<P_V,P_M,Q_V,Q_M> validProblem(env,
                                                    "");

  //******************************************************
  // Step 2 of 4: Set the first validation problem stage
  // There are 6 substeps: 2.1, 2.2, 2.3, 2.4, 2.5 and 2.6
  //******************************************************

  //*****************************************************
  // Substep 2.0: Code specific to this TGA example
  //*****************************************************
  // Open input file on experimental data
  FILE *inp;
  inp = fopen("global4.dat","r");

  // Read kinetic parameters and convert heating rate to K/s
  double tmpA, tmpE, beta1, variance1, criticalMass1;
  fscanf(inp,"%lf %lf %lf %lf %lf",&tmpA,&tmpE,&beta1,&variance1,&criticalMass1);
  beta1 /= 60.;
  
  // Read experimental data
  std::vector<double> Te1(11,0.);
  std::vector<double> Me1(11,0.);

  unsigned int numObservations = 0;
  double tmpTe;
  double tmpMe;
  while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
    UQ_FATAL_TEST_MACRO((numObservations >= Te1.size()),
                        env.rank(),
                        "uqAppl(), in uqTgaEx.h",
                        "input file has too many observations");
    Te1[numObservations] = tmpTe;
    Me1[numObservations] = tmpMe;
    numObservations++;
  }
  UQ_FATAL_TEST_MACRO((numObservations != Te1.size()),
                      env.rank(),
                      "uqAppl(), in uqTgaEx.h",
                      "input file has a smaller number of observations than expected");

  // Close input file on experimental data
  fclose(inp);

  //******************************************************
  // Substep 2.1: Define the prior probability density function object: -2*ln[prior]
  //******************************************************

  // Use the default prior density, provided by the UQ library
  uqProbDensity_BaseClass<P_V,P_M>* s1_priorParamDensityObj = NULL;

  //******************************************************
  // Substep 2.2: Define the likelihood function object: -2*ln[likelihood]
  //******************************************************

  s1_likelihoodRoutine_DataType<P_V,P_M> s1_likelihoodRoutine_Data;
  s1_likelihoodRoutine_Data.beta1     = beta1;
  s1_likelihoodRoutine_Data.variance1 = variance1;
  s1_likelihoodRoutine_Data.Te1       = &Te1; // temperatures
  s1_likelihoodRoutine_Data.Me1       = &Me1; // relative masses
  uqRoutineProbDensity_Class<P_V,P_M> s1_likelihoodFunctionObj(s1_likelihoodRoutine<P_V,P_M>,
                                                               (void *) &s1_likelihoodRoutine_Data,
                                                               true); // the routine computes [-2.*ln(Likelihood)]

  //******************************************************
  // Substep 2.3: Define the proposal density and proposal generator
  //******************************************************

  // Use the default gaussian proposal and default gaussian generator with the default covariance matrix, all provided by the UQ library
  P_M*                                    s1_proposalCovMatrix    = NULL;
  uqProposalDensity_BaseClass<P_V,P_M>*   s1_proposalDensityObj   = NULL;
  uqProposalGenerator_BaseClass<P_V,P_M>* s1_proposalGeneratorObj = NULL;

  //uqParamSpaceClass<P_V,P_M> s1_paramSpace(env,"val_stage_0_");
  //P_M* s1_proposalCovMatrix = s1_paramSpace.newMatrix();
  //(*s1_proposalCovMatrix)(0,0) = 1.65122e+10;
  //(*s1_proposalCovMatrix)(0,1) = 0.;
  //(*s1_proposalCovMatrix)(1,0) = 0.;
  //(*s1_proposalCovMatrix)(1,1) = 9.70225e+04;

  //******************************************************
  // Substep 2.4: Define the propagation input parameters (density and generator)
  //******************************************************

  // Let the 'uqValidProblem' class internally use the parameter generator obtained with the calibration.
  // The 'uqTgaEx2.inp' file should set 'val_stage_0_pro_inputStageId = 0'.
  uqRandomVariableClass<P_V,P_M>* s1_propagParamRV = NULL;

  //******************************************************
  // Substep 2.5: Define the qoi function object
  //******************************************************

  s1_qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> s1_qoiRoutine_Data;
  s1_qoiRoutine_Data.beta1         = beta1;
  s1_qoiRoutine_Data.criticalMass1 = criticalMass1;
  uqVectorFunctionClass<P_V,P_M,Q_V,Q_M> s1_qoiFunctionObj(s1_qoiRoutine<P_V,P_M,Q_V,Q_M>,
                                                           (void *) &s1_qoiRoutine_Data);

  //******************************************************
  // Substep 2.6: Instantiate the first validation problem stage
  //******************************************************

  validProblem.instantiateStage(1,                         // stage id
                                s1_priorParamDensityObj,   // substep 2.1: calibration prior parameter density
                                &s1_likelihoodFunctionObj, // substep 2.2: calibration likelihood function
                                s1_proposalCovMatrix,      // substep 2.3: calibration gaussian proposal (density and generator)
                                s1_proposalDensityObj,     // substep 2.3: calibration proposal density
                                s1_proposalGeneratorObj,   // substep 2.3: calibration proposal generator
                                s1_propagParamRV,          // substep 2.4: propagation input parameter density
                                &s1_qoiFunctionObj);       // substep 2.5: propagation qoi function

  //******************************************************
  // Step 3 of 4: Set the second validation problem stage
  //******************************************************

  // This TGA example deals with only one stage

  //******************************************************
  // Step 4 of 4: solve the validation problem
  //******************************************************
  validProblem.solve();

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqTgaEx4' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_EX4_H__
