/* uq/examples/mcmc/tga/uqTgaEx.h
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

#ifndef __UQ_TGA_EX_H__
#define __UQ_TGA_EX_H__

//#include <uqStateSpace.h>
#include <uqProblem.h>
#include <uqLikelihoodFunction.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
int func(double t, const double Mass[], double f[], void *info)
{
  double* params = (double *)info;
  double A    = params[0];
  double E    = params[1];
  double beta = params[2];

  f[0] = -A*Mass[0]*exp(-E/(R_CONSTANT*t))/beta;

  return GSL_SUCCESS;
}

template<class S_V, class S_M>
struct
slice0_CompleteLikelihoodRoutine_DataType
{
  double               beta1;
  double               variance1;
  std::vector<double>* Te1; // temperatures
  std::vector<double>* Me1; // relative masses
};

template<class P_V,class P_M,class S_V,class S_M,class L_V,class L_M>
void slice0_CompleteLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr, L_V& resultValues)
{
  double A                       = paramValues[0];
  double E                       = paramValues[1];
  double beta1                   =  ((slice0_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->beta1;
  double variance1               =  ((slice0_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->variance1;
  const std::vector<double>& Te1 = *((slice0_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->Te1;
  const std::vector<double>& Me1 = *((slice0_CompleteLikelihoodRoutine_DataType<S_V,S_M> *) functionDataPtr)->Me1;
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
                        "slice0_CompleteLikelihoodRoutine()",
                        "gsl_odeiv_evolve_apply() failed");
    //printf("%6.1lf %10.4lf\n",t,Mass[0]);
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
  resultValues[0] = misfit1/variance1;
	
  //printf("loopSize = %d\n",loopSize);
  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    printf("For A = %g, E = %g, and beta1 = %.3lf, misfit1 = %lf and likelihood1 = %lf.\n",A,E,beta1,misfit1,resultValues[0]);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

//********************************************************
// The MCMC driving routine: called by main()
//********************************************************
template<class P_V,class P_M,class S_V,class S_M,class L_V,class L_M,class Q_V,class Q_M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqTgaEx' example\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //*****************************************************
  // Step 1 of 4: Instantiate the problem
  //*****************************************************
  uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M> problem(env);

  //******************************************************
  // Step 2 of 4: Define first problem slice
  //******************************************************

  // Substep 2.1: Define the prior prob. density function object: -2*ln[prior]
  // This application will pass NULL to 'setSlice()', that is,
  // this application will use the default prior density
  // provided by the MCMC library

  // Substep 2.2: Define the likelihood function object: -2*ln[likelihood]
  // Open input file on experimental data
  FILE *inp;
  inp = fopen("global.dat","r");

  // Read kinetic parameters and convert heating rate to K/s
  double tmpA, tmpE, beta1, variance1;
  fscanf(inp,"%lf %lf %lf %lf",&tmpA,&tmpE,&beta1,&variance1);
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

  slice0_CompleteLikelihoodRoutine_DataType<P_V,P_M> slice0_CompleteLikelihoodRoutine_Data;
  slice0_CompleteLikelihoodRoutine_Data.beta1     = beta1;
  slice0_CompleteLikelihoodRoutine_Data.variance1 = variance1;
  slice0_CompleteLikelihoodRoutine_Data.Te1       = &Te1; // temperatures
  slice0_CompleteLikelihoodRoutine_Data.Me1       = &Me1; // relative masses
  uq_CompleteLikelihoodFunction_Class<P_V,P_M,L_V,L_M> slice0_CompleteLikelihoodFunction_Obj(slice0_CompleteLikelihoodRoutine<P_V,P_M,S_V,S_M,L_V,L_M>,
                                                                                            (void *) &slice0_CompleteLikelihoodRoutine_Data,
                                                                                            true); // the routine computes [-2.*ln(Likelihood)]

  // Substep 2.3: Define the proposal density and proposal generator
  // This application will pass NULL to 'setSlice()', that is,
  // this application will use the default gaussian proposal and default gaussian generator with the default covariance matrix,
  // all provided by the MCMC library
  P_M* proposalCovMatrix = NULL;
  //P_M* proposalCovMatrix = slice0_ParamSpace.newMatrix();
  //(*proposalCovMatrix)(0,0) = 1.65122e+10;
  //(*proposalCovMatrix)(0,1) = 0.;
  //(*proposalCovMatrix)(1,0) = 0.;
  //(*proposalCovMatrix)(1,1) = 9.70225e+04;

  // Substep 2.4: Define the qoi function object

  // Substep 2.5: Set the first problem slice

  problem.instantiateSlice(0,                                      // slice id
                           NULL,                                   // calibration prior parameter density
                           &slice0_CompleteLikelihoodFunction_Obj, // calibration likelihood function
                           proposalCovMatrix,                      // calibration covariance matrix for gaussian proposal
                           NULL,                                   // calibration proposal density
                           NULL,                                   // calibration proposal generator
                           NULL,                                   // qoi input parameter density
                           NULL);                                  // qoi function

  //******************************************************
  // Step 3 of 4: Define second problem slice
  //******************************************************

  //******************************************************
  // Step 4 of 4: Quantify the uncertainty
  //******************************************************
  problem.quantify();

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqTgaEx' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_EX_H__
