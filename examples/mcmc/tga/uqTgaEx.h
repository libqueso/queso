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

double observationsOfTemperatureAndRelativeMass[] = {
  673.034, 0.965855,
  682.003, 0.951549,
  690.985, 0.925048,
  699.979, 0.886353,
  708.989, 0.830585,
  718.02,  0.755306,
  727.089, 0.641003,
  735.96,  0.475474,
  744.904, 0.236777,
  754.062, 0.032234,
  763.049, 0.000855448
};

//#include <uqStateSpace.h>
#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> This application uses the default prior() routine provided by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_M2lPriorRoutine_DataType
{
  const V* nothing1;
  const M* nothing2;
};

template<class V, class M>
double calib_M2lPriorRoutine(const V& paramValues, const void* functionDataPtr)
{
  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_M2lPriorRoutine(), in uqTgaEx.h",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");

  return 0.;
}

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

template<class V, class M>
struct
calib_CompleteLikelihoodRoutine_DataType
{
  double               beta1;
  double               variance1;
  std::vector<double>* Te1; // temperatures
  std::vector<double>* Me1; // relative masses
};

template<class V, class M>
void calib_CompleteLikelihoodRoutine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  double A = paramValues[0];
  double E = paramValues[1];
  double beta1             =  ((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->beta1;
  double variance1         =  ((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->variance1;
  const std::vector<double>& Te1 = *((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->Te1;
  const std::vector<double>& Me1 = *((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->Me1;
  std::vector<double> Mt1(Me1.size(),0.);

  double params[]={A,E,beta1};
      	
  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //gear1;
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
                        "calib_CompleteLikelihoodRoutine()",
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
template<class V, class M>
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
  // Step 1 of 6: Define the finite dimensional linear spaces.
  //*****************************************************
  uqParamSpaceClass<V,M>      calib_ParamSpace     (env,"calib");
  uqObservableSpaceClass<V,M> calib_ObservableSpace(env,"calib");

  //******************************************************
  // Step 2 of 6: Define the prior prob. density function object: -2*ln[prior]
  //******************************************************
  uqDefault_M2lPriorRoutine_DataType<V,M> calib_M2lPriorRoutine_Data; // use default prior() routine
  V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uq_M2lProbDensity_Class<V,M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                             (void *) &calib_M2lPriorRoutine_Data);

  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************

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

  calib_CompleteLikelihoodRoutine_DataType<V,M> calib_CompleteLikelihoodRoutine_Data;
  calib_CompleteLikelihoodRoutine_Data.beta1     = beta1;
  calib_CompleteLikelihoodRoutine_Data.variance1 = variance1;
  calib_CompleteLikelihoodRoutine_Data.Te1       = &Te1; // temperatures
  calib_CompleteLikelihoodRoutine_Data.Me1       = &Me1; // relative masses
  uq_CompleteLikelihoodFunction_Class<V,M> calib_CompleteLikelihoodFunction_Obj(calib_CompleteLikelihoodRoutine<V,M>,
                                                                                (void *) &calib_CompleteLikelihoodRoutine_Data,
                                                                                true); // the routine computes [-2.*ln(Likelihood)]

  //******************************************************
  // Step 4 of 6: Define the Markov chain generator.
  //******************************************************
  uqDRAM_MarkovChainGeneratorClass<V,M> mcg(env,
                                            "calib",
                                            calib_ParamSpace,
                                            calib_ObservableSpace,
                                            calib_M2lPriorProbDensity_Obj,
                                            calib_CompleteLikelihoodFunction_Obj);

  //******************************************************
  // Step 5 of 6: Compute the proposal covariance matrix.
  //******************************************************
  // Proposal covariance matrix will be computed internally
  // by the Markov chain generator.
  M* proposalCovMatrix = NULL;
  //M* proposalCovMatrix = calib_ParamSpace.newMatrix();
  //(*proposalCovMatrix)(0,0) = 1.65122e+10;
  //(*proposalCovMatrix)(0,1) = 0.;
  //(*proposalCovMatrix)(1,0) = 0.;
  //(*proposalCovMatrix)(1,1) = 9.70225e+04;

  //******************************************************
  // Step 6 of 6: Generate chains.
  //              Output data is written (in MATLAB format) to the file
  //              with name specified by the user in the input file.
  //******************************************************
  mcg.generateChains(proposalCovMatrix);

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
