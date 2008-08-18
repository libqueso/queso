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
#include <uqDRAM_MarkovChainGenerator.h>
#include <uqDefaultPrior.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472

//********************************************************
// The prior routine: provided by user and called by the MCMC tool.
// --> The application can use the default prior() routine provided by the MCMC tool.
//********************************************************
template<class V, class M>
struct
calib_M2lPriorRoutine_DataType
{
  const V* maybeYouNeedAVector;
  const M* maybeYouNeedAMatrix;
};

template<class V, class M>
double calib_M2lPriorRoutine(const V& paramValues, const void* functionDataPtr)
{
  const V& v1 = *((calib_M2lPriorRoutine_DataType<V,M> *) functionDataPtr)->maybeYouNeedAVector;
  const M& m1 = *((calib_M2lPriorRoutine_DataType<V,M> *) functionDataPtr)->maybeYouNeedAMatrix;

  UQ_FATAL_TEST_MACRO(true,
                      paramValues.env().rank(),
                      "calib_M2lPriorRoutine(), in uqTgaEx.h",
                      "should not be here, since application is using the default prior() routine provided by PECOS toolkit");

  return (m1 * v1).sumOfComponents();
}

//********************************************************
// The likelihood routine: provided by user and called by the MCMC tool.
//********************************************************
int func(double t, const double Mass[], double f[], void *info)
{
  double* params = (double *)info;
  double A = params[0],E = params[1],beta = params[2];
  f[0] = -A*Mass[0]*exp(-E/(R_CONSTANT*t))/beta;
  return GSL_SUCCESS;
}

template<class V, class M>
struct
calib_CompleteLikelihoodRoutine_DataType
{
  const V* aVectorForInstance;
  const M* aMatrixForInstance;
};

template<class V, class M>
void calib_CompleteLikelihoodRoutine(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
#if 0
  const V& v1 = *((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->aVectorForInstance;
  const M& m1 = *((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->aMatrixForInstance;

  const V& v2 = *((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->aVectorForInstance;
  const M& m2 = *((calib_CompleteLikelihoodRoutine_DataType<V,M> *) functionDataPtr)->aMatrixForInstance;

  V misfit(resultValues);
  misfit[0] = (m1 * v1).norm2Sq();
  misfit[1] = (m2 * v2).norm2Sq();

  V sigma(resultValues);
  sigma[0] = 1.;
  sigma[1] = .34;

  V minus2LogLikelihood(resultValues);
  minus2LogLikelihood = misfit/sigma;

  // The MCMC library will always expect the value [-2. * ln(Likelihood)], so this routine should use this line below
  resultValues = minus2LogLikelihood;

  // If, however, you really need to return the Likelihood(s), then use this line below
  //resultsValues[0] = exp(-.5*minus2LogLikelihood[0]);
  //resultsValues[1] = exp(-.5*minus2LogLikelihood[1]);
#endif
  double E2=0.;
  double A,E,beta,te[11],Me[11],Mt[11];
  int num_data;
  FILE *inp;

  inp = fopen("global.dat","r");
  fscanf(inp,"%lf %lf %lf",&A,&E,&beta);    /* read kinetic parameters */
  A = paramValues[0];
  E = paramValues[1];
	
  beta/=60.;	/* Convert heating rate to K/s */
  
  double params[]={A,E,beta};
      	
  // read experimental data
  int i=0;
  int status;
  while (1){
    status = fscanf(inp,"%lf %lf",&te[i],&Me[i]);
    if (status == EOF) break;
    i++;
  }
  num_data = i;
  fclose(inp);
	
  // integration
  const gsl_odeiv_step_type * T = gsl_odeiv_step_gear1;
	
  gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,1);
  gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-6,0.0);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(1);
	
  gsl_odeiv_system sys = {func, NULL, 1, (void *)params}; 
	
  double t = 0.1, t1 = 800.;
  double h = 1e-3;
  double Mass[1];
         Mass[0]=1.;
  
  i=0;
  double t_old=0., M_old[1];
  M_old[0]=1.;
	
  while (t < t1){
    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, Mass);
		
    if (status != GSL_SUCCESS) break;
			
    //printf("%6.1lf %10.4lf\n",t,Mass[0]);
		
    if ( (t >= te[i]) && (t_old <= te[i]) ) {
      Mt[i] = (te[i]-t_old)*(Mass[0]-M_old[0])/(t-t_old) + M_old[0];
      E2+=(Me[i]-Mt[i])*(Me[i]-Mt[i]);
      //printf("%i %lf %lf %lf %lf\n",i,te[i],Me[i],Mt[i],E2);
      i++;
    }
		
    t_old=t;
    M_old[0]=Mass[0];
  }
  resultValues[0] = E2/1.0e-5;
	
  printf("For A = %g, E = %g, and beta = %.3lf\n",A,E,beta);
  printf("the sum of squared errors is %lf.",E2);
  printf("\n");

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

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
#if 0 // You might need to use your own prior function
  calib_M2lPriorRoutine_DataType<V,M> calib_M2lPriorRoutine_Data;
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.maybeYouNeedAVector = &calib_ParamPriorSigmas;
  calib_M2lPriorRoutine_Data.maybeYouNeedAMatrix = NULL;

  uq_M2lProbDensity_Class<V,M> calib_M2lPriorProbDensity_Obj(calib_M2lPriorRoutine<V,M>,
                                                             (void *) &calib_M2lPriorRoutine_Data);
#else // Or you might want to use the default prior function provided by the MCMC tool
  uqDefault_M2lPriorRoutine_DataType<V,M> calib_M2lPriorRoutine_Data; // use default prior() routine
  V calib_ParamPriorMus   (calib_ParamSpace.priorMuValues   ());
  V calib_ParamPriorSigmas(calib_ParamSpace.priorSigmaValues());
  calib_M2lPriorRoutine_Data.paramPriorMus    = &calib_ParamPriorMus;
  calib_M2lPriorRoutine_Data.paramPriorSigmas = &calib_ParamPriorSigmas;

  uq_M2lProbDensity_Class<V,M> calib_M2lPriorProbDensity_Obj(uqDefault_M2lPriorRoutine<V,M>, // use default prior() routine
                                                             (void *) &calib_M2lPriorRoutine_Data);
#endif

  //******************************************************
  // Step 3 of 6: Define the likelihood prob. density function object: just misfits
  //******************************************************
  calib_CompleteLikelihoodRoutine_DataType<V,M> calib_CompleteLikelihoodRoutine_Data;
  V calib_ParamInitials(calib_ParamSpace.initialValues());
  calib_CompleteLikelihoodRoutine_Data.aVectorForInstance = &calib_ParamInitials;
  calib_CompleteLikelihoodRoutine_Data.aMatrixForInstance = NULL;

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
