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

#include <uqCalibProblem.h>
#include <uqPropagProblem.h>
#include <uqAsciiTable.h>
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
calibLikelihoodRoutine_DataType
{
  double               beta1;
  double               variance1;
  std::vector<double>* Te1; // temperatures
  std::vector<double>* Me1; // relative masses
};

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
calibLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  double A                       = paramValues[0];
  double E                       = paramValues[1];
  double beta1                   =  ((calibLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta1;
  double variance1               =  ((calibLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->variance1;
  const std::vector<double>& Te1 = *((calibLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Te1;
  const std::vector<double>& Me1 = *((calibLikelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Me1;
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
                        "calibLikelihoodRoutine()",
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
    printf("In calibLikelihoodRoutine(), A = %g, E = %g, beta1 = %.3lf: misfit1 = %lf, likelihood1 = %lf.\n",A,E,beta1,misfit1,resultValue);
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
propagQoiRoutine_DataType
{
  double beta1;
  double criticalMass1;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void propagQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  double A             = paramValues[0];
  double E             = paramValues[1];
  double beta1         = ((propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->beta1;
  double criticalMass1 = ((propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->criticalMass1;

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
                        "propagQoiRoutine()",
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
    printf("In propagQoiRoutine(), A = %g, E = %g, beta1 = %.3lf, criticalMass1 = %.3lf: qoi = %lf.\n",A,E,beta1,criticalMass1,qoiValues[0]);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

//********************************************************
// The driving routine "uqAppl()": called by main()
// Step   I.1 of 3: code very user specific to stage I
// Step   I.2 of 3: deal with the calibration problem of stage I 
// Step   I.3 of 3: deal with the propagation problem of stage I
// Step  II.1 of 3: code very user specific to stage II
// Step  II.2 of 3: deal with the calibration problem of stage II
// Step  II.3 of 3: deal with the propagation problem of stage II
// Step III.1 of 1: compare the cdf's of stages I and II
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

  //******************************************************
  // Read Ascii file with important information on both calibration problems.
  //******************************************************
  uqAsciiTableClass<P_V,P_M> calTable(env,
                                      2,    // # of rows
                                      3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                      NULL, // All extra columns are of 'double' type
                                      "cal.tab");

  const EpetraExt::DistArray<std::string>& paramNames = calTable.stringColumn(0);
  P_V                                      s1_minValues    (calTable.doubleColumn(1));
  P_V                                      s1_maxValues    (calTable.doubleColumn(2));
  P_V                                      s1_initialValues(calTable.doubleColumn(3));

  //******************************************************
  // Read Ascii file with important information on both propagation problems.
  //******************************************************
  uqAsciiTableClass<P_V,P_M> proTable(env,
                                      1,    // # of rows
                                      0,    // # of cols after 'parameter name': none
                                      NULL, // All extra columns are of 'double' type
                                      "pro.tab");

  const EpetraExt::DistArray<std::string>& qoiNames = proTable.stringColumn(0);

  //*****************************************************
  // Step I.1 of 3: Code very specific to this TGA example
  //*****************************************************

  // Open input file on experimental data
  FILE *inp;
  inp = fopen("s1_global4.dat","r");

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
  // Usually, spaces are the same throughout different problems.
  // If this is the case, we can instantiate them here, just once.
  //******************************************************
  uqVectorSpaceClass<P_V,P_M> paramSpace(env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         calTable.numRows(),
                                         &paramNames);
  uqVectorSpaceClass<Q_V,Q_M> qoiSpace  (env,
                                         "qoi_",   // Extra prefix before the default "space_" prefix
                                         proTable.numRows(),
                                         &qoiNames);

  //******************************************************
  // Step I.2 of 3: deal with the calibration problem
  //******************************************************

  // Stage I (s1): Prior vector rv
  uqUniformVectorRVClass<P_V,P_M> s1_calibPriorRv("s1_cal_prior_", // Extra prefix before the default "rv_" prefix
                                                   paramSpace,
                                                   s1_minValues,
                                                   s1_maxValues);

  // Stage I (s1): 1Likelihood function object: -2*ln[likelihood]
  calibLikelihoodRoutine_DataType<P_V,P_M> s1_calibLikelihoodRoutine_Data;
  s1_calibLikelihoodRoutine_Data.beta1     = beta1;
  s1_calibLikelihoodRoutine_Data.variance1 = variance1;
  s1_calibLikelihoodRoutine_Data.Te1       = &Te1; // temperatures
  s1_calibLikelihoodRoutine_Data.Me1       = &Me1; // relative masses
  uqGenericVectorPdfClass<P_V,P_M> s1_calibLikelihoodFunctionObj("s1_cal_prior_like_", // Extra prefix before the default "genpd_" prefix
                                                                 paramSpace,
                                                                 calibLikelihoodRoutine<P_V,P_M>,
                                                                 (void *) &s1_calibLikelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(Likelihood)]

  // Stage I (s1): Posterior vector rv
  uqGenericVectorRVClass<P_V,P_M> s1_calibPostRv("s1_cal_post_", // Extra prefix before the default "rv_" prefix
                                                 paramSpace);

  // Stage I (s1): Calibration problem
  uqCalibProblemClass<P_V,P_M> s1_calibProblem("s1_", // No extra prefix before the default "cal_" prefix
                                               s1_calibPriorRv,
                                               s1_calibLikelihoodFunctionObj,
                                               s1_calibPostRv);

  // Stage I (s1): Solve the calibration problem: set 'pdf' and 'realizer' of 's1_calibPostRv'
  P_M* s1_calibProposalCovMatrix = s1_calibPostRv.imageSpace().newGaussianMatrix(s1_calibPriorRv.pdf().domainVarianceValues(),
                                                                                 s1_initialValues);
  s1_calibProblem.solveWithBayesMarkovChain(s1_initialValues,
                                           *s1_calibProposalCovMatrix,
                                            NULL); // use default kernel from library

  //******************************************************
  // Step I.3 of 3: deal with the propagation problem
  //******************************************************

  // Stage I (s1): Input param vector rv for propagation = output posterior vector rv of calibration

  // Stage I (s1): Qoi function object
  propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> s1_propagQoiRoutine_Data;
  s1_propagQoiRoutine_Data.beta1         = beta1;
  s1_propagQoiRoutine_Data.criticalMass1 = criticalMass1;
  uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M> s1_propagQoiFunctionObj("s1_pro_qoi_", // Extra prefix before the default "func_" prefix
                                                                        paramSpace,
                                                                        qoiSpace,
                                                                        propagQoiRoutine<P_V,P_M,Q_V,Q_M>,
                                                                        (void *) &s1_propagQoiRoutine_Data);

  // Stage I (s1): Qoi vector rv
  uqGenericVectorRVClass<Q_V,Q_M> s1_propagQoiRv("s1_pro_qoi_", // Extra prefix before the default "rv_" prefix
                                                 qoiSpace);

  // Stage I (s1): Propagation problem
  uqPropagProblemClass<P_V,P_M,Q_V,Q_M> s1_propagProblem("s1_",          // No extra prefix before the default "pro_" prefix
                                                         s1_calibPostRv, // propagation input = calibration output
                                                         s1_propagQoiFunctionObj,
                                                         s1_propagQoiRv);

  // Stage I (s1): Solve the propagation problem: set 'realizer' and 'cdf' of 's1_propagQoiRv'
  s1_propagProblem.solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  //*****************************************************
  // Step II.1 of 3: Code very specific to this TGA example
  //*****************************************************

  // Open input file on experimental data
  inp = fopen("s2_global4.dat","r");

  // Read kinetic parameters and convert heating rate to K/s
  fscanf(inp,"%lf %lf %lf %lf %lf",&tmpA,&tmpE,&beta1,&variance1,&criticalMass1);
  beta1 /= 60.;
  
  // Read experimental data
  Te1.clear();
  Te1.resize(11,0.);
  Me1.clear();
  Me1.resize(11,0.);

  numObservations = 0;
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
#if 1
  //******************************************************
  // Step II.2 of 3: deal with the calibration problem
  //******************************************************

  // Stage II (s2): Prior vector rv = posterior vector rv of stage I (s1)

  // Stage II (s2): 1Likelihood function object: -2*ln[likelihood]
  calibLikelihoodRoutine_DataType<P_V,P_M> s2_calibLikelihoodRoutine_Data;
  s2_calibLikelihoodRoutine_Data.beta1     = beta1;
  s2_calibLikelihoodRoutine_Data.variance1 = variance1;
  s2_calibLikelihoodRoutine_Data.Te1       = &Te1; // temperatures
  s2_calibLikelihoodRoutine_Data.Me1       = &Me1; // relative masses
  uqGenericVectorPdfClass<P_V,P_M> s2_calibLikelihoodFunctionObj("s2_cal_prior_like_", // Extra prefix before the default "genpd_" prefix
                                                                 paramSpace,
                                                                 calibLikelihoodRoutine<P_V,P_M>,
                                                                 (void *) &s2_calibLikelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(Likelihood)]

  // Stage II (s2): Posterior vector rv
  uqGenericVectorRVClass<P_V,P_M> s2_calibPostRv("s2_cal_post_", // Extra prefix before the default "rv_" prefix
                                                 paramSpace);

  // Stage II (s2): Calibration problem
  uqCalibProblemClass<P_V,P_M> s2_calibProblem("s2_", // No extra prefix before the default "cal_" prefix
                                               s1_calibPostRv, // s2 calibration input = s1 calibration output
                                               s2_calibLikelihoodFunctionObj,
                                               s2_calibPostRv);

  // Stage II (s2): Solve the calibration problem: set 'pdf' and 'realizer' of 's2_calibPostRv'
  P_M* s2_calibProposalCovMatrix = s1_calibPostRv.imageSpace().newGaussianMatrix(s1_calibPostRv.realizer().imageVarianceValues(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                 s1_calibPostRv.realizer().imageExpectedValues()); // Use these values as the initial values
  s2_calibProblem.solveWithBayesMarkovChain(s1_calibPostRv.realizer().imageExpectedValues(),
                                           *s2_calibProposalCovMatrix,
                                            NULL); // use default kernel from library

  //******************************************************
  // Step II.3 of 3: deal with the propagation problem
  //******************************************************

  // Stage II (s2): Input param vector rv for propagation = output posterior vector rv of calibration

  // Stage II (s2): Qoi function object
  propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> s2_propagQoiRoutine_Data;
  s2_propagQoiRoutine_Data.beta1         = beta1;
  s2_propagQoiRoutine_Data.criticalMass1 = criticalMass1;
  uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M> s2_propagQoiFunctionObj("s2_pro_qoi_", // Extra prefix before the default "func_" prefix
                                                                        paramSpace,
                                                                        qoiSpace,
                                                                        propagQoiRoutine<P_V,P_M,Q_V,Q_M>,
                                                                        (void *) &s2_propagQoiRoutine_Data);

  // Stage II (s2): Qoi vector rv: set 'realizer' and 'cdf' of 's2_propagQoiRv'
  uqGenericVectorRVClass<Q_V,Q_M> s2_propagQoiRv("s2_pro_qoi_", // Extra prefix before the default "rv_" prefix
                                                 qoiSpace);

  // Stage II (s2): Propagation problem
  uqPropagProblemClass<P_V,P_M,Q_V,Q_M> s2_propagProblem("s2_",          // No extra prefix before the default "pro_" prefix
                                                         s2_calibPostRv, // propagation input = calibration output
                                                         s2_propagQoiFunctionObj,
                                                         s2_propagQoiRv);

  // Stage II (s2): Solve the propagation problem
  s2_propagProblem.solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************
  delete s2_calibProposalCovMatrix;
#endif
  delete s1_calibProposalCovMatrix;

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqTgaEx4' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_EX4_H__
