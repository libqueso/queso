/* uq/examples/queso/tga/uqTgaEx4.h
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
// Likelihood function object for the first validation problem stage (with prefix "cal_").
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************

// The (user defined) data type for the data needed by the (user defined) likelihood routine
template<class P_V, class P_M>
struct
likelihoodRoutine_DataType
{
  double               beta1;
  double               variance1;
  std::vector<double>* Te1; // temperatures
  std::vector<double>* Me1; // relative masses

  double               beta2;
  double               variance2;
  std::vector<double>* Te2; // temperatures
  std::vector<double>* Me2; // relative masses

  double               beta3;
  double               variance3;
  std::vector<double>* Te3; // temperatures
  std::vector<double>* Me3; // relative masses
};

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
likelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  double resultValue = 0.;

  // Compute likelihood for scenario 1
  double betaTest = ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta1;
  if (betaTest) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta1;
    double variance                =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->variance1;
    const std::vector<double>& Te  = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Te1;
    const std::vector<double>& Me  = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Me1;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {func, NULL, 1, (void *)params}; 

    double t = 0.1, t_final = 1900.;
    double h = 1e-3;
    double Mass[1];
    Mass[0]=1.;
  
    unsigned int i = 0;
    double t_old = 0.;
    double M_old[1];
    M_old[0]=1.;
	
    double misfit=0.;
    //unsigned int loopSize = 0;
    while ((t < t_final) && (i < Me.size())) {
      int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_final, &h, Mass);
      UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                          paramValues.env().rank(),
                          "likelihoodRoutine()",
                          "gsl_odeiv_evolve_apply() failed");
      //printf("t = %6.1lf, mass = %10.4lf\n",t,Mass[0]);
      //loopSize++;
		
      while ( (i < Me.size()) && (t_old <= Te[i]) && (Te[i] <= t) ) {
        Mt[i] = (Te[i]-t_old)*(Mass[0]-M_old[0])/(t-t_old) + M_old[0];
        misfit += (Me[i]-Mt[i])*(Me[i]-Mt[i]);
        //printf("%i %lf %lf %lf %lf\n",i,Te[i],Me[i],Mt[i],misfit);
        i++;
      }
		
      t_old=t;
      M_old[0]=Mass[0];
    }
    resultValue += misfit/variance;
	
    //printf("loopSize = %d\n",loopSize);
    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      printf("In likelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  // Compute likelihood for scenario 2
  betaTest = ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta2;
  if (betaTest > 0.) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta2;
    double variance                =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->variance2;
    const std::vector<double>& Te  = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Te2;
    const std::vector<double>& Me  = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Me2;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {func, NULL, 1, (void *)params}; 

    double t = 0.1, t_final = 1900.;
    double h = 1e-3;
    double Mass[1];
    Mass[0]=1.;
  
    unsigned int i = 0;
    double t_old = 0.;
    double M_old[1];
    M_old[0]=1.;
	
    double misfit=0.;
    //unsigned int loopSize = 0;
    while ((t < t_final) && (i < Me.size())) {
      int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_final, &h, Mass);
      UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                          paramValues.env().rank(),
                          "likelihoodRoutine()",
                          "gsl_odeiv_evolve_apply() failed");
      //printf("t = %6.1lf, mass = %10.4lf\n",t,Mass[0]);
      //loopSize++;
		
      while ( (i < Me.size()) && (t_old <= Te[i]) && (Te[i] <= t) ) {
        Mt[i] = (Te[i]-t_old)*(Mass[0]-M_old[0])/(t-t_old) + M_old[0];
        misfit += (Me[i]-Mt[i])*(Me[i]-Mt[i]);
        //printf("%i %lf %lf %lf %lf\n",i,Te[i],Me[i],Mt[i],misfit);
        i++;
      }
		
      t_old=t;
      M_old[0]=Mass[0];
    }
    resultValue += misfit/variance;
	
    //printf("loopSize = %d\n",loopSize);
    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      printf("In likelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  // Compute likelihood for scenario 3
  betaTest = ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta3;
  if (betaTest > 0.) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->beta3;
    double variance                =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->variance3;
    const std::vector<double>& Te  = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Te3;
    const std::vector<double>& Me  = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->Me3;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {func, NULL, 1, (void *)params}; 

    double t = 0.1, t_final = 1900.;
    double h = 1e-3;
    double Mass[1];
    Mass[0]=1.;
  
    unsigned int i = 0;
    double t_old = 0.;
    double M_old[1];
    M_old[0]=1.;
	
    double misfit=0.;
    //unsigned int loopSize = 0;
    while ((t < t_final) && (i < Me.size())) {
      int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_final, &h, Mass);
      UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                          paramValues.env().rank(),
                          "likelihoodRoutine()",
                          "gsl_odeiv_evolve_apply() failed");
      //printf("t = %6.1lf, mass = %10.4lf\n",t,Mass[0]);
      //loopSize++;
		
      while ( (i < Me.size()) && (t_old <= Te[i]) && (Te[i] <= t) ) {
        Mt[i] = (Te[i]-t_old)*(Mass[0]-M_old[0])/(t-t_old) + M_old[0];
        misfit += (Me[i]-Mt[i])*(Me[i]-Mt[i]);
        //printf("%i %lf %lf %lf %lf\n",i,Te[i],Me[i],Mt[i],misfit);
        i++;
      }
		
      t_old=t;
      M_old[0]=Mass[0];
    }
    resultValue += misfit/variance;
	
    //printf("loopSize = %d\n",loopSize);
    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      printf("In likelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  return resultValue;
}

//********************************************************
// QoI function object for the first validation problem stage (with prefix "cal_").
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data type for the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
propagQoiRoutine_DataType
{
  double beta;
  double criticalMass;
  double criticalTime;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void propagQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  double A             = paramValues[0];
  double E             = paramValues[1];
  double beta          = ((propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->beta;
  double criticalMass  = ((propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->criticalMass;
  double criticalTime  = ((propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->criticalTime;

  double params[]={A,E,beta};
      	
  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
        gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
        gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sys = {func, NULL, 1, (void *)params}; 
	
  double temperature = 0.1;
  double h = 1e-3;
  double Mass[1];
  Mass[0]=1.;
  
  double temperature_old = 0.;
  double M_old[1];
  M_old[0]=1.;
	
  double crossingTemperature = 0.;
  //unsigned int loopSize = 0;
  while ((temperature < criticalTime*beta) &&
         (Mass[0]     > criticalMass     )) {
    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &temperature, criticalTime*beta, &h, Mass);
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        paramValues.env().rank(),
                        "propagQoiRoutine()",
                        "gsl_odeiv_evolve_apply() failed");
    //printf("t = %6.1lf, mass = %10.4lf\n",t,Mass[0]);
    //loopSize++;

    if (Mass[0] <= criticalMass) {
      crossingTemperature = temperature_old + (temperature - temperature_old) * (M_old[0]-criticalMass)/(M_old[0]-Mass[0]);
    }
		
    temperature_old=temperature;
    M_old[0]=Mass[0];
  }

  if (criticalMass > 0.) qoiValues[0] = crossingTemperature/beta; // QoI = time to achieve critical mass
  if (criticalTime > 0.) qoiValues[0] = Mass[0];                  // QoI = mass fraction remaining at critical time
	
  //printf("loopSize = %d\n",loopSize);
  if ((paramValues.env().verbosity() >= 3) && (paramValues.env().rank() == 0)) {
    printf("In propagQoiRoutine(), A = %g, E = %g, beta = %.3lf, criticalTime = %.3lf, criticalMass = %.3lf: qoi = %lf.\n",A,E,beta,criticalTime,criticalMass,qoiValues[0]);
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
  uqAsciiTableClass<P_V,P_M> paramsTable(env,
                                         2,    // # of rows
                                         3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                         NULL, // All extra columns are of 'double' type
                                         "params.tab");

  const EpetraExt::DistArray<std::string>& paramNames = paramsTable.stringColumn(0);
  P_V                                      s1_minValues    (paramsTable.doubleColumn(1));
  P_V                                      s1_maxValues    (paramsTable.doubleColumn(2));
  P_V                                      s1_initialValues(paramsTable.doubleColumn(3));

  //******************************************************
  // Read Ascii file with important information on both propagation problems.
  //******************************************************
  uqAsciiTableClass<P_V,P_M> qoisTable(env,
                                       1,    // # of rows
                                       0,    // # of cols after 'parameter name': none
                                       NULL, // All extra columns are of 'double' type
                                       "qois.tab");

  const EpetraExt::DistArray<std::string>& qoiNames = qoisTable.stringColumn(0);

  double beta_prediction         = 250.;
  double criticalMass_prediction = 0.;
  double criticalTime_prediction = 3.9;

  //*****************************************************
  // Step I.1 of 3: Code very specific to this TGA example
  //*****************************************************

  int iRC;
  struct timeval timevalRef;
  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning stage 1 at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Read experimental data
  double beta_cal1, variance_cal1;
  std::vector<double> Te_cal1(11,0.);
  std::vector<double> Me_cal1(11,0.);

  {
    // Open input file on experimental data
    FILE *inp;
    inp = fopen("scenario_5_K_min.dat","r");

    // Read kinetic parameters and convert heating rate to K/s
    fscanf(inp,"%lf %lf",&beta_cal1,&variance_cal1);
    beta_cal1 /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= Te_cal1.size()),
                          env.rank(),
                          "uqAppl(), in uqTgaEx.h",
                          "input file has too many observations");
      Te_cal1[numObservations] = tmpTe;
      Me_cal1[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != Te_cal1.size()),
                        env.rank(),
                        "uqAppl(), in uqTgaEx.h",
                        "input file has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }

  // Read experimental data
  double beta_cal2, variance_cal2;
  std::vector<double> Te_cal2(11,0.);
  std::vector<double> Me_cal2(11,0.);

  {
    // Open input file on experimental data
    FILE *inp;
    inp = fopen("scenario_25_K_min.dat","r");

    // Read kinetic parameters and convert heating rate to K/s
    fscanf(inp,"%lf %lf",&beta_cal2,&variance_cal2);
    beta_cal2 /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= Te_cal2.size()),
                          env.rank(),
                          "uqAppl(), in uqTgaEx.h",
                          "input file has too many observations");
      Te_cal2[numObservations] = tmpTe;
      Me_cal2[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != Te_cal2.size()),
                        env.rank(),
                        "uqAppl(), in uqTgaEx.h",
                        "input file has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }

  // Read experimental data
  double beta_cal3, variance_cal3;
  std::vector<double> Te_cal3(11,0.);
  std::vector<double> Me_cal3(11,0.);

  {
    // Open input file on experimental data
    FILE *inp;
    inp = fopen("scenario_50_K_min.dat","r");

    // Read kinetic parameters and convert heating rate to K/s
    fscanf(inp,"%lf %lf",&beta_cal3,&variance_cal3);
    beta_cal3 /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= Te_cal3.size()),
                          env.rank(),
                          "uqAppl(), in uqTgaEx.h",
                          "input file has too many observations");
      Te_cal3[numObservations] = tmpTe;
      Me_cal3[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != Te_cal3.size()),
                        env.rank(),
                        "uqAppl(), in uqTgaEx.h",
                        "input file has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }

  //******************************************************
  // Usually, spaces are the same throughout different problems.
  // If this is the case, we can instantiate them here, just once.
  //******************************************************
  uqVectorSpaceClass<P_V,P_M> paramSpace(env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         paramsTable.numRows(),
                                         &paramNames);
  uqVectorSpaceClass<Q_V,Q_M> qoiSpace  (env,
                                         "qoi_",   // Extra prefix before the default "space_" prefix
                                         qoisTable.numRows(),
                                         &qoiNames);

  //******************************************************
  // Step I.2 of 3: deal with the calibration problem
  //******************************************************

  // Stage I (s1): Prior vector rv
  uqUniformVectorRVClass<P_V,P_M> cal_ip_PriorRv("cal_ip_prior_", // Extra prefix before the default "rv_" prefix
                                                 paramSpace,
                                                 s1_minValues,
                                                 s1_maxValues);

  // Stage I (s1): Likelihood function object: -2*ln[likelihood]
  likelihoodRoutine_DataType<P_V,P_M> cal_LikelihoodRoutine_Data;
  cal_LikelihoodRoutine_Data.beta1     = beta_cal1;
  cal_LikelihoodRoutine_Data.variance1 = variance_cal1;
  cal_LikelihoodRoutine_Data.Te1       = &Te_cal1; // temperatures
  cal_LikelihoodRoutine_Data.Me1       = &Me_cal1; // relative masses
  cal_LikelihoodRoutine_Data.beta2     = beta_cal2;
  cal_LikelihoodRoutine_Data.variance2 = variance_cal2;
  cal_LikelihoodRoutine_Data.Te2       = &Te_cal2; // temperatures
  cal_LikelihoodRoutine_Data.Me2       = &Me_cal2; // relative masses
  cal_LikelihoodRoutine_Data.beta3     = beta_cal3;
  cal_LikelihoodRoutine_Data.variance3 = variance_cal3;
  cal_LikelihoodRoutine_Data.Te3       = &Te_cal3; // temperatures
  cal_LikelihoodRoutine_Data.Me3       = &Me_cal3; // relative masses
  uqGenericVectorPdfClass<P_V,P_M> cal_LikelihoodFunctionObj("cal_ip_prior_Like_", // Extra prefix before the default "genpd_" prefix
                                                                paramSpace,
                                                                likelihoodRoutine<P_V,P_M>,
                                                                (void *) &cal_LikelihoodRoutine_Data,
                                                                true); // the routine computes [-2.*ln(Likelihood)]

  // Stage I (s1): Posterior vector rv
  uqGenericVectorRVClass<P_V,P_M> cal_ip_PostRv("cal_ip_post_", // Extra prefix before the default "rv_" prefix
                                                paramSpace);

  // Stage I (s1): Calibration problem
  uqCalibProblemClass<P_V,P_M> cal_ip_Problem("cal_", // No extra prefix before the default "cal_" prefix
                                              cal_ip_PriorRv,
                                              cal_LikelihoodFunctionObj,
                                              cal_ip_PostRv);

  // Stage I (s1): Solve the calibration problem: set 'pdf' and 'realizer' of 'cal_ip_PostRv'
  P_M* cal_ip_ProposalCovMatrix = cal_ip_PostRv.imageSpace().newGaussianMatrix(cal_ip_PriorRv.pdf().domainVarianceValues(),
                                                                                 s1_initialValues);
  cal_ip_Problem.solveWithBayesMarkovChain(s1_initialValues,
                                           *cal_ip_ProposalCovMatrix,
                                            NULL); // use default kernel from library

  //******************************************************
  // Step I.3 of 3: deal with the propagation problem
  //******************************************************

  // Stage I (s1): Input param vector rv for propagation = output posterior vector rv of calibration

  // Stage I (s1): Qoi function object
  propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> cal_fp_QoiRoutine_Data;
  cal_fp_QoiRoutine_Data.beta         = beta_prediction;
  cal_fp_QoiRoutine_Data.criticalMass = criticalMass_prediction;
  cal_fp_QoiRoutine_Data.criticalTime = criticalTime_prediction;
  uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M> cal_fp_QoiFunctionObj("cal_fp_qoi_", // Extra prefix before the default "func_" prefix
                                                                      paramSpace,
                                                                      qoiSpace,
                                                                      propagQoiRoutine<P_V,P_M,Q_V,Q_M>,
                                                                      (void *) &cal_fp_QoiRoutine_Data);

  // Stage I (s1): Qoi vector rv
  uqGenericVectorRVClass<Q_V,Q_M> cal_fp_QoiRv("cal_fp_qoi_", // Extra prefix before the default "rv_" prefix
                                               qoiSpace);

  // Stage I (s1): Propagation problem
  uqPropagProblemClass<P_V,P_M,Q_V,Q_M> cal_fp_Problem("cal_",          // No extra prefix before the default "pro_" prefix
                                                       cal_ip_PostRv, // propagation input = calibration output
                                                       cal_fp_QoiFunctionObj,
                                                       cal_fp_QoiRv);

  // Stage I (s1): Solve the propagation problem: set 'realizer' and 'cdf' of 'cal_fp_QoiRv'
  cal_fp_Problem.solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  struct timeval timevalNow;
  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending stage 1 at " << ctime(&timevalNow.tv_sec)
              << "Total s1 run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //*****************************************************
  // Step II.1 of 3: Code very specific to this TGA example
  //*****************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning stage 2 at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Read experimental data
  double beta_val, variance_val;
  std::vector<double> Te_val(11,0.);
  std::vector<double> Me_val(11,0.);

  {
    // Open input file on experimental data
    FILE *inp = fopen("scenario_100_K_min.dat","r");

    // Read kinetic parameters and convert heating rate to K/s
    fscanf(inp,"%lf %lf",&beta_val,&variance_val);
    beta_val /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= Te_val.size()),
                          env.rank(),
                          "uqAppl(), in uqTgaEx.h",
                          "input file has too many observations");
      Te_val[numObservations] = tmpTe;
      Me_val[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != Te_val.size()),
                        env.rank(),
                        "uqAppl(), in uqTgaEx.h",
                        "input file has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }

  //******************************************************
  // Step II.2 of 3: deal with the calibration problem
  //******************************************************

  // Stage II (s2): Prior vector rv = posterior vector rv of stage I (s1)

  // Stage II (s2): Likelihood function object: -2*ln[likelihood]
  likelihoodRoutine_DataType<P_V,P_M> val_LikelihoodRoutine_Data;
  val_LikelihoodRoutine_Data.beta1     = beta_val;
  val_LikelihoodRoutine_Data.variance1 = variance_val;
  val_LikelihoodRoutine_Data.Te1       = &Te_val; // temperatures
  val_LikelihoodRoutine_Data.Me1       = &Me_val; // relative masses
  val_LikelihoodRoutine_Data.beta2     = 0.;
  val_LikelihoodRoutine_Data.variance2 = 0.;
  val_LikelihoodRoutine_Data.Te2       = NULL; // temperatures
  val_LikelihoodRoutine_Data.Me2       = NULL; // relative masses
  val_LikelihoodRoutine_Data.beta3     = 0.;
  val_LikelihoodRoutine_Data.variance3 = 0.;
  val_LikelihoodRoutine_Data.Te3       = NULL; // temperatures
  val_LikelihoodRoutine_Data.Me3       = NULL; // relative masses
  uqGenericVectorPdfClass<P_V,P_M> val_LikelihoodFunctionObj("s2_cal_prior_Like_", // Extra prefix before the default "genpd_" prefix
                                                                 paramSpace,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &val_LikelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(Likelihood)]

  // Stage II (s2): Posterior vector rv
  uqGenericVectorRVClass<P_V,P_M> val_ip_PostRv("s2_cal_post_", // Extra prefix before the default "rv_" prefix
                                                paramSpace);

  // Stage II (s2): Calibration problem
  uqCalibProblemClass<P_V,P_M> val_ip_Problem("s2_", // No extra prefix before the default "cal_" prefix
                                              cal_ip_PostRv, // s2 calibration input = s1 calibration output
                                              val_LikelihoodFunctionObj,
                                              val_ip_PostRv);

  // Stage II (s2): Solve the calibration problem: set 'pdf' and 'realizer' of 'val_ip_PostRv'
  P_M* val_ip_ProposalCovMatrix = cal_ip_PostRv.imageSpace().newGaussianMatrix(cal_ip_PostRv.realizer().imageVarianceValues(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                 cal_ip_PostRv.realizer().imageExpectedValues()); // Use these values as the initial values
  val_ip_Problem.solveWithBayesMarkovChain(cal_ip_PostRv.realizer().imageExpectedValues(),
                                           *val_ip_ProposalCovMatrix,
                                           NULL); // use default kernel from library

  //******************************************************
  // Step II.3 of 3: deal with the propagation problem
  //******************************************************

  // Stage II (s2): Input param vector rv for propagation = output posterior vector rv of calibration

  // Stage II (s2): Qoi function object
  propagQoiRoutine_DataType<P_V,P_M,Q_V,Q_M> val_fp_QoiRoutine_Data;
  val_fp_QoiRoutine_Data.beta          = beta_prediction;
  val_fp_QoiRoutine_Data.criticalMass  = criticalMass_prediction;
  val_fp_QoiRoutine_Data.criticalTime  = criticalTime_prediction;
  uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M> val_fp_QoiFunctionObj("s2_pro_qoi_", // Extra prefix before the default "func_" prefix
                                                                        paramSpace,
                                                                        qoiSpace,
                                                                        propagQoiRoutine<P_V,P_M,Q_V,Q_M>,
                                                                        (void *) &val_fp_QoiRoutine_Data);

  // Stage II (s2): Qoi vector rv: set 'realizer' and 'cdf' of 'val_fp_QoiRv'
  uqGenericVectorRVClass<Q_V,Q_M> val_fp_QoiRv("s2_pro_qoi_", // Extra prefix before the default "rv_" prefix
                                                 qoiSpace);

  // Stage II (s2): Propagation problem
  uqPropagProblemClass<P_V,P_M,Q_V,Q_M> val_fp_Problem("s2_",          // No extra prefix before the default "pro_" prefix
                                                         val_ip_PostRv, // propagation input = calibration output
                                                         val_fp_QoiFunctionObj,
                                                         val_fp_QoiRv);

  // Stage II (s2): Solve the propagation problem
  val_fp_Problem.solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending stage 2 at " << ctime(&timevalNow.tv_sec)
              << "Total s2 run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //******************************************************
  // Step III.1 of 1: compare the cdf's of stages I and II
  //******************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning stage 3 at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  if (cal_fp_Problem.computeSolutionFlag() &&
      val_fp_Problem.computeSolutionFlag()) {
    Q_V* epsilonVec = cal_fp_QoiRv.imageSpace().newVector(0.02);
    Q_V cdfDistancesVec(cal_fp_QoiRv.imageSpace().zeroVector());
    horizontalDistances(cal_fp_QoiRv.cdf(),
                        val_fp_QoiRv.cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (env.rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(val_fp_QoiRv.cdf(),
                        cal_fp_QoiRv.cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (env.rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec (swithced order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(cal_fp_QoiRv.cdf(),
                        val_fp_QoiRv.cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (env.rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(cal_fp_QoiRv.cdf(),
                        val_fp_QoiRv.cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (env.rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(cal_fp_QoiRv.cdf(),
                        val_fp_QoiRv.cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (env.rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(cal_fp_QoiRv.cdf(),
                        val_fp_QoiRv.cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (env.rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending stage 3 at " << ctime(&timevalNow.tv_sec)
              << "Total s3 run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************
  delete val_ip_ProposalCovMatrix;
  delete cal_ip_ProposalCovMatrix;

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqTgaEx4' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_TGA_EX4_H__
