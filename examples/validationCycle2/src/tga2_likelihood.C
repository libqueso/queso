/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <tga2_func.h>
#include <tga2_likelihood.h>
#include <queso/Environment.h>
#include <gsl/gsl_odeiv.h>

likelihoodRoutine_Data::likelihoodRoutine_Data(
  const QUESO::BaseEnvironment& env,
  const char* inpName1,
  const char* inpName2,
  const char* inpName3)
  :
  m_beta1    (0.),
  m_variance1(0.),
  m_Te1      (0),
  m_Me1      (0),
  m_beta2    (0.),
  m_variance2(0.),
  m_Te2      (0),
  m_Me2      (0),
  m_beta3    (0.),
  m_variance3(0.),
  m_Te3      (0),
  m_Me3      (0),
  m_env      (&env)
{
  // Read experimental data
  if (inpName1) {
    m_Te1.resize(11,0.);
    m_Me1.resize(11,0.);

    // Open input file on experimental data
    FILE *inp;
    inp = fopen(inpName1,"r");

    queso_require_msg(inp, "Unable to open " << inpName1);

    // Read kinetic parameters and convert heating rate to K/s
    int aux1 = fscanf(inp,"%lf %lf",&m_beta1,&m_variance1);
    m_beta1 /= 60.;

    if(aux1) {}; // just to eliminate warnings

    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= m_Te1.size()),
                          env.fullRank(),
                          "uqAppl(), in uqTgaEx4.h",
                          "input file 1 has too many observations");
      m_Te1[numObservations] = tmpTe;
      m_Me1[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != m_Te1.size()),
                        env.fullRank(),
                        "uqAppl(), in uqTgaEx4.h",
                        "input file 1 has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }

  // Read experimental data
  if (inpName2) {
    m_Te2.resize(11,0.);
    m_Me2.resize(11,0.);

    // Open input file on experimental data
    FILE *inp;
    inp = fopen(inpName2,"r");

    queso_require_msg(inp, "Unable to open " << inpName2);

    // Read kinetic parameters and convert heating rate to K/s
    int aux2 = fscanf(inp,"%lf %lf",&m_beta2,&m_variance2);
    m_beta2 /= 60.;

    if(aux2) {}; // just to eliminate warnings

    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= m_Te2.size()),
                          env.fullRank(),
                          "uqAppl(), in uqTgaEx4.h",
                          "input file 2 has too many observations");
      m_Te2[numObservations] = tmpTe;
      m_Me2[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != m_Te2.size()),
                        env.fullRank(),
                        "uqAppl(), in uqTgaEx4.h",
                        "input file 2 has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }

  // Read experimental data
  if (inpName3) {
    m_Te3.resize(11,0.);
    m_Me3.resize(11,0.);

    // Open input file on experimental data
    FILE *inp;
    inp = fopen(inpName3,"r");

    // Read kinetic parameters and convert heating rate to K/s
    int aux3 = fscanf(inp,"%lf %lf",&m_beta3,&m_variance3);
    m_beta3 /= 60.;

    if(aux3) {}; // just to eliminate warnings

    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= m_Te3.size()),
                          env.fullRank(),
                          "uqAppl(), in uqTgaEx4.h",
                          "input file 3 has too many observations");
      m_Te3[numObservations] = tmpTe;
      m_Me3[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != m_Te3.size()),
                        env.fullRank(),
                        "uqAppl(), in uqTgaEx4.h",
                        "input file 3 has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }
}

likelihoodRoutine_Data::~likelihoodRoutine_Data()
{
}

//********************************************************
// The actual (user defined) likelihood routine
//********************************************************
double
likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  double resultValue = 0.;

  const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);

  if (paramDirection  &&
      functionDataPtr &&
      gradVector      &&
      hessianMatrix   &&
      hessianEffect) {
    // Just to eliminate INTEL compiler warnings
  }

  env.subComm().Barrier();
  //env.syncPrintDebugMsg("Entering likelihoodRoutine()",1,env.fullComm());

  // Compute likelihood for scenario 1
  double betaTest = ((likelihoodRoutine_Data*) functionDataPtr)->m_beta1;
  if (betaTest) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    = ((likelihoodRoutine_Data*) functionDataPtr)->m_beta1;
    double variance                = ((likelihoodRoutine_Data*) functionDataPtr)->m_variance1;
    const std::vector<double>& Te  = ((likelihoodRoutine_Data*) functionDataPtr)->m_Te1;
    const std::vector<double>& Me  = ((likelihoodRoutine_Data*) functionDataPtr)->m_Me1;
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
                          paramValues.env().fullRank(),
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
    if ((paramValues.env().displayVerbosity() >= 10) && (paramValues.env().fullRank() == 0)) {
      printf("In likelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  // Compute likelihood for scenario 2
  betaTest = ((likelihoodRoutine_Data*) functionDataPtr)->m_beta2;
  if (betaTest > 0.) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    = ((likelihoodRoutine_Data*) functionDataPtr)->m_beta2;
    double variance                = ((likelihoodRoutine_Data*) functionDataPtr)->m_variance2;
    const std::vector<double>& Te  = ((likelihoodRoutine_Data*) functionDataPtr)->m_Te2;
    const std::vector<double>& Me  = ((likelihoodRoutine_Data*) functionDataPtr)->m_Me2;
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
                          paramValues.env().fullRank(),
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
    if ((paramValues.env().displayVerbosity() >= 10) && (paramValues.env().fullRank() == 0)) {
      printf("In likelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  // Compute likelihood for scenario 3
  betaTest = ((likelihoodRoutine_Data*) functionDataPtr)->m_beta3;
  if (betaTest > 0.) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    = ((likelihoodRoutine_Data*) functionDataPtr)->m_beta3;
    double variance                = ((likelihoodRoutine_Data*) functionDataPtr)->m_variance3;
    const std::vector<double>& Te  = ((likelihoodRoutine_Data*) functionDataPtr)->m_Te3;
    const std::vector<double>& Me  = ((likelihoodRoutine_Data*) functionDataPtr)->m_Me3;
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
                          paramValues.env().fullRank(),
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
    if ((paramValues.env().displayVerbosity() >= 10) && (paramValues.env().fullRank() == 0)) {
      printf("In likelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  env.subComm().Barrier();
  //env.syncPrintDebugMsg("Leaving likelihoodRoutine()",1,env.fullComm());

  return -.5*resultValue;
}

