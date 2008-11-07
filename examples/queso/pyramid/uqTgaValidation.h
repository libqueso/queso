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
#include <uqAsciiTable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472

// The ODE (state dot) function
int tgaOdeFunc(double t, const double Mass[], double f[], void *info)
{
  double* params = (double *)info;
  double A    = params[0];
  double E    = params[1];
  double beta = params[2];

  f[0] = -A*Mass[0]*exp(-E/(R_CONSTANT*t))/beta;

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
  tgaLikelihoodRoutine_DataClass(const uqEnvironmentClass& env,
                                 const char* inpName1,
                                 const char* inpName2,
                                 const char* inpName3);
 ~tgaLikelihoodRoutine_DataClass();

  double              m_beta1;
  double              m_variance1;
  std::vector<double> m_Te1; // temperatures
  std::vector<double> m_Me1; // relative masses

  double              m_beta2;
  double              m_variance2;
  std::vector<double> m_Te2; // temperatures
  std::vector<double> m_Me2; // relative masses

  double              m_beta3;
  double              m_variance3;
  std::vector<double> m_Te3; // temperatures
  std::vector<double> m_Me3; // relative masses
};

template<class P_V, class P_M>
tgaLikelihoodRoutine_DataClass<P_V,P_M>::tgaLikelihoodRoutine_DataClass(
  const uqEnvironmentClass& env,
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
  m_Me3      (0)
{
  // Read experimental data
  if (inpName1) {
    m_Te1.resize(11,0.);
    m_Me1.resize(11,0.);

    // Open input file on experimental data
    FILE *inp;
    inp = fopen(inpName1,"r");

    // Read kinetic parameters and convert heating rate to K/s
    fscanf(inp,"%lf %lf",&m_beta1,&m_variance1);
    m_beta1 /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= m_Te1.size()),
                          env.rank(),
                          "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
                          "input file 1 has too many observations");
      m_Te1[numObservations] = tmpTe;
      m_Me1[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != m_Te1.size()),
                        env.rank(),
                        "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
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

    // Read kinetic parameters and convert heating rate to K/s
    fscanf(inp,"%lf %lf",&m_beta2,&m_variance2);
    m_beta2 /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= m_Te2.size()),
                          env.rank(),
                          "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
                          "input file 2 has too many observations");
      m_Te2[numObservations] = tmpTe;
      m_Me2[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != m_Te2.size()),
                        env.rank(),
                        "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
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
    fscanf(inp,"%lf %lf",&m_beta3,&m_variance3);
    m_beta3 /= 60.;
  
    unsigned int numObservations = 0;
    double tmpTe;
    double tmpMe;
    while (fscanf(inp,"%lf %lf",&tmpTe,&tmpMe) != EOF) {
      UQ_FATAL_TEST_MACRO((numObservations >= m_Te3.size()),
                          env.rank(),
                          "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
                          "input file 3 has too many observations");
      m_Te3[numObservations] = tmpTe;
      m_Me3[numObservations] = tmpMe;
      numObservations++;
    }
    UQ_FATAL_TEST_MACRO((numObservations != m_Te3.size()),
                        env.rank(),
                        "tgaLikelihoodRoutine_DataClass(), in uqTgaValidation.h",
                        "input file 3 has a smaller number of observations than expected");

    // Close input file on experimental data
    fclose(inp);
  }
}

template<class P_V, class P_M>
tgaLikelihoodRoutine_DataClass<P_V,P_M>::~tgaLikelihoodRoutine_DataClass()
{
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
tgaLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  double resultValue = 0.;

  // Compute likelihood for scenario 1
  double betaTest = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_beta1;
  if (betaTest) {
    double A                      = paramValues[0];
    double E                      = paramValues[1];
    double beta                   = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_beta1;
    double variance               = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_variance1;
    const std::vector<double>& Te = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Te1;
    const std::vector<double>& Me = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Me1;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {tgaOdeFunc, NULL, 1, (void *)params}; 

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
                          "tgaLikelihoodRoutine()",
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
      printf("In tgaLikelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  // Compute likelihood for scenario 2
  betaTest = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_beta2;
  if (betaTest > 0.) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_beta2;
    double variance                = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_variance2;
    const std::vector<double>& Te  = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Te2;
    const std::vector<double>& Me  = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Me2;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {tgaOdeFunc, NULL, 1, (void *)params}; 

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
                          "tgaLikelihoodRoutine()",
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
      printf("In tgaLikelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  // Compute likelihood for scenario 3
  betaTest = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_beta3;
  if (betaTest > 0.) {
    double A                       = paramValues[0];
    double E                       = paramValues[1];
    double beta                    = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_beta3;
    double variance                = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_variance3;
    const std::vector<double>& Te  = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Te3;
    const std::vector<double>& Me  = ((tgaLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Me3;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {tgaOdeFunc, NULL, 1, (void *)params}; 

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
                          "tgaLikelihoodRoutine()",
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
      printf("In tgaLikelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
  }

  return resultValue;
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
  double m_beta;
  double m_criticalMass;
  double m_criticalTime;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void tgaQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  double A            = paramValues[0];
  double E            = paramValues[1];
  double beta         = ((tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_beta;
  double criticalMass = ((tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_criticalMass;
  double criticalTime = ((tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_criticalTime;

  double params[]={A,E,beta};
      	
  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
        gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
        gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sys = {tgaOdeFunc, NULL, 1, (void *)params}; 
	
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
                        "tgaQoiRoutine()",
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
    printf("In tgaQoiRoutine(), A = %g, E = %g, beta = %.3lf, criticalTime = %.3lf, criticalMass = %.3lf: qoi = %lf.\n",A,E,beta,criticalTime,criticalMass,qoiValues[0]);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

//********************************************************
// The 'comparison stage' of the validation cycle
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void 
tgaComparisonStage(uqValidationCycleClass<P_V,P_M,Q_V,Q_M>& cycle)
{
  if (cycle.calFP().computeSolutionFlag() &&
      cycle.valFP().computeSolutionFlag()) {
    Q_V* epsilonVec = cycle.calFP().qoiRv().imageSpace().newVector(0.02);
    Q_V cdfDistancesVec(cycle.calFP().qoiRv().imageSpace().zeroVector());
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(cycle.valFP().qoiRv().cdf(),
                        cycle.calFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  return;
}

//********************************************************
// The class related to "TGA validation", instantiated by main()
//********************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
class uqTgaValidationClass : public uqModelValidationClass<P_V,P_M,Q_V,Q_M>
{
public:
  uqTgaValidationClass(const uqEnvironmentClass& env,
                       const char*               prefix);
 ~uqTgaValidationClass();

  void run();

private:
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_cycle;

  uqAsciiTableClass<P_V,P_M>*               m_paramsTable;
  const EpetraExt::DistArray<std::string>*  m_paramNames;
  P_V*                                      m_paramMinValues;
  P_V*                                      m_paramMaxValues;
  P_V*                                      m_paramInitialValues;
  uqVectorSpaceClass<P_V,P_M>*              m_paramSpace;

  uqAsciiTableClass<P_V,P_M>*               m_qoisTable;
  const EpetraExt::DistArray<std::string>*  m_qoiNames;
  uqVectorSpaceClass<Q_V,Q_M>*              m_qoiSpace;

  double                                    m_predBeta;
  double                                    m_predCriticalMass;
  double                                    m_predCriticalTime;

  uqUniformVectorRVClass<P_V,P_M>*          m_calPriorRv;
  tgaLikelihoodRoutine_DataClass<P_V,P_M>*  m_calLikelihoodRoutine_Data;
  tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>* m_calQoiRoutine_Data;

  tgaLikelihoodRoutine_DataClass<P_V,P_M>*  m_valLikelihoodRoutine_Data;
  tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>* m_valQoiRoutine_Data;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::uqTgaValidationClass(
  const uqEnvironmentClass& env,
  const char*               prefix)
  :
  uqModelValidationClass<P_V,P_M,Q_V,Q_M>(env,prefix)
{
  if (m_env.rank() == 0) {
    std::cout << "Beginning run of 'uqTgaValidation' example\n"
              << std::endl;
  }

  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  //******************************************************
  // TGA validation cycle: instantiation
  //******************************************************

  // Read Ascii file with information on parameters
  m_paramsTable = new uqAsciiTableClass<P_V,P_M> (m_env,
                                                  2,    // # of rows
                                                  3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                                  NULL, // All extra columns are of 'double' type
                                                  "params.tab");

  m_paramNames = &(m_paramsTable->stringColumn(0));
  m_paramMinValues     = new P_V(m_paramsTable->doubleColumn(1));
  m_paramMaxValues     = new P_V(m_paramsTable->doubleColumn(2));
  m_paramInitialValues = new P_V(m_paramsTable->doubleColumn(3));

  uqVectorSpaceClass<P_V,P_M> paramSpace(m_env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         m_paramsTable->numRows(),
                                         m_paramNames);

  // Read Ascii file with information on qois
  uqAsciiTableClass<P_V,P_M> qoisTable(m_env,
                                       1,    // # of rows
                                       0,    // # of cols after 'parameter name': none
                                       NULL, // All extra columns are of 'double' type
                                       "qois.tab");

  const EpetraExt::DistArray<std::string>& qoiNames = qoisTable.stringColumn(0);

  uqVectorSpaceClass<Q_V,Q_M> qoiSpace(m_env,
                                       "qoi_", // Extra prefix before the default "space_" prefix
                                       qoisTable.numRows(),
                                       &qoiNames);

  // Instantiate the validation cycle
  uqValidationCycleClass<P_V,P_M,Q_V,Q_M> cycle(m_env,
                                                "", // No extra prefix
                                                paramSpace,
                                                qoiSpace);

  double m_predBeta         = 250.;
  double m_predCriticalMass = 0.;
  double m_predCriticalTime = 3.9;

  //********************************************************
  // TGA validation cycle: calibration stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Beginning 'calibration stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  uqUniformVectorRVClass<P_V,P_M> calPriorRv("cal_prior_", // Extra prefix before the default "rv_" prefix
                                             paramSpace,
                                             *m_paramMinValues,
                                             *m_paramMaxValues);

  tgaLikelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(m_env,
                                                                    "scenario_5_K_min.dat",
                                                                    "scenario_25_K_min.dat",
                                                                    "scenario_50_K_min.dat");

  cycle.setCalIP(calPriorRv,
                 tgaLikelihoodRoutine<P_V,P_M>,
                 (void *) &calLikelihoodRoutine_Data,
                 true); // the likelihood routine computes [-2.*ln(Likelihood)]

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* calProposalCovMatrix = cycle.calIP().postRv().imageSpace().newGaussianMatrix(cycle.calIP().priorRv().pdf().domainVarianceValues(),
                                                                                    *m_paramInitialValues);
  cycle.calIP().solveWithBayesMarkovChain(*m_paramInitialValues,
                                          *calProposalCovMatrix,
                                          NULL); // use default kernel from library
  delete calProposalCovMatrix;

  // Deal with forward problem
  tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> calQoiRoutine_Data;
  calQoiRoutine_Data.m_beta         = m_predBeta;
  calQoiRoutine_Data.m_criticalMass = m_predCriticalMass;
  calQoiRoutine_Data.m_criticalTime = m_predCriticalTime;

  cycle.setCalFP(tgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
                 (void *) &calQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  cycle.calFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Ending 'calibration stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'calibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //********************************************************
  // TGA validation cycle: validation stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Beginning 'validation stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  tgaLikelihoodRoutine_DataClass<P_V,P_M> valLikelihoodRoutine_Data(m_env,
                                                                    "scenario_100_K_min.dat",
                                                                    NULL,
                                                                    NULL);

  cycle.setValIP(tgaLikelihoodRoutine<P_V,P_M>,
                 (void *) &valLikelihoodRoutine_Data,
                 true); // the likelihood routine computes [-2.*ln(Likelihood)]

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* valProposalCovMatrix = cycle.calIP().postRv().imageSpace().newGaussianMatrix(cycle.calIP().postRv().realizer().imageVarianceValues(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                    cycle.calIP().postRv().realizer().imageExpectedValues()); // Use these values as the initial values
  cycle.valIP().solveWithBayesMarkovChain(cycle.calIP().postRv().realizer().imageExpectedValues(),
                                          *valProposalCovMatrix,
                                          NULL); // use default kernel from library
  delete valProposalCovMatrix;

  // Deal with forward problem
  tgaQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> valQoiRoutine_Data;
  valQoiRoutine_Data.m_beta         = m_predBeta;
  valQoiRoutine_Data.m_criticalMass = m_predCriticalMass;
  valQoiRoutine_Data.m_criticalTime = m_predCriticalTime;

  cycle.setValFP(tgaQoiRoutine<P_V,P_M,Q_V,Q_M>,
                 (void *) &valQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  cycle.valFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Ending 'validation stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'validation stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //********************************************************
  // TGA validation cycle: comparison stage
  //********************************************************

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Beginning 'comparison stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  tgaComparisonStage(cycle);

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Ending 'comparison stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'comparison stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (m_env.rank() == 0) {
    std::cout << "Finishing run of 'uqTgaValidation' example"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::~uqTgaValidationClass()
{
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::run()
{
  return;
}
#endif // __UQ_TGA_VALIDATION_H__
