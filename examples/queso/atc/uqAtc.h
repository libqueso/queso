/* uq/examples/queso/pyramid/uqAtcValidation.h
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

#ifndef __UQ_ATC_VALIDATION_H__
#define __UQ_ATC_VALIDATION_H__

#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R_CONSTANT 8.314472

// The ODE (state dot) function
int atcOdeFunc(double t, const double Mass[], double f[], void *info)
{
#if 1 // KENJI
  double* params = (double *)info;
  double A    = params[0];
  double E    = params[1];
  double beta = params[2];

  f[0] = -A*Mass[0]*exp(-E/(R_CONSTANT*t))/beta;
#endif // KENJI
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
atcLikelihoodRoutine_DataClass
{
#if 1 // KENJI
  atcLikelihoodRoutine_DataClass(const uqBaseEnvironmentClass& env,
                                 const char* experimentDescriptionFileName,
                                 const char* reactionFileName,
                                 const char* scenarioFileName1,
                                 const char* experimentalDataFileName1,
                                 const char* scenarioFileName2,
                                 const char* experimentalDataFileName2,
                                 const char* scenarioFileName3,
                                 const char* experimentalDataFileName3);
 ~atcLikelihoodRoutine_DataClass();

  // Description of experiments
  int  NSPECI, NO_STEPS, NO_PARAM;
  double  UREF, TREF, PREF, RHOREF;
  char Species_name[13][100];

  int  reac[19][13],prod[19][13];

 //  std::vector<std::string> Species_name;
  std::vector<double>      comno;
  std::vector<double>      comno_old;
  std::vector<double>      Species_con;
  std::vector<double>      Species_mass;
  std::vector<double>      Species_conm;

  // Scenario parameters
  double ex_temperature_s1;
  std::vector<double> lkc_s1;
  std::vector<double> la_s1;
  std::vector<double> p_s1;
  std::vector<double> s_s1;

  double ex_temperature_s2;
  std::vector<double> lkc_s2;
  std::vector<double> la_s2;
  std::vector<double> p_s2;
  std::vector<double> s_s2;

  // Experimental data
  std::vector<double> Time_ex1;
  std::vector<double> Ocon_ex1;
  std::vector<double> Variance_ex1;
  std::vector<double> Time_ex2;
  std::vector<double> Ocon_ex2;
  std::vector<double> Variance_ex2;

#endif // KENJI
};

template<class P_V, class P_M>
atcLikelihoodRoutine_DataClass<P_V,P_M>::atcLikelihoodRoutine_DataClass(
  const uqBaseEnvironmentClass& env,
  const char* experimentDescriptionFileName,
  const char* reactionFileName,
  const char* scenarioFileName1,
  const char* experimentalDataFileName1,
  const char* scenarioFileName2,
  const char* experimentalDataFileName2,
  const char* scenarioFileName3,
  const char* experimentalDataFileName3)
  :
#if 1 // KENJI
//  Species_name(13,""),
  comno       (13,0.),
  comno_old   (13,0.),
  Species_con (13,0.),
  Species_mass(13,0.),
  Species_conm(13,0.),
  lkc_s1      (18,0.),
  la_s1        (18,0.),
  p_s1        (18,0.),
  s_s1        (18,0.),
  Time_ex1    (4,0.),
  Ocon_ex1    (4,0.),
  Variance_ex1(4,0.),
  lkc_s2      (18,0.),
  la_s2       (18,0.),
  p_s2        (18,0.),
  s_s2        (18,0.),
  Time_ex2    (4,0.),
  Ocon_ex2    (4,0.),
  Variance_ex2(4,0.)
{
  // Read description of experiment
  char dummyNAME1[100],dummyNAME2[100];
  char dummyNAME3[100],dummyNAME4[100];
  double tmp_con, tmp_mass;
  unsigned int tmpCount = 0;

  FILE *inp = fopen(experimentDescriptionFileName,"r");
   fscanf(inp,"%s %s %s",dummyNAME1,dummyNAME2,dummyNAME3);  // NSPECI     NO_STEPS  NO_PARAM
   fscanf(inp,"%d %d %d",&NSPECI,&NO_STEPS,&NO_PARAM);                    
   fscanf(inp,"%s %s %s",dummyNAME1,dummyNAME2,dummyNAME3);  // UREF  TREF_K  PREF_Pa
   fscanf(inp,"%lf %lf %lf",&UREF, &TREF, &PREF);
   fscanf(inp,"%s %s %s",dummyNAME1,dummyNAME2,dummyNAME3);  // SPNAME  CONCO_mass_base  MWT

  while (fscanf(inp,"%s %lf %lf",Species_name[tmpCount],&tmp_con,&tmp_mass) != EOF) {
    Species_con[tmpCount]  = tmp_con;
    Species_mass[tmpCount] =tmp_mass;
    tmpCount++;
  }
  // Close file
  fclose(inp);

  int  tmp_Num,tmpREAC,tmpPROD,i,j;
  inp = fopen(reactionFileName,"r");
  for (j=0; j<NO_STEPS; j++){
       fscanf(inp,"%s",dummyNAME1);
       fscanf(inp,"%s %d %s",dummyNAME2,&tmp_Num,dummyNAME3);
       fscanf(inp,"%s",dummyNAME4);
       for (i=0; i<NSPECI; i++) {
          fscanf(inp,"%d",&tmpREAC);
          reac[j][i] = tmpREAC;
       }
       for (i=0; i<NSPECI; i++) {
          fscanf(inp,"%d",&tmpPROD);
          prod[j][i] = tmpPROD;
       }
  }
  // Close reaction file
  fclose(inp);

  if (scenarioFileName1 && experimentalDataFileName1) {
    // Read scenario parameters
     if (env.rank() == 0) {
     std::cout << "In atcLikelihoodRoutine_DataClass(), reading file '"
               << scenarioFileName1 << "'\n"
               << std::endl;
     }
       inp = fopen(scenarioFileName1,"r");    
       char dummyNAME1[100],dummyNAME2[100];
       char dummyNAME3[100],dummyNAME4[100];
       int i;
       double tmpEQC,tmpA,tmpms,tmpTs;
       fscanf(inp,"%lf",&ex_temperature_s1); 
       for (i=0; i<NO_STEPS; i++){
		fscanf(inp,"%s %s %s %s",dummyNAME1,dummyNAME2,dummyNAME3,dummyNAME4);
                fscanf(inp,"%lf %lf %lf %lf",&tmpEQC,&tmpA,&tmpms,&tmpTs);
                lkc_s1[i] = tmpEQC;
                la_s1[i]   = tmpA;
                p_s1[i]   = tmpms;
                s_s1[i]   = tmpTs;
       }

    // Close file
    fclose(inp);

    // Read experimental data
    if (env.rank() == 0) {
      std::cout << "In atcLikelihoodRoutine_DataClass(), reading file '"
                << experimentalDataFileName1 << "'\n"
                << std::endl;
    }

    inp = fopen(experimentalDataFileName1,"r");
    unsigned int numObservations = 0;
    double tmpTime;
    double tmpOcon;
    double tmpVar;

    while (fscanf(inp,"%lf %lf %lf",&tmpTime,&tmpOcon,&tmpVar) != EOF) {
      Time_ex1[numObservations] = tmpTime;
      Ocon_ex1[numObservations] = tmpOcon;
      Variance_ex1[numObservations] = tmpVar;
      numObservations++;
     printf("%lf %lf", tmpTime, tmpOcon); 
    }

    // Close file
    fclose(inp);
  }

  // Read experimental data
  if (scenarioFileName2 && experimentalDataFileName2) {
    if (env.rank() == 0) {
      std::cout << "In atcLikelihoodRoutine_DataClass(), reading file '"
                << experimentalDataFileName2 << "'\n"
                << std::endl;
    }
       inp = fopen(scenarioFileName1,"r");    
       char dummyNAME1[100],dummyNAME2[100];
       char dummyNAME3[100],dummyNAME4[100];
       int i;
       double tmpEQC,tmpA,tmpms,tmpTs;
       fscanf(inp,"%lf",&ex_temperature_s2); 
       for (i=0; i<NO_STEPS; i++){
		fscanf(inp,"%s %s %s %s",dummyNAME1,dummyNAME2,dummyNAME3,dummyNAME4);
                fscanf(inp,"%lf %lf %lf %lf",&tmpEQC,&tmpA,&tmpms,&tmpTs);
                lkc_s1[i] = tmpEQC;
                la_s2[i]  = tmpA;
                p_s2[i]   = tmpms;
                s_s2[i]   = tmpTs;
       }

    // Close file
    fclose(inp);

    // Read experimental data
    if (env.rank() == 0) {
      std::cout << "In atcLikelihoodRoutine_DataClass(), reading file '"
                << experimentalDataFileName1 << "'\n"
                << std::endl;
    }

    inp = fopen(experimentalDataFileName1,"r");
    unsigned int numObservations = 0;
    double tmpTime;
    double tmpOcon;
    double tmpVar;

    while (fscanf(inp,"%lf %lf %lf",&tmpTime,&tmpOcon,&tmpVar) != EOF) {
      Time_ex2[numObservations] = tmpTime;
      Ocon_ex2[numObservations] = tmpOcon;
      Variance_ex1[numObservations] = tmpVar;
      numObservations++;
     printf("%lf %lf", tmpTime, tmpOcon); 
    }
    // Close file
  }

  // Read experimental data
  if (scenarioFileName3 && experimentalDataFileName3) {
    if (env.rank() == 0) {
      std::cout << "In atcLikelihoodRoutine_DataClass(), reading file '"
                << experimentalDataFileName3 << "'\n"
                << std::endl;
    }
  }

#if 0 //KENJI2
  initial_setup() 
#endif 
   
#endif // KENJI
}

template<class P_V, class P_M>
atcLikelihoodRoutine_DataClass<P_V,P_M>::~atcLikelihoodRoutine_DataClass()
{
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
atcLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  double resultValue = 0.;

  // Compute likelihood for scenario 1
  double la_param  = paramValues[0];    // 19th Reaction parameters  Log_10 (A)
  double p_param   = paramValues[1];    // 19th Reaction parameters  T^p
  double s_param   = paramValues[2];    // 19th Reaction parameters  exp(-s/T)
  double lkc_param = paramValues[3];    // 19th Reaction parameters  K_b=k_f/Kc

#if 0 // KENJI2
  reactionRoutine();
#endif

#if 0 // KENJI
    double variance               = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_variance1;
    const std::vector<double>& Te = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Te1;
    const std::vector<double>& Me = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Me1;
    std::vector<double> Mt(Me.size(),0.);

    double params[]={A,E,beta};
      	
    // integration
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
          gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
          gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
          gsl_odeiv_system     sys = {atcOdeFunc, NULL, 1, (void *)params}; 

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
                          "atcLikelihoodRoutine()",
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
      printf("In atcLikelihoodRoutine(), A = %g, E = %g, beta = %.3lf: misfit = %lf, likelihood = %lf.\n",A,E,beta,misfit,resultValue);
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
#endif // KENJI

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
atcQoiRoutine_DataClass
{
#if 1 // KENJI
  double m_beta;
  double m_criticalMass;
  double m_criticalTime;
#endif // KENJI
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void atcQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
#if 0 // KENJI
  double A            = paramValues[0];
  double E            = paramValues[1];
  double beta         = ((atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_beta;
  double criticalMass = ((atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_criticalMass;
  double criticalTime = ((atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_criticalTime;

  double params[]={A,E,beta};
      	
  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
        gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
        gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sys = {atcOdeFunc, NULL, 1, (void *)params}; 
	
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
                        "atcQoiRoutine()",
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
    printf("In atcQoiRoutine(), A = %g, E = %g, beta = %.3lf, criticalTime = %.3lf, criticalMass = %.3lf: qoi = %lf.\n",A,E,beta,criticalTime,criticalMass,qoiValues[0]);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);
#endif // KENJI
  return;
}

//********************************************************
// Class "uqAtcValidation", instantiated by main()
//********************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
class uqAtcValidationClass : public uqModelValidationClass<P_V,P_M,Q_V,Q_M>
{
public:
  uqAtcValidationClass(const uqBaseEnvironmentClass& env,
                       const char*               prefix);
 ~uqAtcValidationClass();

  void run();

private:
  void  runCalibrationStage();
  void  runValidationStage();
  void  runComparisonStage();

  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_cycle;

  uqAsciiTableClass<P_V,P_M>*               m_paramsTable;
  const EpetraExt::DistArray<std::string>*  m_paramNames;         // instantiated outside this class!!
  P_V*                                      m_paramMinValues;     // instantiated outside this class!!
  P_V*                                      m_paramMaxValues;     // instantiated outside this class!!
  P_V*                                      m_paramInitialValues; // instantiated outside this class!!
  uqVectorSpaceClass<P_V,P_M>*              m_paramSpace;
  uqVectorSetClass<P_V,P_M>*                m_paramDomain;

  uqAsciiTableClass<P_V,P_M>*               m_qoisTable;
  const EpetraExt::DistArray<std::string>*  m_qoiNames; // instantiated outside this class!!
  uqVectorSpaceClass<Q_V,Q_M>*              m_qoiSpace;

#if 1 // KENJI
  double                                    m_predBeta;
  double                                    m_predCriticalMass;
  double                                    m_predCriticalTime;
#endif

  uqBaseVectorRVClass<P_V,P_M>*             m_calPriorRv;
  atcLikelihoodRoutine_DataClass<P_V,P_M>*  m_calLikelihoodRoutine_Data;
  uqBaseScalarFunctionClass<P_V,P_M>*       m_calLikelihoodFunctionObj;
  atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>* m_calQoiRoutine_Data;

  atcLikelihoodRoutine_DataClass<P_V,P_M>*  m_valLikelihoodRoutine_Data;
  uqBaseScalarFunctionClass<P_V,P_M>*       m_valLikelihoodFunctionObj;
  atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>* m_valQoiRoutine_Data;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqAtcValidationClass<P_V,P_M,Q_V,Q_M>::uqAtcValidationClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  uqModelValidationClass<P_V,P_M,Q_V,Q_M>(env,prefix),
  m_paramsTable              (NULL),
  m_paramNames               (NULL),
  m_paramMinValues           (NULL),
  m_paramMaxValues           (NULL),
  m_paramInitialValues       (NULL),
  m_paramSpace               (NULL),
  m_paramDomain              (NULL),
  m_qoisTable                (NULL),
  m_qoiNames                 (NULL),
  m_qoiSpace                 (NULL),
  m_calPriorRv               (NULL),
  m_calLikelihoodRoutine_Data(NULL),
  m_calLikelihoodFunctionObj (NULL),
  m_calQoiRoutine_Data       (NULL),
  m_valLikelihoodRoutine_Data(NULL),
  m_valLikelihoodFunctionObj (NULL),
  m_valQoiRoutine_Data       (NULL)
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqAtcValidation::constructor()\n"
              << std::endl;
  }

  // Read Ascii file with information on parameters
  m_paramsTable = new uqAsciiTableClass<P_V,P_M> (m_env,
                                                  4,    // # of rows
                                                  3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                                  NULL, // All extra columns are of 'double' type
                                                  "params2.tab");

  m_paramNames = &(m_paramsTable->stringColumn(0));
  m_paramMinValues     = new P_V(m_paramsTable->doubleColumn(1));
  m_paramMaxValues     = new P_V(m_paramsTable->doubleColumn(2));
  m_paramInitialValues = new P_V(m_paramsTable->doubleColumn(3));

  m_paramSpace = new uqVectorSpaceClass<P_V,P_M>(m_env,
                                                 "param_", // Extra prefix before the default "space_" prefix
                                                 m_paramsTable->numRows(),
                                                 m_paramNames);

  m_paramDomain = new uqBoxSubsetClass<P_V,P_M>("param_",
                                                *m_paramSpace,
                                                *m_paramMinValues,
                                                *m_paramMaxValues);

  // Read Ascii file with information on qois
  m_qoisTable = new uqAsciiTableClass<P_V,P_M>(m_env,
                                               1,    // # of rows
                                               0,    // # of cols after 'parameter name': none
                                               NULL, // All extra columns are of 'double' type
                                               "qois.tab");

  m_qoiNames = &(m_qoisTable->stringColumn(0));

  m_qoiSpace = new uqVectorSpaceClass<Q_V,Q_M>(m_env,
                                               "qoi_", // Extra prefix before the default "space_" prefix
                                               m_qoisTable->numRows(),
                                               m_qoiNames);

  // Instantiate the validation cycle
  m_cycle = new uqValidationCycleClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                        m_prefix.c_str(), // Use the prefix passed above
                                                        *m_paramSpace,
                                                        *m_qoiSpace);

#if 1 // KENJI
  m_predBeta         = 250.;
  m_predCriticalMass = 0.;
  m_predCriticalTime = 3.9;
#endif // KENJI

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqAtcValidation::constructor()\n"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqAtcValidationClass<P_V,P_M,Q_V,Q_M>::~uqAtcValidationClass()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqAtcValidation::destructor()"
              << std::endl;
  }

  if (m_valQoiRoutine_Data)        delete m_valQoiRoutine_Data;
  if (m_valLikelihoodFunctionObj)  delete m_valLikelihoodFunctionObj;
  if (m_valLikelihoodRoutine_Data) delete m_valLikelihoodRoutine_Data;
  if (m_calQoiRoutine_Data)        delete m_calQoiRoutine_Data;
  if (m_calLikelihoodFunctionObj)  delete m_calLikelihoodFunctionObj;
  if (m_calLikelihoodRoutine_Data) delete m_calLikelihoodRoutine_Data;
  if (m_calPriorRv)                delete m_calPriorRv;
  if (m_qoiSpace)                  delete m_qoiSpace;
//if (m_qoiNames)                  delete m_qoiNames; // instantiated outside this class!!
  if (m_qoisTable)                 delete m_qoisTable;
  if (m_paramDomain)               delete m_paramDomain;
  if (m_paramSpace)                delete m_paramSpace;
//if (m_paramInitialValues)        delete m_paramInitialValues; // instantiated outside this class!!
//if (m_paramMaxValues)            delete m_paramMaxValues;     // instantiated outside this class!!
//if (m_paramMinValues)            delete m_paramMinValues;     // instantiated outside this class!!
//if (m_paramNames)                delete m_paramNames;         // instantiated outside this class!!
  if (m_paramsTable)               delete m_paramsTable;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqAtcValidation::destructor()"
              << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqAtcValidationClass<P_V,P_M,Q_V,Q_M>::run()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqAtcValidation::run()"
              << std::endl;
  }
  
  printf("Hello Kenji3\n");
  runCalibrationStage();
  runValidationStage();
  //runComparisonStage();

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqAtcValidation::run()"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAtcValidationClass<P_V,P_M,Q_V,Q_M>::runCalibrationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqAtcValidation::runCalibrationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_calPriorRv = new uqUniformVectorRVClass<P_V,P_M> ("cal_prior_", // Extra prefix before the default "rv_" prefix
                                                      *m_paramDomain);

  m_calLikelihoodRoutine_Data = new atcLikelihoodRoutine_DataClass<P_V,P_M>(m_env,
                                                                            "input.dat",
                                                                            "reaction.dat", 
                                                                            "scenario_001.dat",
                                                                            "Experimental_001.dat",
                                                                            NULL,
                                                                            NULL,
                                                                            NULL,
                                                                            NULL);

  m_calLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("cal_like_",
                                                                         *m_paramDomain,
                                                                         atcLikelihoodRoutine<P_V,P_M>,
                                                                         NULL,
                                                                         NULL,
                                                                         (void *) m_calLikelihoodRoutine_Data,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setCalIP(*m_calPriorRv,
                    *m_calLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* calProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().priorRv().pdf().domainVarVector(),
                                                                                                  *m_paramInitialValues);
  m_cycle->calIP().solveWithBayesMarkovChain(*m_paramInitialValues,
                                             calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Deal with forward problem
  m_calQoiRoutine_Data = new atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>();
  m_calQoiRoutine_Data->m_beta         = m_predBeta;
  m_calQoiRoutine_Data->m_criticalMass = m_predCriticalMass;
  m_calQoiRoutine_Data->m_criticalTime = m_predCriticalTime;

  m_cycle->setCalFP(atcQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_calQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->calFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqAtcValidation::runCalibrationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqAtcValidation::runCalibrationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAtcValidationClass<P_V,P_M,Q_V,Q_M>::runValidationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqAtcValidation::runValidationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_valLikelihoodRoutine_Data = new atcLikelihoodRoutine_DataClass<P_V,P_M>(m_env,
                                                                            "input.dat",
                                                                            "reaction.dat",
                                                                            "scenario_002.dat",
                                                                            "Experimental_002.dat",
                                                                            NULL,
                                                                            NULL,
                                                                            NULL,
                                                                            NULL);

  m_valLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("val_like_",
                                                                         *m_paramDomain,
                                                                         atcLikelihoodRoutine<P_V,P_M>,
                                                                         NULL,
                                                                         NULL,
                                                                         (void *) m_valLikelihoodRoutine_Data,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setValIP(*m_valLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* valProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().postRv().realizer().imageVarVector(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                                   m_cycle->calIP().postRv().realizer().imageExpVector()); // Use these values as the initial values
  m_cycle->valIP().solveWithBayesMarkovChain(m_cycle->calIP().postRv().realizer().imageExpVector(),
                                             valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Deal with forward problem
  m_valQoiRoutine_Data = new atcQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>();
  m_valQoiRoutine_Data->m_beta         = m_predBeta;
  m_valQoiRoutine_Data->m_criticalMass = m_predCriticalMass;
  m_valQoiRoutine_Data->m_criticalTime = m_predCriticalTime;

  m_cycle->setValFP(atcQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_valQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->valFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqAtcValidation::runValidationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqAtcValidation::runValidationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAtcValidationClass<P_V,P_M,Q_V,Q_M>::runComparisonStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqAtcValidation::runComparisonStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  if (m_cycle->calFP().computeSolutionFlag() &&
      m_cycle->valFP().computeSolutionFlag()) {
    Q_V* epsilonVec = m_cycle->calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    Q_V cdfDistancesVec(m_cycle->calFP().qoiRv().imageSet().vectorSpace().zeroVector());
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(m_cycle->valFP().qoiRv().cdf(),
                        m_cycle->calFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqAtcValidation::runComparisonStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqAtcValidation::runComparisonStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}
#endif // __UQ_ATC_VALIDATION_H__
