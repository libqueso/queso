//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

using namespace std;

#include <basic_classes.h>
#include <basic_int.h>
#include <grvy.h>

namespace QUESO_Basic_API {

  char log_prefix[] = "queso_basic";

  //------------------
  // Member Functions
  //------------------

  QUESO_Basic_Class::QUESO_Basic_Class()
  {
    m_initialized           = 0;		
    m_silent                = 0; 		
    m_paramSpace            = NULL;
    m_paramDomain           = NULL;
    m_user_likelihood_func  = NULL;
    m_env                   = NULL; // prudenci 2010-05-17
  }

  QUESO_Basic_Class::~QUESO_Basic_Class()
  {
    if (m_env) delete m_env; // prudenci 2010-05-17
  }

  void QUESO_Basic_Class::Initialize(const char *inputfile)
  {

    grvy_timer_init("QUESO");

    // Define new QUESO environment and store inputfile information

    m_env         = new uqFullEnvironmentClass(MPI_COMM_WORLD,inputfile,"",NULL);
    m_inputfile   = new string(inputfile);

    m_initialized = 1;

}

  void QUESO_Basic_Class:: DefineParameterSpace()
  {
    double *param_min;		// min value of each parameter
    double *param_max;		// max value of each parameter
    double *param_ini;		// initial value of each parameter

    int     ierr = 1;

    VerifyInit();

    // Verify presence of required input file

    ierr *= grvy_input_fopen(m_inputfile->c_str());

    if(ierr == 0)
      QUESO_fatal("Unable to access QUESO input file");

    // Derive the UQ parameters from input file and create QUESO parameter vector

    ierr *= grvy_input_fread_int("queso/parameters/num_params",&m_num_params);

    if(ierr == 0)
      QUESO_fatal("Unable to read num_params from input file.");

    grvy_printf(GRVY_INFO,"%s: Num params defined\n",log_prefix,m_num_params);
    
    printf("--> Setup parameter space variables/ranges...\n");
    m_paramSpace  = new uqVectorSpaceClass <basicV,basicM> (*m_env,"queso_basic_",m_num_params,NULL);
    
    param_min = (double *) calloc(m_num_params,sizeof(double));
    param_max = (double *) calloc(m_num_params,sizeof(double));
    param_ini = (double *) calloc(m_num_params,sizeof(double));
  
    if(param_min == NULL || param_max == NULL || param_ini == NULL)
      QUESO_fatal("Unable to allocate memory for desired parameter space");
  
    ierr *= grvy_input_fread_double_vec("queso/parameters/param_mins",  param_min, m_num_params);
    ierr *= grvy_input_fread_double_vec("queso/parameters/param_maxs",  param_max, m_num_params);
    ierr *= grvy_input_fread_double_vec("queso/parameters/param_inits", param_ini, m_num_params);

    printf("max = %f\n",param_max[0]);

    if(ierr == 0)
      QUESO_fatal("Unable to read parameter min/max/init values");

    for(int i = 0;i<m_num_params;i++)
      {
	grvy_printf(GRVY_INFO,"%s: min  value = \n",log_prefix,param_min[i]);
	grvy_printf(GRVY_INFO,"%s: max  value = \n",log_prefix,param_max[i]);
	grvy_printf(GRVY_INFO,"%s: init value = \n",log_prefix,param_ini[i]);
      }

    m_queso_var_min = new basicV ( m_paramSpace->zeroVector() );
    m_queso_var_max = new basicV ( m_paramSpace->zeroVector() );
    m_queso_var_ini = new basicV ( m_paramSpace->zeroVector() );
    
    for(int i = 0;i<m_num_params;i++)
      {
	(*m_queso_var_min)[i] = param_min[i];
	(*m_queso_var_max)[i] = param_max[i];
	(*m_queso_var_ini)[i] = param_ini[i];
      }
    
    printf("--> Setup parameter search box...\n");
    m_paramDomain = new uqBoxSubsetClass<basicV, basicM> ("queso_basic_",*m_paramSpace,
							  *m_queso_var_min,*m_queso_var_max);
    // un poquito de clean up.

    free(param_min);
    free(param_max);
    free(param_ini);

    grvy_input_fclose();

  }

  void QUESO_Basic_Class:: VerifyInit()
  {
    if(!m_initialized)
      QUESO_fatal("QUESO not initialized prior to use");
  }

  void QUESO_Basic_Class::Likelihood_Register(double (*fp)(double *) )
  {

    VerifyInit();

    printf("--> Setting prior and post vectors...\n");

    m_priorRV = new uqUniformVectorRVClass <basicV, basicM> ("prior_",*m_paramDomain);
    m_postRV  = new uqGenericVectorRVClass <basicV, basicM> ("post_" ,*m_paramSpace );
    
    m_likelihoodObj = new uqGenericScalarFunctionClass 
      <basicV,basicM> ("like_",*m_paramDomain,Likelihood_Wrapper,NULL,true);
    
    m_user_likelihood_func = fp;
    
    // define the inverse (calibration) problem
    
    printf("--> Defining inverse problem...\n");
    
    m_ip = new uqStatisticalInverseProblemClass <basicV,basicM> ("",NULL,*m_priorRV,*m_likelihoodObj,*m_postRV);
    
    // Default covariance matrix for now - default assumption assumes
    // 6-sigma distribution range falls over 1/3 of the max parameter range.
    //
    // koomie TODO: store relevant constants regarding this default
    // assumption in input file and register as default values, but
    // allow the savvy user to override
    
    double cov_param = 1/(3.*6.);	// 6-sigma over 1/3 of the range
    double param_range = 0.0;
    
    printf("--> Defining default covariance matrix...\n");
    
    m_CovMatrix = m_postRV->imageSet().vectorSpace().newProposalMatrix(NULL,
								       m_queso_var_ini);

    for(int i=0;i<m_num_params;i++)
      {
	param_range = (*m_queso_var_max)[i] - (*m_queso_var_min)[i];
	(*m_CovMatrix)(i,i) = (cov_param*param_range)*(cov_param*param_range);
      }
  }

  void QUESO_Basic_Class::SolveInverseProblem()
  {
    VerifyInit();  

    // Launch the Markov Chain 

    m_ip->solveWithBayesMetropolisHastings(NULL,*m_queso_var_ini,m_CovMatrix);
  }

  //----------------------------------------------
  // Fatal error(): example only - to be replaced.
  //----------------------------------------------

  void QUESO_fatal(const char *message)
  {
    printf("%s\n",message);
    exit(1);

    // koomie note: likely need mpi_abort and better exception handling aqui.
    // koomie TODO: switch to queso error macro
  }

  //---------------------------------------------------------------------
  // Wrapper for user likelihood routine: culls data from QUESO uqVector
  // and passes to the user routine as an array of doubles.
  //---------------------------------------------------------------------

  double Likelihood_Wrapper(const basicV &paramValue,
			    const basicV *paramDirection,
			    const void *Data,basicV *gradV,basicM *hessianM,
			    basicV *hessianE)
  {
    // Logic just to avoid warnings from INTEL compiler: added by prudenci on 2009/Sep/06
    const uqGslVectorClass* aux1 = paramDirection;
    if (aux1) {};
    aux1 = gradV;
    aux1 = hessianE;
    uqGslMatrixClass* aux2 = hessianM;
    if (aux2) {};
    const void* aux3 = Data;
    if (aux3) {};

    // Actual code
    static int first_entry = 1;
    double *uqParams = NULL;
    int num_params = paramValue.sizeGlobal();

    if(first_entry)
      {
	if(num_params < 1)
	  QUESO_fatal("Invalid number of parameters");

	if(_QUESO_Basic->m_user_likelihood_func  == NULL )
	  QUESO_fatal("Invalid user-supplied likelihood function");

	uqParams = (double *)calloc(num_params,sizeof(double));
	if(uqParams == NULL)
	  QUESO_fatal("Unable to allocate emmory for uqParams");
      }
	
    for(int i=0;i<num_params;i++)
      {
	uqParams[i] = paramValue[i];
	grvy_printf(GRVY_DEBUG,"%s: sending param to likelihood\n",log_prefix,paramValue[i]);
      }

    grvy_timer_begin("Likelihood Routine");

    double lhood_return = _QUESO_Basic->m_user_likelihood_func(uqParams);

    grvy_timer_end("Likelihood Routine");
    return(lhood_return);

    //      return( _QUESO_Basic->m_user_likelihood_func(uqParams) );
      
  }

}   //  QUESO_Basic_API namespace
