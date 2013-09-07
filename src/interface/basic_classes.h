//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#ifndef __BASIC_CLASSES_H__
#define __BASIC_CLASSES_H__

using namespace std;

#include <uqStatisticalInverseProblem.h>
#include <uqGslMatrix.h>
#include <string.h>

namespace QUESO_Basic_API {

#define basicV GslVector
#define basicM GslMatrix

  void   QUESO_fatal        (const char *message);
  double Likelihood_Wrapper (const basicV &,const basicV *,const void *,basicV *,basicM *,basicV *);

  // Note: the Basic QUESO API uses GSL for the base vector/matrix class

  class QUESO_Basic_ {
  private:
    short int m_initialized;	                                       // basic API initialized?
    short int m_silent; 		                               // silence error messages?
    
    //------------------------				           
    // General QUESO variables				           
    //------------------------				           

    string *m_inputfile;                                               // QUESO ascii input file
    BaseEnvironment *m_env;                                     // QUESO environment
    int     m_num_params;	                                       // number of defined UQ parameter variables
    basicV *m_queso_var_min;                                           // min UQ paramater values
    basicV *m_queso_var_max;                                           // max UQ paramater values
    basicV *m_queso_var_ini;                                           // initial UQ paramater values
    VectorSpace <basicV,basicM> *m_paramSpace;                  // QUESO parameter space
    
    //--------------------------------------
    // Inverse/calibration related variables
    //--------------------------------------
    
    BoxSubset                 <basicV,basicM> *m_paramDomain;   // QUESO parameter domain
    UniformVectorRV           <basicV,basicM> *m_priorRV;       // QUESO prior vector
    GenericVectorRV           <basicV,basicM> *m_postRV;        // QUESO post vector
    StatisticalInverseProblem <basicV,basicM> *m_ip;            // QUESO inverse problem
    GenericScalarFunction     <basicV,basicM> *m_likelihoodObj; // QUESO likelihood object
    basicM *m_CovMatrix;                                               // QUESO covariance matrix
    double (*m_user_likelihood_func) (double *);	               // User supplied likelihood function 
    
  public:
    QUESO_Basic_                 ();
   ~QUESO_Basic_                 (); // prudenci 2010-05-17
    void Initialize                   (const char *inputfile);
    void VerifyInit                   ();
    void DefineParameterSpace         ();
    
    //--------------------------------------
    // Statistical Inverse Problem Functions
    //--------------------------------------
    
    void Likelihood_Register         (double (*fp)(double *) );
    void SolveInverseProblem         ();
    friend double Likelihood_Wrapper (const basicV &,const basicV *,const void *,basicV *,basicM *,basicV *);
    
  };
  
}

#endif //  __BASIC_CLASSES_H__
