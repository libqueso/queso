/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008,2009 The PECOS Development Team
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
 *
 * $Id$
 *
 * Basic API: Class definitions for basic API (uses GSL)
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __BASIC_CLASSES_H__
#define __BASIC_CLASSES_H__

using namespace std;

#include <uqStatisticalInverseProblem.h>
#include <uqGslMatrix.h>
#include <string.h>

namespace QUESO_Basic_API {

#define basicV uqGslVectorClass
#define basicM uqGslMatrixClass

  void   QUESO_fatal        (const char *message);
  double Likelihood_Wrapper (const basicV &,const basicV *,const void *,basicV *,basicM *,basicV *);

  // Note: the Basic QUESO API uses GSL for the base vector/matrix class

  class QUESO_Basic_Class {
  private:
    short int m_initialized;	                                       // basic API initialized?
    short int m_silent; 		                               // silence error messages?
    
    //------------------------				           
    // General QUESO variables				           
    //------------------------				           

    string *m_inputfile;                                               // QUESO ascii input file
    uqBaseEnvironmentClass *m_env;                                     // QUESO environment
    int     m_num_params;	                                       // number of defined UQ parameter variables
    basicV *m_queso_var_min;                                           // min UQ paramater values
    basicV *m_queso_var_max;                                           // max UQ paramater values
    basicV *m_queso_var_ini;                                           // initial UQ paramater values
    uqVectorSpaceClass <basicV,basicM> *m_paramSpace;                  // QUESO parameter space
    
    //--------------------------------------
    // Inverse/calibration related variables
    //--------------------------------------
    
    uqBoxSubsetClass                 <basicV,basicM> *m_paramDomain;   // QUESO parameter domain
    uqUniformVectorRVClass           <basicV,basicM> *m_priorRV;       // QUESO prior vector
    uqGenericVectorRVClass           <basicV,basicM> *m_postRV;        // QUESO post vector
    uqStatisticalInverseProblemClass <basicV,basicM> *m_ip;            // QUESO inverse problem
    uqGenericScalarFunctionClass     <basicV,basicM> *m_likelihoodObj; // QUESO likelihood object
    basicM *m_CovMatrix;                                               // QUESO covariance matrix
    double (*m_user_likelihood_func) (double *);	               // User supplied likelihood function 
    
  public:
    QUESO_Basic_Class                 ();
   ~QUESO_Basic_Class                 (); // prudenci 2010-05-17
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
