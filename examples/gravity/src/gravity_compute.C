/*-------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
 *-------------------------------------------------------------------
 *
 * $Id$
 */
 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file is divided in two parts:
 * - the first one handles the statistical inverse problem (SIP) 
 * for estimating the magnitude 'g' of gravity acceleration; and
 * - the second part handles the statistical forward problem (SFP) 
 * for predicting the maximum distance traveled by a projectile.
 *
 * The SIP definition requires a user defined likelihood function. 
 * See files 'gravity_likelihood.h' and 'gravity_likelihood.C'.
 *
 * The SFP definition requires a user defined qoi function. 
 * See files 'gravity_qoi.h' and 'gravity_qoi.C'.
 *-----------------------------------------------------------------*/

#include <gravity_compute.h>
#include <gravity_likelihood.h>
#include <gravity_qoi.h>
#include <uqGslMatrix.h>
#include <uqStatisticalInverseProblem.h>
#include <uqStatisticalForwardProblem.h>
#include <sys/time.h>
#include <cmath>

//================================================================
// If PRIOR_IS_GAUSSIAN is defined, then:
// --> Gaussian prior for gravity.
// Otherwise:
// --> Uniform prior for gravity.
//================================================================

//#define PRIOR_IS_GAUSSIAN

void computeGravityAndTraveledDistance(const QUESO::FullEnvironment& env) {
  struct timeval timevalNow;
  
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\nBeginning run of 'Gravity + Projectile motion' example at "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  // Just examples of possible calls
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Beginning run of 'Gravity + Projectile motion' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  env.fullComm().Barrier();
  env.subComm().Barrier();  // Just an example of a possible call
  
  //================================================================
  // Statistical inverse problem (SIP): find posterior PDF for 'g'
  //================================================================
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> Gravity estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", 1, NULL);

  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());
  
  paramMinValues[0] = 8.;
  paramMaxValues[0] = 11.;
  
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",
      paramSpace,
      paramMinValues,
      paramMaxValues);
  
  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function 
  // object to be used by QUESO.
  //------------------------------------------------------
  likelihoodRoutine_Data likelihoodRoutine_Data(env);
  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("like_",
			  paramDomain,
			  likelihoodRoutine,
			  (void *) &likelihoodRoutine_Data,
			  true); // the routine computes [ln(function)]
    
  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  
#ifdef PRIOR_IS_GAUSSIAN
  QUESO::GslVector meanVector( paramSpace.zeroVector() );
  meanVector[0] = 9;
 
  QUESO::GslMatrix covMatrix = QUESO::GslMatrix(paramSpace.zeroVector());
  covMatrix(0,0) = 1.; 
  
  // Create a Gaussian prior RV
  QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix> priorRv("prior_",paramDomain,meanVector,covMatrix);
  
#else
  // Create an uniform prior RV
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix> priorRv("prior_",paramDomain);
  
#endif
  
  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_",  // Extra prefix before the default "rv_" prefix
           paramSpace);
        
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("",          // No extra prefix before the default "ip_" prefix
       NULL, 
       priorRv, 
       likelihoodFunctionObj, 
       postRv); 

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV
  //------------------------------------------------------
  std::cout << "Solving the SIP with Metropolis Hastings" 
	    << std::endl << std::endl;  

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  priorRv.realizer().realization(paramInitials);

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = std::pow( fabs(paramInitials[0])/20. , 2. );

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  //================================================================
  // Statistical forward problem (SFP): find the max distance 
  // traveled by an object in projectile motion; input pdf for 'g' 
  // is the solution of the SIP above.
  //================================================================
  gettimeofday(&timevalNow, NULL);
  std::cout << "Beginning 'SFP -> Projectile motion' at " 
            << ctime(&timevalNow.tv_sec)
            << std::endl;
	    
  //------------------------------------------------------
  // SFP Step 1 of 6: Instantiate the parameter *and* qoi spaces. 
  // SFP input RV = FIP posterior RV, so SFP parameter space
  // has been already defined.
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> qoiSpace(env, "qoi_", 1, NULL);
  
  //------------------------------------------------------
  // SFP Step 2 of 6: Instantiate the parameter domain 
  //------------------------------------------------------
  
  // Not necessary because input RV of the SFP = output RV of SIP. 
  // Thus, the parameter domain has been already defined.
  
  //------------------------------------------------------ 
  // SFP Step 3 of 6: Instantiate the qoi function object 
  // to be used by QUESO.
  //------------------------------------------------------
  qoiRoutine_Data qoiRoutine_Data;
  qoiRoutine_Data.m_angle          = M_PI/4.0; //45 degrees (radians)
  qoiRoutine_Data.m_initialVelocity= 5.;      //initial speed (m/s) 
  qoiRoutine_Data.m_initialHeight  = 0.;       //initial height (m)
  
  QUESO::GenericVectorFunction<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    qoiFunctionObj("qoi_",
                   paramDomain,
                   qoiSpace,
                   qoiRoutine,
                   (void *) &qoiRoutine_Data);
      
  //------------------------------------------------------
  // SFP Step 4 of 6: Define the input RV
  //------------------------------------------------------
  
  // Not necessary because input RV of SFP = output RV of SIP 
  // (postRv).
      
  //------------------------------------------------------
  // SFP Step 5 of 6: Instantiate the forward problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> qoiRv("qoi_", qoiSpace);
  
  QUESO::StatisticalForwardProblem<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    fp("",
       NULL,
       postRv,
       qoiFunctionObj,
       qoiRv);

  //------------------------------------------------------
  // SFP Step 6 of 6: Solve the forward problem
  //------------------------------------------------------
  std::cout << "Solving the SFP with Monte Carlo" 
            << std::endl << std::endl;  
  fp.solveWithMonteCarlo(NULL);

  //------------------------------------------------------
  gettimeofday(&timevalNow, NULL);
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Ending run of 'Gravity + Projectile motion' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  if (env.fullRank() == 0) {
    std::cout << "Ending run of 'Gravity + Projectile motion' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  return;
}
