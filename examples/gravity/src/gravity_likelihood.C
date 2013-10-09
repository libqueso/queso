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
 * This file contains the code for the user defined likelihood data 
 * class and the user defined likelihood routine.
 *-----------------------------------------------------------------*/

#include <gravity_likelihood.h>
#include <cmath>
#include <stdio.h>
#include <fstream>

// Construtor
likelihoodRoutine_Data::likelihoodRoutine_Data(const QUESO::BaseEnvironment& env)
  :
  m_heights(0),
  m_times  (0),
  m_stdDevs(0),
  m_env    (&env)
{
  // Data available in /inputData/data02.dat
  double const heights[] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140};
  double const times  [] = {1.41,2.14,2.49,2.87,3.22,3.49,3.81,4.07,4.32,4.47,
      4.75,4.99,5.16,5.26};
  double const stdDevs[] = {0.020,0.120,0.020,0.010,0.030,0.010,0.030,0.030,
      0.030,0.050,0.010,0.040,0.010,0.09};
  
  std::size_t const n = sizeof(heights)/sizeof(*heights); 
  m_heights.assign(heights, heights + n);
  m_times.assign  (times,   times   + n);
  m_stdDevs.assign(stdDevs, stdDevs + n);
}

// Destructor
likelihoodRoutine_Data::~likelihoodRoutine_Data()
{
}

//------------------------------------------------------
// The user defined likelihood routine
//------------------------------------------------------
double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    
  if (paramDirection && functionDataPtr && gradVector && hessianMatrix && hessianEffect) 
  {
    // Just to eliminate INTEL compiler warnings
  }
  
  env.subComm().Barrier(); 
  
  // The user, at the application level, should have set
  // the vector 'paramValues' to have size 1.
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 1,
                      env.fullRank(),
                      "likelihoodRoutine()",
                      "paramValues vector does not have size 1");
  
  // Compute likelihood 
  double g = paramValues[0];  
  const std::vector<double>& heights=((likelihoodRoutine_Data*) functionDataPtr)->m_heights;
  const std::vector<double>& times  =((likelihoodRoutine_Data*) functionDataPtr)->m_times;
  const std::vector<double>& stdDevs=((likelihoodRoutine_Data*) functionDataPtr)->m_stdDevs;
  
  double misfitValue = 0.;
  for (unsigned int i = 0; i < heights.size(); ++i) {
    double modelTime = sqrt(2.0 * heights[i]/g);
    double ratio = (modelTime - times[i])/stdDevs[i];
    misfitValue += ratio*ratio;
  }  
  return (-0.5*misfitValue);  
}
