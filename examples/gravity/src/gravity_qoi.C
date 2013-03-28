/*------------------------------------------------------------------------
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
 *------------------------------------------------------------------------
 *
 * $Id$
 */
 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the user defined qoi routine.
 *-----------------------------------------------------------------*/

#include <gravity_qoi.h>
#include <cmath>

//------------------------------------------------------
/// The actual (user-defined) qoi routine
//------------------------------------------------------
void
qoiRoutine(
  const uqGslVectorClass&                    paramValues,
  const uqGslVectorClass*                    paramDirection,
  const void*                                functionDataPtr,
        uqGslVectorClass&                    qoiValues,
        uqDistArrayClass<uqGslVectorClass*>* gradVectors,
        uqDistArrayClass<uqGslMatrixClass*>* hessianMatrices,
        uqDistArrayClass<uqGslVectorClass*>* hessianEffects)
{
  const uqBaseEnvironmentClass& env = paramValues.env();

  if (paramDirection && 
      gradVectors    &&
      hessianEffects &&
      hessianMatrices) {
    // Logic just to avoid warnings from INTEL compiler
  }

  // The user, at the application level, should have set
  // the vector 'paramValues' to have size 1 and
  // the vector 'qoiValues' to have size 1.
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 1,
                      env.fullRank(),
                      "qoiRoutine()",
                      "paramValues vector does not have size 1");

  UQ_FATAL_TEST_MACRO(qoiValues.sizeGlobal() != 1,
                      env.fullRank(),
                      "qoiRoutine()",
                      "qoiValues vector does not have size 1");
  
  // Compute qoi(s)
  double g = paramValues[0]; // Sample of the RV 'gravity acceleration'
  double distanceTraveled = 0.;
  if (env.subRank() == 0) {
    double velocity = ((qoiRoutine_DataClass *) functionDataPtr)->m_initialVelocity;
    double heights  = ((qoiRoutine_DataClass *) functionDataPtr)->m_initialHeight;
    double alpha    = ((qoiRoutine_DataClass *) functionDataPtr)->m_angle;
    
    double aux       = velocity * sin(alpha);
    distanceTraveled = (velocity * cos(alpha) / g) * ( aux + sqrt(pow(aux,2) + 2.*g*heights) );
  }

  qoiValues[0] = distanceTraveled;
    
  return;  
}
