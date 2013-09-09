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
 * This is the header file for gravity_likelihood.C. 
 *-----------------------------------------------------------------*/

#ifndef __GRAVITY_LIKELIHOOD_H__
#define __GRAVITY_LIKELIHOOD_H__

#include <uqGslMatrix.h>

struct likelihoodRoutine_DataClass // user defined class
{
  likelihoodRoutine_DataClass(const QUESO::BaseEnvironmentClass& env);
 ~likelihoodRoutine_DataClass();

  std::vector<double> m_heights; // heights
  std::vector<double> m_times;   // times
  std::vector<double> m_stdDevs; // account for uncertainties in 
				 // time measurement: sigmas

  const QUESO::BaseEnvironmentClass* m_env;
};

double likelihoodRoutine( // user defined routine
  const QUESO::GslVectorClass& paramValues,
  const QUESO::GslVectorClass* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVectorClass*       gradVector,
  QUESO::GslMatrixClass*       hessianMatrix,
  QUESO::GslVectorClass*       hessianEffect);

#endif
