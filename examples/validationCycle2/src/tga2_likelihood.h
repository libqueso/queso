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

#ifndef __TGA2_LIKELIHOOD_H__
#define __TGA2_LIKELIHOOD_H__

#include <queso/GslMatrix.h>

//********************************************************
// The (user defined) data class that carries the data
// needed by the (user defined) likelihood routine
//********************************************************
struct
likelihoodRoutine_Data
{
  likelihoodRoutine_Data(const QUESO::BaseEnvironment& env,
                              const char* inpName1,
                              const char* inpName2,
                              const char* inpName3);
 ~likelihoodRoutine_Data();

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

  const QUESO::BaseEnvironment* m_env;
};

double
likelihoodRoutine(
  const QUESO::GslVector&  paramValues,
  const QUESO::GslVector*  paramDirection,
  const void*              functionDataPtr,
  QUESO::GslVector*        gradVector,
  QUESO::GslMatrix*        hessianMatrix,
  QUESO::GslVector*        hessianEffect);

#endif // __TGA2_LIKELIHOOD_H__
