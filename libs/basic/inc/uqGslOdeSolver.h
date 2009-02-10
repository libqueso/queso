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
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_GSL_ODE_SOLVER_H__
#define __UQ_GSL_ODE_SOLVER_H__

#include <uqGslVector.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

class uqGslOdeSolverClass
{
public:
  uqGslOdeSolverClass(const std::string& solverType,
                      unsigned int       sizeOfStateVector);
 ~uqGslOdeSolverClass();

  void solveODE(int (*stateDot   )(double t, const double currentState[], double answer[],             void* infoForComputingStateDot),
                int (*stateDotJac)(double t, const double currentState[], double* dfdy, double dfdt[], void* infoForComputingStateDot),
                const std::vector<double>&      instants,
                const uqGslVectorClass&         state0,
                void*                           infoForComputingStateDot,
                double                          suggestedDeltaT,
                std::vector<uqGslVectorClass*>& states);
protected:
  unsigned int               m_sizeOfStateVector;
  const gsl_odeiv_step_type* m_solver;
  gsl_odeiv_step*            m_s;
  gsl_odeiv_control*         m_c;
  gsl_odeiv_evolve*          m_e;
};

#endif // __UQ_GSL_ODE_SOLVER_H__
