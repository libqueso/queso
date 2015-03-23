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

#include <tga2_func.h>
#include <tga2_qoi.h>
#include <gsl/gsl_odeiv.h>

// The actual (user defined) qoi routine
void qoiRoutine(const QUESO::GslVector&                    paramValues,
                const QUESO::GslVector*                    paramDirection,
                const void*                                functionDataPtr,
                      QUESO::GslVector&                    qoiValues,
                      QUESO::DistArray<QUESO::GslVector*>* gradVectors,
                      QUESO::DistArray<QUESO::GslMatrix*>* hessianMatrices,
                      QUESO::DistArray<QUESO::GslVector*>* hessianEffects)
{
  if (paramDirection  &&
      functionDataPtr &&
      gradVectors     &&
      hessianMatrices &&
      hessianEffects) {
    // Just to eliminate INTEL compiler warnings
  }

  double A             = paramValues[0];
  double E             = paramValues[1];
  double beta          = ((qoiRoutine_Data *) functionDataPtr)->m_beta;
  double criticalMass  = ((qoiRoutine_Data *) functionDataPtr)->m_criticalMass;
  double criticalTime  = ((qoiRoutine_Data *) functionDataPtr)->m_criticalTime;

  double params[]={A,E,beta};

  // integration
  const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s   = gsl_odeiv_step_alloc(T,1);
        gsl_odeiv_control   *c   = gsl_odeiv_control_y_new(1e-6,0.0);
        gsl_odeiv_evolve    *e   = gsl_odeiv_evolve_alloc(1);
        gsl_odeiv_system     sys = {func, NULL, 1, (void *)params};

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
                        paramValues.env().fullRank(),
                        "qoiRoutine()",
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
  if ((paramValues.env().displayVerbosity() >= 3) && (paramValues.env().fullRank() == 0)) {
    printf("In qoiRoutine(), A = %g, E = %g, beta = %.3lf, criticalTime = %.3lf, criticalMass = %.3lf: qoi = %lf.\n",A,E,beta,criticalTime,criticalMass,qoiValues[0]);
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}
