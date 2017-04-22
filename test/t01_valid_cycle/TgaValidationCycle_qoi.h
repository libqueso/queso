//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef EX_TGA_VALIDATION_CYCLE_QOI_H
#define EX_TGA_VALIDATION_CYCLE_QOI_H

#include <queso/DistArray.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

//********************************************************
// The (user defined) data class that carries the data
// needed by the (user defined) qoi routine
//********************************************************
template<class P_V,class P_M,class Q_V, class Q_M>
struct
qoiRoutine_Data
{
  double m_beta;
  double m_criticalMass;
  double m_criticalTime;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void qoiRoutine(const P_V&                    paramValues,
                const P_V*                    paramDirection,
                const void*                   functionDataPtr,
                      Q_V&                    qoiValues,
                      QUESO::DistArray<P_V*>* gradVectors,
                      QUESO::DistArray<P_M*>* hessianMatrices,
                      QUESO::DistArray<P_V*>* hessianEffects)
{
  double A             = paramValues[0];
  double E             = paramValues[1];
  double beta          = ((qoiRoutine_Data<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_beta;
  double criticalMass  = ((qoiRoutine_Data<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_criticalMass;
  double criticalTime  = ((qoiRoutine_Data<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_criticalTime;

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
                        paramValues.env().worldRank(),
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

#endif // EX_TGA_VALIDATION_CYCLE_QOI_H
