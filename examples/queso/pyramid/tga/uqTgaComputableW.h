/* uq/examples/queso/pyramid/uqTgaComputableW.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_COMPUTABLE_W_H__
#define __UQ_TGA_COMPUTABLE_W_H__

#include <uqVectorSpace.h>
#include <uq1D1DFunction.h>
#include <uqTgaDefines.h>
#include <uqTgaMeasuredW.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

// The "W dot" function
typedef struct
{
  double                         A;
  double                         E;
  const uqBase1D1DFunctionClass* temperatureFunctionObj;
  bool                           computeGradAlso;
} uqTgaWDotInfoStruct;

int uqTgaWDotWrtTimeRoutine(double time, const double w[], double f[], void *voidInfo)
{
  std::cout << "Should not call ode-time routine(), case 2" << std::endl;
  exit(1);

  const uqTgaWDotInfoStruct& info = *((uqTgaWDotInfoStruct *)voidInfo);
  double A    = info.A;
  double E    = info.E;
  double temp = info.temperatureFunctionObj->value(time);

  if (info.computeGradAlso) {
    f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp));
    f[1] = 0.;
    f[2] = 0.;
  }
  else {
    f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp));
  }

  return GSL_SUCCESS;
}

int uqTgaWDotWrtTempRoutine(double temp, const double w[], double f[], void *voidInfo)
{
  const uqTgaWDotInfoStruct& info = *((uqTgaWDotInfoStruct *)voidInfo);
  double A     = info.A;
  double E     = info.E;
  double deriv = info.temperatureFunctionObj->deriv(0.); // COMPATIBILITY WITH OLD VERSION (valid only for linear temperature function)

  if (info.computeGradAlso) {
    f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp))/deriv; // COMPATIBILITY WITH OLD VERSION
    f[1] = 0.;
    f[2] = 0.;
  }
  else {
    f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp))/deriv;
  }

  return GSL_SUCCESS;
}

template<class P_V, class P_M>
class
uqTgaComputableWClass
{
public:
  uqTgaComputableWClass(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                        const uqBase1D1DFunctionClass&     temperatureFunctionObj,
                        bool useOdeWithDerivativeWrtTime);
 ~uqTgaComputableWClass();

        void                 compute(const P_V&                          params,
                                     bool                                computeGradAlso,
                                     const uqTgaMeasuredWClass<P_V,P_M>* referenceW,
                                     double                              maximumTemp); // COMPATIBILITY WITH OLD VERSION
        double               w      (double time) const;
        double               wAtTemp(double temp) const;
  const P_V&                 grad   (double time) const;
  const P_V&                 params () const;
  const std::vector<double>& times  () const;
  const std::vector<double>& temps  () const;
  const std::vector<double>& ws     () const;
  const std::vector<P_V*  >& grads  () const;
  const std::vector<double>& diffsForMisfit() const;

protected:
        void                 resetInternalValues();

  const uqBaseEnvironmentClass&  m_env;
  const uqBase1D1DFunctionClass& m_temperatureFunctionObj;
        bool                     m_useOdeWithDerivativeWrtTime;

        P_V                      m_params;
        std::vector<double>      m_times;
        std::vector<double>      m_temps;
        std::vector<double>      m_ws;
        std::vector<P_V*  >      m_grads;
	std::vector<double>      m_diffsForMisfit;
};

template<class P_V, class P_M>
uqTgaComputableWClass<P_V,P_M>::uqTgaComputableWClass(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const uqBase1D1DFunctionClass&     temperatureFunctionObj,
  bool useOdeWithDerivativeWrtTime)
  :
  m_env                        (paramSpace.env()),
  m_temperatureFunctionObj     (temperatureFunctionObj),
  m_useOdeWithDerivativeWrtTime(useOdeWithDerivativeWrtTime),
  m_params(paramSpace.zeroVector()),
  m_times(0),
  m_temps(0),
  m_ws   (0),
  m_grads(0),
  m_diffsForMisfit(0)
{
  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "Entering uqTgaComputableWClass::constructor()"
              << std::endl;
  }

  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqTgaComputableWClass::constructor()"
              << std::endl;
  }
}

template<class P_V, class P_M>
uqTgaComputableWClass<P_V,P_M>::~uqTgaComputableWClass()
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
}

template<class P_V, class P_M>
void
uqTgaComputableWClass<P_V,P_M>::resetInternalValues()
{
  m_params *= 0.;
  m_times.clear();
  m_temps.clear();
  m_ws.clear();
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
  m_grads.clear();
  m_diffsForMisfit.clear();
}

template<class P_V, class P_M>
void
uqTgaComputableWClass<P_V,P_M>::compute(
  const P_V&                          params,
  bool                                computeGradAlso,
  const uqTgaMeasuredWClass<P_V,P_M>* referenceW,
  double                              maximumTemp) // COMPATIBILITY WITH OLD VERSION
{
  this->resetInternalValues();
  m_params = params;
  m_times.resize(1000,0.  );
  m_temps.resize(1000,0.  );
  m_ws.resize   (1000,0.  );
  m_grads.resize(1000,NULL);
  if (referenceW) m_diffsForMisfit.resize(referenceW->times().size(),0.);

  uqTgaWDotInfoStruct wDotWrtTimeInfo = {params[0],params[1],&m_temperatureFunctionObj,computeGradAlso};

  unsigned int numWComponents = 1;
  if (computeGradAlso) numWComponents = 3;

  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,numWComponents);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(numWComponents);
        gsl_odeiv_system     sysTime = {uqTgaWDotWrtTimeRoutine, NULL, numWComponents, (void *)&wDotWrtTimeInfo};
        gsl_odeiv_system     sysTemp = {uqTgaWDotWrtTempRoutine, NULL, numWComponents, (void *)&wDotWrtTimeInfo}; 

  double currentTime = 0.;
  double deltaTime   = 1e-3; // COMPATIBILITY WITH OLD VERSION

  double currentTemp = m_temperatureFunctionObj.value(currentTime);
  double deltaTemp   = 1e-3; // COMPATIBILITY WITH OLD VERSION

  double currentW[numWComponents];
  currentW[0]=1.;
  if (computeGradAlso) {
    currentW[1]=0.;
    currentW[2]=0.;
  }

  unsigned int loopId = 0;
  m_times[loopId] = currentTime;
  m_temps[loopId] = currentTemp;
  m_ws   [loopId] = currentW[0];
  if (computeGradAlso) {
    m_grads[loopId] = new P_V(params);
    (*(m_grads[loopId]))[0] = currentW[1];
    (*(m_grads[loopId]))[1] = currentW[2];
  }

  double previousTemp = 0.;
  double previousW[1];
  previousW[0]=1.;
  unsigned int misfitId = 0;

  bool continueOnWhile = true;
  if (referenceW) {
    continueOnWhile = (currentTemp < maximumTemp) && (misfitId < referenceW->ws().size());
  }
  else if (maximumTemp > 0) {
    continueOnWhile = (currentTemp < maximumTemp); // COMPATIBILITY WITH OLD VERSION
  }
  else {
    continueOnWhile = (0. < currentW[0]);
  }
  while (continueOnWhile) {
    int status = 0;
    if (m_useOdeWithDerivativeWrtTime) {
      //deltaTime = .1;
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, maximumTemp/*currentTime+deltaTime*/, &deltaTime, currentW);
      currentTemp = m_temperatureFunctionObj.value(currentTime);
    }
    else {
      //deltaTemp = .1;
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, maximumTemp/*currentTemp+deltaTemp*/, &deltaTemp, currentW); // COMPATIBILITY WITH OLD VERSION
      currentTime = m_temperatureFunctionObj.inverseValue(currentTemp);
    }
    if (currentW[0] < 1.e-6) currentW[0] = 0.;

    loopId++;
    if (loopId >= m_times.size()) {
      m_times.resize(m_times.size()+1000,0.);
      m_temps.resize(m_temps.size()+1000,0.);
      m_ws.resize   (m_ws.size()   +1000,0.);
      m_grads.resize(m_grads.size()+1000,NULL);
    }

    m_times[loopId] = currentTime;
    m_temps[loopId] = currentTemp;
    m_ws   [loopId] = currentW[0];
    if (computeGradAlso) {
      m_grads[loopId] = new P_V(params);
      (*(m_grads[loopId]))[0] = currentW[1];
      (*(m_grads[loopId]))[1] = currentW[2];
    }

    if (referenceW) {
      while ( (misfitId < referenceW->ws().size()) && (previousTemp <= referenceW->temps()[misfitId]) && (referenceW->temps()[misfitId] <= currentTemp) ) {
        double tmpValue = (referenceW->temps()[misfitId]-previousTemp)*(currentW[0]-previousW[0])/(currentTemp-previousTemp) + previousW[0];
        m_diffsForMisfit[misfitId] = tmpValue - referenceW->ws()[misfitId];
        misfitId++;
      }
    }

    previousTemp = currentTemp;
    previousW[0] = currentW[0];

    if (referenceW) {
      continueOnWhile = (currentTemp < maximumTemp) && (misfitId < referenceW->ws().size());
    }
    else if (maximumTemp > 0) {
      continueOnWhile = (currentTemp < maximumTemp); // COMPATIBILITY WITH OLD VERSION
    }
    else {
      continueOnWhile = (0. < currentW[0]);
    }
  }

  m_times.resize(loopId+1);
  m_temps.resize(loopId+1);
  m_ws.resize   (loopId+1);
  m_grads.resize(loopId+1);

  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "uqTgaComputableWClass<P_V,P_M>::compute()"
              << ", with A = "       << params[0]
              << ", E = "            << params[1]
              << ": computed w has " << m_times.size()
              << " samples"
              << std::endl;
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

template<class P_V, class P_M>
double
uqTgaComputableWClass<P_V,P_M>::w(double time) const
{
  double value = 0.;

  return value;
}

template<class P_V, class P_M>
double
uqTgaComputableWClass<P_V,P_M>::wAtTemp(double temp) const
{
  double value = 0.;

  return value;
}

template<class P_V, class P_M>
const P_V&
uqTgaComputableWClass<P_V,P_M>::grad(double time) const
{
  return;
}

template<class P_V, class P_M>
const P_V&
uqTgaComputableWClass<P_V,P_M>::params() const
{
  return m_params;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaComputableWClass<P_V,P_M>::times() const
{
  return m_times;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaComputableWClass<P_V,P_M>::temps() const
{
  return m_temps;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaComputableWClass<P_V,P_M>::ws() const
{
  return m_ws;
}

template<class P_V, class P_M>
const std::vector<P_V*>&
uqTgaComputableWClass<P_V,P_M>::grads() const
{
  return m_grads;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaComputableWClass<P_V,P_M>::diffsForMisfit() const
{
  return m_diffsForMisfit;
}
#endif // __UQ_TGA_COMPUTABLE_W_H__
