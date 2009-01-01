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
#include <uqTgaStorage.h>
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
  //std::cout << "Should not call ode-time routine(), case 2" << std::endl;
  //exit(1);

  const uqTgaWDotInfoStruct& info = *((uqTgaWDotInfoStruct *)voidInfo);
  double A    = info.A;
  double E    = info.E;
  double temp = info.temperatureFunctionObj->value(time);
  double expTerm = exp(-E/(R_CONSTANT*temp));

  if (info.computeGradAlso) {
    f[0] =  -A*w[0]                            *expTerm;
    f[1] = (-A*w[1] -   w[0]                  )*expTerm;
    f[2] = (-A*w[2] + A*w[0]/(R_CONSTANT*temp))*expTerm;
  }
  else {
    f[0] = -A*w[0]*expTerm;
  }

  return GSL_SUCCESS;
}

int uqTgaWDotWrtTempRoutine(double temp, const double w[], double f[], void *voidInfo)
{
  const uqTgaWDotInfoStruct& info = *((uqTgaWDotInfoStruct *)voidInfo);
  double A     = info.A;
  double E     = info.E;
  double deriv = info.temperatureFunctionObj->deriv(0.); // COMPATIBILITY WITH OLD VERSION (valid only for linear temperature function)

  UQ_FATAL_TEST_MACRO(info.computeGradAlso,
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaWDotWrtTempRoutine()",
                      "gradient computation was requested but is not supported in this routine");

  f[0] = -A*w[0]*exp(-E/(R_CONSTANT*temp))/deriv; // COMPATIBILITY WITH OLD VERSION

  return GSL_SUCCESS;
}

template<class P_V, class P_M>
class
uqTgaComputableWClass
{
public:
  uqTgaComputableWClass(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                        const uqBase1D1DFunctionClass&     temperatureFunctionObj);
 ~uqTgaComputableWClass();

        void                    computeUsingTime(const P_V&                        params,
                                                 bool                              computeGradAlso,
                                                 const uqTgaStorageClass<P_V,P_M>* referenceW,
                                                 double*                           weigthedMisfitSum,
                                                 uqTgaStorageClass<P_V,P_M>*       weigthedMisfitData);
        void                    computeUsingTemp(const P_V&                        params,
                                                 double                            maximumTemp, // COMPATIBILITY WITH OLD VERSION
                                                 const uqTgaStorageClass<P_V,P_M>* referenceW,
                                                 double*                           weigthedMisfitSum);

        void                    interpolate(double        time,
                                            unsigned int& startingTimeId,
                                            double*       wValue,
                                            P_V*          wGrad) const;
  const std::vector<double>&    times() const;
  const std::vector<double>&    temps() const;
  const std::vector<double>&    ws   () const;
  const std::vector<P_V*  >&    grads() const;
  const uqBaseEnvironmentClass& env  () const;

  void  printForMatlab(std::ofstream& ofs, const std::string& prefixName) const;

protected:
        void                    resetInternalValues();

  const uqBaseEnvironmentClass&  m_env;
  const uqBase1D1DFunctionClass& m_temperatureFunctionObj;

        std::vector<double>      m_times;
        std::vector<double>      m_temps;
        std::vector<double>      m_ws;
        std::vector<P_V*  >      m_grads;
};

template<class P_V, class P_M>
uqTgaComputableWClass<P_V,P_M>::uqTgaComputableWClass(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const uqBase1D1DFunctionClass&     temperatureFunctionObj)
  :
  m_env                   (paramSpace.env()),
  m_temperatureFunctionObj(temperatureFunctionObj),
  m_times                 (0),
  m_temps                 (0),
  m_ws                    (0),
  m_grads                 (0)
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
  m_times.clear();
  m_temps.clear();
  m_ws.clear();
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
  m_grads.clear();
}

template<class P_V, class P_M>
void
uqTgaComputableWClass<P_V,P_M>::computeUsingTime(
  const P_V&                        params,
  bool                              computeGradAlso,
  const uqTgaStorageClass<P_V,P_M>* referenceW,
  double*                           weigthedMisfitSum,
  uqTgaStorageClass<P_V,P_M>*       weigthedMisfitData)
{
  UQ_FATAL_TEST_MACRO((weigthedMisfitSum != NULL) && (referenceW == NULL),
                      m_env.rank(),
                      "uqTgaComputableWClass<P_V,P_M>::computeUsingTime()",
                      "weigthedMisfitSum is being requested but not referenceW is supplied");
  UQ_FATAL_TEST_MACRO((weigthedMisfitData != NULL) && (referenceW == NULL),
                      m_env.rank(),
                      "uqTgaComputableWClass<P_V,P_M>::computeUsingTime()",
                      "weigthedMisfitData is being requested but not referenceW is supplied");

  this->resetInternalValues();
  m_times.resize(1000,0.  );
  m_temps.resize(1000,0.  );
  m_ws.resize   (1000,0.  );
  m_grads.resize(1000,NULL);
  if (weigthedMisfitData) weigthedMisfitData->resizeData(referenceW->times().size());

  uqTgaWDotInfoStruct wDotWrtTimeInfo = {params[0],params[1],&m_temperatureFunctionObj,computeGradAlso};

  unsigned int numWComponents = 1;
  if (computeGradAlso) numWComponents = 3;

  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,numWComponents);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(numWComponents);
        gsl_odeiv_system     sysTime = {uqTgaWDotWrtTimeRoutine, NULL, numWComponents, (void *)&wDotWrtTimeInfo};

  double currentTime = 0.;
  double deltaTime   = .1;
  double currentTemp = m_temperatureFunctionObj.value(currentTime);

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

  double previousTime = 0.;
  double previousW[1];
  previousW[0]=1.;
  unsigned int misfitId = 0;

  double maximumTime = 1.e+9;
  bool continueOnWhile = true;
  if (referenceW) {
    continueOnWhile = (misfitId < referenceW->values().size());
    unsigned int tmpSize = referenceW->times().size();
    maximumTime = referenceW->times()[tmpSize-1];
  }
  else {
    continueOnWhile = (0. < currentW[0]);
  }
  while (continueOnWhile) {
    int status = 0;
    double nextTime = std::min(maximumTime,currentTime+deltaTime);
    status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, nextTime, &deltaTime, currentW);
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        params.env().rank(),
                        "uqTgaComputableWClass<P_V,P_M>::computeUsingTime()",
                        "gsl_odeiv_evolve_apply() failed");
    deltaTime = std::min(.1,deltaTime);
    currentTemp = m_temperatureFunctionObj.value(currentTime);
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
      while ((misfitId < referenceW->values().size()       ) &&
             (previousTime <= referenceW->times()[misfitId]) &&
             (referenceW->times()[misfitId] <= currentTime )) {
        double tmpValue = (referenceW->times()[misfitId]-previousTime)*(currentW[0]-previousW[0])/(currentTime-previousTime) + previousW[0];
        double diff = tmpValue - referenceW->values()[misfitId];

        if (weigthedMisfitData) weigthedMisfitData->setInstantData(misfitId,
                                                                   referenceW->times()[misfitId],
                                                                   referenceW->temps()[misfitId],
                                                                   diff/referenceW->variances()[misfitId],
                                                                   1.);
        if (weigthedMisfitSum) {
          if (referenceW->dataIsContinuousWithTime()) {
            // Properly scale in order to compute integral correctly
            double tmpTimeInterval = 0.;
            if (misfitId == 0) tmpTimeInterval = referenceW->times()[misfitId];
            else               tmpTimeInterval = referenceW->times()[misfitId]-referenceW->times()[misfitId-1];
            *weigthedMisfitSum += diff*diff*tmpTimeInterval/referenceW->variances()[misfitId];
          }
          else {
            *weigthedMisfitSum += diff*diff/referenceW->variances()[misfitId];
          }
        }

        if ((m_env.verbosity() >= 99) && (m_env.rank() == 0)) {
          std::cout << "In uqTgaComputableWClass<P_V,P_M>::computeUsingTime()"
                    << ", currentTime = "  << currentTime
                    << ", misfitId = "     << misfitId
                    << ", measuredTime = " << referenceW->times()[misfitId]
                    << ", measuredTemp = " << referenceW->temps()[misfitId]
                    << ", measuredW = "    << referenceW->values()[misfitId]
                    << ", variance = "     << referenceW->variances()[misfitId]
                    << ": computedW = "    << tmpValue
                    << ", diffValue = "    << diff
                    << std::endl;
        }

        misfitId++;
      }
    }

    previousTime = currentTime;
    previousW[0] = currentW[0];

    if (referenceW) {
      continueOnWhile = (misfitId < referenceW->values().size());
    }
    else {
      continueOnWhile = (0. < currentW[0]);
    }
  }

  m_times.resize(loopId+1);
  m_temps.resize(loopId+1);
  m_ws.resize   (loopId+1);
  m_grads.resize(loopId+1);

  if ((m_env.verbosity() >= 0) && (m_env.rank() == 0)) {
    char stringA[64];
    char stringE[64];
    sprintf(stringA,"%12.6e",params[0]);
    sprintf(stringE,"%12.6e",params[1]);
    std::cout << "In uqTgaComputableWClass<P_V,P_M>::computeUsingTime()"
              << ", A = "                          << stringA
              << ", E = "                          << stringE
              << ", beta = "                       << m_temperatureFunctionObj.deriv(0.)
              << ": finished ode loop after "      << m_times.size()
              << " iterations, with final time = " << m_times[m_times.size()-1]
              << ", final temp = "                 << m_temps[m_temps.size()-1]
              << ", final w = "                    << m_ws   [m_ws.size()   -1]
              << std::endl;
    if (weigthedMisfitSum) {
      std::cout << " and with referenceW->times()[max] = " << referenceW->times()[referenceW->times().size()-1]
                << ", weigthedMisfitSum = "                << *weigthedMisfitSum
                << std::endl;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

template<class P_V, class P_M>
void
uqTgaComputableWClass<P_V,P_M>::computeUsingTemp(
  const P_V&                        params,
  double                            maximumTemp, // COMPATIBILITY WITH OLD VERSION
  const uqTgaStorageClass<P_V,P_M>* referenceW,
  double*                           weigthedMisfitSum)
{
  this->resetInternalValues();
  m_times.resize(1000,0.);
  m_temps.resize(1000,0.);
  m_ws.resize   (1000,0.);
  m_grads.clear();

  uqTgaWDotInfoStruct wDotInfo = {params[0],params[1],&m_temperatureFunctionObj,false};

  unsigned int numWComponents = 1;

  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,numWComponents);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(numWComponents);
        gsl_odeiv_system     sysTemp = {uqTgaWDotWrtTempRoutine, NULL, numWComponents, (void *)&wDotInfo}; 

  double currentTime = 0.;
  double currentTemp = m_temperatureFunctionObj.value(currentTime);
  double deltaTemp   = 1.e-3; // COMPATIBILITY WITH OLD VERSION

  double currentW[numWComponents];
  currentW[0]=1.;

  unsigned int loopId = 0;
  m_times[loopId] = currentTime;
  m_temps[loopId] = currentTemp;
  m_ws   [loopId] = currentW[0];

  double previousTemp = 0.;
  double previousW[1];
  previousW[0]=1.;
  unsigned int misfitId = 0;

  bool continueOnWhile = true;
  if (referenceW) {
    continueOnWhile = (currentTemp < maximumTemp) && (misfitId < referenceW->values().size());
  }
  else if (maximumTemp > 0) {
    continueOnWhile = (currentTemp < maximumTemp); // COMPATIBILITY WITH OLD VERSION
  }
  else {
    continueOnWhile = (currentW[0] > 0.);
  }
  while (continueOnWhile) {
    int status = 0;
    status = gsl_odeiv_evolve_apply(e, c, s, &sysTemp, &currentTemp, maximumTemp, &deltaTemp, currentW); // COMPATIBILITY WITH OLD VERSION
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        params.env().rank(),
                        "uqTgaComputableWClass<P_V,P_M>::computeUsingTemp()",
                        "gsl_odeiv_evolve_apply() failed");
    currentTime = m_temperatureFunctionObj.inverseValue(currentTemp);
    if (currentW[0] < 1.e-6) currentW[0] = 0.;

    loopId++;
    if (loopId >= m_times.size()) {
      m_times.resize(m_times.size()+1000,0.);
      m_temps.resize(m_temps.size()+1000,0.);
      m_ws.resize   (m_ws.size()   +1000,0.);
      m_grads.clear();
    }

    m_times[loopId] = currentTime;
    m_temps[loopId] = currentTemp;
    m_ws   [loopId] = currentW[0];

    if (referenceW) {
      while ((misfitId < referenceW->values().size()       ) &&
             (previousTemp <= referenceW->temps()[misfitId]) &&
             (referenceW->temps()[misfitId] <= currentTemp )) {
        double tmpValue = (referenceW->temps()[misfitId]-previousTemp)*(currentW[0]-previousW[0])/(currentTemp-previousTemp) + previousW[0];
        double diff = tmpValue - referenceW->values()[misfitId];

        if (weigthedMisfitSum) *weigthedMisfitSum += diff*diff/referenceW->variances()[misfitId]; //COMPATIBILITY WITH OLD VERSION

        if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
          std::cout << "In uqTgaComputableWClass<P_V,P_M>::computeUsingTemp()"
                    << ", misfitId = "     << misfitId
                    << ", measuredTemp = " << referenceW->temps()[misfitId]
                    << ", measuredW = "    << referenceW->values()[misfitId]
                    << ": computedW = "    << tmpValue
                    << ", diffValue = "    << diff
                    << std::endl;
        }

        misfitId++;
      }
    }

    previousTemp = currentTemp;
    previousW[0] = currentW[0];

    if (referenceW) {
      continueOnWhile = (currentTemp < maximumTemp) && (misfitId < referenceW->values().size());
    }
    else if (maximumTemp > 0) {
      continueOnWhile = (currentTemp < maximumTemp); // COMPATIBILITY WITH OLD VERSION
    }
    else {
      continueOnWhile = (currentW[0] > 0.);
    }
  }

  m_times.resize(loopId+1);
  m_temps.resize(loopId+1);
  m_ws.resize   (loopId+1);
  m_grads.clear();

  if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
    char stringA[64];
    char stringE[64];
    sprintf(stringA,"%12.6e",params[0]);
    sprintf(stringE,"%12.6e",params[1]);
    std::cout << "In uqTgaComputableWClass<P_V,P_M>::computeUsingTemp()"
              << ", A = "                          << stringA
              << ", E = "                          << stringE
              << ", beta = "                       << m_temperatureFunctionObj.deriv(0.)
              << ": finished ode loop after "      << m_times.size()
              << " iterations, with final time = " << m_times[m_times.size()-1]
              << ", final temp = "                 << m_temps[m_temps.size()-1]
              << ", final w = "                    << m_ws   [m_ws.size()   -1]
              << std::endl;
    if (weigthedMisfitSum) {
      std::cout << " and with referenceW->times()[max] = " << referenceW->times()[referenceW->times().size()-1]
                << ", weigthedMisfitSum = "                << *weigthedMisfitSum
                << std::endl;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

template<class P_V, class P_M>
void
uqTgaComputableWClass<P_V,P_M>::interpolate(
  double        time,
  unsigned int& startingTimeId, // input and output
  double*       wValue,
  P_V*          wGrad) const
{
  unsigned int tmpSize = m_times.size(); // Yes, 'm_grads'
  //std::cout << "In uqTgaComputableWClass<P_V,P_M>::grad()"
  //          << ": time = "           << time
  //          << ", m_times.size() = " << tmpSize
  //          << ", m_times[0] = "     << m_times[0]
  //          << ", m_times[max] = "   << m_times[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      m_env.rank(),
                      "uqTgaComputableW<P_V,P_M>::grad()",
                      "m_times.size() = 0");

  UQ_FATAL_TEST_MACRO(startingTimeId >= tmpSize,
                      m_env.rank(),
                      "uqTgaComputableW<P_V,P_M>::grad()",
                      "startingTimeId is too big");

  UQ_FATAL_TEST_MACRO(time < m_times[0],
                      m_env.rank(),
                      "uqTgaComputableW<P_V,P_M>::grad()",
                      "time < m_times[0]");

  UQ_FATAL_TEST_MACRO(m_times[tmpSize-1] < time,
                      m_env.rank(),
                      "uqTgaComputableW<P_V,P_M>::grad()",
                      "m_times[max] < time");

  UQ_FATAL_TEST_MACRO(wGrad && (m_grads[0] == NULL),
                      m_env.rank(),
                      "uqTgaComputableW<P_V,P_M>::grad()",
                      "m_grads[0] == NULL");

  unsigned int i = 0;
  for (i = startingTimeId; i < tmpSize; ++i) {
    if (time <= m_times[i]) break;
  }
  startingTimeId = i;

  if (time == m_times[i]) {
    if (wValue) *wValue = m_ws[i];
    if (wGrad)  *wGrad  = *(m_grads[i]);
  }
  else {
    //if ((9130.0 < time) && (time < 9131.0)) {
    //  std::cout << "time = " << time
    //            << "i = " << i
    //            << "time[i-1] = " << m_times[i-1]
    //            << "value[i-1] = " << m_values[i-1]
    //            << "time[i] = " << m_times[i]
    //            << "value[i] = " << m_values[i]
    //            << std::endl;
    //}
    double ratio = (time - m_times[i-1])/(m_times[i]-m_times[i-1]);
    if (wValue) *wValue =      m_ws[i-1]  + ratio * (      m_ws[i]  -      m_ws[i-1]  );
    if (wGrad)  *wGrad  = *(m_grads[i-1]) + ratio * ( *(m_grads[i]) - *(m_grads[i-1]) );
  }

  return;
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
const uqBaseEnvironmentClass&
uqTgaComputableWClass<P_V,P_M>::env() const
{
  return m_env;
}

template<class P_V, class P_M>
void
uqTgaComputableWClass<P_V,P_M>::printForMatlab(
  std::ofstream&     ofs,
  const std::string& prefixName) const
{
  unsigned int tmpSize = m_times.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofs << "\n" << prefixName << "Time = zeros(" << tmpSize << ",1);"
        << "\n" << prefixName << "Temp = zeros(" << tmpSize << ",1);"
        << "\n" << prefixName << "W = zeros("    << tmpSize << ",1);";
  }
  else {
    ofs << "\n" << prefixName << "Time = zeros(" << tmpSize << ",1);"
        << "\n" << prefixName << "Temp = zeros(" << tmpSize << ",1);"
        << "\n" << prefixName << "W = zeros("    << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofs << "\n" << prefixName << "Time(" << i+1 << ",1) = " << m_times[i] << ";"
          << "\n" << prefixName << "Temp(" << i+1 << ",1) = " << m_temps[i] << ";"
          << "\n" << prefixName << "W("    << i+1 << ",1) = " << m_ws   [i] << ";";
    }
  }

  return;
}

#endif // __UQ_TGA_COMPUTABLE_W_H__
