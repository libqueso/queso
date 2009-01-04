/* uq/examples/queso/pyramid/uqTgaW.h
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
struct
uqTgaWDotInfoStruct
{
  double                         A;
  double                         E;
  const uqBase1D1DFunctionClass* temperatureFunctionObj;
  bool                           computeGradAlso;
};

int uqTgaWDotWrtTimeRoutine(double time, const double w[], double f[], void *voidInfo)
{
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

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
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
#endif

template<class P_V, class P_M>
class
uqTgaWClass
{
public:
  uqTgaWClass(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
              const uqBase1D1DFunctionClass&     temperatureFunctionObj);
 ~uqTgaWClass();

        void                    compute(const P_V&                     params,
                                        double                         reqMaxTime,
                                        double                         maxTimeStep,
                                        bool                           computeGradAlso,
                                        const uqBase1D1DFunctionClass* referenceW,
                                        const uqBase1D1DFunctionClass* weightFunction,
                                        double*                        misfitValue,
                                        uqSampled1D1DFunctionClass*    diffFunction);
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
        void                    computeUsingTemp(const P_V&                        params,
                                                 double                            maximumTemp, // COMPATIBILITY WITH OLD VERSION
                                                 const uqTgaStorageClass<P_V,P_M>* referenceW,
                                                 double*                           misfitValue);
#endif

        void                    interpolate(double        time,
                                            unsigned int& startingTimeId,
                                            double*       wValue,
                                            P_V*          wGrad,
                                            bool*         timeWasMatchedExactly) const;
  const std::vector<double>&    times() const;
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
	std::vector<double>      m_diffs;
};

template<class P_V, class P_M>
uqTgaWClass<P_V,P_M>::uqTgaWClass(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const uqBase1D1DFunctionClass&     temperatureFunctionObj)
  :
  m_env                   (paramSpace.env()),
  m_temperatureFunctionObj(temperatureFunctionObj),
  m_times                 (0),
  m_temps                 (0),
  m_ws                    (0),
  m_grads                 (0),
  m_diffs                 (0)
{
  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "Entering uqTgaWClass::constructor()"
              << std::endl;
  }

  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqTgaWClass::constructor()"
              << std::endl;
  }
}

template<class P_V, class P_M>
uqTgaWClass<P_V,P_M>::~uqTgaWClass()
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
}

template<class P_V, class P_M>
void
uqTgaWClass<P_V,P_M>::resetInternalValues()
{
  m_times.clear();
  m_temps.clear();
  m_ws.clear();
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
  m_grads.clear();
  m_diffs.clear();
}

template<class P_V, class P_M>
void
uqTgaWClass<P_V,P_M>::compute(
  const P_V&                     params,          // input
  double                         reqMaxTime,      // input
  double                         maxTimeStep,     // input
  bool                           computeGradAlso, // input
  const uqBase1D1DFunctionClass* referenceW,      // input (possible)
  const uqBase1D1DFunctionClass* weightFunction,  // input (possible)
  double*                        misfitValue,     // output 
  uqSampled1D1DFunctionClass*    diffFunction)    // output
{
  UQ_FATAL_TEST_MACRO((misfitValue != NULL) && (referenceW == NULL),
                      m_env.rank(),
                      "uqTgaWClass<P_V,P_M>::compute()",
                      "misfitValue is being requested but not referenceW is supplied");
  UQ_FATAL_TEST_MACRO((diffFunction != NULL) && (referenceW == NULL),
                      m_env.rank(),
                      "uqTgaWClass<P_V,P_M>::compute()",
                      "diffFunction is being requested but not referenceW is supplied");
  UQ_FATAL_TEST_MACRO((referenceW != NULL) && (misfitValue == NULL) && (diffFunction == NULL),
                      m_env.rank(),
                      "uqTgaWClass<P_V,P_M>::compute()",
                      "refenceW is being supplied for no purpose");

  if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
    std::cout << "Entering uqTgaWClass<P_V,P_M>::compute()"
              << ": params = "          << params
              << ", reqMaxTime = "      << reqMaxTime
              << ", maxTimeStep = "     << maxTimeStep
              << ", computeGradAlso = " << computeGradAlso
              << ", referenceW = "      << referenceW
              << ", weightFunction = "  << weightFunction
              << ", misfitValue = "     << misfitValue
              << ", diffFunction = "    << diffFunction
              << std::endl;
  }

  const std::vector<double> dummyVec(0);

  // Initialize variables related to the delta sequence function
  const uqDeltaSeq1D1DFunctionClass* deltaSeqFunction = NULL;
  unsigned int deltaSeqSize = 0;
  const std::vector<double>* deltaSeqTimesPtr(&dummyVec);
  const std::vector<double>* deltaSeqIntegratedValuesPtr(&dummyVec);
  if (weightFunction) {
    deltaSeqFunction = dynamic_cast< const uqDeltaSeq1D1DFunctionClass* >(weightFunction);
    if (deltaSeqFunction != NULL) {
      deltaSeqSize = deltaSeqFunction->domainValues().size();
      deltaSeqTimesPtr = &(deltaSeqFunction->domainValues());
      deltaSeqIntegratedValuesPtr = &(deltaSeqFunction->integratedValues());
    }
  }
  const std::vector<double>& deltaSeqTimes           (*deltaSeqTimesPtr);
  const std::vector<double>& deltaSeqIntegratedValues(*deltaSeqIntegratedValuesPtr);
  unsigned int deltaSeqId = 0;

  // Initialize other variables
  this->resetInternalValues();
  m_times.resize(1000,0.  );
  m_ws.resize   (1000,0.  );
  if (computeGradAlso) m_grads.resize(1000,NULL);
  if (diffFunction)    m_diffs.resize(1000,0.  );

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
  double timeStep   = .1;
  if (maxTimeStep > 0) timeStep = maxTimeStep;
  double currentTemp = m_temperatureFunctionObj.value(currentTime);

  double currentW[numWComponents];
  currentW[0]=1.;
  if (computeGradAlso) {
    currentW[1]=0.;
    currentW[2]=0.;
  }

  unsigned int loopId = 0;
  m_times[loopId] = currentTime;
  m_ws   [loopId] = currentW[0];
  if (computeGradAlso) {
    m_grads[loopId] = new P_V(params);
    (*(m_grads[loopId]))[0] = currentW[1];
    (*(m_grads[loopId]))[1] = currentW[2];
  }
  if (diffFunction) m_diffs[loopId] = currentW[0] - referenceW->value(currentTime);

  double previousTime = 0.;
  double previousW[1];
  previousW[0]=1.;

  double maximumTime = 1.e+9;
  if (reqMaxTime > 0.) maximumTime = reqMaxTime;
  if (referenceW     ) maximumTime = referenceW->maxDomainValue();     // FIX ME
  if (weightFunction ) maximumTime = weightFunction->maxDomainValue(); // FIX ME

  bool continueOnWhile = true;
#if 1
  if ((reqMaxTime > 0.) || referenceW || weightFunction) {
    continueOnWhile = (currentTime < maximumTime);
  }
#else
  if (weightFunctionIsDeltaSeq) {
    continueOnWhile = (deltaSeqId < deltaSeqSize);
  }
  else if (referenceW) {
    continueOnWhile = (currentTime < maximumTime);
  }
  else if (reqMaxTime > 0) {
    continueOnWhile = (currentTime < maximumTime);
  }
#endif
  else {
    continueOnWhile = (0. < currentW[0]);
  }
  while (continueOnWhile) {
    int status = 0;
    double nextTime = maximumTime;
    if (maxTimeStep > 0) nextTime = std::min(nextTime,currentTime+maxTimeStep);
    status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTime, nextTime, &timeStep, currentW);
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        params.env().rank(),
                        "uqTgaWClass<P_V,P_M>::compute()",
                        "gsl_odeiv_evolve_apply() failed");
    if (currentW[0] < 1.e-6) currentW[0] = 0.;
    if (maxTimeStep > 0) timeStep = std::min(maxTimeStep,timeStep);
    currentTemp = m_temperatureFunctionObj.value(currentTime);

    double diff = 0.;
    if (referenceW) diff = currentW[0] - referenceW->value(currentTime);

    if ((m_env.verbosity() >= 99) && (m_env.rank() == 0)) {
      std::cout << "In uqTgaWClass<P_V,P_M>::compute()"
                << ": loopId = "      << loopId
                << ", currentTime = " << currentTime
                << ", currentW[0] = " << currentW[0]
                << ", maximumTime-currentTime = " << maximumTime-currentTime
                << ", deltaSeqFunction = "        << deltaSeqFunction
                << std::endl;
    }

    loopId++;
    if (loopId >= m_times.size()) {
      m_times.resize(m_times.size()+1000,0.);
      m_ws.resize   (m_ws.size()   +1000,0.);
      if (computeGradAlso) m_grads.resize(m_grads.size()+1000,NULL);
      if (diffFunction)    m_diffs.resize(m_diffs.size()+1000,0.  );
    }

    m_times[loopId] = currentTime;
    m_ws   [loopId] = currentW[0];
    if (computeGradAlso) {
      m_grads[loopId] = new P_V(params);
      (*(m_grads[loopId]))[0] = currentW[1];
      (*(m_grads[loopId]))[1] = currentW[2];
    }
    if (diffFunction) m_diffs[loopId] = diff;

    if (misfitValue) {
      if (deltaSeqFunction != NULL) {
        while ((deltaSeqId                < deltaSeqSize              ) &&
               (previousTime              <= deltaSeqTimes[deltaSeqId]) &&
               (deltaSeqTimes[deltaSeqId] <= currentTime              )) {
          double tmpValue = (deltaSeqTimes[deltaSeqId]-previousTime)*(currentW[0]-previousW[0])/(currentTime-previousTime) + previousW[0];
          double diffForDelta = tmpValue - referenceW->value(deltaSeqTimes[deltaSeqId]);
          *misfitValue += diffForDelta*diffForDelta*deltaSeqIntegratedValues[deltaSeqId]; // Consider delta as delta, i.e., integral = 1;
#if 0
          if ((m_env.verbosity() >= 99) && (m_env.rank() == 0)) {
            std::cout << "In uqTgaWClass<P_V,P_M>::compute()"
                      << ", currentTime = "   << currentTime
                      << ", deltaSeqId = "    << deltaSeqId
                      << ", measuredTime = "  << deltaSeqTimes[deltaSeqId]
                      << ", measuredW = "     << refValues[deltaSeqId]
                      << ", variance = "      << refVariances[deltaSeqId]
                      << ": computedW = "     << tmpValue
                      << ", diffForDelta = "  << diffForDelta
                      << std::endl;
          }
#endif
          deltaSeqId++;
        }
      }
      else {
        // Multiply by [time interval] in order to compute integral correctly
        double weightScale = 1.;
        if (weightFunction) weightScale *= weightFunction->value(currentTime);
        *misfitValue += diff*diff*weightScale*(currentTime-previousTime);
      }
    }

    previousTime = currentTime;
    previousW[0] = currentW[0];

#if 1
    if ((reqMaxTime > 0.) || referenceW || weightFunction) {
      continueOnWhile = (currentTime < maximumTime);
    }
#else
    if (weightFunctionIsDeltaSeq) {
      continueOnWhile = (deltaSeqId < deltaSeqSize);
    }
    else if (referenceW) {
      continueOnWhile = (currentTime < maximumTime);
    }
    else if (reqMaxTime > 0) {
      continueOnWhile = (currentTime < maximumTime);
    }
#endif
    else {
      continueOnWhile = (0. < currentW[0]);
    }
  }

  m_times.resize(loopId+1);
  m_ws.resize   (loopId+1);
  if (computeGradAlso) m_grads.resize(loopId+1);
  if (diffFunction)    m_diffs.resize(loopId+1);

  if (diffFunction) diffFunction->set(m_times,m_diffs);

  if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
    char stringA[64];
    char stringE[64];
    sprintf(stringA,"%12.6e",params[0]);
    sprintf(stringE,"%12.6e",params[1]);
    std::cout << "In uqTgaWClass<P_V,P_M>::compute()"
              << ", A = "                          << stringA
              << ", E = "                          << stringE
              << ", beta = "                       << m_temperatureFunctionObj.deriv(0.)
              << ": finished ode loop after "      << m_times.size()
              << " iterations, with final time = " << m_times[m_times.size()-1]
              << ", final w = "                    << m_ws   [m_ws.size()   -1]
              << std::endl;
    if (misfitValue) {
      std::cout << ", and with misfitValue = " << *misfitValue
                << std::endl;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
template<class P_V, class P_M>
void
uqTgaWClass<P_V,P_M>::computeUsingTemp(
  const P_V&                        params,
  double                            maximumTemp, // COMPATIBILITY WITH OLD VERSION
  const uqTgaStorageClass<P_V,P_M>* referenceW,
  double*                           misfitValue)
{
  ////UQ_FATAL_TEST_MACRO((misfitValue != NULL) && (referenceW == NULL),
  ////                    m_env.rank(),
  ////                    "uqTgaWClass<P_V,P_M>::computeUsingTemp()",
  ////                    "misfitValue is being requested but not referenceW is supplied");

  const std::vector<double> dummyVec(0);
  const std::vector<double>* refTimesPtr    (&dummyVec);
  const std::vector<double>* refTempsPtr    (&dummyVec);
  const std::vector<double>* refValuesPtr   (&dummyVec);
  const std::vector<double>* refVariancesPtr(&dummyVec);
  if (referenceW) {
    refTimesPtr     = &(referenceW->times    ());
    refTempsPtr     = &(referenceW->temps    ());
    refValuesPtr    = &(referenceW->values   ());
    refVariancesPtr = &(referenceW->variances());
  }
  const std::vector<double>& refTimes    (*refTimesPtr    );
  const std::vector<double>& refTemps    (*refTempsPtr    );
  const std::vector<double>& refValues   (*refValuesPtr   );
  const std::vector<double>& refVariances(*refVariancesPtr);
  unsigned int refSize = refTimes.size();

  unsigned int sizeForResize = 100;
  this->resetInternalValues();
  ////m_times.resize(sizeForResize,0.);
  ////m_temps.resize(sizeForResize,0.);
  m_ws.resize   (sizeForResize,0.);
  ////m_grads.clear();

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
  ////m_times[loopId] = currentTime;
  ////m_temps[loopId] = currentTemp;
  m_ws   [loopId] = currentW[0];

  double previousTemp = 0.;
  double previousW[1];
  previousW[0]=1.;
  unsigned int misfitId = 0;

  double sumValue = 0.;
  bool continueOnWhile = true;
  if (referenceW) {
    continueOnWhile = (currentTemp < maximumTemp) && (misfitId < refSize);
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
                        "uqTgaWClass<P_V,P_M>::computeUsingTemp()",
                        "gsl_odeiv_evolve_apply() failed");
    ////currentTime = m_temperatureFunctionObj.inverseValue(currentTemp); // Comment for performance reasons
    if (currentW[0] < 1.e-6) currentW[0] = 0.;

    loopId++;
    if (loopId >= m_ws.size()) {
      ////m_times.resize(m_times.size()+sizeForResize,0.);
      ////m_temps.resize(m_temps.size()+sizeForResize,0.);
      m_ws.resize   (m_ws.size()   +sizeForResize,0.);
      ////m_grads.clear();
    }

    ////m_times[loopId] = currentTime;
    ////m_temps[loopId] = currentTemp;
    m_ws   [loopId] = currentW[0];

    if (referenceW) {
      while ((misfitId           < refSize            ) &&
             (previousTemp       <= refTemps[misfitId]) &&
             (refTemps[misfitId] <= currentTemp       )) {
        double tmpValue = (refTemps[misfitId]-previousTemp)*(currentW[0]-previousW[0])/(currentTemp-previousTemp) + previousW[0];
        double diff = tmpValue - refValues[misfitId];

        sumValue += diff*diff/refVariances[misfitId]; //COMPATIBILITY WITH OLD VERSION
#if 0
        if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
          std::cout << "In uqTgaWClass<P_V,P_M>::computeUsingTemp()"
                    << ", misfitId = "     << misfitId
                    << ", measuredTemp = " << refTemps[misfitId]
                    << ", measuredW = "    << refValues[misfitId]
                    << ": computedW = "    << tmpValue
                    << ", diffValue = "    << diff
                    << std::endl;
        }
#endif
        misfitId++;
      }
    }

    previousTemp = currentTemp;
    previousW[0] = currentW[0];

    if (referenceW) {
      continueOnWhile = (currentTemp < maximumTemp) && (misfitId < refSize);
    }
    else if (maximumTemp > 0) {
      continueOnWhile = (currentTemp < maximumTemp); // COMPATIBILITY WITH OLD VERSION
    }
    else {
      continueOnWhile = (currentW[0] > 0.);
    }
  }
  if (misfitValue) *misfitValue = sumValue;

  ////m_times.resize(loopId+1);
  ////m_temps.resize(loopId+1);
  m_ws.resize   (loopId+1);
  ////m_grads.clear();
#if 0
  if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
    char stringA[64];
    char stringE[64];
    sprintf(stringA,"%12.6e",params[0]);
    sprintf(stringE,"%12.6e",params[1]);
    std::cout << "In uqTgaWClass<P_V,P_M>::computeUsingTemp()"
              << ", A = "                          << stringA
              << ", E = "                          << stringE
              << ", beta = "                       << m_temperatureFunctionObj.deriv(0.)
              << ": finished ode loop after "      << m_times.size()
              << " iterations, with final time = " << m_times[m_times.size()-1]
              << ", final temp = "                 << m_temps[m_temps.size()-1]
              << ", final w = "                    << m_ws   [m_ws.size()   -1]
              << std::endl;
    if (misfitValue) {
      std::cout << " and with refTimes[max] = "  << refTimes[refSize-1]
                << ", misfitValue = " << *misfitValue
                << std::endl;
    }
  }
#endif
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}
#endif

template<class P_V, class P_M>
void
uqTgaWClass<P_V,P_M>::interpolate(
  double        time,
  unsigned int& startingTimeId, // input and output
  double*       wValue,
  P_V*          wGrad,
  bool*         timeWasMatchedExactly) const
{
  unsigned int tmpSize = m_times.size(); // Yes, 'm_grads'
  //std::cout << "In uqTgaWClass<P_V,P_M>::grad()"
  //          << ": time = "           << time
  //          << ", m_times.size() = " << tmpSize
  //          << ", m_times[0] = "     << m_times[0]
  //          << ", m_times[max] = "   << m_times[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      m_env.rank(),
                      "uqTgaW<P_V,P_M>::grad()",
                      "m_times.size() = 0");

  UQ_FATAL_TEST_MACRO(startingTimeId >= tmpSize,
                      m_env.rank(),
                      "uqTgaW<P_V,P_M>::grad()",
                      "startingTimeId is too big");

  UQ_FATAL_TEST_MACRO(time < m_times[0],
                      m_env.rank(),
                      "uqTgaW<P_V,P_M>::grad()",
                      "time < m_times[0]");

  UQ_FATAL_TEST_MACRO(m_times[tmpSize-1] < time,
                      m_env.rank(),
                      "uqTgaW<P_V,P_M>::grad()",
                      "m_times[max] < time");

  UQ_FATAL_TEST_MACRO(wGrad && (m_grads[0] == NULL),
                      m_env.rank(),
                      "uqTgaW<P_V,P_M>::grad()",
                      "m_grads[0] == NULL");

  unsigned int i = 0;
  for (i = startingTimeId; i < tmpSize; ++i) {
    if (time <= m_times[i]) break;
  }
  startingTimeId = i;

  if (time == m_times[i]) {
    if (timeWasMatchedExactly) *timeWasMatchedExactly = true;
    if (wValue) *wValue = m_ws[i];
    if (wGrad)  *wGrad  = *(m_grads[i]);
  }
  else {
    if (timeWasMatchedExactly) *timeWasMatchedExactly = false;
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
uqTgaWClass<P_V,P_M>::times() const
{
  return m_times;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaWClass<P_V,P_M>::ws() const
{
  return m_ws;
}

template<class P_V, class P_M>
const std::vector<P_V*>&
uqTgaWClass<P_V,P_M>::grads() const
{
  return m_grads;
}

template<class P_V, class P_M>
const uqBaseEnvironmentClass&
uqTgaWClass<P_V,P_M>::env() const
{
  return m_env;
}

template<class P_V, class P_M>
void
uqTgaWClass<P_V,P_M>::printForMatlab(
  std::ofstream&     ofs,
  const std::string& prefixName) const
{
  unsigned int tmpSize = m_times.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofs << "\n" << prefixName << "Time = zeros(" << tmpSize << ",1);"
        << "\n" << prefixName << "W = zeros("    << tmpSize << ",1);";
  }
  else {
    ofs << "\n" << prefixName << "Time = zeros(" << tmpSize << ",1);"
        << "\n" << prefixName << "W = zeros("    << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofs << "\n" << prefixName << "Time(" << i+1 << ",1) = " << m_times[i] << ";"
          << "\n" << prefixName << "W("    << i+1 << ",1) = " << m_ws   [i] << ";";
    }
  }

  return;
}

#endif // __UQ_TGA_COMPUTABLE_W_H__
