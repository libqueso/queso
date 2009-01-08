/* uq/examples/queso/pyramid/uqTgaLambda.h
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

#ifndef __UQ_TGA_LAMBDA_W_H__
#define __UQ_TGA_LAMBDA_W_H__

#include <uqTgaW.h>
#include <uqTgaDefines.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

// The "Lambda dot" function
template<class P_V,class P_M>
struct
uqTgaLambdaInfoStruct
{
  const P_V*                     params;
  const P_V*                     paramDirection;
  const uqBase1D1DFunctionClass* temperatureFunctionObj;
  bool                           computeGradAlso;
  const uqBase1D1DFunctionClass* diffFunction;
  const uqBase1D1DFunctionClass* tildeWeightFunction;
  const uqTgaWClass<P_V,P_M>*    wObj;
  unsigned int                   suggestedWTimeId;
};

template<class P_V,class P_M>
int uqTgaLambdaTildeDotWrtTimeRoutine(double tildeTime, const double lambdaTilde[], double f[], void *voidInfo)
{
  const uqTgaLambdaInfoStruct<P_V,P_M>& info = *((uqTgaLambdaInfoStruct<P_V,P_M> *)voidInfo);
  double A = (*info.params)[0];
  double E = (*info.params)[1];

  double maxTildeTime   = info.diffFunction->maxDomainValue(); // FIX ME
  if (info.tildeWeightFunction) maxTildeTime = info.tildeWeightFunction->maxDomainValue(); // FIX ME

  const uqDeltaSeq1D1DFunctionClass* tildeDeltaSeqFunction = NULL;
  if (info.tildeWeightFunction) {
    tildeDeltaSeqFunction = dynamic_cast< const uqDeltaSeq1D1DFunctionClass* >(info.tildeWeightFunction);
  }

  double equivalentTime = maxTildeTime - tildeTime;
  double temp           = info.temperatureFunctionObj->value(equivalentTime);
  double expTerm        = exp(-E/(R_CONSTANT*temp));
  double tildeDiff      = 0.;
  if (tildeDeltaSeqFunction == NULL) tildeDiff = info.diffFunction->value(equivalentTime); // Might be slow (non delta seq case)

  double weightValue = 1.;
  if (info.tildeWeightFunction) {
    if (tildeDeltaSeqFunction == NULL) weightValue = info.tildeWeightFunction->value(tildeTime); // Might be slow (non delta seq case)
    else                               weightValue = 0.;
  }

  if ((info.wObj->env().verbosity() >= 99) && (info.wObj->env().rank() == 0)) {
    std::cout << "In uqTgaLambdaTildeDotWrtTimeRoutine()"
              << ": tildeTime = "      << tildeTime
              << ", equivalentTime = " << equivalentTime
              << ", tildeDiff = "      << tildeDiff
              << ", weightValue = "    << weightValue
              << std::endl;
  }

  f[0] = -A*lambdaTilde[0]*expTerm - 2*tildeDiff*weightValue;
  if (info.computeGradAlso || info.paramDirection) {
    f[1] = 0.;
    f[2] = 0.;
    f[3] = 0.;
  }

  double wA = 0.;
  double wE = 0.;
  P_V  wGrad(*info.params);
  wGrad *= 0.;
  P_V* wGradPtr   = NULL;

  double wDir = 0.;
  double* wDirPtr = NULL;

  if ((tildeDeltaSeqFunction == NULL) && (info.computeGradAlso || info.paramDirection)) {
    if (info.computeGradAlso) wGradPtr = &wGrad;
    if (info.paramDirection ) wDirPtr  = &wDir;

    // Always reset to this value because 'gsl' algorithm might backward in tildeTime
    unsigned int suggestedTimeId = info.suggestedWTimeId; // For performance on w->interpolate()
    info.wObj->interpolate(equivalentTime,
                           suggestedTimeId,
                           -1., // For performance on w->interpolate()
                           NULL,
                           wGradPtr,
                           wDirPtr,
                           NULL);

    if (info.computeGradAlso) {
      wA = wGrad[0];
      wE = wGrad[1];
      //std::cout << "HERE, equivalentTime = " << equivalentTime
      //          << ", suggId = " << suggestedTimeId
      //          << ", wA = " << wA
      //          << ", wE = " << wE
      //          << std::endl;
    }
  }

  if (info.computeGradAlso) {
    f[1] = (-A*lambdaTilde[1] -   lambdaTilde[0]                  )*expTerm - 2*wA*weightValue;
    f[2] = (-A*lambdaTilde[2] + A*lambdaTilde[0]/(R_CONSTANT*temp))*expTerm - 2*wE*weightValue;
  }

  if (info.paramDirection) {
    double q1 = (*info.paramDirection)[0];
    double q2 = (*info.paramDirection)[1];
    double aux = -q1 + A*q2/(R_CONSTANT*temp);
    f[3] = (-A*lambdaTilde[3] + aux*lambdaTilde[0])*expTerm - 2*wDir*weightValue;
  }

  return GSL_SUCCESS;
}

template<class P_V, class P_M>
class
uqTgaLambdaClass
{
public:
  uqTgaLambdaClass(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                   const uqBase1D1DFunctionClass&     temperatureFunctionObj);
 ~uqTgaLambdaClass();

        void                 compute(const P_V&                     params,
                                     const P_V*                     paramDirection,
                                     double                         maxTimeStep,
                                     bool                           computeGradAlso,
                                     const uqBase1D1DFunctionClass& diffFunction,
                                     const uqBase1D1DFunctionClass* tildeWeightFunction,
                                     const uqTgaWClass<P_V,P_M>&    wObj);
        void                 interpolate(double        time,
                                         unsigned int& startingTimeId,
                                         double*       lambdaValue,
                                         P_V*          lambdaGrad,
                                         double*       lambdaDir,
                                         bool*         timeWasMatchedExactly) const;
  const std::vector<double>& times  () const;
  const std::vector<double>& lambdas() const;
  const std::vector<P_V*  >& grads  () const;

protected:
        void                 resetInternalValues();

  const uqBaseEnvironmentClass&      m_env;
  const uqVectorSpaceClass<P_V,P_M>& m_paramSpace;
  const uqBase1D1DFunctionClass&     m_temperatureFunctionObj;

        std::vector<double>          m_times;
        std::vector<double>          m_lambdas;
        std::vector<P_V*  >          m_grads;
	std::vector<double>          m_lambdaDirs;
};

template<class P_V, class P_M>
uqTgaLambdaClass<P_V,P_M>::uqTgaLambdaClass(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const uqBase1D1DFunctionClass&     temperatureFunctionObj)
  :
  m_env                   (paramSpace.env()),
  m_paramSpace            (paramSpace),
  m_temperatureFunctionObj(temperatureFunctionObj),
  m_times                 (0),
  m_lambdas               (0),
  m_grads                 (0),
  m_lambdaDirs            (0)
{
  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "Entering uqTgaLambdaClass::constructor()"
              << std::endl;
  }

  if ((m_env.verbosity() >= 30) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqTgaLambdaClass::constructor()"
              << std::endl;
  }
}

template<class P_V, class P_M>
uqTgaLambdaClass<P_V,P_M>::~uqTgaLambdaClass()
{
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
}

template<class P_V, class P_M>
void
uqTgaLambdaClass<P_V,P_M>::resetInternalValues()
{
  m_times.clear();
  m_lambdas.clear();
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
  m_grads.clear();
  m_lambdaDirs.clear();
}

template<class P_V, class P_M>
void
uqTgaLambdaClass<P_V,P_M>::compute(
  const P_V&                     params,
  const P_V*                     paramDirection,
  double                         maxTimeStep,
  bool                           computeGradAlso,
  const uqBase1D1DFunctionClass& diffFunction,
  const uqBase1D1DFunctionClass* tildeWeightFunction,
  const uqTgaWClass<P_V,P_M>&    wObj)
{
  if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
    std::cout << "Entering uqTgaLambdaClass<P_V,P_M>::compute()"
              << ": params = "          << params
              << ", paramDirection = "  << paramDirection
              << ", maxTimeStep = "     << maxTimeStep
              << ", computeGradAlso = " << computeGradAlso
              << std::endl;
  }

  // Initialize variables related to the weight function
  const uqDeltaSeq1D1DFunctionClass* tildeDeltaSeqFunction = NULL;
  if (tildeWeightFunction) {
    tildeDeltaSeqFunction = dynamic_cast< const uqDeltaSeq1D1DFunctionClass* >(tildeWeightFunction);
  }
  unsigned int tildeDeltaWeightId = 1; // Yes, '1', not '0'

  // Initialize other variables
  this->resetInternalValues();
  std::vector<double> tildeTimes  (1000,0.);
  std::vector<double> tildeLambdas(1000,0.);
  std::vector<P_V*  > tildeGrads(1000,(P_V*) NULL);
  std::vector<double> tildeLambdaDirs(1000,0.);

  uqTgaLambdaInfoStruct<P_V,P_M> lambdaDotWrtTimeInfo = {&params,
                                                         paramDirection,
                                                         &m_temperatureFunctionObj,
                                                         computeGradAlso,
                                                         &diffFunction,
                                                         tildeWeightFunction,
                                                         &wObj,
                                                         wObj.times().size()-1}; // For performance on w->interpolate()

  unsigned int numLambdaComponents = 1;
  if (computeGradAlso || paramDirection) numLambdaComponents = 4;

  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,numLambdaComponents);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(numLambdaComponents);
        gsl_odeiv_system     sysTime = {uqTgaLambdaTildeDotWrtTimeRoutine<P_V,P_M>, NULL, numLambdaComponents, (void *)&lambdaDotWrtTimeInfo};

  double maxTildeTime     = diffFunction.maxDomainValue(); // FIX ME
  if (tildeWeightFunction) maxTildeTime = tildeWeightFunction->maxDomainValue(); // FIX ME
  double currentTildeTime = 0.;
  double equivalentTime   = maxTildeTime - currentTildeTime;
  double tildeTimeStep    = .1;
  if (maxTimeStep > 0) tildeTimeStep = maxTimeStep;

  double currentLambdaTilde[numLambdaComponents];
  currentLambdaTilde[0]=0.;
  if (computeGradAlso || paramDirection) {
    currentLambdaTilde[1]=0.;
    currentLambdaTilde[2]=0.;
    currentLambdaTilde[3]=0.;
  }

  unsigned int loopId = 0;
  tildeTimes  [loopId] = currentTildeTime;
  tildeLambdas[loopId] = currentLambdaTilde[0];
  if (computeGradAlso) {
    tildeGrads[loopId] = new P_V(m_paramSpace.zeroVector());
    (*(tildeGrads[loopId]))[0] = currentLambdaTilde[1];
    (*(tildeGrads[loopId]))[1] = currentLambdaTilde[2];
  }
  if (paramDirection) tildeLambdaDirs[loopId] = currentLambdaTilde[3];

  unsigned int suggestedWTimeId = wObj.times().size()-1; // For performance on w->interpolate()
  while (currentTildeTime < maxTildeTime) {
    int status = 0;
    double nextTildeTime = maxTildeTime;
    if (maxTimeStep > 0) nextTildeTime = std::min(nextTildeTime,currentTildeTime+maxTimeStep);
    if (tildeDeltaSeqFunction != NULL) {
      nextTildeTime = std::min(nextTildeTime,tildeDeltaSeqFunction->domainValues()[tildeDeltaWeightId]);
      if (currentTildeTime == tildeDeltaSeqFunction->domainValues()[tildeDeltaWeightId-1]) { // Yes, '[...-1]'
        currentLambdaTilde[0] -= 2.*diffFunction.value(equivalentTime); // Might be slow (delta seq case)
        tildeLambdas[loopId] = currentLambdaTilde[0];

        if (computeGradAlso || paramDirection) {
          P_V  wGrad(m_paramSpace.zeroVector());
          P_V* wGradPtr = NULL;

          double wDir = 0.;
          double* wDirPtr = NULL;

          if (computeGradAlso) wGradPtr = &wGrad;
          if (paramDirection ) wDirPtr  = &wDir;

          wObj.interpolate(equivalentTime,
                           suggestedWTimeId,
                           -1., // For performance on w->interpolate()
                           NULL,
                           wGradPtr,
                           wDirPtr,
                           NULL);

          if (computeGradAlso) {
            currentLambdaTilde[1] -= 2.*wGrad[0];
            currentLambdaTilde[2] -= 2.*wGrad[1];
            (*(tildeGrads[loopId]))[0] = currentLambdaTilde[1];
            (*(tildeGrads[loopId]))[1] = currentLambdaTilde[2];
          }

          if (paramDirection) { 
            currentLambdaTilde[3] -= 2.*wDir;
            tildeLambdaDirs[loopId] = currentLambdaTilde[3];
          }
        }
      }
    }

    lambdaDotWrtTimeInfo.suggestedWTimeId = suggestedWTimeId; // For performance on w->interpolate()
    if ((m_env.verbosity() >= 99) && (m_env.rank() == 0)) {
      std::cout << "In uqTgaLambdaClass::compute(), before"
                << ": loopId = "                << loopId
                << ", currentTildeTime = "      << currentTildeTime
                << ", currentLambdaTilde[0] = " << currentLambdaTilde[0]
                << ", equivalentTime = "        << equivalentTime
                << ", nextTildeTime = "         << nextTildeTime
                << ", tildeDeltaWeightId = "    << tildeDeltaWeightId
                << std::endl;
    }
    status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTildeTime, nextTildeTime, &tildeTimeStep, currentLambdaTilde);
    UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                        params.env().rank(),
                        "uqTgaWClass<P_V,P_M>::compute()",
                        "gsl_odeiv_evolve_apply() failed");
    if (maxTimeStep > 0) tildeTimeStep = std::min(maxTimeStep,tildeTimeStep);
    if (tildeDeltaSeqFunction != NULL &&
        (currentTildeTime == tildeDeltaSeqFunction->domainValues()[tildeDeltaWeightId])) {
      tildeDeltaWeightId++;
    }
    equivalentTime = maxTildeTime - currentTildeTime;
    if ((m_env.verbosity() >= 99) && (m_env.rank() == 0)) {
      std::cout << "In uqTgaLambdaClass::compute(), after"
                << ": loopId = "                << loopId
                << ", currentTildeTime = "      << currentTildeTime
                << ", currentLambdaTilde[0] = " << currentLambdaTilde[0]
                << ", equivalentTime = "        << equivalentTime
                << ", nextTildeTime = "         << nextTildeTime
                << ", tildeDeltaWeightId = "    << tildeDeltaWeightId
                << std::endl;
    }

    loopId++;
    if (loopId >= tildeTimes.size()) {
      tildeTimes.resize  (tildeTimes.size()  +1000,0.  );
      tildeLambdas.resize(tildeLambdas.size()+1000,0.  );
      if (computeGradAlso) tildeGrads.resize     (tildeGrads.size()     +1000,NULL);
      if (paramDirection ) tildeLambdaDirs.resize(tildeLambdaDirs.size()+1000,0.  );
    }

    tildeTimes  [loopId] = currentTildeTime;
    tildeLambdas[loopId] = currentLambdaTilde[0];
    if (computeGradAlso) {
      tildeGrads[loopId] = new P_V(m_paramSpace.zeroVector());
      (*(tildeGrads[loopId]))[0] = currentLambdaTilde[1];
      (*(tildeGrads[loopId]))[1] = currentLambdaTilde[2];
    }
    if (paramDirection) tildeLambdaDirs[loopId] = currentLambdaTilde[3];
  }

  tildeTimes.resize  (loopId+1);
  tildeLambdas.resize(loopId+1);
  if (computeGradAlso) tildeGrads.resize     (loopId+1);
  if (paramDirection ) tildeLambdaDirs.resize(loopId+1);

  m_times.resize  (loopId+1,0.);
  m_lambdas.resize(loopId+1,0.);
  if (computeGradAlso) m_grads.resize(loopId+1,NULL);
  if (paramDirection)  m_lambdaDirs.resize(loopId+1,NULL);

  for (unsigned int i = 0; i <= loopId; ++i) {
    unsigned int tildeI = loopId-i;
    m_times  [i] = maxTildeTime - tildeTimes[tildeI];
    m_lambdas[i] = tildeLambdas[tildeI];
    if (computeGradAlso) {
      m_grads[i] = new P_V(*(tildeGrads[tildeI]));
      delete tildeGrads[tildeI];
    }
    if (paramDirection) m_lambdaDirs[i] = tildeLambdaDirs[tildeI];
  }

  if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
    char stringA[64];
    char stringE[64];
    sprintf(stringA,"%12.6e",params[0]);
    sprintf(stringE,"%12.6e",params[1]);
    std::cout << "In uqTgaLambdaClass<P_V,P_M>::compute()"
              << ", A = "                          << stringA
              << ", E = "                          << stringE
              << ", beta = "                       << m_temperatureFunctionObj.deriv(0.)
              << ": finished ode loop after "      << m_times.size()
              << " iterations, with final time = " << m_times  [m_times.size()  -1]
              << ", lambda[0] = "                  << m_lambdas[0]
              << ", final lambda = "               << m_lambdas[m_lambdas.size()-1]
              << std::endl;
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

template<class P_V, class P_M>
void
uqTgaLambdaClass<P_V,P_M>::interpolate(
  double        time,
  unsigned int& startingTimeId, // input and output
  double*       lambdaValue,
  P_V*          lambdaGrad,
  double*       lambdaDir,
  bool*         timeWasMatchedExactly) const
{
  unsigned int tmpSize = m_times.size();
  //std::cout << "In uqTgaLambdaClass<P_V,P_M>::interpolate()"
  //          << ": time = "           << time
  //          << ", m_times.size() = " << tmpSize
  //          << ", m_times[0] = "     << m_times[0]
  //          << ", m_times[max] = "   << m_times[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "m_times.size() = 0");

  UQ_FATAL_TEST_MACRO(startingTimeId >= tmpSize,
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "startingTimeId is too big");

  UQ_FATAL_TEST_MACRO(time < m_times[0],
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "time < m_times[0]");

  UQ_FATAL_TEST_MACRO(m_times[tmpSize-1] < time,
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "m_times[max] < time");

  UQ_FATAL_TEST_MACRO(lambdaGrad && (m_grads.size() == 0),
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "m_grads.size() == 0");

  UQ_FATAL_TEST_MACRO(lambdaGrad && (m_grads[0] == NULL),
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "m_grads[0] == NULL");

  UQ_FATAL_TEST_MACRO(lambdaDir && (m_lambdaDirs.size() == 0),
                      m_env.rank(),
                      "uqTgaLambda<P_V,P_M>::interpolate()",
                      "m_lambdaDirs.size() == 0");

  unsigned int i = 0;
  for (i = startingTimeId; i < tmpSize; ++i) {
    if (time <= m_times[i]) break;
  }
  startingTimeId = i;

  if (time == m_times[i]) {
    if (timeWasMatchedExactly) *timeWasMatchedExactly = true;
    if (lambdaValue) *lambdaValue = m_lambdas   [i];
    if (lambdaGrad)  *lambdaGrad  = *(m_grads   [i]);
    if (lambdaDir)   *lambdaDir   = m_lambdaDirs[i];
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
    if (lambdaValue) *lambdaValue = m_lambdas   [i-1]  + ratio * ( m_lambdas   [i]  - m_lambdas   [i-1]  );
    if (lambdaGrad)  *lambdaGrad  = *(m_grads   [i-1]) + ratio * ( *(m_grads   [i]) - *(m_grads   [i-1]) );
    if (lambdaDir) *lambdaDir     = m_lambdaDirs[i-1]  + ratio * ( m_lambdaDirs[i]  - m_lambdaDirs[i-1]  );
  }

  return;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaLambdaClass<P_V,P_M>::times() const
{
  return m_times;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaLambdaClass<P_V,P_M>::lambdas() const
{
  return m_lambdas;
}

template<class P_V, class P_M>
const std::vector<P_V*>&
uqTgaLambdaClass<P_V,P_M>::grads() const
{
  return m_grads;
}
#endif // __UQ_TGA_LAMBDA_W_H__
