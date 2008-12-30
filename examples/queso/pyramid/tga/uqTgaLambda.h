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

#include <uqTgaComputableW.h>
#include <uqTgaDefines.h>
#include <uqTgaStorage.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

// The "Lambda dot" function
template<class P_V,class P_M>
struct
uqTgaLambdaInfoStruct
{
  double                                A;
  double                                E;
  const uqBase1D1DFunctionClass*        temperatureFunctionObj;
  bool                                  computeGradAlso;
  const uqTgaStorageClass<P_V,P_M>*     weigthedMisfitData;
  const uqTgaComputableWClass<P_V,P_M>* wObj;
};

template<class P_V,class P_M>
int uqTgaLambdaTildeDotWrtTimeRoutine(double timeTilde, const double lambdaTilde[], double f[], void *voidInfo)
{
  const uqTgaLambdaInfoStruct<P_V,P_M>& info = *((uqTgaLambdaInfoStruct<P_V,P_M> *)voidInfo);
  double A = info.A;
  double E = info.E;

  if (info.weigthedMisfitData->dataIsContinuousWithTime()) {
    unsigned int tmpSize    = info.weigthedMisfitData->times().size();
    double maximumTimeTilde = info.weigthedMisfitData->times()[tmpSize-1];
    double equivalentTime   = maximumTimeTilde - timeTilde;
    double temp             = info.temperatureFunctionObj->value(equivalentTime);
    double expTerm          = exp(-E/(R_CONSTANT*temp));
    double weigthedMisfit   = info.weigthedMisfitData->value(equivalentTime);

    if ((info.wObj->env().verbosity() >= 99) && (info.wObj->env().rank() == 0)) {
      std::cout << "In uqTgaLambdaTildeDotWrtTimeRoutine()"
                << ", continuous case"
                << ": timeTilde = "      << timeTilde
                << ", equivalentTime = " << equivalentTime
                << ", weigthedMisfit = " << weigthedMisfit
               << std::endl;
    }

    if (info.computeGradAlso) {
      P_V wGrad(*(info.wObj->grads()[0]));
      info.wObj->grad(equivalentTime,wGrad);
      double wA = wGrad[0];
      double wE = wGrad[1];

      f[0] =  -A*lambdaTilde[0]                                      *expTerm - 2*weigthedMisfit;
      f[1] = (-A*lambdaTilde[1] -   lambdaTilde[0]                  )*expTerm - 2*wA;
      f[2] = (-A*lambdaTilde[2] + A*lambdaTilde[0]/(R_CONSTANT*temp))*expTerm - 2*wE;
    }
    else {
      f[0] = (-A*lambdaTilde[0]*exp(-E/(R_CONSTANT*temp)))-2*weigthedMisfit;
    }
  }
  else {
    // AQUI
    UQ_FATAL_TEST_MACRO(true,
                        UQ_UNAVAILABLE_RANK,
                        "uqTgaLambdaTildeDotWrtTimeRoutine(), discrete case",
                        "INCOMPLETE CODE");
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

        void                 compute(const P_V&                            params,
                                     bool                                  computeGradAlso,
                                     const uqTgaStorageClass<P_V,P_M>&     weigthedMisfitData,
                                     const uqTgaComputableWClass<P_V,P_M>& wObj);
        double               lambda (double time) const;
  const P_V&                 grad   (double time) const;
  const std::vector<double>& times  () const;
  const std::vector<double>& temps  () const;
  const std::vector<double>& lambdas() const;
  const std::vector<P_V*  >& grads  () const;

protected:
        void                 resetInternalValues();

  const uqBaseEnvironmentClass&  m_env;
  const uqBase1D1DFunctionClass& m_temperatureFunctionObj;

        std::vector<double>      m_times;
        std::vector<double>      m_temps;
        std::vector<double>      m_lambdas;
        std::vector<P_V*  >      m_grads;
};

template<class P_V, class P_M>
uqTgaLambdaClass<P_V,P_M>::uqTgaLambdaClass(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const uqBase1D1DFunctionClass&     temperatureFunctionObj)
  :
  m_env                   (paramSpace.env()),
  m_temperatureFunctionObj(temperatureFunctionObj),
  m_times                 (0),
  m_temps                 (0),
  m_lambdas               (0),
  m_grads                 (0)
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
  m_temps.clear();
  m_lambdas.clear();
  for (unsigned int i = 0; i < m_grads.size(); ++i) {
    delete m_grads[i];
  }
  m_grads.clear();
}

template<class P_V, class P_M>
void
uqTgaLambdaClass<P_V,P_M>::compute(
  const P_V&                            params,
  bool                                  computeGradAlso,
  const uqTgaStorageClass<P_V,P_M>&     weigthedMisfitData,
  const uqTgaComputableWClass<P_V,P_M>& wObj)
{
  this->resetInternalValues();

  std::vector<double> timesTilde  (1000,0.  );
  std::vector<double> tempsTilde  (1000,0.  );
  std::vector<double> lambdasTilde(1000,0.  );
  std::vector<P_V*  > gradsTilde  (1000,NULL);

  uqTgaLambdaInfoStruct<P_V,P_M> lambdaDotWrtTimeInfo = {params[0],
                                                         params[1],
                                                         &m_temperatureFunctionObj,
                                                         computeGradAlso,
                                                         &weigthedMisfitData,
                                                         &wObj};

  unsigned int numLambdaComponents = 1;
  if (computeGradAlso) numLambdaComponents = 3;

  // Integration
  const gsl_odeiv_step_type *st = gsl_odeiv_step_rkf45; //rkf45; //gear1;
        gsl_odeiv_step      *s  = gsl_odeiv_step_alloc(st,numLambdaComponents);
        gsl_odeiv_control   *c  = gsl_odeiv_control_y_new(GSL_ODE_CONTROL_ABS_PRECISION_FOR_QUESO,0.0);
        gsl_odeiv_evolve    *e  = gsl_odeiv_evolve_alloc(numLambdaComponents);
        gsl_odeiv_system     sysTime = {uqTgaLambdaTildeDotWrtTimeRoutine<P_V,P_M>, NULL, numLambdaComponents, (void *)&lambdaDotWrtTimeInfo};

  unsigned int tmpSize    = weigthedMisfitData.times().size();
  double maximumTimeTilde = weigthedMisfitData.times()[tmpSize-1];
  double currentTimeTilde = 0.;
  double equivalentTime   = maximumTimeTilde - currentTimeTilde;
  double deltaTimeTilde   = 5.;

  double currentLambdaTilde[numLambdaComponents];
  currentLambdaTilde[0]=0.;
  if (computeGradAlso) {
    currentLambdaTilde[1]=0.;
    currentLambdaTilde[2]=0.;
  }

  unsigned int loopId = 0;
  timesTilde  [loopId] = currentTimeTilde;
  tempsTilde  [loopId] = m_temperatureFunctionObj.value(equivalentTime);
  lambdasTilde[loopId] = currentLambdaTilde[0];
  if (computeGradAlso) {
    gradsTilde[loopId] = new P_V(params);
    (*(gradsTilde[loopId]))[0] = currentLambdaTilde[1];
    (*(gradsTilde[loopId]))[1] = currentLambdaTilde[2];
  }

  unsigned int misfitId = weigthedMisfitData.times().size()-1;
  while (currentTimeTilde < maximumTimeTilde) {
    int status = 0;
    if (weigthedMisfitData.dataIsContinuousWithTime()) {
      double nextTimeTilde = std::min(maximumTimeTilde,currentTimeTilde+deltaTimeTilde);
      status = gsl_odeiv_evolve_apply(e, c, s, &sysTime, &currentTimeTilde, nextTimeTilde, &deltaTimeTilde, currentLambdaTilde);
      UQ_FATAL_TEST_MACRO((status != GSL_SUCCESS),
                          params.env().rank(),
                          "uqTgaComputableWClass<P_V,P_M>::compute()",
                          "gsl_odeiv_evolve_apply() failed");
      deltaTimeTilde = std::min(5.,deltaTimeTilde);
      if ((m_env.verbosity() >= 99) && (m_env.rank() == 0)) {
        std::cout << "In uqTgaLambdaClass::compute()"
                  << ": currentTimeTilde = "      << currentTimeTilde
                  << ", currentLambdaTilde[0] = " << currentLambdaTilde[0]
                  << std::endl;
      }
    }
    else {
      // AQUI
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqTgaLambdaClass<P_V,P_M>::compute()",
                          "INCOMPLETE CODE");
      misfitId--;
    }
    equivalentTime = maximumTimeTilde - currentTimeTilde;

    loopId++;
    if (loopId >= timesTilde.size()) {
      timesTilde.resize  (timesTilde.size()  +1000,0.  );
      tempsTilde.resize  (tempsTilde.size()  +1000,0.  );
      lambdasTilde.resize(lambdasTilde.size()+1000,0.  );
      gradsTilde.resize  (gradsTilde.size()  +1000,NULL);
    }

    timesTilde  [loopId] = currentTimeTilde;
    tempsTilde  [loopId] = m_temperatureFunctionObj.value(equivalentTime);
    lambdasTilde[loopId] = currentLambdaTilde[0];
    if (computeGradAlso) {
      gradsTilde[loopId] = new P_V(params);
      (*(gradsTilde[loopId]))[0] = currentLambdaTilde[1];
      (*(gradsTilde[loopId]))[1] = currentLambdaTilde[2];
    }
  }

  timesTilde.resize  (loopId+1);
  tempsTilde.resize  (loopId+1);
  lambdasTilde.resize(loopId+1);
  gradsTilde.resize  (loopId+1);

  m_times.resize  (loopId+1,0.  );
  m_temps.resize  (loopId+1,0.  );
  m_lambdas.resize(loopId+1,0.  );
  m_grads.resize  (loopId+1,NULL);
  for (unsigned int i = 0; i <= loopId; ++i) {
    unsigned int iTilde = loopId-i;
    m_times  [i] = maximumTimeTilde - timesTilde[iTilde];
    m_temps  [i] = tempsTilde  [iTilde];
    m_lambdas[i] = lambdasTilde[iTilde];
    if (computeGradAlso) {
      m_grads  [i] = new P_V(*(gradsTilde[iTilde]));
      delete gradsTilde[iTilde];
    }
  }

  if ((m_env.verbosity() >= 0) && (m_env.rank() == 0)) {
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
              << ", final temp = "                 << m_temps  [m_temps.size()  -1]
              << ", final lambdas = "              << m_lambdas[m_lambdas.size()-1]
              << std::endl;
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free   (s);

  return;
}

template<class P_V, class P_M>
double
uqTgaLambdaClass<P_V,P_M>::lambda(double time) const
{
  double value = 0.;

  // AQUI
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqTgaLambdaClass<P_V,P_M>::lambda()",
                      "INCOMPLETE CODE");

  return value;
}

template<class P_V, class P_M>
const P_V&
uqTgaLambdaClass<P_V,P_M>::grad(double time) const
{
  // AQUI
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqTgaLambdaClass<P_V,P_M>::grad()",
                      "INCOMPLETE CODE");

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
uqTgaLambdaClass<P_V,P_M>::temps() const
{
  return m_temps;
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
