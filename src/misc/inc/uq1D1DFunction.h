//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_1D_1D_FUNCTION_H__
#define __UQ_1D_1D_FUNCTION_H__

#include <uqScalarSequence.h>
#include <uq1DQuadrature.h>
#include <uqEnvironment.h>
#include <uqDefines.h>
#include <vector>
#include <math.h>
#include <fstream>

//*****************************************************
// Base 1D->1D class
//*****************************************************
class
uqBase1D1DFunctionClass {
public:
           uqBase1D1DFunctionClass(double minDomainValue,
                                   double maxDomainValue);
  virtual ~uqBase1D1DFunctionClass();

           double minDomainValue() const;
           double maxDomainValue() const;
  virtual  double value         (double domainValue) const = 0;
  virtual  double deriv         (double domainValue) const = 0;
  virtual  double multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const;

protected:
  double m_minDomainValue;
  double m_maxDomainValue;
};

//*****************************************************
// Generic 1D->1D class
//*****************************************************
class uqGeneric1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqGeneric1D1DFunctionClass(double minDomainValue,
                             double maxDomainValue,
                             double (*valueRoutinePtr)(double domainValue, const void* routinesDataPtr),
                             double (*derivRoutinePtr)(double domainValue, const void* routinesDataPtr),
                             const void* routinesDataPtr);
 ~uqGeneric1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double (*m_valueRoutinePtr)(double domainValue, const void* routinesDataPtr);
  double (*m_derivRoutinePtr)(double domainValue, const void* routinesDataPtr);
  const void* m_routinesDataPtr;
};

//*****************************************************
// Constant 1D->1D class
//*****************************************************
class uqConstant1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqConstant1D1DFunctionClass(double minDomainValue,
                              double maxDomainValue,
                              double constantValue);
 ~uqConstant1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double m_constantValue;
};

//*****************************************************
// Linear 1D->1D class
//*****************************************************
class uqLinear1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqLinear1D1DFunctionClass(double minDomainValue,
                            double maxDomainValue,
                            double referenceDomainValue,
                            double referenceImageValue,
                            double rateValue);
 ~uqLinear1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double m_referenceDomainValue;
  double m_referenceImageValue;
  double m_rateValue;
};

//*****************************************************
// PiecewiseLinear 1D->1D class
//*****************************************************
class uqPiecewiseLinear1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqPiecewiseLinear1D1DFunctionClass(double minDomainValue,
                                     double maxDomainValue,
                                     const std::vector<double>& referenceDomainValues,
                                     double referenceImageValue0,
                                     const std::vector<double>& rateValues);
 ~uqPiecewiseLinear1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  unsigned int        m_numRefValues;
  std::vector<double> m_referenceDomainValues;
  std::vector<double> m_referenceImageValues;
  std::vector<double> m_rateValues;
};

//*****************************************************
// Quadratic 1D->1D class
//*****************************************************
class uqQuadratic1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqQuadratic1D1DFunctionClass(double minDomainValue,
                               double maxDomainValue,
                               double a,
                               double b,
                               double c);
 ~uqQuadratic1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double m_a;
  double m_b;
  double m_c;
};

//*****************************************************
// Sampled 1D->1D class
//*****************************************************
class uqSampled1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqSampled1D1DFunctionClass();
  uqSampled1D1DFunctionClass(const std::vector<double>& domainValues,
                             const std::vector<double>& imageValues);
  virtual ~uqSampled1D1DFunctionClass();

  virtual double               value(double domainValue) const;
          double               deriv(double domainValue) const;

  const   std::vector<double>& domainValues() const;
  const   std::vector<double>& imageValues () const;
          bool                 domainValueMatchesExactly(double domainValue) const;
          void                 set(const std::vector<double>& domainValues,
                                   const std::vector<double>& imageValues);

  virtual void                 printForMatlab(const uqBaseEnvironmentClass& env, std::ofstream& ofsvar, const std::string& prefixName) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_domainValues;
  std::vector<double> m_imageValues;
};

//*****************************************************
// Delta Sequence 1D->1D class
//*****************************************************
class uqDeltaSeq1D1DFunctionClass : public uqSampled1D1DFunctionClass {
public:
  uqDeltaSeq1D1DFunctionClass();
  uqDeltaSeq1D1DFunctionClass(const std::vector<double>& domainValues,
                              const std::vector<double>& imageValues,
                              const std::vector<double>& integratedValues);
 ~uqDeltaSeq1D1DFunctionClass();

  double                     value(double domainValue) const;
  const std::vector<double>& integratedValues() const;
  void                       set(const std::vector<double>& domainValues,
                                 const std::vector<double>& imageValues,
                                 const std::vector<double>& integratedValues);
  void                       printForMatlab(const uqBaseEnvironmentClass& env, std::ofstream& ofsvar, const std::string& prefixName) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;
  using uqSampled1D1DFunctionClass::m_domainValues;
  using uqSampled1D1DFunctionClass::m_imageValues;

  std::vector<double> m_integratedValues;
};

//*****************************************************
// Gaussian 1d Kde 1D->1D class
//*****************************************************
class uqGaussian1dKde1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqGaussian1dKde1D1DFunctionClass(const uqScalarSequenceClass<double>* chain,
                                   double chainMin,
                                   double chainMax,
                                   double gaussian1DScale);
 ~uqGaussian1dKde1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;
  double multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqScalarSequenceClass<double>* m_chain;
  double m_chainMin;
  double m_chainMax;
  double m_gaussian1DScale;

  mutable int                 m_lastInputK;
  mutable unsigned int        m_lastQuadratureOrder;
  mutable std::vector<double> m_level_0;
  mutable std::vector<double> m_level_m1;
};

//*****************************************************
// 'ScalarTimesFunc' 1D->1D class
//*****************************************************
class uqScalarTimesFunc1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqScalarTimesFunc1D1DFunctionClass(double scalar,
                                     const uqBase1D1DFunctionClass& func);
 ~uqScalarTimesFunc1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double m_scalar;
  const uqBase1D1DFunctionClass& m_func;
};

//*****************************************************
// 'FuncTimesFunc' 1D->1D class
//*****************************************************
class uqFuncTimesFunc1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqFuncTimesFunc1D1DFunctionClass(const uqBase1D1DFunctionClass& func1,
                                   const uqBase1D1DFunctionClass& func2);
 ~uqFuncTimesFunc1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func1;
  const uqBase1D1DFunctionClass& m_func2;
};

//*****************************************************
// 'FuncPlusFunc' 1D->1D class
//*****************************************************
class uqFuncPlusFunc1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqFuncPlusFunc1D1DFunctionClass(const uqBase1D1DFunctionClass& func1,
                                  const uqBase1D1DFunctionClass& func2);
 ~uqFuncPlusFunc1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func1;
  const uqBase1D1DFunctionClass& m_func2;
};

//*****************************************************
// 'QuadDenominator' 1D->1D class
//*****************************************************
class uqQuadDenominator1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqQuadDenominator1D1DFunctionClass(const uqBase1D1DFunctionClass& func);
 ~uqQuadDenominator1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;
  const uqBase1D1DFunctionClass& constitutiveFunction() const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func;
};

//*****************************************************
// 'QuadNumerator' 1D->1D class
//*****************************************************
class uqQuadNumerator1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqQuadNumerator1D1DFunctionClass(const uqBase1D1DFunctionClass& func);
 ~uqQuadNumerator1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func;
};

//*****************************************************
// 'QuadCRecursion' 1D->1D class
//*****************************************************
class uqQuadCRecursion1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqQuadCRecursion1D1DFunctionClass(int k);
  uqQuadCRecursion1D1DFunctionClass(int                        k,
                                    const std::vector<double>& alpha,
                                    const std::vector<double>& beta);
 ~uqQuadCRecursion1D1DFunctionClass();

        double               value(double domainValue) const;
        double               deriv(double domainValue) const;
        int                  k    ()                   const;
  const std::vector<double>& alpha()                   const;
  const std::vector<double>& beta ()                   const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;
  int                 m_k;
  std::vector<double> m_alpha;
  std::vector<double> m_beta;
};

//*****************************************************
// Templated isolated function
//*****************************************************
template<class V>
void
computeAlphasAndBetasWithRiemann(
  const uqBase1D1DFunctionClass& function1D1D,
  unsigned int numIntegrationPoints,
  V&           alpha,
  V&           beta)
{
  unsigned int n = alpha.sizeLocal();
  UQ_FATAL_TEST_MACRO((n < 1),
                      UQ_UNAVAILABLE_RANK,
                      "computeAlphasAndBetasWithRiemann()",
                      "invalid n");

  std::vector<double> samplePositions(numIntegrationPoints,0.);
  std::vector<double> sampleValues   (numIntegrationPoints,0.);

  double delta = (function1D1D.maxDomainValue()-function1D1D.minDomainValue())/((double) numIntegrationPoints);
  for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
    samplePositions[i] = function1D1D.minDomainValue() + (.5 + ((double) i))*delta;
    sampleValues   [i] = function1D1D.value(samplePositions[i]);
    UQ_FATAL_TEST_MACRO((sampleValues[i] < 0.),
                        UQ_UNAVAILABLE_RANK,
                        "computeAlphasAndBetasWithRiemann<V>()",
                        "sampleValue is negative");
  }

  std::vector<double> pi_m1(numIntegrationPoints,0.);
  std::vector<double> pi_0 (numIntegrationPoints,1.); // Yes, '1'
  std::vector<double> pi_p1(numIntegrationPoints,0.);
  double pi_pi_m1 = 0.;

  for (unsigned int k = 0; k < n; ++k) {
    double pi_pi_0   = 0.;
    double t_pi_pi_0 = 0.;
    for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
      pi_pi_0   += delta                      * pi_0[i] * pi_0[i] * sampleValues[i];
      t_pi_pi_0 += delta * samplePositions[i] * pi_0[i] * pi_0[i] * sampleValues[i];
    }

    // Alpha and beta
    alpha[k] = t_pi_pi_0/pi_pi_0;
    if (k == 0) {
      beta[k] = 0.;
      for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
        beta[k] += delta * sampleValues[i];
      }
    }
    else {
      beta[k] = pi_pi_0/pi_pi_m1;
    }
    UQ_FATAL_TEST_MACRO((beta[k] < 0.),
                        UQ_UNAVAILABLE_RANK,
                        "computeAlphasAndBetasWithRiemann<V>()",
                        "beta is negative");

    // Prepare for next k
    if (k < (n-1)) {
      for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
        pi_p1[i] = (samplePositions[i] - alpha[k])*pi_0[i] - beta[k]*pi_m1[i];
      }
      pi_m1    = pi_0;
      pi_0     = pi_p1;
      pi_pi_m1 = pi_pi_0;
    }
  }

  return;
}

//*****************************************************
// Isolated function
//*****************************************************
void
alphaBetaCLoop(
  int                                      k,
  double                                   pi_pi_m1,
  const uqQuadCRecursion1D1DFunctionClass& pi_m1,
  const uqQuadCRecursion1D1DFunctionClass& pi_0,
  const uqBase1D1DFunctionClass&           rho,
  unsigned int                             integrationOrder,
  std::vector<double>&                     alpha,
  std::vector<double>&                     beta);

//*****************************************************
// Isolated templated function
//*****************************************************
template<class V>
void
computeAlphasAndBetasWithCQuadrature(
  const uqBase1D1DFunctionClass& function1D1D,
  unsigned int                   integrationOrder,
  V&                             alpha,
  V&                             beta)
{
  uqQuadCRecursion1D1DFunctionClass pi_m1(-1);
  uqQuadCRecursion1D1DFunctionClass pi_0 (0);

  std::vector<double> stdAlpha(alpha.sizeLocal(),0.);
  std::vector<double> stdBeta (beta.sizeLocal(), 0.);
  alphaBetaCLoop(0, // Yes, '0'
                 1.,
                 pi_m1,
                 pi_0,
                 function1D1D,
                 integrationOrder,
                 stdAlpha,
                 stdBeta);

  for (unsigned int i = 0; i < alpha.sizeLocal(); ++i) {
    alpha[i] = stdAlpha[i];
    beta [i] = stdBeta [i];
  }

  return;
}

//#include <uq1D1DFunction2.h>

//*****************************************************
// Isolated templated function
//*****************************************************
template<class V, class M>
void
computeQuadPtsAndWeights(const uqBase1D1DFunctionClass& function1D1D,
                         unsigned int integrationMethod,
                         unsigned int integrationOrder,
                         double       exactValueToForceOnSumOfWeights,
                         V&           quadPositions,
                         V&           quadWeights)
{
  const uqBaseEnvironmentClass& env = quadPositions.env();
  unsigned int n = quadPositions.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      UQ_UNAVAILABLE_RANK,
                      "computeQuadPtsAndWeights<V,M>()",
                      "invalid input vector size");

  UQ_FATAL_TEST_MACRO((quadPositions.sizeLocal() != quadWeights.sizeLocal()),
                      UQ_UNAVAILABLE_RANK,
                      "computeQuadPtsAndWeights<V,M>()",
                      "different input vector sizes");

  // Compute alphas and betas
  V alpha(quadPositions);
  V beta (quadPositions);
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()"
                          << ": computing alphas and betas ..."
                          << std::endl;
  }

  if (integrationMethod == 0) { // Riemann
    computeAlphasAndBetasWithRiemann<V>(function1D1D,
                                        integrationOrder,
                                        alpha,
                                        beta);
  }
  else if (integrationMethod == 1) { // Quadrature
    //computeAlphasAndBetasWithQuadrature<V>   (function1D1D, // IMPORTANT
    //computeAlphasAndBetasWithPolQuadrature<V>(function1D1D, // IMPORTANT
    computeAlphasAndBetasWithCQuadrature<V>(function1D1D, // IMPORTANT
                                            integrationOrder,
                                            alpha,
                                            beta);
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        UQ_UNAVAILABLE_RANK,
                        "computeQuadPtsAndWeights<V,M>()",
                        "invalid 'quadIntegrationMethod'");
  }

  if ((env.displayVerbosity() >= 2) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()";
    for (unsigned int k = 0; k < n; ++k) {
      *env.subDisplayFile() << "\n alpha[" << k << "] = " << alpha[k]
                            << ", beta["   << k << "] = " << beta[k];
    }
    *env.subDisplayFile() << std::endl;
  }

  // Form mJ
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()"
                          << ": forming J..."
                          << std::endl;
  }
  V zeroVector(quadPositions);
  zeroVector.cwSet(0.);
  M mJ(zeroVector);
  for (unsigned int k = 0; k < n; ++k) {
    mJ(k,k) = alpha[k];
    if (mJ.numRowsGlobal() > 1) {
      if (k < (mJ.numRowsGlobal()-1)) {      
        mJ(k,k+1) = sqrt(beta[k+1]);
        mJ(k+1,k) = sqrt(beta[k+1]);
      }
    }
  } // end for 'k'
  if ((env.displayVerbosity() >= 2) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>():"
                          << "\n  mJ = " << mJ
                          << std::endl;
  }

  // Compute eigen values and vectors of mJ
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()"
                          << ": computing eigen stuff of J..."
                          << std::endl;
  }
  V eigenValues(quadPositions);
  M eigenVectors(zeroVector);
  mJ.eigen(eigenValues,&eigenVectors);

  // Prepare output information
  for (unsigned int k = 0; k < n; ++k) {
    quadPositions[k] = eigenValues[k];
    quadWeights  [k] = beta[0] * eigenVectors(0,k) * eigenVectors(0,k);
  } 

  double wSum = 0.;
  for (unsigned int k = 0; k < n; ++k) {
    wSum += quadWeights[k];
  } 

  if ((env.displayVerbosity() >= 2) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "IMPORTANT: In computeQuadPtsAndWeights<V,M>()"
                          << ", before consulting 'exactValueToForceOnSumOfWeights'"
                          << ": exactValueToForceOnSumOfWeights = " << exactValueToForceOnSumOfWeights
                          << ", beta[0] = "                         << beta[0]
                          << ", wSum = "                            << wSum
                          << "\n  eigenValues = "                   << eigenValues
                          << "\n  eigenVectors = "                  << eigenVectors
                          << std::endl;
  }

  if (exactValueToForceOnSumOfWeights > 0.) {
    for (unsigned int k = 0; k < n; ++k) {
      quadWeights[k] *= (exactValueToForceOnSumOfWeights/wSum);
    }
  }

  return;
}

//*****************************************************
// Lagrange Polynomial 1D->1D class
//*****************************************************
class uqLagrangePolynomial1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqLagrangePolynomial1D1DFunctionClass(const std::vector<double>& positionValues,
                                        const std::vector<double>* functionValues);
 ~uqLagrangePolynomial1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_positionValues;
  std::vector<double> m_functionValues;
};

//*****************************************************
// Lagrange Basis 1D->1D class
//*****************************************************
class uqLagrangeBasis1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqLagrangeBasis1D1DFunctionClass(const std::vector<double>& positionValues,
                                   unsigned int basisIndex);
 ~uqLagrangeBasis1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_positionValues;
  unsigned int        m_basisIndex;
};

template <class T>
double
uqSubF1F2Gaussian2dKdeIntegral(const uqScalarSequenceClass<T>& scalarSeq1,
                               const uqScalarSequenceClass<T>& scalarSeq2,
                               unsigned int                    initialPos,
                               double                          scaleValue1,
                               double                          scaleValue2,
                               const uqBase1D1DFunctionClass&  func1,
                               const uqBase1D1DFunctionClass&  func2,
                               unsigned int                    quadratureOrder)
{
  double resultValue = 0.;

  UQ_FATAL_TEST_MACRO(initialPos != 0,
                      scalarSeq1.env().worldRank(),
                      "uqSubF1F2Gaussian2dKdeIntegral()",
                      "not implemented yet for initialPos != 0");
  UQ_FATAL_TEST_MACRO(scalarSeq1.subSequenceSize() != scalarSeq2.subSequenceSize(),
                      scalarSeq1.env().worldRank(),
                      "uqSubF1F2Gaussian2dKdeIntegral()",
                      "different sizes");

  uqGaussianHermite1DQuadratureClass quadObj(0.,1.,quadratureOrder);
  const std::vector<double>& quadPositions = quadObj.positions();
  const std::vector<double>& quadWeights   = quadObj.weights  ();
  UQ_FATAL_TEST_MACRO(quadPositions.size() != quadWeights.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqSubF1F2Gaussian2dKdeIntegral()",
                      "quadObj has invalid state");

  unsigned int numQuadraturePositions = quadPositions.size();
  unsigned int dataSize = scalarSeq1.subSequenceSize();
  for (unsigned int k = 0; k < dataSize; ++k) {
    double value1 = 0.;
    double value2 = 0.;
    double x1k = scalarSeq1[k];
    double x2k = scalarSeq2[k];
    for (unsigned int j = 0; j < numQuadraturePositions; ++j) {
      value1 += func1.value(scaleValue1*quadPositions[j]+x1k)*quadWeights[j];
      value2 += func2.value(scaleValue2*quadPositions[j]+x2k)*quadWeights[j];
    }
    resultValue += value1*value2;
  }
  resultValue *= 1./(2.*M_PI)/((double) dataSize);

  return resultValue;
}
#endif // __UQ_1D_1D_FUNCTION_H__

