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

#ifndef __UQ_1D_1D_FUNCTION_H__
#define __UQ_1D_1D_FUNCTION_H__

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

           template<class V, class M>
	     void quadPtsWeights(unsigned int numIntervals, bool forceWeightsToSum1, V& quadPositions, V& quadWeights) const;

protected:
  double m_minDomainValue;
  double m_maxDomainValue;
};

template<class V, class M>
void
uqBase1D1DFunctionClass::quadPtsWeights(
  unsigned int numIntervals,
  bool         forceWeightsToSum1,
  V&           quadPositions,
  V&           quadWeights) const
{
  const uqBaseEnvironmentClass& env = quadPositions.env();
  unsigned int n = quadPositions.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      UQ_UNAVAILABLE_RANK,
                      "uqBase1D1DFunctionClass::quadPtsWeights()",
                      "invalid input vector size");

  UQ_FATAL_TEST_MACRO((quadPositions.sizeLocal() != quadWeights.sizeLocal()),
                      UQ_UNAVAILABLE_RANK,
                      "uqBase1D1DFunctionClass::quadPtsWeights()",
                      "different input vector sizes");

  double delta = (m_maxDomainValue-m_minDomainValue)/((double) numIntervals);
  std::vector<double> samplePositions(numIntervals,0.);
  std::vector<double> sampleValues   (numIntervals,0.);
  for (unsigned int i = 0; i < numIntervals; ++i) {
    samplePositions[i] = m_minDomainValue + (.5 + ((double) i))*delta;
    sampleValues   [i] = this->value(samplePositions[i]);
    UQ_FATAL_TEST_MACRO((sampleValues[i] < 0.),
                        UQ_UNAVAILABLE_RANK,
                        "uqBase1D1DFunctionClass::quadPtsWeights()",
                        "sampleValue is negative");
  }

  // Compute alphas and betas
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In uqBase1D1DFunctionClass::quadPtsWeights()"
                          << ": computing alphas and betas ..."
                          << std::endl;
  }
  std::vector<double> pi_m1(numIntervals,0.);
  std::vector<double> pi_0 (numIntervals,1.); // Yes, '1'
  std::vector<double> pi_p1(numIntervals,0.);
  double pi_pi_m1 = 0.;

  V alpha(quadPositions);
  V beta (quadPositions);
  for (unsigned int k = 0; k < n; ++k) {
    double pi_pi_0   = 0.;
    double t_pi_pi_0 = 0.;
    for (unsigned int i = 0; i < numIntervals; ++i) {
      pi_pi_0   += delta                      * pi_0[i] * pi_0[i] * sampleValues[i];
      t_pi_pi_0 += delta * samplePositions[i] * pi_0[i] * pi_0[i] * sampleValues[i];
    }

    // Alpha and beta
    alpha[k] = t_pi_pi_0/pi_pi_0;
    if (k == 0) {
      beta[k] = 0.;
      for (unsigned int i = 0; i < numIntervals; ++i) {
        beta[k] += delta * sampleValues[i];
      }
    }
    else {
      beta[k] = pi_pi_0/pi_pi_m1;
    }
    UQ_FATAL_TEST_MACRO((beta[k] < 0.),
                        UQ_UNAVAILABLE_RANK,
                        "uqBase1D1DFunctionClass::quadPtsWeights()",
                        "beta is negative");

    // Prepare for next k
    for (unsigned int i = 0; i < numIntervals; ++i) {
      pi_p1[i] = (samplePositions[i] - alpha[k])*pi_0[i] - beta[k]*pi_m1[i];
    }
    pi_m1    = pi_0;
    pi_0     = pi_p1;
    pi_pi_m1 = pi_pi_0;
  }
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In uqBase1D1DFunctionClass::quadPtsWeights()";
    for (unsigned int k = 0; k < n; ++k) {
      *env.subDisplayFile() << "\n alpha[" << k << "] = " << alpha[k]
                            << ", beta["   << k << "] = " << beta[k];
    }
    *env.subDisplayFile() << std::endl;
  }

  // Form mJ
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In uqBase1D1DFunctionClass::quadPtsWeights()"
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
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In uqBase1D1DFunctionClass::quadPtsWeights():"
                          << "\n  mJ = " << mJ
                          << std::endl;
  }

  // Compute eigen values and vectors of mJ
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In uqBase1D1DFunctionClass::quadPtsWeights()"
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
    *env.subDisplayFile() << "IMPORTANT: In uqBase1D1DFunctionClass::quadPtsWeights()"
                          << ", before consulting 'forceWeightsToSum1'"
                          << ": forceWeightsToSum1 = " << forceWeightsToSum1
                          << ", beta[0] = "            << beta[0]
                          << ", wSum = "               << wSum
                          << "\n  eigenValues = "      << eigenValues
                          << "\n  eigenVectors = "     << eigenVectors
                          << std::endl;
  }

  if (forceWeightsToSum1) {
    for (unsigned int k = 0; k < n; ++k) {
      quadWeights[k] /= wSum;
    }
  }

  return;
}

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

#endif // __UQ_1D_1D_FUNCTION_H__

