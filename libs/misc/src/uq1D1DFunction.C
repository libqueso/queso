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

#include <uq1D1DFunction.h>
#include <uq1DQuadrature.h>

//*****************************************************
// Base 1D->1D class
//*****************************************************
uqBase1D1DFunctionClass::uqBase1D1DFunctionClass(
  double minDomainValue,
  double maxDomainValue)
  :
  m_minDomainValue(minDomainValue),
  m_maxDomainValue(maxDomainValue)
{
  UQ_FATAL_TEST_MACRO(m_minDomainValue >= m_maxDomainValue,
                      UQ_UNAVAILABLE_RANK,
                      "uqBase1D1DFunctionClass::constructor()",
                      "min >= max");
}

uqBase1D1DFunctionClass::~uqBase1D1DFunctionClass()
{
}

double
uqBase1D1DFunctionClass::minDomainValue() const
{
  return m_minDomainValue;
}

double
uqBase1D1DFunctionClass::maxDomainValue() const
{
  return m_maxDomainValue;
}

double
uqBase1D1DFunctionClass::multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqBase1D1DFunctionClass::multiplyAndIntegrate()",
                      "not implemented yet");

  return value;
}

//*****************************************************
// Generic 1D->1D class
//*****************************************************
uqGeneric1D1DFunctionClass::uqGeneric1D1DFunctionClass(
  double minDomainValue,
  double maxDomainValue,
  double (*valueRoutinePtr)(double domainValue, const void* routinesDataPtr),
  double (*derivRoutinePtr)(double domainValue, const void* routinesDataPtr),
  const void* routinesDataPtr)
  :
  uqBase1D1DFunctionClass(minDomainValue,maxDomainValue),
  m_valueRoutinePtr      (valueRoutinePtr),
  m_derivRoutinePtr      (derivRoutinePtr),
  m_routinesDataPtr      (routinesDataPtr)
{
}

uqGeneric1D1DFunctionClass::~uqGeneric1D1DFunctionClass()
{
}

double
uqGeneric1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqGeneric1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqGeneric1D1DFunctionClass::value()",
                      "x out of range");

  return (*m_valueRoutinePtr)(domainValue,m_routinesDataPtr);
}

double
uqGeneric1D1DFunctionClass::deriv(double domainValue) const
{
  // FIX ME: check range of domainValue
  return (*m_derivRoutinePtr)(domainValue,m_routinesDataPtr);
}

//*****************************************************
// Constant 1D->1D class
//*****************************************************
uqConstant1D1DFunctionClass::uqConstant1D1DFunctionClass(
  double minDomainValue,
  double maxDomainValue,
  double constantValue)
  :
  uqBase1D1DFunctionClass(minDomainValue,maxDomainValue),
  m_constantValue        (constantValue)
{
}

uqConstant1D1DFunctionClass::~uqConstant1D1DFunctionClass()
{
}

double
uqConstant1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqConstant1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqConstant1D1DFunctionClass::value()",
                      "x out of range");

  return m_constantValue;
}

double
uqConstant1D1DFunctionClass::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqConstant1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqConstant1D1DFunctionClass::deriv()",
                      "x out of range");

  return 0.;
}

//*****************************************************
// Linear 1D->1D class
//*****************************************************
uqLinear1D1DFunctionClass::uqLinear1D1DFunctionClass(
  double minDomainValue,
  double maxDomainValue,
  double referenceDomainValue,
  double referenceImageValue,
  double rateValue)
  :
  uqBase1D1DFunctionClass(minDomainValue,maxDomainValue),
  m_referenceDomainValue (referenceDomainValue),
  m_referenceImageValue  (referenceImageValue),
  m_rateValue            (rateValue)
{
}

uqLinear1D1DFunctionClass::~uqLinear1D1DFunctionClass()
{
}

double
uqLinear1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqLinear1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqLinear1D1DFunctionClass::value()",
                      "x out of range");

  double imageValue = m_referenceImageValue + m_rateValue*(domainValue - m_referenceDomainValue);

  return imageValue;
}

double
uqLinear1D1DFunctionClass::deriv(double domainValue) const
{
  // FIX ME: check range of domainValue
  return m_rateValue;
}

//*****************************************************
// Quadratic 1D->1D class
//*****************************************************
uqQuadratic1D1DFunctionClass::uqQuadratic1D1DFunctionClass(
  double minDomainValue,
  double maxDomainValue,
  double a,
  double b,
  double c)
  :
  uqBase1D1DFunctionClass(minDomainValue,maxDomainValue),
  m_a                    (a),
  m_b                    (b),
  m_c                    (c)
{
}

uqQuadratic1D1DFunctionClass::~uqQuadratic1D1DFunctionClass()
{
}

double
uqQuadratic1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqQuadratic1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadratic1D1DFunctionClass::value()",
                      "x out of range");

  double imageValue = m_a*domainValue*domainValue + m_b*domainValue + m_c;

  return imageValue;
}

double
uqQuadratic1D1DFunctionClass::deriv(double domainValue) const
{
  // FIX ME: check range of domainValue
  return 2.*m_a*domainValue + m_b;
}

//*****************************************************
// Sampled 1D->1D class
//*****************************************************
uqSampled1D1DFunctionClass::uqSampled1D1DFunctionClass()
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY)
{
}

uqSampled1D1DFunctionClass::uqSampled1D1DFunctionClass(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues)
  :
  uqBase1D1DFunctionClass(domainValues[0],domainValues[domainValues.size()-1]),
  m_domainValues    (domainValues.size(),0.),
  m_imageValues     (imageValues.size(), 0.)
{
  unsigned int tmpSize = m_domainValues.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_domainValues[i] = domainValues[i];
    m_imageValues [i] = imageValues [i];
  }
}

uqSampled1D1DFunctionClass::~uqSampled1D1DFunctionClass()
{
}

double
uqSampled1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqSampled1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::value()",
                      "x out of range");

  double returnValue = 0.;

  unsigned int tmpSize = m_domainValues.size();
  //std::cout << "In uqSampled1D1DFunctionClass::value()"
  //          << ": domainValue = "         << domainValue
  //          << ", tmpSize = "             << tmpSize
  //          << ", m_domainValues[0] = "   << m_domainValues[0]
  //          << ", m_domainValues[max] = " << m_domainValues[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::value()",
                      "m_domainValues.size() = 0");

  UQ_FATAL_TEST_MACRO(domainValue < m_domainValues[0],
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::value()",
                      "domainValue < m_domainValues[0]");

  UQ_FATAL_TEST_MACRO(m_domainValues[tmpSize-1] < domainValue,
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::value()",
                      "m_domainValues[max] < domainValue");

  unsigned int i = 0;
  for (i = 0; i < tmpSize; ++i) {
    if (domainValue <= m_domainValues[i]) break;
  }

  if (domainValue == m_domainValues[i]) {
    //if (domainValueWasMatchedExactly) *domainValueWasMatchedExactly = true;
    returnValue = m_imageValues[i];
  }
  else {
    //if (domainValueWasMatchedExactly) *domainValueWasMatchedExactly = false;
    double ratio = (domainValue - m_domainValues[i-1])/(m_domainValues[i]-m_domainValues[i-1]);
    returnValue = m_imageValues[i-1] + ratio * (m_imageValues[i]-m_imageValues[i-1]);
  }

  return returnValue;
}

double
uqSampled1D1DFunctionClass::deriv(double domainValue) const
{
  // FIX ME: check range of domainValue
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1d1DFunctionClass::deriv()",
                      "this funciton makes no sense for this class");
  return 0.;
}

const std::vector<double>&
uqSampled1D1DFunctionClass::domainValues() const
{
  return m_domainValues;
}

const std::vector<double>&
uqSampled1D1DFunctionClass::imageValues() const
{
  return m_imageValues;
}

bool
uqSampled1D1DFunctionClass::domainValueMatchesExactly(double domainValue) const
{
  bool result = false;

  unsigned int tmpSize = m_domainValues.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    if (domainValue <= m_domainValues[i]) {
      result = (domainValue == m_domainValues[i]);
      break;
    }
  }

  return result;
}

void
uqSampled1D1DFunctionClass::set(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues)
{
  m_domainValues.clear();
  m_imageValues.clear();

  unsigned int tmpSize = domainValues.size();
  m_minDomainValue = domainValues[0];
  m_maxDomainValue = domainValues[tmpSize-1];

  m_domainValues.resize(tmpSize,0.);
  m_imageValues.resize(tmpSize,0.);
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_domainValues[i] = domainValues[i];
    m_imageValues [i] = imageValues [i];
  }

  return;
}

void
uqSampled1D1DFunctionClass::printForMatlab(
  const uqBaseEnvironmentClass& env,
  std::ofstream&                ofsvar,
  const std::string&            prefixName) const
{
  unsigned int tmpSize = m_domainValues.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofsvar << "\n" << prefixName << "Time_sub" << env.subIdString() << " = zeros("  << tmpSize << ",1);"
           << "\n" << prefixName << "Value_sub" << env.subIdString() << " = zeros(" << tmpSize << ",1);";
  }
  else {
    ofsvar << "\n" << prefixName << "Time_sub" << env.subIdString() << " = zeros("  << tmpSize << ",1);"
           << "\n" << prefixName << "Value_sub" << env.subIdString() << " = zeros(" << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofsvar << "\n" << prefixName << "Time_sub" << env.subIdString() << "("  << i+1 << ",1) = " << m_domainValues[i] << ";"
             << "\n" << prefixName << "Value_sub" << env.subIdString() << "(" << i+1 << ",1) = " << m_imageValues[i]  << ";";
    }
  }

  return;
}

//*****************************************************
// Delta Sequence 1D->1D class
//*****************************************************
uqDeltaSeq1D1DFunctionClass::uqDeltaSeq1D1DFunctionClass()
  :
  uqSampled1D1DFunctionClass()
{
}

uqDeltaSeq1D1DFunctionClass::uqDeltaSeq1D1DFunctionClass(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues,
  const std::vector<double>& integratedValues)
  :
  uqSampled1D1DFunctionClass(domainValues,imageValues),
  m_integratedValues        (integratedValues.size(),0.)
{
  unsigned int tmpSize = m_integratedValues.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_integratedValues[i] = integratedValues[i];
  }
}

uqDeltaSeq1D1DFunctionClass::~uqDeltaSeq1D1DFunctionClass()
{
}

double
uqDeltaSeq1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqDeltaSeq1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqDeltaSeq1D1DFunctionClass::value()",
                      "x out of range");

  double returnValue = 0.;

  unsigned int tmpSize = m_domainValues.size();
  //std::cout << "In uqDeltaSeq1D1DFunctionClass::value()"
  //          << ": domainValue = "         << domainValue
  //          << ", tmpSize = "             << tmpSize
  //          << ", m_domainValues[0] = "   << m_domainValues[0]
  //          << ", m_domainValues[max] = " << m_domainValues[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      UQ_UNAVAILABLE_RANK,
                      "uqDeltaSeq1D1DFunctionClass::value()",
                      "m_domainValues.size() = 0");

  UQ_FATAL_TEST_MACRO(domainValue < m_domainValues[0],
                      UQ_UNAVAILABLE_RANK,
                      "uqDeltaSeq1D1DFunctionClass::value()",
                      "domainValue < m_domainValues[0]");

  UQ_FATAL_TEST_MACRO(m_domainValues[tmpSize-1] < domainValue,
                      UQ_UNAVAILABLE_RANK,
                      "uqDeltaSeq1D1DFunctionClass::value()",
                      "m_domainValues[max] < domainValue");

  unsigned int i = 0;
  for (i = 0; i < tmpSize; ++i) {
    if (domainValue <= m_domainValues[i]) break;
  }

  if (domainValue == m_domainValues[i]) {
    //if (domainValueWasMatchedExactly) *domainValueWasMatchedExactly = true;
    returnValue = m_imageValues[i];
  }
  else {
    // Leave returnValue = 0.;
  }

  return returnValue;
}

const std::vector<double>&
uqDeltaSeq1D1DFunctionClass::integratedValues() const
{
  return m_integratedValues;
}

void
uqDeltaSeq1D1DFunctionClass::set(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues,
  const std::vector<double>& integratedValues)
{
  uqSampled1D1DFunctionClass::set(domainValues,imageValues);
  m_integratedValues.clear();

  unsigned int tmpSize = integratedValues.size();
  m_integratedValues.resize(tmpSize,0.);
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_integratedValues[i] = integratedValues[i];
  }

  return;
}

void
uqDeltaSeq1D1DFunctionClass::printForMatlab(
  const uqBaseEnvironmentClass& env,
  std::ofstream&                ofsvar,
  const std::string&            prefixName) const
{
  uqSampled1D1DFunctionClass::printForMatlab(env,ofsvar,prefixName);

  unsigned int tmpSize = m_integratedValues.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofsvar << "\n" << prefixName << "intValue_sub" << env.subIdString() << " = zeros("  << tmpSize << ",1);";
  }
  else {
    ofsvar << "\n" << prefixName << "intValue_sub" << env.subIdString() << " = zeros("  << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofsvar << "\n" << prefixName << "intValue_sub" << env.subIdString() << "(" << i+1 << ",1) = " << m_integratedValues[i]  << ";";
    }
  }

  return;
}

//*****************************************************
// Gaussian 1d Kde 1D->1D class
//*****************************************************
uqGaussian1dKde1D1DFunctionClass::uqGaussian1dKde1D1DFunctionClass(
  const uqScalarSequenceClass<double>* chain,
  double chainMin,
  double chainMax,
  double gaussian1DScale)
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY), //chainMin-4.*gaussian1DScale,chainMax+4.*gaussian1DScale),
  m_chain              (chain),
  m_chainMin           (chainMin),
  m_chainMax           (chainMax),
  m_gaussian1DScale    (gaussian1DScale),
  m_lastInputK         (-1),
  m_lastQuadratureOrder(0),
  m_level_0            (0),
  m_level_m1           (0)
{
}

uqGaussian1dKde1D1DFunctionClass::~uqGaussian1dKde1D1DFunctionClass()
{
}

double                     
uqGaussian1dKde1D1DFunctionClass::value(double domainValue) const
{
  double value = 0.;

  double scaleInv = 1./m_gaussian1DScale;
  unsigned int dataSize = m_chain->subSequenceSize();
  for (unsigned int k = 0; k < dataSize; ++k) {
    double xk = (*m_chain)[k];
    value += uqMiscGaussianDensity((domainValue-xk)*scaleInv,0.,1.);
  }
  value = scaleInv * (value/(double) dataSize);

  return value;
}

double                     
uqGaussian1dKde1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqGaussian1dKde1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

double
uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const
{
  //std::cout << "Entering uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()" << std::endl;
  double value = 0.;
  double valueWithT = 0.;

  uqGaussianHermite1DQuadratureClass quadObj(0.,1.,quadratureOrder);
  const std::vector<double>& quadPositions = quadObj.positions();
  const std::vector<double>& quadWeights   = quadObj.weights  ();
  UQ_FATAL_TEST_MACRO(quadPositions.size() != quadWeights.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()",
                      "quadObj has invalid state");
  unsigned int numQuadraturePositions = quadPositions.size();

  unsigned int dataSize = m_chain->subSequenceSize();

  if (resultWithMultiplicationByTAsWell) {
    bool useFasterMethod = false;
    const uqQuadDenominator1D1DFunctionClass* functionPtr     = dynamic_cast<const uqQuadDenominator1D1DFunctionClass* >(&func);
    const uqQuadCRecursion1D1DFunctionClass*  constitutivePtr = NULL;
    if (functionPtr != NULL) {
      constitutivePtr = dynamic_cast<const uqQuadCRecursion1D1DFunctionClass* >(&(functionPtr->constitutiveFunction()));
    }

    int currentK = 0;
    if (constitutivePtr != NULL) {
      currentK = constitutivePtr->k();
      if (currentK == 0) {
        useFasterMethod       = true;
        m_lastInputK          = currentK;
        m_lastQuadratureOrder = quadratureOrder;
      }
      else if (currentK == (m_lastInputK+1)) {
        if (m_lastQuadratureOrder == quadratureOrder) {
          useFasterMethod       = true;
          m_lastInputK          = currentK;
        }
      }

      if (useFasterMethod == false) {
        std::cerr << "In uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()"
                  << ": weired situation"
                  << std::endl;
        m_lastInputK          = -1;
        m_lastQuadratureOrder = 0;
        m_level_0.clear();
        m_level_m1.clear();
      }
    }

    if (useFasterMethod) {
      if (currentK == 0) {
        m_level_0.resize (dataSize*numQuadraturePositions,0.); // Yes, '0..'
        m_level_m1.resize(dataSize*numQuadraturePositions,0.); // Yes, '0.'

        //std::cout << "In uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()"
        //          << ": currentK = " << currentK
        //          << std::endl;
        unsigned int auxK = 0;
        for (unsigned int k = 0; k < dataSize; ++k) {
          double xk = (*m_chain)[k];
          for (unsigned int j = 0; j < numQuadraturePositions; ++j) {
            double auxX = m_gaussian1DScale*quadPositions[j]+xk;
            m_level_0[auxK] = 1.; // Yes, '1.'
            double auxValue = m_level_0[auxK]*quadWeights[j];
            value      +=      auxValue;
            valueWithT += auxX*auxValue;
            auxK++;
          }
        }
      }
      else {
        double alpha = constitutivePtr->alpha()[currentK-1];
        double beta  = constitutivePtr->beta ()[currentK-1];

        //std::cout << "In uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()"
        //          << ": currentK = " << currentK
        //          << ", alpha = "    << alpha
        //          << ", beta = "     << beta
        //          << std::endl;
        unsigned int auxK = 0;
        for (unsigned int k = 0; k < dataSize; ++k) {
          double xk = (*m_chain)[k];
          for (unsigned int j = 0; j < numQuadraturePositions; ++j) {
            double auxX = m_gaussian1DScale*quadPositions[j]+xk;

            UQ_FATAL_TEST_MACRO((currentK == 1) && (m_level_0[auxK] != 1.),
                                UQ_UNAVAILABLE_RANK,
                                "uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()",
                                "invalid currentK=1a");

            UQ_FATAL_TEST_MACRO((currentK == 1) && (m_level_m1[auxK] != 0.),
                                UQ_UNAVAILABLE_RANK,
                                "uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()",
                                "invalid currentK=1b");

            UQ_FATAL_TEST_MACRO((currentK == 2) && (m_level_m1[auxK] != 1.),
                                UQ_UNAVAILABLE_RANK,
                                "uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate()",
                                "invalid currentK=2");

            //double tmpValue = (auxX - alpha)*(auxX - alpha)*m_level_0[auxK]
            //                + beta*beta*m_level_m1[auxK]
            //                - 2.*(auxX - alpha)*sqrt(m_level_0[auxK])*beta*sqrt(m_level_m1[auxK]);
            double tmpValue = (auxX - alpha)*m_level_0[auxK] - beta*m_level_m1[auxK];
            m_level_m1[auxK] = m_level_0[auxK];
            m_level_0 [auxK] = tmpValue;
            double auxValue = tmpValue*tmpValue*quadWeights[j];
            value      +=      auxValue;
            valueWithT += auxX*auxValue;
            auxK++;
          }
        }
      }
    }
    else {
      for (unsigned int k = 0; k < dataSize; ++k) {
        double xk = (*m_chain)[k];
        for (unsigned int j = 0; j < numQuadraturePositions; ++j) {
          double auxX = m_gaussian1DScale*quadPositions[j]+xk;
          double auxValue = func.value(auxX)*quadWeights[j];
          value      +=      auxValue;
          valueWithT += auxX*auxValue;
        }
      }
    }

    value      *= (1./sqrt(2.*M_PI)/((double) dataSize));
    valueWithT *= (1./sqrt(2.*M_PI)/((double) dataSize));
    *resultWithMultiplicationByTAsWell = valueWithT;
  }
  else {
    for (unsigned int k = 0; k < dataSize; ++k) {
      double xk = (*m_chain)[k];
      for (unsigned int j = 0; j < numQuadraturePositions; ++j) {
        value += func.value(m_gaussian1DScale*quadPositions[j]+xk)*quadWeights[j];
      }
    }
    value *= (1./sqrt(2.*M_PI)/((double) dataSize));
  }

  //std::cout << "Leaving uqGaussian1dKde1D1DFunctionClass::multiplyAndIntegrate(), value = " << value << std::endl;

  return value;
}

//*****************************************************
// 'ScalarTimesFunc' 1D->1D class
//*****************************************************
uqScalarTimesFunc1D1DFunctionClass::uqScalarTimesFunc1D1DFunctionClass(
  double scalar,
  const uqBase1D1DFunctionClass& func)
  :
  uqBase1D1DFunctionClass(func.minDomainValue(),func.maxDomainValue()),
  m_scalar(scalar),
  m_func  (func)
{
}

uqScalarTimesFunc1D1DFunctionClass::~uqScalarTimesFunc1D1DFunctionClass()
{
}

double
uqScalarTimesFunc1D1DFunctionClass::value(double domainValue) const
{
  double value = 0.;

  value += m_scalar*m_func.value(domainValue);

  return value;
}

double
uqScalarTimesFunc1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarTimesFunc1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

//*****************************************************
// 'FuncTimesFunc' 1D->1D class
//*****************************************************
uqFuncTimesFunc1D1DFunctionClass::uqFuncTimesFunc1D1DFunctionClass(
  const uqBase1D1DFunctionClass& func1,
  const uqBase1D1DFunctionClass& func2)
  :
  uqBase1D1DFunctionClass(std::max(func1.minDomainValue(),func2.minDomainValue()),std::min(func1.maxDomainValue(),func2.maxDomainValue())),
  m_func1(func1),
  m_func2(func2)
{
}

uqFuncTimesFunc1D1DFunctionClass::~uqFuncTimesFunc1D1DFunctionClass()
{
}

double
uqFuncTimesFunc1D1DFunctionClass::value(double domainValue) const
{
  double value = 0.;

  value += m_func1.value(domainValue)*m_func2.value(domainValue);

  return value;
}

double
uqFuncTimesFunc1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqFuncTimesFunc1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

//*****************************************************
// 'FuncPlusFunc' 1D->1D class
//*****************************************************
uqFuncPlusFunc1D1DFunctionClass::uqFuncPlusFunc1D1DFunctionClass(
  const uqBase1D1DFunctionClass& func1,
  const uqBase1D1DFunctionClass& func2)
  :
  uqBase1D1DFunctionClass(std::max(func1.minDomainValue(),func2.minDomainValue()),std::min(func1.maxDomainValue(),func2.maxDomainValue())),
  m_func1(func1),
  m_func2(func2)
{
}

uqFuncPlusFunc1D1DFunctionClass::~uqFuncPlusFunc1D1DFunctionClass()
{
}

double
uqFuncPlusFunc1D1DFunctionClass::value(double domainValue) const
{
  double value = 0.;

  value += m_func1.value(domainValue) + m_func2.value(domainValue);

  return value;
}

double
uqFuncPlusFunc1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqFuncPlusFunc1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

//*****************************************************
// 'QuadDenominator' 1D->1D class
//*****************************************************
uqQuadDenominator1D1DFunctionClass::uqQuadDenominator1D1DFunctionClass(
  const uqBase1D1DFunctionClass& func)
  :
  uqBase1D1DFunctionClass(func.minDomainValue(),func.maxDomainValue()),
  m_func(func)
{
}

uqQuadDenominator1D1DFunctionClass::~uqQuadDenominator1D1DFunctionClass()
{
}

double
uqQuadDenominator1D1DFunctionClass::value(double domainValue) const
{
  double value = m_func.value(domainValue);
  return value*value;
}

double
uqQuadDenominator1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadDenominator1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

const uqBase1D1DFunctionClass&
uqQuadDenominator1D1DFunctionClass::constitutiveFunction() const
{
  return m_func;
}

//*****************************************************
// 'QuadNumerator' 1D->1D class
//*****************************************************
uqQuadNumerator1D1DFunctionClass::uqQuadNumerator1D1DFunctionClass(
  const uqBase1D1DFunctionClass& func)
  :
  uqBase1D1DFunctionClass(func.minDomainValue(),func.maxDomainValue()),
  m_func(func)
{
}

uqQuadNumerator1D1DFunctionClass::~uqQuadNumerator1D1DFunctionClass()
{
}

double
uqQuadNumerator1D1DFunctionClass::value(double domainValue) const
{
  double value = m_func.value(domainValue);
  return domainValue*value*value;
}

double
uqQuadNumerator1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadNumerator1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
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
  std::vector<double>&                     beta)
{
  UQ_FATAL_TEST_MACRO(k < 0,
                      UQ_UNAVAILABLE_RANK,
                      "alphaBetaCLoop()",
                      "k is negative");

  unsigned int n = alpha.size();

  UQ_FATAL_TEST_MACRO((n < 1),
                      UQ_UNAVAILABLE_RANK,
                      "alphaBetaCLoop()",
                      "invalid n");

  uqQuadDenominator1D1DFunctionClass pi_pi_0_func  (pi_0);
//uqQuadNumerator1D1DFunctionClass   t_pi_pi_0_func(pi_0);

  double t_pi_pi_0 = 0.;
  double pi_pi_0   = rho.multiplyAndIntegrate(pi_pi_0_func,integrationOrder,&t_pi_pi_0);//,1+((unsigned int) k) );
//double t_pi_pi_0 = rho.multiplyAndIntegrate(t_pi_pi_0_func,integrationOrder,NULL);    //,1+((unsigned int) k) );

  // Alpha and beta
  alpha[k] = t_pi_pi_0/pi_pi_0;
  beta[k]  = pi_pi_0/pi_pi_m1;
  UQ_FATAL_TEST_MACRO((beta[k] < 0.),
                      UQ_UNAVAILABLE_RANK,
                      "alphaBetaCLoop()",
                      "beta is negative");

  std::cout << "In alphaBetaCLoop()"
            << ": k = "         << k
    //<< ", pi_0(0.) = "  << pi_0.value(0.)
    //<< ", pi_pi_m1 = "  << pi_pi_m1
    //<< ", pi_pi_0 = "   << pi_pi_0
    //<< ", t_pi_pi_0 = " << t_pi_pi_0
    //<< ", alpha[k] = "  << alpha[k]
    //<< ", beta[k] = "   << beta[k]
            << std::endl;

  if (k < (int) (n-1)) {
    uqQuadCRecursion1D1DFunctionClass pi_p1(k+1,
                                            alpha,
                                            beta);

    alphaBetaCLoop(k+1,
                   pi_pi_0,
                   pi_0,
                   pi_p1,
                   rho,
                   integrationOrder,
                   alpha,
                   beta);
  }
  //std::cout << "Leaving alphaBetaCLoop(), k = " << k << ", beta[0] = " << beta[0] << std::endl;

  return;
}

//*****************************************************
// 'QuadCRecursion' 1D->1D class
//*****************************************************
uqQuadCRecursion1D1DFunctionClass::uqQuadCRecursion1D1DFunctionClass(int k)
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY),
  m_k    (0), // Yes, '0' instead of 'k'
  m_alpha(0),
  m_beta (0)
{
  UQ_FATAL_TEST_MACRO(k < -1, // Yes, '-1'
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadCRecursion1D1DFunctionClass::constructor(1)",
                      "too small k");

  UQ_FATAL_TEST_MACRO(k > 0, // Yes, '0'
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadCRecursion1D1DFunctionClass::constructor(1)",
                      "too big k");
  return;
}

uqQuadCRecursion1D1DFunctionClass::uqQuadCRecursion1D1DFunctionClass(
  int                        k,
  const std::vector<double>& alpha,
  const std::vector<double>& beta)
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY),
  m_k    (k),
  m_alpha(k,0.),
  m_beta (k,0.)
{
  UQ_FATAL_TEST_MACRO(m_k < 0, // Yes, '0'
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadCRecursion1D1DFunctionClass::constructor(2)",
                      "invalid k");

  UQ_FATAL_TEST_MACRO(m_k > (int) alpha.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadCRecursion1D1DFunctionClass::constructor(2)",
                      "invalid k and alpha");

  UQ_FATAL_TEST_MACRO(alpha.size() != beta.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadCRecursion1D1DFunctionClass::constructor(2)",
                      "invalid alpha and beta");

  for (unsigned int i = 0; i < (unsigned int) m_k; ++i) {
    m_alpha[i] = alpha[i];
    m_beta [i] = beta [i];
  }

  //std::cout << "In uqQuadCRecursion1D1DFunctionClass constructor: m_k = " << m_k << ", m_alpha = ";
  //for (int i = 0; i < m_k; ++i) {
  //  std::cout << " " << m_alpha[i];
  //}
  //std::cout << std::endl;

  //std::cout << "In uqQuadCRecursion1D1DFunctionClass constructor: m_k = " << m_k << ", m_beta = ";
  // for (int i = 0; i < m_k; ++i) {
  //  std::cout << " " << m_beta[i];
  //}
  //std::cout << std::endl;
}

uqQuadCRecursion1D1DFunctionClass::~uqQuadCRecursion1D1DFunctionClass()
{
}

double
uqQuadCRecursion1D1DFunctionClass::value(double domainValue) const
{
  double value = 0.;
  if (m_k < -1) {
    UQ_FATAL_TEST_MACRO(true,
                        UQ_UNAVAILABLE_RANK,
                        "uqQuadCRecursion1D1DFunctionClass::value()",
                        "invalid m_k");
  }
  else if (m_k == -1) {
    value = 0.;
  }
  else if (m_k == 0) {
    value = 1.;
  }
  else {
    double val_m1 = 0.;
    double val_0  = 1.;
    for (unsigned int i = 0; i < (unsigned int) m_k; ++i) {
      value = (domainValue - m_alpha[i])*val_0 - m_beta[i]*val_m1;
      val_m1 = val_0;
      val_0 = value;
    }
  }

  return value;
}

double
uqQuadCRecursion1D1DFunctionClass::deriv(double domainValue) const
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadCRecursion1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return 0.;
}

int
uqQuadCRecursion1D1DFunctionClass::k() const
{
  return m_k;
}

const std::vector<double>&
uqQuadCRecursion1D1DFunctionClass::alpha() const
{
  return m_alpha;
}

const std::vector<double>&
uqQuadCRecursion1D1DFunctionClass::beta() const
{
  return m_beta;
}
