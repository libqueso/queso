//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// $Id:$
//
//--------------------------------------------------------------------------

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

  if (resultWithMultiplicationByTAsWell) { // Just to eliminate INTEL compiler warnings
    func.value((double) quadratureOrder);
  }

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
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqGeneric1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqGeneric1D1DFunctionClass::deriv()",
                      "x out of range");

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
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqLinear1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqLinear1D1DFunctionClass::deriv()",
                      "x out of range");

  return m_rateValue;
}

//*****************************************************
// PiecewiseLinear 1D->1D class
//*****************************************************
uqPiecewiseLinear1D1DFunctionClass::uqPiecewiseLinear1D1DFunctionClass(
  double                     minDomainValue,
  double                     maxDomainValue,
  const std::vector<double>& referenceDomainValues,
  double                     referenceImageValue0,
  const std::vector<double>& rateValues)
  :
  uqBase1D1DFunctionClass(minDomainValue,maxDomainValue),
  m_numRefValues         (referenceDomainValues.size()),
  m_referenceDomainValues(referenceDomainValues),
  m_rateValues           (rateValues)
{
  UQ_FATAL_TEST_MACRO(m_numRefValues == 0,
                      UQ_UNAVAILABLE_RANK,
                      "uqPiecewiseLinear1D1DFunctionClass::constructor()",
                      "num ref values = 0");

  UQ_FATAL_TEST_MACRO(m_numRefValues != rateValues.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqPiecewiseLinear1D1DFunctionClass::constructor()",
                      "num rate values is inconsistent");

  for (unsigned int i = 1; i < m_numRefValues; ++i) { // Yes, from '1'
    UQ_FATAL_TEST_MACRO(m_referenceDomainValues[i] <= m_referenceDomainValues[i-1],
                        UQ_UNAVAILABLE_RANK,
                        "uqPiecewiseLinear1D1DFunctionClass::constructor()",
                        "reference domain values are inconsistent");
  }

  m_referenceImageValues.clear();
  m_referenceImageValues.resize(m_numRefValues,0.);
  m_referenceImageValues[0] = referenceImageValue0;
  for (unsigned int i = 1; i < m_numRefValues; ++i) { // Yes, from '1'
    m_referenceImageValues[i] = m_referenceImageValues[i-1] + m_rateValues[i-1]*(m_referenceDomainValues[i] - m_referenceDomainValues[i-1]);
  }

  if (false) { // For debug only
    std::cout << "In uqPiecewiseLinear1D1DFunctionClass::constructor():"
              << std::endl;
    for (unsigned int i = 0; i < m_numRefValues; ++i) {
      std::cout << "i = " << i
                << ", m_referenceDomainValues[i] = " << m_referenceDomainValues[i]
                << ", m_referenceImageValues[i] = "  << m_referenceImageValues[i]
                << ", m_rateValues[i] = "            << m_rateValues[i]
                << std::endl;
    }
  }
}

uqPiecewiseLinear1D1DFunctionClass::~uqPiecewiseLinear1D1DFunctionClass()
{
  m_rateValues.clear();
  m_referenceImageValues.clear();
  m_referenceDomainValues.clear();
}

double
uqPiecewiseLinear1D1DFunctionClass::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqPiecewiseLinear1D1DFunctionClass::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqPiecewiseLinear1D1DFunctionClass::value()",
                      "x out of range");

  unsigned int i = 0;
  if (m_numRefValues == 1) {
    // Nothing else to do
  }
  else {
    bool referenceDomainValueFound = false;
    while (referenceDomainValueFound == false) {
      if (domainValue < m_referenceDomainValues[i+1]) {
        referenceDomainValueFound = true;
      }
      else {
        ++i;
        if (i == (m_numRefValues-1)) {
          referenceDomainValueFound = true;
        }
      }
    }
  }
  double imageValue = m_referenceImageValues[i] + m_rateValues[i]*(domainValue - m_referenceDomainValues[i]);
  if (false) { // For debug only
    std::cout << "In uqPiecewiseLinear1D1DFunctionClass::value()"
              << ": domainValue = "                << domainValue
              << ", i = "                          << i
              << ", m_referenceDomainValues[i] = " << m_referenceDomainValues[i]
              << ", m_referenceImageValues[i] = "  << m_referenceImageValues[i]
              << ", m_rateValues[i] = "            << m_rateValues[i]
              << ", imageValue = "                 << imageValue
              << std::endl;
  }

  return imageValue;
}

double
uqPiecewiseLinear1D1DFunctionClass::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqPiecewiseLinear1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqPiecewiseLinear1D1DFunctionClass::deriv()",
                      "x out of range");

  unsigned int i = 0;
  if (m_numRefValues == 1) {
    // Nothing else to do
  }
  else {
    bool referenceDomainValueFound = false;
    while (referenceDomainValueFound == false) {
      if (domainValue < m_referenceDomainValues[i+1]) {
        referenceDomainValueFound = true;
      }
      else {
        ++i;
        UQ_FATAL_TEST_MACRO(i > m_numRefValues,
                            UQ_UNAVAILABLE_RANK,
                           "uqPiecewiseLinear1D1DFunctionClass::deriv()",
                           "too big 'i'");
      }
    }
  }

  return m_rateValues[i];
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
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqQuadratic1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqQuadratic1D1DFunctionClass::deriv()",
                      "x out of range");

  return 2.*m_a*domainValue + m_b;
}

//*****************************************************
// Sampled 1D->1D class
//*****************************************************
uqSampled1D1DFunctionClass::uqSampled1D1DFunctionClass()
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY)
{
  //UQ_FATAL_TEST_MACRO(true,
  //                    UQ_UNAVAILABLE_RANK,
  //                    "uqSampledD1DFunctionClass::deriv()",
  //                    "invalid constructor");
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
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqSampled1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::deriv()",
                      "x out of range");

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1d1DFunctionClass::deriv()",
                      "this function makes no sense for this class");
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

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqScalarTimes1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqScalarTimes1D1DFunctionClass::deriv()",
                      "x out of range");

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

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqFuncTimes1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqFuncTimes1D1DFunctionClass::deriv()",
                      "x out of range");

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

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqFuncPlus1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqFuncPlus1D1DFunctionClass::deriv()",
                      "x out of range");

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqFuncPlusFunc1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

//*****************************************************
// Lagrange Polynomial 1D->1D class
//*****************************************************
uqLagrangePolynomial1D1DFunctionClass::uqLagrangePolynomial1D1DFunctionClass(
  const std::vector<double>& positionValues,
  const std::vector<double>* functionValues)
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY),
  m_positionValues(positionValues),
  m_functionValues(positionValues.size(),1.)
{
  if (functionValues) {
    UQ_FATAL_TEST_MACRO(m_positionValues.size() != functionValues->size(),
                        UQ_UNAVAILABLE_RANK,
                        "uqLagrangePolynomial1D1DFunctionClass::constructor()",
                        "invalid input");
    m_functionValues = *functionValues;
  }
}

uqLagrangePolynomial1D1DFunctionClass::~uqLagrangePolynomial1D1DFunctionClass()
{
}

double                     
uqLagrangePolynomial1D1DFunctionClass::value(double domainValue) const
{
  double value = 0.;

  for (unsigned int k = 0; k < m_positionValues.size(); ++k) {
    double scaleFactor = 1.;
    double posK = m_positionValues[k];
    for (unsigned int j = 0; j < m_positionValues.size(); ++j) {
      if (j != k) {
        double posJ = m_positionValues[j];
        scaleFactor *= (domainValue-posJ)/(posK-posJ);
      }
    }

    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In sddTGClass<K_V,K_M>::lagrange()"
    //                          << ": k = " << k
    //                          << ", scaleFactor = " << scaleFactor
    //                          << std::endl;
    //}

    value += scaleFactor * m_functionValues[k];

    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In sddTGClass<K_V,K_M>::lagrange()"
    //                          << ": k = " << k
    //                          << ", value = " << value
    //                          << std::endl;
    //}
  }

  return value;
}

double                     
uqLagrangePolynomial1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqLagrangePolynomial1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqLagrangePolynomial1D1DFunctionClass::deriv()",
                      "x out of range");

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqLagrangePolynomial1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

//*****************************************************
// Lagrange Basis 1D->1D class
//*****************************************************
uqLagrangeBasis1D1DFunctionClass::uqLagrangeBasis1D1DFunctionClass(
  const std::vector<double>& positionValues,
  unsigned int               basisIndex)
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY),
  m_positionValues(positionValues),
  m_basisIndex    (basisIndex)
{
  UQ_FATAL_TEST_MACRO(m_basisIndex >= m_positionValues.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqLagrangeBasis1D1DFunctionClass::constructor()",
                      "invalid input");
}

uqLagrangeBasis1D1DFunctionClass::~uqLagrangeBasis1D1DFunctionClass()
{
}

double                     
uqLagrangeBasis1D1DFunctionClass::value(double domainValue) const
{
  double scaleFactor = 1.;

  unsigned int k = m_basisIndex;
  double posK = m_positionValues[k];
  for (unsigned int j = 0; j < m_positionValues.size(); ++j) {
    if (j != k) {
      double posJ = m_positionValues[j];
      scaleFactor *= (domainValue-posJ)/(posK-posJ);
    }
  }

  return scaleFactor;
}

double                     
uqLagrangeBasis1D1DFunctionClass::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In uqLagrangeBasis1D1DFunctionClass::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)),
                      UQ_UNAVAILABLE_RANK,
                      "uqLagrangeBasis1D1DFunctionClass::deriv()",
                      "x out of range");

  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqLagrangeBasis1D1DFunctionClass::deriv()",
                      "not implemented yet");

  return value;
}

