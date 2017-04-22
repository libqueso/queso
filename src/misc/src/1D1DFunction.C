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

#include <queso/1D1DFunction.h>
#include <queso/1DQuadrature.h>

namespace QUESO {

//*****************************************************
// Base 1D->1D class
//*****************************************************
Base1D1DFunction::Base1D1DFunction(
  double minDomainValue,
  double maxDomainValue)
  :
  m_minDomainValue(minDomainValue),
  m_maxDomainValue(maxDomainValue)
{
  queso_require_less_msg(m_minDomainValue, m_maxDomainValue, "min >= max");
}

Base1D1DFunction::~Base1D1DFunction()
{
}

double
Base1D1DFunction::minDomainValue() const
{
  return m_minDomainValue;
}

double
Base1D1DFunction::maxDomainValue() const
{
  return m_maxDomainValue;
}

double
Base1D1DFunction::multiplyAndIntegrate(const Base1D1DFunction& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const
{
  double value = 0.;

  queso_not_implemented();

  if (resultWithMultiplicationByTAsWell) { // Just to eliminate INTEL compiler warnings
    func.value((double) quadratureOrder);
  }

  return value;
}

//*****************************************************
// Generic 1D->1D class
//*****************************************************
Generic1D1DFunction::Generic1D1DFunction(
  double minDomainValue,
  double maxDomainValue,
  double (*valueRoutinePtr)(double domainValue, const void* routinesDataPtr),
  double (*derivRoutinePtr)(double domainValue, const void* routinesDataPtr),
  const void* routinesDataPtr)
  :
  Base1D1DFunction(minDomainValue,maxDomainValue),
  m_valueRoutinePtr      (valueRoutinePtr),
  m_derivRoutinePtr      (derivRoutinePtr),
  m_routinesDataPtr      (routinesDataPtr)
{
}

Generic1D1DFunction::~Generic1D1DFunction()
{
}

double
Generic1D1DFunction::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Generic1D1DFunction::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  return (*m_valueRoutinePtr)(domainValue,m_routinesDataPtr);
}

double
Generic1D1DFunction::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Generic1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  return (*m_derivRoutinePtr)(domainValue,m_routinesDataPtr);
}

//*****************************************************
// Constant 1D->1D class
//*****************************************************
Constant1D1DFunction::Constant1D1DFunction(
  double minDomainValue,
  double maxDomainValue,
  double constantValue)
  :
  Base1D1DFunction(minDomainValue,maxDomainValue),
  m_constantValue        (constantValue)
{
}

Constant1D1DFunction::~Constant1D1DFunction()
{
}

double
Constant1D1DFunction::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Constant1D1DFunction::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  return m_constantValue;
}

double
Constant1D1DFunction::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Constant1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  return 0.;
}

//*****************************************************
// Linear 1D->1D class
//*****************************************************
Linear1D1DFunction::Linear1D1DFunction(
  double minDomainValue,
  double maxDomainValue,
  double referenceDomainValue,
  double referenceImageValue,
  double rateValue)
  :
  Base1D1DFunction(minDomainValue,maxDomainValue),
  m_referenceDomainValue (referenceDomainValue),
  m_referenceImageValue  (referenceImageValue),
  m_rateValue            (rateValue)
{
}

Linear1D1DFunction::~Linear1D1DFunction()
{
}

double
Linear1D1DFunction::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Linear1D1DFunction::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  double imageValue = m_referenceImageValue + m_rateValue*(domainValue - m_referenceDomainValue);

  return imageValue;
}

double
Linear1D1DFunction::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Linear1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  return m_rateValue;
}

//*****************************************************
// PiecewiseLinear 1D->1D class
//*****************************************************
PiecewiseLinear1D1DFunction::PiecewiseLinear1D1DFunction(
  double                     minDomainValue,
  double                     maxDomainValue,
  const std::vector<double>& referenceDomainValues,
  double                     referenceImageValue0,
  const std::vector<double>& rateValues)
  :
  Base1D1DFunction(minDomainValue,maxDomainValue),
  m_numRefValues         (referenceDomainValues.size()),
  m_referenceDomainValues(referenceDomainValues),
  m_rateValues           (rateValues)
{
  queso_require_not_equal_to_msg(m_numRefValues, 0, "num ref values = 0");

  queso_require_equal_to_msg(m_numRefValues, rateValues.size(), "num rate values is inconsistent");

  for (unsigned int i = 1; i < m_numRefValues; ++i) { // Yes, from '1'
    queso_require_greater_msg(m_referenceDomainValues[i], m_referenceDomainValues[i-1], "reference domain values are inconsistent");
  }

  m_referenceImageValues.clear();
  m_referenceImageValues.resize(m_numRefValues,0.);
  m_referenceImageValues[0] = referenceImageValue0;
  for (unsigned int i = 1; i < m_numRefValues; ++i) { // Yes, from '1'
    m_referenceImageValues[i] = m_referenceImageValues[i-1] + m_rateValues[i-1]*(m_referenceDomainValues[i] - m_referenceDomainValues[i-1]);
  }

  if (false) { // For debug only
    std::cout << "In PiecewiseLinear1D1DFunction::constructor():"
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

PiecewiseLinear1D1DFunction::~PiecewiseLinear1D1DFunction()
{
  m_rateValues.clear();
  m_referenceImageValues.clear();
  m_referenceDomainValues.clear();
}

double
PiecewiseLinear1D1DFunction::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In PiecewiseLinear1D1DFunction::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

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
    std::cout << "In PiecewiseLinear1D1DFunction::value()"
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
PiecewiseLinear1D1DFunction::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In PiecewiseLinear1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

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
        queso_require_less_equal_msg(i, m_numRefValues, "too big 'i'");
      }
    }
  }

  return m_rateValues[i];
}

//*****************************************************
// Quadratic 1D->1D class
//*****************************************************
Quadratic1D1DFunction::Quadratic1D1DFunction(
  double minDomainValue,
  double maxDomainValue,
  double a,
  double b,
  double c)
  :
  Base1D1DFunction(minDomainValue,maxDomainValue),
  m_a                    (a),
  m_b                    (b),
  m_c                    (c)
{
}

Quadratic1D1DFunction::~Quadratic1D1DFunction()
{
}

double
Quadratic1D1DFunction::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Quadratic1D1DFunction::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  double imageValue = m_a*domainValue*domainValue + m_b*domainValue + m_c;

  return imageValue;
}

double
Quadratic1D1DFunction::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Quadratic1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  return 2.*m_a*domainValue + m_b;
}

//*****************************************************
// Sampled 1D->1D class
//*****************************************************
Sampled1D1DFunction::Sampled1D1DFunction()
  :
  Base1D1DFunction(-INFINITY,INFINITY)
{

  //                    UQ_UNAVAILABLE_RANK,
  //                    "SampledD1DFunction::deriv()",
  //                    "invalid constructor");
}

Sampled1D1DFunction::Sampled1D1DFunction(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues)
  :
  Base1D1DFunction(domainValues[0],domainValues[domainValues.size()-1]),
  m_domainValues    (domainValues.size(),0.),
  m_imageValues     (imageValues.size(), 0.)
{
  unsigned int tmpSize = m_domainValues.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_domainValues[i] = domainValues[i];
    m_imageValues [i] = imageValues [i];
  }
}

Sampled1D1DFunction::~Sampled1D1DFunction()
{
}

double
Sampled1D1DFunction::value(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Sampled1D1DFunction::value()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  double returnValue = 0.;

  unsigned int tmpSize = m_domainValues.size();
  //std::cout << "In Sampled1D1DFunction::value()"
  //          << ": domainValue = "         << domainValue
  //          << ", tmpSize = "             << tmpSize
  //          << ", m_domainValues[0] = "   << m_domainValues[0]
  //          << ", m_domainValues[max] = " << m_domainValues[tmpSize-1]
  //          << std::endl;

  queso_require_not_equal_to_msg(tmpSize, 0, "m_domainValues.size() = 0");

  queso_require_greater_equal_msg(domainValue, m_domainValues[0], "domainValue < m_domainValues[0]");

  queso_require_greater_equal_msg(m_domainValues[tmpSize-1], domainValue, "m_domainValues[max] < domainValue");

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
Sampled1D1DFunction::deriv(double domainValue) const
{
  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In Sampled1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  queso_error_msg("this function makes no sense for this class");
  return 0.;
}

const std::vector<double>&
Sampled1D1DFunction::domainValues() const
{
  return m_domainValues;
}

const std::vector<double>&
Sampled1D1DFunction::imageValues() const
{
  return m_imageValues;
}

bool
Sampled1D1DFunction::domainValueMatchesExactly(double domainValue) const
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
Sampled1D1DFunction::set(
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
Sampled1D1DFunction::printForMatlab(
  const BaseEnvironment& env,
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
ScalarTimesFunc1D1DFunction::ScalarTimesFunc1D1DFunction(
  double scalar,
  const Base1D1DFunction& func)
  :
  Base1D1DFunction(func.minDomainValue(),func.maxDomainValue()),
  m_scalar(scalar),
  m_func  (func)
{
}

ScalarTimesFunc1D1DFunction::~ScalarTimesFunc1D1DFunction()
{
}

double
ScalarTimesFunc1D1DFunction::value(double domainValue) const
{
  double value = 0.;

  value += m_scalar*m_func.value(domainValue);

  return value;
}

double
ScalarTimesFunc1D1DFunction::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In ScalarTimes1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  queso_not_implemented();

  return value;
}

//*****************************************************
// 'FuncTimesFunc' 1D->1D class
//*****************************************************
FuncTimesFunc1D1DFunction::FuncTimesFunc1D1DFunction(
  const Base1D1DFunction& func1,
  const Base1D1DFunction& func2)
  :
  Base1D1DFunction(std::max(func1.minDomainValue(),func2.minDomainValue()),std::min(func1.maxDomainValue(),func2.maxDomainValue())),
  m_func1(func1),
  m_func2(func2)
{
}

FuncTimesFunc1D1DFunction::~FuncTimesFunc1D1DFunction()
{
}

double
FuncTimesFunc1D1DFunction::value(double domainValue) const
{
  double value = 0.;

  value += m_func1.value(domainValue)*m_func2.value(domainValue);

  return value;
}

double
FuncTimesFunc1D1DFunction::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In FuncTimes1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  queso_not_implemented();

  return value;
}

//*****************************************************
// 'FuncPlusFunc' 1D->1D class
//*****************************************************
FuncPlusFunc1D1DFunction::FuncPlusFunc1D1DFunction(
  const Base1D1DFunction& func1,
  const Base1D1DFunction& func2)
  :
  Base1D1DFunction(std::max(func1.minDomainValue(),func2.minDomainValue()),std::min(func1.maxDomainValue(),func2.maxDomainValue())),
  m_func1(func1),
  m_func2(func2)
{
}

FuncPlusFunc1D1DFunction::~FuncPlusFunc1D1DFunction()
{
}

double
FuncPlusFunc1D1DFunction::value(double domainValue) const
{
  double value = 0.;

  value += m_func1.value(domainValue) + m_func2.value(domainValue);

  return value;
}

double
FuncPlusFunc1D1DFunction::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In FuncPlus1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  queso_not_implemented();

  return value;
}

//*****************************************************
// Lagrange Polynomial 1D->1D class
//*****************************************************
LagrangePolynomial1D1DFunction::LagrangePolynomial1D1DFunction(
  const std::vector<double>& positionValues,
  const std::vector<double>* functionValues)
  :
  Base1D1DFunction(-INFINITY,INFINITY),
  m_positionValues(positionValues),
  m_functionValues(positionValues.size(),1.)
{
  if (functionValues) {
    queso_require_equal_to_msg(m_positionValues.size(), functionValues->size(), "invalid input");
    m_functionValues = *functionValues;
  }
}

LagrangePolynomial1D1DFunction::~LagrangePolynomial1D1DFunction()
{
}

double
LagrangePolynomial1D1DFunction::value(double domainValue) const
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
    //  *m_env.subDisplayFile() << "In sddTG<K_V,K_M>::lagrange()"
    //                          << ": k = " << k
    //                          << ", scaleFactor = " << scaleFactor
    //                          << std::endl;
    //}

    value += scaleFactor * m_functionValues[k];

    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In sddTG<K_V,K_M>::lagrange()"
    //                          << ": k = " << k
    //                          << ", value = " << value
    //                          << std::endl;
    //}
  }

  return value;
}

double
LagrangePolynomial1D1DFunction::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In LagrangePolynomial1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  queso_not_implemented();

  return value;
}

//*****************************************************
// Lagrange Basis 1D->1D class
//*****************************************************
LagrangeBasis1D1DFunction::LagrangeBasis1D1DFunction(
  const std::vector<double>& positionValues,
  unsigned int               basisIndex)
  :
  Base1D1DFunction(-INFINITY,INFINITY),
  m_positionValues(positionValues),
  m_basisIndex    (basisIndex)
{
  queso_require_less_msg(m_basisIndex, m_positionValues.size(), "invalid input");
}

LagrangeBasis1D1DFunction::~LagrangeBasis1D1DFunction()
{
}

double
LagrangeBasis1D1DFunction::value(double domainValue) const
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
LagrangeBasis1D1DFunction::deriv(double domainValue) const
{
  double value = 0.;

  if ((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)) {
    std::cerr << "In LagrangeBasis1D1DFunction::deriv()"
              << ": requested x ("            << domainValue
              << ") is out of the interval (" << m_minDomainValue
              << ", "                         << m_maxDomainValue
              << ")"
              << std::endl;
  }

  queso_require_msg(!((domainValue < m_minDomainValue) || (domainValue > m_maxDomainValue)), "x out of range");

  queso_not_implemented();

  return value;
}

template <class T>
double
SubF1F2Gaussian2dKdeIntegral(const ScalarSequence<T>& scalarSeq1,
                               const ScalarSequence<T>& scalarSeq2,
                               unsigned int                    initialPos,
                               double                          scaleValue1,
                               double                          scaleValue2,
                               const Base1D1DFunction&  func1,
                               const Base1D1DFunction&  func2,
                               unsigned int                    quadratureOrder)
{
  double resultValue = 0.;

  queso_require_equal_to_msg(initialPos, 0, "not implemented yet for initialPos != 0");
  queso_require_equal_to_msg(scalarSeq1.subSequenceSize(), scalarSeq2.subSequenceSize(), "different sizes");

  GaussianHermite1DQuadrature quadObj(0.,1.,quadratureOrder);
  const std::vector<double>& quadPositions = quadObj.positions();
  const std::vector<double>& quadWeights   = quadObj.weights  ();
  queso_require_equal_to_msg(quadPositions.size(), quadWeights.size(), "quadObj has invalid state");

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

}  // End namespace QUESO
