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

protected:
  double m_minDomainValue;
  double m_maxDomainValue;
};

uqBase1D1DFunctionClass::uqBase1D1DFunctionClass(
  double minDomainValue,
  double maxDomainValue)
  :
  m_minDomainValue(minDomainValue),
  m_maxDomainValue(maxDomainValue)
{
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
  // FIX ME: check range of domainValue
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
              << ": requested x ("             << domainValue
              << ") is out of the intervarl (" << m_minDomainValue
              << ", "                          << m_maxDomainValue
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
              << ": requested x ("             << domainValue
              << ") is out of the intervarl (" << m_minDomainValue
              << ", "                          << m_maxDomainValue
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
  // FIX ME: check range of domainValue
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
  // FIX ME: check range of domainValue
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

  virtual void                 printForMatlab(std::ofstream& ofs, const std::string& prefixName) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_domainValues;
  std::vector<double> m_imageValues;
};

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
  // FIX ME: check range of domainValue
  double returnValue = 0.;

  unsigned int tmpSize = m_domainValues.size();
#if 0
  std::cout << "In uqSampled1D1DFunctionClass::value()"
            << ": domainValue = "         << domainValue
            << ", tmpSize = "             << tmpSize
            << ", m_domainValues[0] = "   << m_domainValues[0]
            << ", m_domainValues[max] = " << m_domainValues[tmpSize-1]
            << std::endl;
#endif

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
  std::ofstream&     ofs,
  const std::string& prefixName) const
{
  unsigned int tmpSize = m_domainValues.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofs << "\n" << prefixName << "Time = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Value = zeros(" << tmpSize << ",1);";
  }
  else {
    ofs << "\n" << prefixName << "Time = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Value = zeros(" << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofs << "\n" << prefixName << "Time("  << i+1 << ",1) = " << m_domainValues[i]  << ";"
          << "\n" << prefixName << "Value(" << i+1 << ",1) = " << m_imageValues[i] << ";";
    }
  }

  return;
}

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
  void                       printForMatlab(std::ofstream& ofs, const std::string& prefixName) const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;
  using uqSampled1D1DFunctionClass::m_domainValues;
  using uqSampled1D1DFunctionClass::m_imageValues;

  std::vector<double> m_integratedValues;
};

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
  // FIX ME: check range of domainValue
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
  std::ofstream&     ofs,
  const std::string& prefixName) const
{
  uqSampled1D1DFunctionClass::printForMatlab(ofs,prefixName);

  unsigned int tmpSize = m_integratedValues.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofs << "\n" << prefixName << "intValue = zeros("  << tmpSize << ",1);";
  }
  else {
    ofs << "\n" << prefixName << "intValue = zeros("  << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofs << "\n" << prefixName << "intValue(" << i+1 << ",1) = " << m_integratedValues[i]  << ";";
    }
  }

  return;
}
#endif // __UQ_1D_1D_FUNCTION_H__

