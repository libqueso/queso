/* uq/examples/queso/pyramid/uq1D1DFunction.h
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

#ifndef __UQ_1D_1D_FUNCTION_H__
#define __UQ_1D_1D_FUNCTION_H__

#include <uqDefines.h>
#include <vector>
#include <math.h>

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
  // FIX ME: check range of domainValue
  return m_constantValue;
}

double
uqConstant1D1DFunctionClass::deriv(double domainValue) const
{
  // FIX ME: check range of domainValue
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
// Sampled 1D->1D class
//*****************************************************
class uqSampled1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqSampled1D1DFunctionClass();
  uqSampled1D1DFunctionClass(const std::vector<double>& domainValues,
                             const std::vector<double>& imageValues,
                             bool                       dataIsContinuous);
 ~uqSampled1D1DFunctionClass();

  double value(double domainValue) const;
  double deriv(double domainValue) const;

  const std::vector<double>& domainValues    () const;
  const std::vector<double>& imageValues     () const;
  const bool                 dataIsContinuous() const;

        void                 set(const std::vector<double>& domainValues,
                                 const std::vector<double>& imageValues,
                                 bool                       dataIsContinuous);

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_domainValues;
  std::vector<double> m_imageValues;
  bool                m_dataIsContinuous;
};

uqSampled1D1DFunctionClass::uqSampled1D1DFunctionClass()
  :
  uqBase1D1DFunctionClass(-INFINITY,INFINITY)
{
}

uqSampled1D1DFunctionClass::uqSampled1D1DFunctionClass(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues,
  bool                       dataIsContinuous)
  :
  uqBase1D1DFunctionClass(domainValues[0],domainValues[domainValues.size()-1]),
  m_domainValues    (domainValues.size(),0.),
  m_imageValues     (imageValues.size(), 0.),
  m_dataIsContinuous(dataIsContinuous)
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
  //std::cout << "In uqSampled1D1DFunctionClass::getValue()"
  //          << ": domainValue = "         << domainValue
  //          << ", tmpSize = "             << tmpSize
  //          << ", m_domainValues[0] = "   << m_domainValues[0]
  //          << ", m_domainValues[max] = " << m_domainValues[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::getValue()",
                      "m_domainValues.size() = 0");

  UQ_FATAL_TEST_MACRO(domainValue < m_domainValues[0],
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::getValue()",
                      "domainValue < m_domainValues[0]");

  UQ_FATAL_TEST_MACRO(m_domainValues[tmpSize-1] < domainValue,
                      UQ_UNAVAILABLE_RANK,
                      "uqSampled1D1DFunctionClass::getValue()",
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
    if (m_dataIsContinuous) {
      double ratio = (domainValue - m_domainValues[i-1])/(m_domainValues[i]-m_domainValues[i-1]);
      returnValue = m_imageValues[i-1] + ratio * (m_imageValues[i]-m_imageValues[i-1]);
    }
    else {
      // Leave returnValue = 0.;
    }
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

const bool
uqSampled1D1DFunctionClass::dataIsContinuous() const
{
  return m_dataIsContinuous;
}

void
uqSampled1D1DFunctionClass::set(
  const std::vector<double>& domainValues,
  const std::vector<double>& imageValues,
  bool                       dataIsContinuous)
{
  m_domainValues.clear();
  m_imageValues.clear();

  unsigned int tmpSize = m_domainValues.size();
  m_minDomainValue = domainValues[0];
  m_maxDomainValue = domainValues[tmpSize-1];

  m_dataIsContinuous = dataIsContinuous;
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_domainValues[i] = domainValues[i];
    m_imageValues [i] = imageValues [i];
  }

  return;
}
#endif // __UQ_1D_1D_FUNCTION_H__

