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
  virtual  double inverseValue  (double imageValue ) const = 0;
  virtual  double inverseDeriv  (double imageValue ) const = 0;

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

  double value       (double domainValue) const;
  double deriv       (double domainValue) const;
  double inverseValue(double imageValue ) const;
  double inverseDeriv(double imageValue ) const;

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
  m_referenceDomainValue      (referenceDomainValue),
  m_referenceImageValue       (referenceImageValue),
  m_rateValue                 (rateValue)
{
}

uqLinear1D1DFunctionClass::~uqLinear1D1DFunctionClass()
{
}

double
uqLinear1D1DFunctionClass::value(double domainValue) const
{
  double imageValue = m_referenceImageValue + m_rateValue*(domainValue - m_referenceDomainValue);

  return imageValue;
}

double
uqLinear1D1DFunctionClass::deriv(double domainValue) const
{
  return m_rateValue;
}

double
uqLinear1D1DFunctionClass::inverseValue(double imageValue) const
{
  double domainValue = m_referenceDomainValue + (1./m_rateValue)*(imageValue - m_referenceImageValue);

  return domainValue;
}

double
uqLinear1D1DFunctionClass::inverseDeriv(double imageValue) const
{
  return 1./m_rateValue;
}
#endif // __UQ_1D_1D_FUNCTION_H__

