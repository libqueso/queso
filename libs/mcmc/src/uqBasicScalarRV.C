/* uq/libs/mcmc/src/uqBasicScalarRV.C
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqBasicScalarRV.h>

uqBasicScalarRVClass::uqBasicScalarRVClass(
  const std::string& name,
  double             initialValue,
  double             minValue,
  double             maxValue,
  double             priorMu,
  double             priorSigma)
  :
  m_name        (name),
  m_initialValue(initialValue),
  m_minValue    (minValue),
  m_maxValue    (maxValue),
  m_priorMu     (priorMu),
  m_priorSigma  (priorSigma)
{
}

uqBasicScalarRVClass::~uqBasicScalarRVClass()
{
}

std::string
uqBasicScalarRVClass::name() const
{
  return m_name;
}

double
uqBasicScalarRVClass::initialValue() const
{
  return m_initialValue;
}

double
uqBasicScalarRVClass::minValue() const
{
  return m_minValue;
}

double
uqBasicScalarRVClass::maxValue() const
{
  return m_maxValue;
}

double
uqBasicScalarRVClass::priorMu() const
{
  return m_priorMu;
}

double
uqBasicScalarRVClass::priorSigma() const
{
  return m_priorSigma;
}

void
uqBasicScalarRVClass::setName(const std::string& name)
{
  m_name = name;
  return;
}

void
uqBasicScalarRVClass::setInitialValue(double initialValue)
{
  m_initialValue = initialValue;
  return;
}

void
uqBasicScalarRVClass::setMinValue(double minValue)
{
  m_minValue = minValue;
  return;
}

void
uqBasicScalarRVClass::setMaxValue(double maxValue)
{
  m_maxValue = maxValue;
  return;
}

void
uqBasicScalarRVClass::setPriorMu(double priorMu)
{
  m_priorMu = priorMu;
  return;
}

void
uqBasicScalarRVClass::setPriorSigma(double priorSigma)
{
  m_priorSigma = priorSigma;
  return;
}

void
uqBasicScalarRVClass::print(std::ostream& os) const
{
  os << this->m_name         << " "
     << this->m_initialValue << " "
     << this->m_minValue     << " "
     << this->m_maxValue     << " "
     << this->m_priorMu      << " "
     << this->m_priorSigma;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqBasicScalarRVClass& obj)
{
  obj.print(os);

  return os;
}
