/* uq/libs/mcmc/src/uqScalarRV.C
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

#include <uqScalarRV.h>

uqBaseScalarRVClass::uqBaseScalarRVClass(
  double minValue,
  double maxValue,
  double expectValue,
  double stdDevValue)
  :
  m_minValue   (minValue),
  m_maxValue   (maxValue),
  m_expectValue(expectValue),
  m_stdDevValue(stdDevValue)
{
}

uqBaseScalarRVClass::~uqBaseScalarRVClass()
{
}

double
uqBaseScalarRVClass::minValue() const
{
  return m_minValue;
}

double
uqBaseScalarRVClass::maxValue() const
{
  return m_maxValue;
}

double
uqBaseScalarRVClass::expectValue() const
{
  return m_expectValue;
}

double
uqBaseScalarRVClass::stdDevValue() const
{
  return m_stdDevValue;
}

void
uqBaseScalarRVClass::setMinValue(double minValue)
{
  m_minValue = minValue;
  return;
}

void
uqBaseScalarRVClass::setMaxValue(double maxValue)
{
  m_maxValue = maxValue;
  return;
}

void
uqBaseScalarRVClass::setExpectValue(double expectValue)
{
  m_expectValue = expectValue;
  return;
}

void
uqBaseScalarRVClass::setStdDevValue(double stdDevValue)
{
  m_stdDevValue = stdDevValue;
  return;
}

void
uqBaseScalarRVClass::print(std::ostream& os) const
{
  os << this->m_minValue    << " "
     << this->m_maxValue    << " "
     << this->m_expectValue << " "
     << this->m_stdDevValue;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqBaseScalarRVClass& obj)
{
  obj.print(os);

  return os;
}
