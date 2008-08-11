/* uq/libs/mcmc/src/uqParameter.C
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

#include <uqParameter.h>

uqParameterClass::uqParameterClass(
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

uqParameterClass::~uqParameterClass()
{
}

std::string
uqParameterClass::name() const
{
  return m_name;
}

double
uqParameterClass::initialValue() const
{
  return m_initialValue;
}

double
uqParameterClass::minValue() const
{
  return m_minValue;
}

double
uqParameterClass::maxValue() const
{
  return m_maxValue;
}

double
uqParameterClass::priorMu() const
{
  return m_priorMu;
}

double
uqParameterClass::priorSigma() const
{
  return m_priorSigma;
}

void
uqParameterClass::setName(const std::string& name)
{
  m_name = name;
  return;
}

void
uqParameterClass::setInitialValue(double initialValue)
{
  m_initialValue = initialValue;
  return;
}

void
uqParameterClass::setMinValue(double minValue)
{
  m_minValue = minValue;
  return;
}

void
uqParameterClass::setMaxValue(double maxValue)
{
  m_maxValue = maxValue;
  return;
}

void
uqParameterClass::setPriorMu(double priorMu)
{
  m_priorMu = priorMu;
  return;
}

void
uqParameterClass::setPriorSigma(double priorSigma)
{
  m_priorSigma = priorSigma;
  return;
}

void
uqParameterClass::print(std::ostream& os) const
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
operator<<(std::ostream& os, const uqParameterClass& param)
{
  param.print(os);

  return os;
}
