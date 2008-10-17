/* libs/fp/src/uq1DNode.C
 * 
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
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

#include <uq1DNode.h>

uq1DNodeClass::uq1DNodeClass(
  unsigned int       globalId,
  double             x,
  uqNodePositionEnum nodePositionType,
  uqBCEnum           bcType,
  double             bcValue)
  :
  m_myGlobalId    (globalId),
  m_myX           (x),
  m_myPositionType(nodePositionType),
  m_myBCType      (bcType),
  m_myBCValue     (bcValue)
{
}

uq1DNodeClass::uq1DNodeClass(const uq1DNodeClass& obj)
{
  this->copy(obj);
}

uq1DNodeClass&
uq1DNodeClass::operator=(const uq1DNodeClass& rhs)
{
  this->copy(rhs);

  return *this;
}

void
uq1DNodeClass::copy(const uq1DNodeClass& src)
{
  m_myGlobalId     = src.m_myGlobalId;
  m_myX            = src.m_myX;
  m_myPositionType = src.m_myPositionType;
  m_myBCType       = src.m_myBCType;
  m_myBCValue      = src.m_myBCValue;

  return;
}

uq1DNodeClass::~uq1DNodeClass()
{
}

unsigned int
uq1DNodeClass::globalId() const
{
  return m_myGlobalId;
}

double
uq1DNodeClass::x() const
{
  return m_myX;
}

void
uq1DNodeClass::print(std::ostream& os) const
{
  return;
}
