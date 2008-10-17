/* libs/fp/src/uq1DLocalDof.C
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

#include <uq1DLocalDof.h>

uq1DLocalDofClass::uq1DLocalDofClass(
  double                 bcc,
  uqLocalDofPositionEnum localPositionType,
  unsigned int           globalIdOfRespectiveNode)
  :
  m_bcc                     (bcc),
  m_myLocalPositionType     (localPositionType),
  m_globalIdOfRespectiveNode(globalIdOfRespectiveNode),
  m_myGlobalId              (UQ_INVALID_DOF_ID)
{
}

uq1DLocalDofClass::uq1DLocalDofClass(const uq1DLocalDofClass& obj)
{
  this->copy(obj);
}

uq1DLocalDofClass&
uq1DLocalDofClass::operator=(const uq1DLocalDofClass& rhs)
{
  this->copy(rhs);

  return *this;
}

void
uq1DLocalDofClass::copy(const uq1DLocalDofClass& src)
{
  m_bcc                      = src.m_bcc;
  m_myLocalPositionType      = src.m_myLocalPositionType;
  m_globalIdOfRespectiveNode = src.m_globalIdOfRespectiveNode;
  m_myGlobalId               = src.m_myGlobalId;

  return;
}

uq1DLocalDofClass::~uq1DLocalDofClass()
{
}

double
uq1DLocalDofClass::bcc() const
{
  return m_bcc;
}

unsigned int
uq1DLocalDofClass::globalId() const
{
  UQ_FATAL_TEST_MACRO((m_myGlobalId == UQ_INVALID_DOF_ID),
                      UQ_UNAVAILABLE_RANK,
                      "uq1DLocalDofClass::globalId()",
                      "m_myGlobalId was not set yet");

  return m_myGlobalId;
}

void
uq1DLocalDofClass::setGlobalId(unsigned int globalId)
{
  UQ_FATAL_TEST_MACRO((m_myGlobalId != UQ_INVALID_DOF_ID),
                      UQ_UNAVAILABLE_RANK,
                      "uq1DLocalDofClass::globalId()",
                      "m_myGlobalId was already set");

  m_myGlobalId = globalId;

  return;
}

unsigned int
uq1DLocalDofClass::globalIdOfRespectiveNode() const
{
  return m_globalIdOfRespectiveNode;
}

void
uq1DLocalDofClass::print(std::ostream& os) const
{
  return;
}
