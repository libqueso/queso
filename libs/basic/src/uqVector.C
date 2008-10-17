/* libs/basic/src/uqVector.C
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

#include <uqVector.h>
#include <uqDefines.h>

uqVectorClass::uqVectorClass()
  :
  m_env(*(new uqEnvironmentClass())                                    ),
  m_map(*(new Epetra_Map( 1,0,*(new Epetra_MpiComm(MPI_COMM_WORLD)) ) ))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqVectorClass::constructor(), default",
                      "should not be used by user");
}

uqVectorClass::uqVectorClass(const uqEnvironmentClass& env, const Epetra_Map& map)
  :
  m_env              (env),
  m_map              (map),
  m_printHorizontally(true)
{
}

uqVectorClass::uqVectorClass(const uqVectorClass& rhs)
  :
  m_env(rhs.m_env),
  m_map(rhs.m_map)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::constructor(), copy",
                      "code should not execute through here");
}

uqVectorClass::~uqVectorClass()
{
}

void
uqVectorClass::copy(const uqVectorClass& src)
{
  //m_env               = src.env;
  //m_map               = src.map;
  m_printHorizontally = src.m_printHorizontally;

  return;
}

uqVectorClass&
uqVectorClass::operator= (const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator*=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator/=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator/=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator+=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator+=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator-=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator-=()",
                      "code should not execute through here");
  return *this;
}

const uqEnvironmentClass&
uqVectorClass::env() const
{
  return m_env;
}

const Epetra_Map&
uqVectorClass::map() const
{
  return m_map;
}

void
uqVectorClass::setPrintHorizontally(bool value) const
{
  m_printHorizontally = value;
  return;
}

bool
uqVectorClass::getPrintHorizontally() const
{
  return m_printHorizontally;
}
