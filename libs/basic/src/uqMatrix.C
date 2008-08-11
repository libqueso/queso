/* libs/basic/src/uqMatrix.C
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

#include <uqMatrix.h>
#include <uqDefines.h>

uqMatrixClass::uqMatrixClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqMatrixClass::constructor(), default",
                      "should not be used by user");
}

uqMatrixClass::uqMatrixClass(const uqEnvironmentClass& env)
  :
  m_env(env)
{
}

uqMatrixClass::uqMatrixClass(const uqMatrixClass& rhs)
  :
  m_env(rhs.m_env)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::constructor(), copy",
                      "code should not execute through here");
}

uqMatrixClass::~uqMatrixClass()
{
}

void
uqMatrixClass::copy(const uqMatrixClass& src)
{
  //m_env = env;

  return;
}

uqMatrixClass&
uqMatrixClass::operator= (const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator*=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator+=(const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator+=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator-=(const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator-=()",
                      "code should not execute through here");
  return *this;
}

const uqEnvironmentClass&
uqMatrixClass::env() const
{
  return m_env;
}
