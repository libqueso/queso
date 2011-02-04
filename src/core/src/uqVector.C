//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqVector.h>
#include <uqDefines.h>

uqVectorClass::uqVectorClass()
  :
  m_env(*(new uqEmptyEnvironmentClass())                                    ),
  m_map(*(new uqMapClass( 1,0,*(new uqMpiCommClass(m_env,MPI_COMM_SELF)) ) )) // avoid using MPI_COMM_WORLD
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqVectorClass::constructor(), default",
                      "should not be used by user");
}

uqVectorClass::uqVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map)
  :
  m_env              (env),
  m_map              (map),
  m_printHorizontally(true),
  m_printScientific  (false) // for compatibility with previous regression tests
{
  //std::cout << "Entering uqVectorClass::constructor(env,map)" << std::endl;

  //std::cout << "Leaving uqVectorClass::constructor(env,map)" << std::endl;
}

uqVectorClass::uqVectorClass(const uqVectorClass& rhs)
  :
  m_env(rhs.m_env),
  m_map(rhs.m_map)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
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
  //m_map               = src.map; // prudenci 2010-06-17
  m_printHorizontally = src.m_printHorizontally;
  m_printScientific   = src.m_printScientific;

  return;
}

uqVectorClass&
uqVectorClass::operator=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "uqVectorClass::operator=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "uqVectorClass::operator*=()",
                      "code should not execute through here");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqVectorClass&
uqVectorClass::operator/=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "uqVectorClass::operator/=()",
                      "code should not execute through here");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqVectorClass&
uqVectorClass::operator+=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "uqVectorClass::operator+=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator-=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "uqVectorClass::operator-=()",
                      "code should not execute through here");
  return *this;
}

const uqBaseEnvironmentClass&
uqVectorClass::env() const
{
  return m_env;
}

const uqMapClass&
uqVectorClass::map() const
{
  return m_map;
}

unsigned int
uqVectorClass::numOfProcsForStorage() const
{
  return m_map.Comm().NumProc();
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

void
uqVectorClass::setPrintScientific(bool value) const
{
  m_printScientific = value;
  return;
}

bool
uqVectorClass::getPrintScientific() const
{
  return m_printScientific;
}
