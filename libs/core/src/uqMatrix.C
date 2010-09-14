/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <uqMatrix.h>
#include <uqDefines.h>

uqMatrixClass::uqMatrixClass()
  :
  m_env(*(new uqEmptyEnvironmentClass())                               ),
  m_map(*(new Epetra_Map( 1,0,*(new Epetra_MpiComm(MPI_COMM_WORLD)) ) ))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqMatrixClass::constructor(), default",
                      "should not be used by user");
}

uqMatrixClass::uqMatrixClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map)
  :
  m_env              (env),
  m_map              (map),
  m_printHorizontally(true),
  m_inDebugMode      (false)
{
}

uqMatrixClass::uqMatrixClass(const uqMatrixClass& rhs)
  :
  m_env(rhs.m_env),
  m_map(rhs.m_map)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "uqMatrixClass::constructor(), copy",
                      "code should not execute through here");
}

uqMatrixClass::~uqMatrixClass()
{
}

void
uqMatrixClass::copy(const uqMatrixClass& src)
{
  //m_env = src.env;
  //m_map = src.map;
  m_printHorizontally = src.m_printHorizontally;
  m_inDebugMode       = src.m_inDebugMode;

  return;
}

uqMatrixClass&
uqMatrixClass::operator= (const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "uqMatrixClass::operator=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "uqMatrixClass::operator*=()",
                      "code should not execute through here");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator+=(const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "uqMatrixClass::operator+=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator-=(const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "uqMatrixClass::operator-=()",
                      "code should not execute through here");
  return *this;
}

const uqBaseEnvironmentClass&
uqMatrixClass::env() const
{
  return m_env;
}

const Epetra_Map&
uqMatrixClass::map() const
{
  return m_map;
  //return (const Epetra_Map&) (m_mat->Map());
}

unsigned int
uqMatrixClass::numOfProcsForStorage() const
{
  return m_map.Comm().NumProc();
}

void
uqMatrixClass::setPrintHorizontally(bool value) const
{
  m_printHorizontally = value;
  return;
}

bool
uqMatrixClass::getPrintHorizontally() const
{
  return m_printHorizontally;
}

void
uqMatrixClass::setInDebugMode(bool value) const
{
  m_inDebugMode = value;
  return;
}

bool
uqMatrixClass::getInDebugMode() const
{
  return m_inDebugMode;
}
