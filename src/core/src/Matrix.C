//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <queso/Matrix.h>
#include <queso/Defines.h>

namespace QUESO {

// Default constructor ------------------------------
Matrix::Matrix()
  :
  m_env(*(new EmptyEnvironment())                                               ),
  m_map(*(new Map( 1,0,*(new MpiComm(m_env,RawValue_MPI_COMM_SELF)) ) )) // avoid using MPI_COMM_WORLD
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "Matrix::constructor(), default",
                      "should not be used by user");
}

// Shaped constructor --------------------------------
Matrix::Matrix(const BaseEnvironment& env, const Map& map)
  :
  m_env              (env),
  m_map              (map),
  m_printHorizontally(true),
  m_inDebugMode      (false)
{
}

// Copy constructor -----------------------------------
Matrix::Matrix(const Matrix& rhs)
  :
  m_env(rhs.m_env),
  m_map(rhs.m_map)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "Matrix::constructor(), copy",
                      "code should not execute through here");
}

// Destructor ---------------------------------------
Matrix::~Matrix()
{
}

// Set methods --------------------------------------
Matrix&
Matrix::operator= (const Matrix& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "Matrix::operator=()",
                      "code should not execute through here");
  return *this;
}

// --------------------------------------------------
Matrix&
Matrix::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "Matrix::operator*=()",
                      "code should not execute through here");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

// --------------------------------------------------
Matrix&
Matrix::operator+=(const Matrix& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "Matrix::operator+=()",
                      "code should not execute through here");
  return *this;
}

// --------------------------------------------------
Matrix&
Matrix::operator-=(const Matrix& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "Matrix::operator-=()",
                      "code should not execute through here");
  return *this;
}

// Environment and Map methods ----------------------
const BaseEnvironment&
Matrix::env() const
{
  return m_env;
}

// --------------------------------------------------
const Map&
Matrix::map() const
{
  return m_map;
  //return (const Map&) (m_mat->Map());
}

// --------------------------------------------------
unsigned int
Matrix::numOfProcsForStorage() const
{
  return m_map.Comm().NumProc();
}

// I/O and Miscellaneous methods --------------------
void
Matrix::setPrintHorizontally(bool value) const
{
  m_printHorizontally = value;
  return;
}

// --------------------------------------------------
bool
Matrix::getPrintHorizontally() const
{
  return m_printHorizontally;
}

// --------------------------------------------------
void
Matrix::setInDebugMode(bool value) const
{
  m_inDebugMode = value;
  return;
}

// --------------------------------------------------
bool
Matrix::getInDebugMode() const
{
  return m_inDebugMode;
}

// --------------------------------------------------
void
Matrix::copy(const Matrix& src)
{
  //m_env = src.env;
  //m_map = src.map;
  m_printHorizontally = src.m_printHorizontally;
  m_inDebugMode       = src.m_inDebugMode;

  return;
}

}  // End namespace QUESO
