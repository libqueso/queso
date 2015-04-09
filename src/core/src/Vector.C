//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/Vector.h>
#include <queso/Defines.h>

namespace QUESO {

// Default constructor ------------------------------

Vector::Vector()
  :
  m_env(*(new EmptyEnvironment())                                               ),
  m_map(*(new Map( 1,0,*(new MpiComm(m_env,RawValue_MPI_COMM_SELF)) ) )) // avoid using MPI_COMM_WORLD
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "Vector::constructor(), default",
                      "should not be used by user");
}

// Shaped constructor --------------------------------
Vector::Vector(const BaseEnvironment& env, const Map& map)
  :
  m_env              (env),
  m_map              (map),
  m_printHorizontally(true),
  m_printScientific  (false) // for compatibility with previous regression tests
{
  //std::cout << "Entering Vector::constructor(env,map)" << std::endl;

  //std::cout << "In Vector::constructor(env,map)" // mox
  //          << ", m_map.NumMyElements() = " << m_map.NumMyElements()
  //          << std::endl;

  //std::cout << "Leaving Vector::constructor(env,map)" << std::endl;
}

// Copy constructor -----------------------------------
Vector::Vector(const Vector& rhs)
  :
  m_env(rhs.m_env),
  m_map(rhs.m_map)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "Vector::constructor(), copy",
                      "code should not execute through here");
}

// Destructor ---------------------------------------
Vector::~Vector()
{
}

// Set methods --------------------------------------
Vector&
Vector::operator=(const Vector& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "Vector::operator=()",
                      "code should not execute through here");
  return *this;
}

// --------------------------------------------------
Vector&
Vector::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "Vector::operator*=()",
                      "code should not execute through here");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

// --------------------------------------------------
Vector&
Vector::operator/=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.worldRank(),
                      "Vector::operator/=()",
                      "code should not execute through here");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

// --------------------------------------------------
Vector&
Vector::operator+=(const Vector& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "Vector::operator+=()",
                      "code should not execute through here");
  return *this;
}

// --------------------------------------------------
Vector&
Vector::operator-=(const Vector& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.m_env.worldRank(),
                      "Vector::operator-=()",
                      "code should not execute through here");
  return *this;
}

// Environment and Map methods ----------------------
const BaseEnvironment&
Vector::env() const
{
  return m_env;
}

// --------------------------------------------------
const Map&
Vector::map() const
{
  return m_map;
}

// --------------------------------------------------
unsigned int
Vector::numOfProcsForStorage() const
{
  return m_map.Comm().NumProc();
}

// I/O and Miscellaneous methods --------------------
void
Vector::setPrintHorizontally(bool value) const
{
  m_printHorizontally = value;
  return;
}
// --------------------------------------------------
bool
Vector::getPrintHorizontally() const
{
  return m_printHorizontally;
}
// --------------------------------------------------
void
Vector::setPrintScientific(bool value) const
{
  m_printScientific = value;
  return;
}
// --------------------------------------------------
bool
Vector::getPrintScientific() const
{
  return m_printScientific;
}
// --------------------------------------------------
void
Vector::copy(const Vector& src)
{
  //m_env               = src.env;
  //m_map               = src.map; // prudenci 2010-06-17
  m_printHorizontally = src.m_printHorizontally;
  m_printScientific   = src.m_printScientific;

  return;
}

}  // End namespace QUESO
