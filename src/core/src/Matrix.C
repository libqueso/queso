//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <queso/Matrix.h>
#include <queso/Defines.h>

namespace QUESO {

// Shaped constructor --------------------------------
Matrix::Matrix(const BaseEnvironment& env, const Map& map)
  :
  m_env              (env),
  m_map              (map),
  m_printHorizontally(true),
  m_inDebugMode      (false)
{
}

// Destructor ---------------------------------------
Matrix::~Matrix()
{
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
Matrix::base_copy(const Matrix& src)
{
  //m_env = src.env;
  //m_map = src.map;
  m_printHorizontally = src.m_printHorizontally;
  m_inDebugMode       = src.m_inDebugMode;

  return;
}

}  // End namespace QUESO
