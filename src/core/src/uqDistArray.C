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

#include <uqDistArray.h>
#ifdef QUESO_HAS_TRILINOS

// 'uqDistArrayClass<T>::type' is just an alias to the 'EpetraExt::DistArray<T>' class of Trilinos

#else // QUESO_HAS_TRILINOS

uqDistArrayClass::uqDistArrayClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      UQ_UNAVAILABLE_RANK,
                      "uqDistArrayClass::constructor()"
                      "should not be called");
}

uqDistArrayClass::uqDistArrayClass(
  unsigned int          numGlobalElements,
  unsigned int          numNotUsed,
  const uqMpiCommClass& comm)
  :
  m_comm             (comm),
  m_numGlobalElements(numGlobalElements),
  m_numMyElements    (numGlobalElements)
{
  if (numNotUsed > 1) { // just to avoid compiler warning
    // Do nothing
  }
}

uqDistArrayClass::uqDistArrayClass(const uqDistArrayClass& src)
{
  this->copy(src);
}

uqDistArrayClass&
uqDistArrayClass::operator=(const uqDistArrayClass& rhs)
{
  this->copy(rhs);
  return *this;
}

uqDistArrayClass::~uqDistArrayClass()
{
  // Nothing to do
}

void
uqDistArrayClass::copy(const uqDistArrayClass& src)
{
  m_comm              = src.m_comm;
  m_numGlobalElements = src.m_numGlobalElements;
  m_numMyElements     = src.m_numMyElements;

  return;
}

const uqMpiCommClass&
uqDistArrayClass::Comm() const
{
  return m_comm;
}

unsigned int
uqDistArrayClass::NumGlobalElements() const
{
  return m_numGlobalElements;
}

unsigned int
uqDistArrayClass::NumMyElements() const
{
  return m_numMyElements;
}

#endif // QUESO_HAS_TRILINOS
