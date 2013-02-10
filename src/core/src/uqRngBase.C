//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#include <uqRngBase.h>
#include <mpi.h>

uqRngBaseClass::uqRngBaseClass()
  :
  m_seed     (0),
  m_worldRank(UQ_UNAVAILABLE_RANK)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "uqRngBaseClass::constructor(), default",
                      "should not be used by user");
}

uqRngBaseClass::uqRngBaseClass(int seed, int worldRank)
  :
  m_seed     (seed),
  m_worldRank(worldRank)
{
  privateResetSeed();
}

uqRngBaseClass::~uqRngBaseClass()
{
}

int
uqRngBaseClass::seed() const
{
  return m_seed;
}

void
uqRngBaseClass::resetSeed(int newSeed)
{
  m_seed = newSeed;
  privateResetSeed();
  return;
}

void
uqRngBaseClass::privateResetSeed()
{
  if (m_seed >= 0) {
    // Do nothing
  }
  else if (m_seed < 0) {
    m_seed = (-m_seed+m_worldRank);
  }
  //else {
  //  struct timeval timevalNow;
  //  /*iRC = */gettimeofday(&timevalNow, NULL);
  //  gsl_rng_default_seed = (unsigned long int) timevalNow.tv_usec;
  //}

  return;
}
