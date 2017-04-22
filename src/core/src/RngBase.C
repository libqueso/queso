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

#include <queso/RngBase.h>

namespace QUESO {

RngBase::RngBase(int seed, int worldRank)
  :
  m_seed     (seed),
  m_worldRank(worldRank)
{
  privateResetSeed();
}

RngBase::~RngBase()
{
}

int
RngBase::seed() const
{
  return m_seed;
}

void
RngBase::resetSeed(int newSeed)
{
  m_seed = newSeed;
  privateResetSeed();
  return;
}

void
RngBase::privateResetSeed()
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
  //  m_seed = (int) timevalNow.tv_usec;
  //}

  return;
}

}  // End namespace QUESO
