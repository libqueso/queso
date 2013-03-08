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

#include <uqRngBoost.h>
#include <mpi.h>

uqRngBoostClass::uqRngBoostClass()
  :
  uqRngBaseClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "uqRngBoostClass::constructor(), default",
                      "should not be used by user");
}

uqRngBoostClass::uqRngBoostClass(int seed, int worldRank)
  :
  uqRngBaseClass(seed,worldRank)
{
  // Kemelli todo
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "uqRngBoosClass::constructor()",
                      "Kemelli todo: boost rng");
}

uqRngBoostClass::~uqRngBoostClass()
{
  // Kemelli todo
}

void
uqRngBoostClass::resetSeed(int newSeed)
{
  uqRngBaseClass::resetSeed(newSeed);

  // Kemelli todo

  return;
}

double
uqRngBoostClass::uniformSample() const
{
  // Kemelli todo
  return 0.;
}

double
uqRngBoostClass::gaussianSample(double stdDev) const
{
  // Kemelli todo
  return 0.;
}

double
uqRngBoostClass::betaSample(double alpha, double beta) const
{
  // Kemelli todo
  return 0.;
}

double
uqRngBoostClass::gammaSample(double a, double b) const
{
  // Kemelli todo
  return 0.;
}
