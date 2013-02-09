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

#include <uqRngGsl.h>

uqRngGslClass::uqRngGslClass()
  :
  uqRngBaseClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "uqRngGslClass::constructor(), default",
                      "should not be used by user");
}

uqRngGslClass::uqRngGslClass(int seed, int worldRank)
  :
  uqRngBaseClass(seed,worldRank),
  m_rng         (NULL)
{
  if (m_seed >= 0) {
    gsl_rng_default_seed = (unsigned long int) m_seed;
  }
  else if (m_seed < 0) {
    gsl_rng_default_seed = (unsigned long int) (-m_seed+m_worldRank);
  }
  //else {
  //  struct timeval timevalNow;
  //  /*iRC = */gettimeofday(&timevalNow, NULL);
  //  gsl_rng_default_seed = (unsigned long int) timevalNow.tv_usec;
  //}

  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  UQ_FATAL_TEST_MACRO((m_rng == NULL),
                      m_worldRank,
                      "uqRngGslClass::constructor()",
                      "null m_rng");

  //gsl_rng_set(m_rng, gsl_rng_default_seed);

  if ((m_subDisplayFile)/* && (this->displayVerbosity() > 0)*/) {
    *m_subDisplayFile << "In uqRngGslClass::constructor():"
                      << "\n  m_seed = "                                              << m_seed
                      << "\n  internal seed = "                                       << gsl_rng_default_seed
                    //<< "\n  first generated sample from uniform distribution = "    << gsl_rng_uniform(m_rng)
                    //<< "\n  first generated sample from std normal distribution = " << gsl_ran_gaussian(m_rng,1.)
                      << std::endl;
  }
}

uqRngGslClass::~uqRngGslClass()
{
  if (m_rng) gsl_rng_free(m_rng);
}

void
uqRngGslClass::resetSeed(int newSeed)
{
  uqRngBaseClass::resetSeed(int newSeed);

  gsl_rng_free(m_rng);

  if (newSeed >= 0) {
    gsl_rng_default_seed = (unsigned long int) newSeed;
  }
  else if (newSeed < 0) {
    gsl_rng_default_seed = (unsigned long int) (-newSeed+m_worldRank);
  }
  //else {
  //  struct timeval timevalNow;
  //  /*iRC = */gettimeofday(&timevalNow, NULL);
  //  gsl_rng_default_seed = (unsigned long int) timevalNow.tv_usec;
  //}

  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  UQ_FATAL_TEST_MACRO((m_rng == NULL),
                      m_worldRank,
                      "uqRngGslClass::resetGslSeed()",
                      "null m_rng");


double
uqRngGslClass::uniformSample() const
{
}

double
uqRngGslClass::gaussainSample() const
{
}

double
uqRngGslClass::betaSample() const
{
}

double
uqRngGslClass::gammaSample() const
{
}
