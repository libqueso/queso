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

#include <queso/Defines.h>
#include <queso/RngGsl.h>
#include <gsl/gsl_randist.h>

namespace QUESO {

//! Constructor with seed ---------------------------
RngGsl::RngGsl(int seed, int worldRank)
  :
  RngBase(seed,worldRank),
  m_rng         (NULL)
{
  gsl_rng_default_seed = (unsigned long int) m_seed;
  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  queso_require_msg(m_rng, "null m_rng");

//gsl_rng_set(m_rng, gsl_rng_default_seed);
#if 0
  if (m_worldRank == 0) {
    std::cout << "In RngGsl::constructor():"
              << "\n  m_seed = "                                              << m_seed
              << "\n  internal seed = "                                       << gsl_rng_default_seed
            //<< "\n  first generated sample from uniform distribution = "    << gsl_rng_uniform(m_rng)
            //<< "\n  first generated sample from std normal distribution = " << gsl_ran_gaussian(m_rng,1.)
              << std::endl;
  }
#endif
}

// Destructor ---------------------------------------
RngGsl::~RngGsl()
{
  if (m_rng) gsl_rng_free(m_rng);
}

// Sampling methods ---------------------------------
void
RngGsl::resetSeed(int newSeed)
{
  RngBase::resetSeed(newSeed);
  gsl_rng_free(m_rng);

  gsl_rng_default_seed = (unsigned long int) m_seed;
  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  queso_require_msg(m_rng, "null m_rng");
  return;
}

// --------------------------------------------------
double
RngGsl::uniformSample() const
{
  return gsl_rng_uniform(m_rng);
}

// --------------------------------------------------
double
RngGsl::gaussianSample(double stdDev) const
{
  return gsl_ran_gaussian(m_rng,stdDev);
}

// --------------------------------------------------
double
RngGsl::betaSample(double alpha, double beta) const
{
  return gsl_ran_beta(m_rng,alpha,beta);
}

// --------------------------------------------------
double
RngGsl::gammaSample(double a, double b) const
{
  return gsl_ran_gamma(m_rng,a,b);
}

const gsl_rng*
RngGsl::rng() const
{
  return m_rng;
}

}  // End namespace QUESO
