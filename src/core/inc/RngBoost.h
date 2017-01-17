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

#ifndef UQ_RNG_BOOST_H
#define UQ_RNG_BOOST_H

#include <queso/RngBase.h>
#include <boost/random/mersenne_twister.hpp>


/*! \file RngBoost.h
    \brief Boost Random Number Generation class.
*/

/*! \class RngBoost
    \brief Class for random number generation using Boost library.

    This class is  class of random number generator, using Boost library. The Boost Random Number
    Generator Library provides a framework for random number generators with well-defined properties
    so that the generators can be used in the demanding numerics and security domain.
*/

namespace QUESO {

class RngBoost : public RngBase
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor with seed.
  RngBoost(int seed, int worldRank);

  //! Destructor
 ~RngBoost();
  //@}

  //! @name Sampling methods
  //@{
  //! Resets the seed with value \c newSeed.
  void     resetSeed     (int newSeed);


  //! Samples a value from a uniform distribution. Support: [0,1] or [a,b].
  /*! This function samples from continuous uniform distribution on the range [0,1). It is
   * possible to scale this distribution so the support is defined by the two parameters,
   * a and b, which are its minimum and maximum values. Support: -infinity < a < x< b< infinity.
   * Uses boost::uniform_01<boost::mt19937> zeroone(m_rng).*/
  double   uniformSample ()                          const;


  //! Samples a value from a Gaussian distribution with standard deviation given by \c stdDev.
  //! Support:  (-infinity, infinity).
  /*! The parameter mu (mean or expectation of the distribution) in this Gaussian sample isn't
   * present, and thus, needs to be provided. e.g., in the form of a sum. The parameter stdDev is
   * its standard deviation; its variance is therefore stdDev^2. A random variable with a Gaussian
   * distribution is said to be normally distributed and is called a normal deviate. Uses
   * boost::math::normal_distribution<double>  gaussian_dist(mean, stdDev) Support (domain):
   * (-infinity, infinity). */
  double   gaussianSample(double stdDev)             const;

  //! Samples a value from a Beta distribution. Support: [0,1]
  /*! The Beta Distribution is a continuous probability distribution; it has two free
   * parameters, which are labeled alpha and beta. The beta distribution is used as a
   * prior distribution for binomial proportions in Bayesian analysis  (Evans et al.
   * 2000, p. 34). Uses boost::math::beta_distribution<double> beta_dist(alpha, beta).
   * Support (domain): [0,1].*/
  double   betaSample    (double alpha, double beta) const;

  //! Samples a value from a Gamma distribution. Support: [0,infinity).
   /*! The Gamma Distribution is a continuous probability distribution; it has two free
   * parameters, which may be labeled: a shape parameter a and an inverse scale parameter
   * b, called a rate parameter. The parameterization with a and b is more common in
   * Bayesian statistics, where the gamma distribution is used as a conjugate prior distribution
   * for various types of inverse scale (aka rate) parameters. Uses
   * boost::math::gamma_distribution<double>  gamma_dist(a,b). Support (domain): [0,infinity).*/
  double   gammaSample   (double a, double b)        const;

private:
  //! Default Constructor: it should not be used.
  RngBoost();

  //! Random number generator from class boost::mt19937.
  /*! mt19937 are models for a pseudo-random number generator. Here it is cannot be static,
   *  as it has not been initialized yet. mt19937 has length cycle of 2^(19937)-1, requires
   * approximately 625*sizeof(uint32_t) of memory, has relatively high speed (93% of the
   * fastest available in Boost library), and provides good uniform distribution in up to
   * 623 dimensions. */
  boost::mt19937 m_rng; // it cannot be static, as it is not initialized yet
};

}  // End namespace QUESO

#endif // UQ_RNG_BOOST_H
