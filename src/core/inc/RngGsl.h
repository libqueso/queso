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

#ifndef UQ_RNG_GSL_H
#define UQ_RNG_GSL_H

#include <queso/RngBase.h>

#include <queso/Defines.h>
#include <gsl/gsl_rng.h>


/*! \file RngGsl.h
    \brief GSL Random Number Generation class.
*/

/*! \class RngGsl
    \brief Class for random number generation using GSL library.

    This class is  class of random number generator, using GSL library. The GSL library
    provides a large collection of random number generators which can be accessed through
    a uniform interface.
*/


extern unsigned long int gsl_rng_default_seed;

namespace QUESO {

class RngGsl : public RngBase
{
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Constructor with seed.
  /* Uses GSL with the generator of type: gsl_rng_ranlxd2. This is the second, of a total of three
   luxury random numbers, all an extension of a second-generation version of the ranlux algorithm of
   Lüscher. The period of the generator is about 10^171. The algorithm has mathematically proven
   properties and can provide truly decorrelated numbers at a known level of randomness. The higher
   luxury levels provide increased decorrelation between samples as an additional safety margin.*/
  RngGsl(int seed, int worldRank);

  //! Destructor
 ~RngGsl();
 //@}

   //! @name Sampling methods
  //@{
  //! Resets the seed with value \c newSeed.
  void     resetSeed     (int newSeed);

  //! Samples a value from a uniform distribution. Support: [0,1] or [a,b].
  /*! This function samples from continuous uniform distribution on the range [0,1). It is
   * possible to scale this distribution so the support is defined by the two parameters,
   * a and b, which are its minimum and maximum values. Support: -infinity < a < x< b< infinity.
   * Uses gsl_rng_uniform(m_rng).*/
  double   uniformSample ()                          const;

  //! Samples a value from a Gaussian distribution with standard deviation given by \c stdDev.
  //! Support:  (-infinity, infinity).
  /*! The parameter mu (mean or expectation of the distribution) in this Gaussian sample is
   * set to zero, and thus, needs to be provided in an alternative way (e.g., in the form of
   * a sum. The parameter stdDev is its standard deviation; its variance is therefore stdDev^2.
   * A random variable with a Gaussian distribution is said to be normally distributed and is
   * called a normal deviate. Uses gsl_ran_gaussian(). Support: (-infinity, infinity). */
  double   gaussianSample(double stdDev)             const;

    //! Samples a value from a Beta distribution. Support: [0,1]
  /*! The Beta Distribution is a continuous probability distribution; it has two free
   * parameters, which are labeled alpha and beta. The beta distribution is used as a
   * prior distribution for binomial proportions in Bayesian analysis  (Evans et al.
   * 2000, p. 34). Uses gsl_ran_beta(). Support (domain): [0,1].*/
  double   betaSample    (double alpha, double beta) const;

  //! Samples a value from a Gamma distribution. Support: [0,infinity).
  /*! The Gamma Distribution is a continuous probability distribution; it has two
   * free parameters, which may be labeled: a shape parameter a and an inverse
   * scale parameter b, called a rate parameter. Uses gsl_ran_gamma(). Support
   * (domain): [0,infinity).*/
  double   gammaSample   (double a, double b)        const;

  //! GSL random number generator.
  const gsl_rng* rng           () const;

protected:
  //! Default Constructor: it should not be used.
  RngGsl();

  //! GSL random number generator.
  /*! It is chosen, in the constructor, to be of type gsl_rng_ranlxd2. */
  gsl_rng* m_rng;
};

}  // End namespace QUESO

#endif // UQ_RNG_GSL_H
