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

#ifndef UQ_SCALAR_GAUSSIAN_RANDOM_FIELD_H
#define UQ_SCALAR_GAUSSIAN_RANDOM_FIELD_H

#include <queso/ScalarCovarianceFunction.h>
#include <queso/ScalarFunction.h>
#include <queso/GaussianVectorRV.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \file ScalarGaussianRandomField.h
 * \brief A class for handling Gaussian random fields (GRF).
 *
 * \class ScalarGaussianRandomField
 * \brief A class for handling scalar Gaussian random fields (GRF).
 *
 * This class implements a scalar Gaussian random field (GRF); i.e. a random field involving
 * Gaussian probability density functions (PDFs) of the variables. A one-dimensional GRF is
 * also called a Gaussian process.*/

template <class V = GslVector, class M = GslMatrix>
class ScalarGaussianRandomField
{
 public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! Constructs a new object, given a prefix, an index set, and both a mean and a
   * covariance function. This method deletes the previous saved positions. */
  ScalarGaussianRandomField(const char*                                     prefix,
                                   const VectorSet<V,M>&                    indexSet,
                                   const BaseScalarFunction<V,M>&           meanFunction,
                                   const BaseScalarCovarianceFunction<V,M>& covarianceFunction);

  //! TODO: Copy constructor.
  /*! \todo: implement me!*/
  ScalarGaussianRandomField(const ScalarGaussianRandomField&         obj);

  //! Destructor.
  ~ScalarGaussianRandomField();
  //@}

  //! @name Set methods
  //@{
  //! TODO: Assignment operator; it copies \c rhs to \c this.
  /*! \todo: implement me!*/
  ScalarGaussianRandomField& operator=(const ScalarGaussianRandomField& rhs);
  //@}

  //! @name Math methods
  //@{
  //! Index set; access to protected attribute m_indexSet.
  const VectorSet<V,M>&           indexSet          () const;

  //! Mean function; access to protected attribute m_meanFunction.
  const BaseScalarFunction<V,M>&  meanFunction      () const;

  //! Covariance function; access to protected attribute m_covarianceFunction.
  const BaseScalarCovarianceFunction<V,M>& covarianceFunction() const;

  //! Function that samples from a Gaussian PDF.
  /*! Given the field positions, this method performs a number of tests, calculates the mean vector,
   * the covariance matrix and then it samples from a Gaussian random vector as
   * many positions as required.*/
  void                                  sampleFunction(const std::vector<V*>& fieldPositions, V& sampleValues);
  //@}

protected:
  //! Copy method.
  void                                  copy          (const ScalarGaussianRandomField& src);

  //! Environment.
  const BaseEnvironment&         m_env;

  //! Prefix.
  std::string                           m_prefix;

  //! Index set.
  const VectorSet<V,M>&          m_indexSet;

  //! Mean function.
  const BaseScalarFunction<V,M>& m_meanFunction;

  //! Covariance function.
  const BaseScalarCovarianceFunction<V,M>& m_covarianceFunction;

  //! Saved positions.
  std::vector<V*>                       m_savedPositions;

  //! Image set of the RV.
  VectorSpace<V,M>*              m_savedRvImageSpace;

  //! Vector of the mean value of the RV.
  V*                                    m_savedRvLawExpVector;

  //! Covariance matrix of the RV.
  M*                                    m_savedRvLawCovMatrix;

  //! My RV.
  GaussianVectorRV<V,M>*         m_savedRv;
};

}  // End namespace QUESO

#endif // UQ_SCALAR_GAUSSIAN_RANDOM_FIELD_H
