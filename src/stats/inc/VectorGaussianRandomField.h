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

#ifndef UQ_VECTOR_GAUSSIAN_RANDOM_FIELD_H
#define UQ_VECTOR_GAUSSIAN_RANDOM_FIELD_H

#include <queso/MatrixCovarianceFunction.h>
#include <queso/VectorFunction.h>
#include <queso/GaussianVectorRV.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \file VectorGaussianRandomField.h
 * \brief A class for handling Gaussian random fields (GRF).
 *
 * \class VectorGaussianRandomField
 * \brief A class for handling vector Gaussian random fields (GRF).
 *
 * This class implements a vector Gaussian random field (GRF); i.e. a random field involving
 * vector Gaussian probability density functions (PDFs) of the variables. */

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class VectorGaussianRandomField
{
 public:
   //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! Constructs a new object, given a prefix, an index set, and both a mean and a
   * covariance function. This method deletes the previous saved positions. */
  VectorGaussianRandomField(const char*                                                 prefix,
                                   const VectorSet<P_V,P_M>&                            indexSet,
                                   const VectorSet<Q_V,Q_M>&                            imageSetPerIndex,
                                   const BaseVectorFunction<P_V,P_M,Q_V,Q_M>&           meanFunction,
                                   const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& covarianceFunction);

  //! TODO: Copy constructor.
  /*! \todo: implement me!*/
  VectorGaussianRandomField(const VectorGaussianRandomField&                     obj);

  //! Destructor.
  ~VectorGaussianRandomField();
  //@}

  //! @name Set methods
  //@{
  //! TODO: Assignment operator; it copies \c rhs to \c this.
  /*! \todo: implement me!*/
  VectorGaussianRandomField& operator=(const VectorGaussianRandomField& rhs);
  //@}

  //! @name Math methods
  //@{
  //! Index set; access to protected attribute m_indexSet.
  const VectorSet<P_V,P_M>&                   indexSet          () const;

  //! Mean function; access to protected attribute m_meanFunction.
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>&  meanFunction      () const;

  //! Covariance function; access to protected attribute m_covarianceFunction.
  const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& covarianceFunction() const;

  //! Function that samples from a Gaussian PDF.
  /*! Given the field positions, this method performs a number of tests, calculates the mean vector,
   * the covariance matrix and then it samples from a Gaussian random vector as
   * many positions as required.*/
  void                                              sampleFunction(const std::vector<P_V*>& fieldPositions, Q_V& sampleValues);
  //@}
protected:
  //! Copy method.
  void                                              copy          (const VectorGaussianRandomField& src);

  //! Environment.
  const BaseEnvironment&                     m_env;

  //! Prefix.
  std::string                                       m_prefix;

  //! Index set.
  const VectorSet<P_V,P_M>&                  m_indexSet;

  //! Image set of the RV, per index.
  const VectorSet<Q_V,Q_M>&                  m_imageSetPerIndex;

  //! Mean function.
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& m_meanFunction;

  //! Covariance function.
  const BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>& m_covarianceFunction;

  //! Saved positions.
  std::vector<P_V*>                                 m_savedPositions;

  //! Image set of the RV.
  VectorSpace<Q_V,Q_M>*                      m_savedRvImageSpace;

  //! Vector of the mean value of the RV.
  Q_V*                                              m_savedRvLawExpVector;

  //! Covariance matrix of the RV.
  Q_M*                                              m_savedRvLawCovMatrix;

  //! My RV.
  GaussianVectorRV<Q_V,Q_M>*                 m_savedRv;
};

}  // End namespace QUESO

#endif // UQ_VECTOR_GAUSSIAN_RANDOM_FIELD_H
