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

#ifndef UQ_INV_LOGIT_GAUSSIAN_REALIZER_H
#define UQ_INV_LOGIT_GAUSSIAN_REALIZER_H

#include <queso/VectorRealizer.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class InvLogitGaussianVectorRealizer
 * \brief A class for handling sampling from (transformed) Gaussian probability density distributions with bounds
 *
 * To get a realization from a transformed Gaussian (i.e. with bounds), an
 * inverse \c logit transform is applied.  I.e., samples are \f$ f(X) \f$ where
 * \f$ X \f$ is drawn from a Gaussian and where
 *
 * \f[
 *   f(x) = \frac{b \exp(x) + a}{1 + \exp(x)}.
 * \f]
 *
 * This will produce a sample in the closed interval [a, b].
 */

template <class V = GslVector, class M = GslMatrix>
class InvLogitGaussianVectorRealizer : public BaseVectorRealizer<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*!
   * Constructs a new object, given a prefix and the image set of the vector
   * realizer, a vector of mean values, \c lawExpVector (of the Gaussian, not
   * the transformed Gaussian), and a lower triangular matrix resulting from
   * Cholesky decomposition of the covariance matrix, \c lowerCholLawCovMatrix.
   */
  InvLogitGaussianVectorRealizer(const char * prefix,
      const BoxSubset<V, M> & unifiedImageBoxSubset, const V & lawExpVector,
      const M & lowerCholLawCovMatrix);

  //! Constructor
  /*!
   * Constructs a new object, given a prefix and the image set of the vector
   * realizer, a vector of mean values, \c lawExpVector (of the Gaussian, not
   * the transformed Gaussian), and a set of two matrices and one vector
   * resulting from the Single Value Decomposition of the covariance matrix,
   * \c matU, \c vecSsqrt and \c matVt.
   */
  InvLogitGaussianVectorRealizer(const char * prefix,
      const BoxSubset<V, M> & unifiedImageBoxSubset, const V & lawExpVector,
      const M & matU, const V & vecSsqrt, const M & matVt);

  //! Destructor
  ~InvLogitGaussianVectorRealizer();
  //@}

  //! @name Realization-related methods
  //@{
  //! Access to the vector of mean values of the Gaussian and private attribute:  m_unifiedLawExpVector.
  const V & unifiedLawExpVector() const;

  //! Access to the vector of variance values and private attribute:  m_unifiedLawVarVector.
  const V & unifiedLawVarVector() const;

  //! Draws a realization.
  /*!
   * This function draws a realization of a (transformed) Gaussian
   * distribution with mean \c m_unifiedLawExpVector and variance
   * \c m_unifiedLawVarVector and saves it in \c nextValues.
   */
  void realization(V & nextValues) const;

  //! Updates the mean of the Gaussian with the new value \c newLawExpVector.
  void updateLawExpVector(const V & newLawExpVector);

  //! Updates the lower triangular matrix from Cholesky decomposition of the covariance matrix to the new value \c newLowerCholLawCovMatrix.
  /*! The lower triangular matrix results resulting from a Cholesky
   * decomposition of the covariance matrix. This routine deletes old expected
   * values: m_lowerCholLawCovMatrix; m_matU, m_vecSsqrt, m_matVt.
   */
  void updateLowerCholLawCovMatrix(const M & newLowerCholLawCovMatrix);

  //! Updates the SVD matrices from SVD decomposition of the covariance matrix to the new values: \c matU, \c vecSsqrt, and \c matVt.
  /*! The lower triangular matrix results resulting from a Cholesky
   * decomposition of the covariance matrix. This routine deletes old expected
   * values: m_lowerCholLawCovMatrix; m_matU, m_vecSsqrt, m_matVt.
   */
  void updateLowerCholLawCovMatrix(const M & matU, const V & vecSsqrt,
      const M & matVt);
  //@}

private:
  V * m_unifiedLawExpVector;
  V * m_unifiedLawVarVector;
  M * m_lowerCholLawCovMatrix;
  M * m_matU;
  V * m_vecSsqrt;
  M * m_matVt;

  using BaseVectorRealizer<V, M>::m_env;
  using BaseVectorRealizer<V, M>::m_prefix;
  using BaseVectorRealizer<V, M>::m_unifiedImageSet;
  using BaseVectorRealizer<V, M>::m_subPeriod;

  // For easy access to the bounds of the domain
  const BoxSubset<V, M> & m_unifiedImageBoxSubset;
};

}  // End namespace QUESO

#endif // UQ_INV_LOGIT_GAUSSIAN_REALIZER_H
