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

#ifndef UQ_TRANSFORMED_SCALEDCOV_TK_GROUP_H
#define UQ_TRANSFORMED_SCALEDCOV_TK_GROUP_H

#include <queso/TKGroup.h>
#include <queso/VectorRV.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/InvLogitGaussianVectorRV.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class TransformedScaledCovMatrixTKGroup
 * \brief This class represents a transition kernel with a scaled covariance matrix on hybrid bounded/unbounded state spaces.
 *
 * The unbounded directions utilise a standard Gaussian proposal.  The bounded
 * or half-bounded directions utilise a transformed Gaussian proposal, so that
 * no realizations are generated outside of the state space.
 */

template <class V = GslVector, class M = GslMatrix>
class TransformedScaledCovMatrixTKGroup : public BaseTKGroup<V, M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  TransformedScaledCovMatrixTKGroup(const char * prefix,
      const BoxSubset<V, M> & boxSubset, const std::vector<double> & scales,
      const M & covMatrix);

  //! Destructor.
  ~TransformedScaledCovMatrixTKGroup();
  //@}

  //! @name Statistical/Mathematical methods
  //@{
  //! Whether or not the matrix is symmetric.  Always 'false'.
  bool symmetric() const;

  //! InvLogitGaussian increment property to construct a transition kernel.
  const InvLogitGaussianVectorRV<V, M> & rv(unsigned int stageId) const;

  //! InvLogitGaussian increment property to construct a transition kernel.
  const InvLogitGaussianVectorRV<V, M> & rv(
      const std::vector<unsigned int> & stageIds);

  virtual const InvLogitGaussianVectorRV<V, M> & rv(const V & position) const;

  //! Scales the covariance matrix of the underlying Gaussian distribution.
  /*! The covariance matrix is scaled by a factor of \f$ 1/scales^2 \f$.*/
  void updateLawCovMatrix(const M & covMatrix);
  //@}

  //! @name Misc methods
  //@{
  //! Sets the pre-computing positions \c m_preComputingPositions[stageId] with a new vector of size \c position.
  /*!
   * The vector \c position is in *physical* space.  This is then transformed
   * using transformToGaussianSpace to map to a point in Gaussian space where
   * we can, for example, update the mean of the underlying Gaussian RV.
   */
  bool setPreComputingPosition(const V & position, unsigned int stageId);

  //! Clears the pre-computing positions \c m_preComputingPositions[stageId]
  void clearPreComputingPositions();

  virtual unsigned int set_dr_stage(unsigned int stageId);

  virtual bool covMatrixIsDirty() { return false; }
  virtual void cleanCovMatrix() { }
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the transition kernel.
  /*! \todo: implement me!*/
  void print(std::ostream & os) const;
  //@}

  // Convenience method that transforms a point in physical (user) space to a
  // space with no boundaries and where the proposal covariance matrix is
  // Gaussian
  void transformToGaussianSpace(const V & physicalPoint,
      V & transformedPoint) const;

private:
  //! Sets the mean of the underlying Gaussian RVs to zero.
  void setRVsWithZeroMean();

  using BaseTKGroup<V, M>::m_env;
  using BaseTKGroup<V, M>::m_prefix;
  using BaseTKGroup<V, M>::m_vectorSpace;
  using BaseTKGroup<V, M>::m_scales;
  using BaseTKGroup<V, M>::m_preComputingPositions;
  using BaseTKGroup<V, M>::m_rvs;

  const BoxSubset<V, M> & m_boxSubset;
  M m_originalCovMatrix;

};

}  // End namespace QUESO

#endif // UQ_TRANSFORMED_SCALEDCOV_TK_GROUP_H
