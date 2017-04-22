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

#ifndef UQ_MALA_TK_H
#define UQ_MALA_TK_H

#include <queso/TKGroup.h>

namespace QUESO {

class GslVector;
class GslMatrix;
template <class V, class M> class BayseianJointPdf;

/*!
 * \class MetropolisAdjustedLangevinTK
 *
 * \brief This class allows the representation of the MALA transition kernel
 * with a scaled covariance matrix for the purposes of delayed rejection.
 */
template <class V = GslVector, class M = GslMatrix>
class MetropolisAdjustedLangevinTK : public BaseTKGroup<V, M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  MetropolisAdjustedLangevinTK(const char * prefix,
                                const BayesianJointPdf<V, M> & targetPdf,
                                const std::vector<double> & scales,
                                const M & covMatrix);
  //! Destructor.
  ~MetropolisAdjustedLangevinTK();
  //@}

  //! @name Statistical/Mathematical methods
  //@{
  //! Whether or not the matrix is symmetric. Always 'true'.
  /*! \todo: It only returns 'true', thus a test for its symmetricity must be included.*/
  bool symmetric() const;

  //! Gaussian increment property to construct a transition kernel.
  const GaussianVectorRV<V, M> & rv(unsigned int stageId) const;

  //! Gaussian increment property to construct a transition kernel.
  const GaussianVectorRV<V, M> & rv(const std::vector<unsigned int> & stageIds);

  //! Constructs transition kernel pdf based on internal \c m_stageId variable
  /*
   * This uses the formula for the transition kernel in Roberts & Tweedie 2001.
   *
   * http://projecteuclid.org/download/pdf_1/euclid.bj/1178291835
   */
  virtual const GaussianVectorRV<V, M> & rv(const V & position) const;

  //! Scales the covariance matrix.
  /*! The covariance matrix is scaled by a factor of \f$ 1/scales^2 \f$.*/
  void updateLawCovMatrix(const M & covMatrix);
  //@}

  //! @name Misc methods
  //@{
  //! Sets the pre-computing positions \c m_preComputingPositions[stageId] with a new vector of size \c position.
  bool setPreComputingPosition(const V & position, unsigned int stageId);

  //! Clears the pre-computing positions \c m_preComputingPositions[stageId]
  void clearPreComputingPositions();

  virtual bool covMatrixIsDirty() { return false; }
  virtual void cleanCovMatrix() { }
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the transition kernel.
  /*! \todo: implement me!*/
  void print(std::ostream & os) const;
  //@}
private:
  //! Sets the mean of the RVs to zero.
  void setRVsWithZeroMean();
  using BaseTKGroup<V,M>::m_env;
  using BaseTKGroup<V,M>::m_prefix;
  using BaseTKGroup<V,M>::m_vectorSpace;
  using BaseTKGroup<V,M>::m_scales;
  using BaseTKGroup<V,M>::m_preComputingPositions;
  using BaseTKGroup<V,M>::m_rvs;

  M m_originalCovMatrix;

  const BayesianJointPdf<V, M> & m_targetPdf;

  double m_time_step;
};

}  // End namespace QUESO

#endif  // UQ_MALA_TK_H
