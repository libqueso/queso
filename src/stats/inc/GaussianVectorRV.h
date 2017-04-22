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

#ifndef UQ_GAUSSIAN_VECTOR_RV_H
#define UQ_GAUSSIAN_VECTOR_RV_H

#include <queso/VectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/JointPdf.h>
#include <queso/VectorRealizer.h>
#include <queso/VectorCdf.h>
#include <queso/VectorMdf.h>
#include <queso/SequenceOfVectors.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Gaussian class [RV-03]
//*****************************************************

/*!
 * \class GaussianVectorRV
 * \brief A class representing a Gaussian vector RV.
 *
 * This class allows the user to compute the value of a Gaussian PDF and to generate realizations
 * (samples) from it.
 *
 * In probability theory, the normal (or Gaussian) distribution is a continuous probability
 * distribution, defined by the formula:
 * \f[    f(x| \mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}} e^{ -\frac{(x-\mu)^2}{2\sigma^2} }. \f]
 *
 * The parameter \f$ \mu \f$  in this formula is the mean or expectation of the distribution (and also
 * its median and mode). The parameter \f$ \sigma \f$  is its standard deviation; its variance is therefore
 * \f$ \sigma^2 \f$ . */

template <class V = GslVector, class M = GslMatrix>
class GaussianVectorRV : public BaseVectorRV<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Construct a Gaussian vector RV with mean \c lawExpVector and diagonal covariance matrix
   * \c lawVarVector whose variates live in \c imageSet.*/
  GaussianVectorRV(const char*                  prefix,
                          const VectorSet<V,M>& imageSet,
                          const V&                     lawExpVector,
                          const V&                     lawVarVector);

  //! Constructor
  /*! Construct a Gaussian vector RV with mean \c lawExpVector and covariance matrix
   * \c lawCovMatrix whose variates live in \c imageSet.*/
  GaussianVectorRV(const char*                  prefix,
                          const VectorSet<V,M>& imageSet,
                          const V&                     lawExpVector,
                          const M&                     lawCovMatrix);

  //! Virtual destructor
  virtual ~GaussianVectorRV();
  //@}

  //! @name Statistical methods
  //@{
  //! Updates the vector that contains the mean values.
  void updateLawExpVector(const V& newLawExpVector);

  //! Updates the covariance matrix.
  /*! This method tries to use Cholesky decomposition; and if it fails, the method then
   *  calls a SVD decomposition.*/
  void updateLawCovMatrix(const M& newLawCovMatrix);
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
 //@}

private:
  using BaseVectorRV<V,M>::m_env;
  using BaseVectorRV<V,M>::m_prefix;
  using BaseVectorRV<V,M>::m_imageSet;
  using BaseVectorRV<V,M>::m_pdf;
  using BaseVectorRV<V,M>::m_realizer;
  using BaseVectorRV<V,M>::m_subCdf;
  using BaseVectorRV<V,M>::m_unifiedCdf;
  using BaseVectorRV<V,M>::m_mdf;
};

//---------------------------------------------------
// Method declared outside class definition ---------
//---------------------------------------------------
template <class V, class M>
void
ComputeConditionalGaussianVectorRV(
  const V& muVec1,
  const V& muVec2,
  const M& sigmaMat11,
  const M& sigmaMat12,
  const M& sigmaMat21,
  const M& sigmaMat22,
  const V& sampleVec2,
        V& muVec1_cond_on_2,
        M& sigmaMat11_cond_on_2);

}  // End namespace QUESO

#endif // UQ_GAUSSIAN_VECTOR_RV_H
