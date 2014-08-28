//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef UQ_INVLOGIT_GAUSSIAN_VECTOR_RV_H
#define UQ_INVLOGIT_GAUSSIAN_VECTOR_RV_H

#include <queso/VectorRV.h>
#include <queso/BoxSubset.h>

namespace QUESO {

/*!
 * \class InvLogitGaussianVectorRV
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

template<class V, class M>
class InvLogitGaussianVectorRV : public BaseVectorRV<V, M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor  
  /*! Construct a Gaussian vector RV with mean \c lawExpVector and diagonal covariance matrix
   * \c lawVarVector whose variates live in \c imageSet.*/
  InvLogitGaussianVectorRV(const char * prefix,
      const BoxSubset<V, M> & imageBoxSubset, const V & lawExpVector,
      const V & lawVarVector);
  
  //! Constructor  
  /*! Construct a Gaussian vector RV with mean \c lawExpVector and covariance matrix
   * \c lawCovMatrix whose variates live in \c imageSet.*/
  InvLogitGaussianVectorRV(const char * prefix,
      const BoxSubset<V, M> & imageBoxSubset, const V & lawExpVector,
      const M & lawCovMatrix);
  
  //! Virtual destructor
  virtual ~InvLogitGaussianVectorRV();
  //@}

  //! @name Statistical methods
  //@{
  //! Updates the vector that contains the mean values.
  void updateLawExpVector(const V & newLawExpVector);
  
  //! Updates the covariance matrix.
  /*! This method tries to use Cholesky decomposition; and if it fails, the method then 
   *  calls a SVD decomposition.*/
  void updateLawCovMatrix(const M & newLawCovMatrix);
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream & os) const;
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

}  // End namespace QUESO

#endif // UQ_INVLOGIT_GAUSSIAN_VECTOR_RV_H
