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

#ifndef UQ_GAUSSIAN_JOINT_PROB_DENSITY_H
#define UQ_GAUSSIAN_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Gaussian probability density class [PDF-03]
//*****************************************************
/*!
 * \class GaussianJointPdf
 * \brief A class for handling Gaussian joint PDFs.
 *
 * This class allows the mathematical definition of a Gaussian Joint PDF.*/

template <class V = GslVector, class M = GslMatrix>
class GaussianJointPdf : public BaseJointPdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix and the domain of the PDF, a vector of mean
   * values, \c lawExpVector, and a vector of covariance values \c lawVarVector (an
   * alternative representation for a diagonal covariance matrix).  */
  GaussianJointPdf(const char*                  prefix,
                          const VectorSet<V,M>& domainSet,
                          const V&                     lawExpVector,
                          const V&                     lawVarVector);
  //! Constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer, a
   * vector of mean values, \c lawExpVector, and a covariance matrix, \c lawCovMatrix. */
  GaussianJointPdf(const char*                  prefix,
                          const VectorSet<V,M>& domainSet,
                          const V&                     lawExpVector,
                          const M&                     lawCovMatrix);
  //! Destructor
 ~GaussianJointPdf();
 //@}

  //! @name Math methods
  //@{

  //! Actual value of the Gaussian PDF:
  /*! This method calls lnValue() and applies the exponential to it.*/
 double   actualValue       (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Logarithm of the value of the Gaussian PDF (scalar function).
 /*! The ln(value) comes from a summation of the Gaussian density:
  * \f[ lnValue =- \sum_i \frac{1}{\sqrt{|covMatrix|} \sqrt{2 \pi}} exp(-\frac{(domainVector_i - lawExpVector_i)* covMatrix^{-1}* (domainVector_i - lawExpVector_i) }{2},  \f]
  * where the \f$ covMatrix \f$ may recovered via \c this->lawVarVector(), in case of diagonal
  * matrices or via \c this->m_lawCovMatrix, otherwise.*/
  double   lnValue           (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Mean value of the underlying random variable.
  virtual void   distributionMean (V & meanVector) const;

  //! Covariance matrix of the underlying random variable.
  virtual void   distributionVariance (M & covMatrix) const;

  //! Computes the logarithm of the normalization factor.
  /*! This routine calls BaseJointPdf::commonComputeLogOfNormalizationFactor().*/
  double   computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

  //! Updates the mean with the new value \c newLawExpVector.
  /*! This method deletes old expected values (allocated at construction or last call to this method).*/
  void     updateLawExpVector(const V& newLawExpVector);

  //! Updates the lower triangular matrix from Cholesky decomposition of the covariance matrix to the new value \c newLowerCholLawCovMatrix.
  /*! This method deletes old expected values (allocated at construction or last call to this method).*/
  void     updateLawCovMatrix(const M& newLawCovMatrix);

  //! Returns the covariance matrix; access to protected attribute m_lawCovMatrix.
  const M& lawCovMatrix      () const;

  //! Access to the vector of mean values and private attribute:  m_lawExpVector.
  const V& lawExpVector() const;

  //! Access to the vector of variance values and private attribute:  m_lawVarVector.
  const V& lawVarVector() const;
  //@}

  //! Print method for informational and logging purposes
  virtual void print(std::ostream & os) const;

protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_normalizationStyle;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;
  V*       m_lawExpVector;
  V*       m_lawVarVector;
  bool     m_diagonalCovMatrix;
  const M* m_lawCovMatrix;
};

}  // End namespace QUESO

#endif // UQ_GAUSSIAN_JOINT_PROB_DENSITY_H
