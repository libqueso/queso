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

#ifndef UQ_LOGNORM_JOINT_PROB_DENSITY_H
#define UQ_LOGNORM_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// LogNormal probability density class [PDF-10]
//*****************************************************
/*!
 * \class LogNormalJointPdf
 * \brief A class for handling Log-Normal joint PDFs.
 *
 * This class allows the mathematical definition of a Log-Normal Joint PDF.*/

template <class V = GslVector, class M = GslMatrix>
class LogNormalJointPdf : public BaseJointPdf<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object of the class, given a prefix and the domain of the PDF,
   * a vector of mean values, \c lawExpVector, and a vector of covariance values
   * \c lawVarVector (an alternative representation for a diagonal covariance matrix).*/
  LogNormalJointPdf(const char*                  prefix,
                           const VectorSet<V,M>& domainSet,
                           const V&                     lawExpVector,
                           const V&                     lawVarVector);
  //! Destructor
 ~LogNormalJointPdf();
 //@}

   //! @name Math methods
  //@{
  //! Actual value of the Log-Normal PDF (scalar function).
  /*! This method calls lnValue() and applies the exponential to its result.*/
  double   actualValue (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Logarithm of the value of the Log-Normal PDF (scalar function).
   /*! The logarithm of the value of the Log-Normal density of diagonal covariance matrix (sigma^2) comes from the summation:
  * \f[ lnValue =- \sum_i \frac{1}{domainVector_i * \sqrt{2 \pi * lawVarVector_i}} exp(-\frac{(\ln( domainVector_i) - lawExpVector_i)^2}{2 lawVarVector_i}) \f] as long as
  * \f$ domainVector_i > 0 \f$ for all \f$ i \f$.*/
  double   lnValue     (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Mean value of the underlying random variable.
  virtual void   distributionMean (V & meanVector) const;

  //! Covariance matrix of the underlying random variable.
  virtual void   distributionVariance (M & covMatrix) const;

  //! Computes the logarithm of the normalization factor.
  /*! This routine calls BaseJointPdf::commonComputeLogOfNormalizationFactor().*/
  double   computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

  //! Access to the vector of mean values and private attribute:  m_lawExpVector.
  const V& lawExpVector() const;

  //! Access to the vector of variance values and private attribute:  m_lawVarVector.
  const V& lawVarVector() const;
//@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_normalizationStyle;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;
  V*   m_lawExpVector;
  V*   m_lawVarVector;
  bool m_diagonalCovMatrix;
};

}  // End namespace QUESO

#endif // UQ_LOGNORM_JOINT_PROB_DENSITY_H
