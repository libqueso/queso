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

#ifndef UQ_BETA_JOINT_PROB_DENSITY_H
#define UQ_BETA_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class BetaJointPdf
 * \brief A class for handling Beta joint PDFs.
 *
 * This class allows the mathematical definition of a Beta Joint PDF.*/

template <class V = GslVector, class M = GslMatrix>
class BetaJointPdf : public BaseJointPdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object of the class, given a prefix, the domain set, and the parameters
   * \c alpha and \c beta of the Beta PDF.  */
  BetaJointPdf(const char*                  prefix,
                      const VectorSet<V,M>& domainSet,
                      const V&                     alpha,
                      const V&                     beta);
  //! Destructor
 ~BetaJointPdf();
 //@}

    //! @name Math methods
  //@{
  //! Actual value of the Beta PDF.
  /*! This routine calls method lnValue() and returns the exponent of the returning value of such method.*/
  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Logarithm of the value of the Beta PDF.
  /*! If the normalization style (m_normalizationStyle) is zero, then this routine calls a environment method
   * which handles basic PDFs, e.g. basicPdfs()->betaPdfActualValue() and adds the log of the normalization
   * factor (m_logOfNormalizationFactor) to it; otherwise the method uses the formula: \f$ lnValue =
   * \sum[ (alpha_i-1)*log(domainVector_i) + (beta_i-1)*log(1-domainVector_i)] + m_logOfNormalizationFactor \f$. */
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Mean value of the underlying random variable.
  virtual void   distributionMean (V & meanVector) const;

  //! Covariance matrix of the underlying random variable.
  virtual void   distributionVariance (M & covMatrix) const;

  //! Computes the logarithm of the normalization factor.
  /*! This routine calls BaseJointPdf::commonComputeLogOfNormalizationFactor().*/
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_normalizationStyle;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;

  V m_alpha;
  V m_beta;
};

}  // End namespace QUESO

#endif // UQ_BETA_JOINT_PROB_DENSITY_H
