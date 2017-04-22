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

#ifndef UQ_GAMMA_JOINT_PROB_DENSITY_H
#define UQ_GAMMA_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class GammaJointPdf
 * \brief A class for handling Gamma joint PDFs.
 *
 * This class allows the mathematical definition of a Gamma Joint PDF.
 */
template <class V = GslVector, class M = GslMatrix>
class GammaJointPdf : public BaseJointPdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor for a shape-scale parameterisation of a Gamma(a,b) pdf
  //! defined on \c domainSet.
  /*!
   * Constructs a new object of the class, given a prefix, the domain set, and
   * the parameters \c a and \c b of the Gamma PDF.
   *
   * The parameters \c a and \c b are the shape and scale parameters,
   * respectively.  That is, the pdf is
   *
   * \f[
   *   f(x) = C x^{a-1} \exp \left( x / b \right),
   * \f]
   *
   * where \c C is a normalising constant.
   *
   * Indeed, \c a and \c b are vectors and in multiple dimensions the pdf is
   * just the product of Gamma(a_i, b_i) pdfs distributions in each dimension.
   * That is, they are independent.
   *
   * Note: the parameters \c a and \c b correspond to the shape and scale
   * parameters \c k and \f$\theta\f$ in the
   * <a href="https://en.wikipedia.org/wiki/Gamma_distribution">Wikipedia entry for the Gamma distribution</a>.
   */
  GammaJointPdf(const char*                  prefix,
                       const VectorSet<V,M>& domainSet,
                       const V&                     a,
                       const V&                     b);
  //! Destructor
  ~GammaJointPdf();
  //@}

  //! @name Math methods
  //@{
  //! Actual value of the Gamma PDF.
  /*! This routine calls method lnValue() and returns the exponent of the returning value of such method.*/
  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Logarithm of the value of the Gamma PDF.
  /*! If the normalization style (m_normalizationStyle) is zero, then this routine calls a environment method
   * which handles basic PDFs, e.g. basicPdfs()->gammaPdfActualValue() and adds the log of the normalization
   * factor (m_logOfNormalizationFactor) to it; otherwise the method uses the formula: \f$ lnValue =
   * \sum[ (a_i-1)*log(domainVector_i) -domainVector_i/b_i + m_logOfNormalizationFactor \f$, where a and b
   * are the parameters of the Gamma PDF. */
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

  V m_a;
  V m_b;
};

}  // End namespace QUESO

#endif // UQ_GAMMA_JOINT_PROB_DENSITY_H
