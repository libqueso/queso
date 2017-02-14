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

#ifndef UQ_GAUSSIAN_LIKELIHOOD_SCALAR_COV_H
#define UQ_GAUSSIAN_LIKELIHOOD_SCALAR_COV_H

#include <queso/LikelihoodBase.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \file GaussianLikelihoodScalarCovariance.h
 *
 * \class GaussianLikelihoodScalarCovariance
 * \brief A class that represents a Gaussian likelihood with scalar covariance
 */

template <class V = GslVector, class M = GslMatrix>
class GaussianLikelihoodScalarCovariance : public LikelihoodBase<V, M> {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default constructor.
  /*!
   * Instantiates a Gaussian likelihood function, given a prefix, its domain,
   * a set of observations and a scalar covariance matrix.  The scalar
   * 'covariance matrix' is just passed as a \c double.
   */
  GaussianLikelihoodScalarCovariance(const char * prefix,
      const VectorSet<V, M> & domainSet, const V & observations,
      double covariance);

  //! Constructor for likelihood that includes marginalization
  /*!
   * If the likelihood requires marginalization, the user can provide the pdf of the
   * marginal parameter(s) and the integration to be used. Additionally, the user will
   * be required to have provided an implementation of
   * evaluateModel(const V & domainVector, const V & marginalVector, V & modelOutput).
   *
   * Mathematically, this likelihood evaluation will be
   * \f[ \pi(d|m) = \int \pi(d|m,q) \pi(q)\; dq \approx \sum_{i=1}^{N} \pi(d|m,q_i) \pi(q_i) w_i\f]
   * where \f$ N \f$ is the number of quadrature points. However, the PDF for the
   * marginal parameter(s) may be such that it is convenient to interpret it as a
   * weighting function for Gaussian quadrature. In that case, then,
   * \f[ \int \pi(d|m,q) \pi(q)\; dq \approx \sum_{i=1}^{N} \pi(d|m,q_i) w_i \f]
   * If this is the case, the user should set the argument marg_pdf_is_weight_func = true.
   * If it is set to false, then the former quadrature equation will be used.
   * For example, if the marginal parameter(s) pdf is Gaussian, a Gauss-Hermite quadrature
   * rule could make sense (GaussianHermite1DQuadrature).
   */
  GaussianLikelihoodScalarCovariance(const char * prefix,
      const VectorSet<V, M> & domainSet, const V & observations,
      double covariance,
      typename SharedPtr<BaseVectorRV<V,M> >::Type & marg_param_pdf,
      typename SharedPtr<MultiDQuadratureBase<V,M> >::Type & marg_integration,
      bool marg_pdf_is_weight_func);

  //! Destructor
  virtual ~GaussianLikelihoodScalarCovariance();
  //@}

protected:

  //! Logarithm of the value of the likelihood function.
  virtual double lnLikelihood(const V & domainVector, V & modelOutput) const;

private:
  double m_covariance;
};

}  // End namespace QUESO

#endif  // UQ_GAUSSIAN_LIKELIHOOD_SCALAR_COV_H
