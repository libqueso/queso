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

#ifndef UQ_GAUSSIAN_LLHD_H
#define UQ_GAUSSIAN_LLHD_H

#include <vector>
#include <cmath>
#include <queso/ScalarFunction.h>
#include <queso/SharedPtr.h>

namespace QUESO {

template<class V, class M>
class BaseVectorRV;

template<class V, class M>
class MultiDQuadratureBase;

class GslVector;
class GslMatrix;

/*!
 * \file LikelihoodBase.h
 *
 * \class LikelihoodBase
 * \brief Base class for canned Gaussian likelihoods
 *
 * This class is an abstract base class for 'canned' Gaussian likelihoods.  All
 * this class does is add a pure virtual function called \c evaluateModel that
 * the user will implement to interact with the forward code.
 */

template <class V = GslVector, class M = GslMatrix>
class LikelihoodBase : public BaseScalarFunction<V, M> {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default constructor.
  /*!
   * The vector of observations must be passed.  This will be used when
   * evaluating the likelihood functional
   */
  LikelihoodBase(const char * prefix,
      const VectorSet<V, M> & domainSet,
      const V & observations);

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
  LikelihoodBase(const char * prefix,
                 const VectorSet<V, M> & domainSet,
                 const V & observations,
                 typename SharedPtr<BaseVectorRV<V,M> >::Type & marg_param_pdf,
                 typename SharedPtr<MultiDQuadratureBase<V,M> >::Type & marg_integration,
                 bool marg_pdf_is_weight_func);

  //! Destructor, pure to make this class abstract
  virtual ~LikelihoodBase() =0;
  
  //@}

  //! Deprecated. Evaluates the user's model at the point \c domainVector
  /*!
   * This method is deprecated. The user instead should subclass
   * evaluateModel(const V & domainVector, V & modelOutput).
   */
  virtual void evaluateModel(const V & domainVector, const V * domainDirection,
                             V & modelOutput, V * gradVector, M * hessianMatrix,
                             V * hessianEffect) const;

  //! Evaluates the user's model at the point \c domainVector
  /*!
   * Subclass implementations will fill up the \c modelOutput vector with output
   * from the model. This represents a vector of synthetic observations that will
   * be compared to actual observations when computing the likelihood functional.
   *
   * The first \c n components of \c domainVector are the model parameters.
   * The rest of \c domainVector contains the hyperparameters, if any. Mathematically, this
   * function is \f$ f(m) \f$. For
   * example, in \c GaussianLikelihoodFullCovarainceRandomCoefficient, the last
   * component of \c domainVector contains the multiplicative coefficient of
   * the observational covariance matrix.  In this case, the user need not
   * concern themselves with this parameter as it is handled not in the model
   * evaluation but by the likelihood evaluation.
   *
   * By default, to ensure backward compatibility, this method will call the
   * deprecated
   * evaluateModel( const V & domainVector, const V * domainDirection,
   * V & modelOutput, V * gradVector, M * hessianMatrix, V * hessianEffect ).
   */
  virtual void evaluateModel(const V & domainVector, V & modelOutput) const
  { this->evaluateModel(domainVector,NULL,modelOutput,NULL,NULL,NULL); }

  //! Evaluates the user's model at the point \c (domainVector,marginalVector)
  /*!
   * Subclass implementations will fill up the \c modelOutput vector with output
   * from the model. This function will be passed the current value of the parameters in domainVector
   * and the current values of the marginalization parameters in marginalVector. Mathematically, this
   * function is \f$ f(m,q) \f$. This function
   * is only called if the likelihood includes marginalization. Otherwise,
   * evaluateModel(const V & domainVector, V & modelOutput) will be called by lnValue().
   *
   * Not every user will do marginalization, so this is not a pure function, but there is
   * no default, so we error if marginalization is attempted without overriding this function.
   */
  virtual void evaluateModel(const V & domainVector, const V & marginalVector, V & modelOutput) const;

  //! Actual value of the scalar function.
  virtual double actualValue(const V & domainVector, const V * /*domainDirection*/,
                             V * /*gradVector*/, M * /*hessianMatrix*/, V * /*hessianEffect*/) const
  { return std::exp(this->lnValue(domainVector)); }

  //! Logarithm of the value of the scalar function.
  /*!
   * This method well evaluate the users model (evaluateModel()) and then
   * call lnLikelihood() and, if there is marginalization, will handle the
   * numerical integration over the marginal parameter space.
   */
  virtual double lnValue(const V & domainVector) const;

protected:

  const V & m_observations;

  typename SharedPtr<const BaseVectorRV<V,M> >::Type m_marg_param_pdf;

  typename SharedPtr<MultiDQuadratureBase<V,M> >::Type m_marg_integration;

  bool m_marg_pdf_is_weight_func;

  //! Compute log-likelihood value given the current parameters and model output
  /*!
   * Subclasses should override this method to compute the likelihood value, given
   * the current parameters in domainVector and the output from the model evaluated
   * at those parameter values. Note that this function may alter the values of
   * modelOutput to reduce copying.
   */
  virtual double lnLikelihood(const V & domainVector, V & modelOutput) const =0;
};

}  // End namespace QUESO

#endif  // UQ_GAUSSIAN_LLHD_H
