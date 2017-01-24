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

namespace QUESO {

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
   * The rest of \c domainVector contains the hyperparameters, if any.  For
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

  //! Actual value of the scalar function.
  virtual double actualValue(const V & domainVector, const V * /*domainDirection*/,
                             V * /*gradVector*/, M * /*hessianMatrix*/, V * /*hessianEffect*/) const
  { return std::exp(this->lnValue(domainVector)); }

protected:
  const V & m_observations;
};

}  // End namespace QUESO

#endif  // UQ_GAUSSIAN_LLHD_H
