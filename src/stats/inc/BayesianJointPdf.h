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

#ifndef UQ_BAYESIAN_JOINT_PROB_DENSITY_H
#define UQ_BAYESIAN_JOINT_PROB_DENSITY_H

#include <cmath>
#include <queso/math_macros.h>
#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class BayesianJointPdf
 * \brief A class for handling Bayesian joint PDFs.
 *
 * This class allows the mathematical definition of a Bayesian Joint PDF.*/

template <class V = GslVector, class M = GslMatrix>
class BayesianJointPdf : public BaseJointPdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of this class given a prefix and a scalar function.
   * The domain of the scalar function is assigned to the protected attribute m_domainSet,
   * and the scalar fiction is also itself copied to the protected attribute m_scalarFunction.*/
  BayesianJointPdf(const char*                           prefix,
                          const BaseJointPdf      <V,M>& priorDensity,
                          const BaseScalarFunction<V,M>& likelihoodFunction,
                          double                                likelihoodExponent,
                          const VectorSet         <V,M>& intersectionDomain);
  //! Destructor
  ~BayesianJointPdf();


  //! @name Math methods
  //@{
  //! Actual value of the PDF (scalar function).
  /*! If the exponent of the likelihood function (likelihoodExponent) is zero, i.e. the likelihood is
   * constant and unitary, then the actual value is the value of the prior PDF; otherwise, the actual
   * value is scaled (multiplied) by a power of the value of the likelihood function.*/
  double actualValue              (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Computes the logarithm of the value of the function.
  /*! Analogously to the method actualValue(), if the exponent of the likelihood function
   * (likelihoodExponent) is zero then the Logarithm of the value of the function is the logarithm of
   * the value of the prior PDF; otherwise, the value is scaled (added) by a power of the value of the
   * likelihood function.*/
  virtual double lnValue(const V & domainVector) const;
  virtual double lnValue(const V & domainVector, V & gradVector) const;

  //! Mean value of the underlying random variable.
  virtual void   distributionMean (V & meanVector) const { queso_not_implemented(); }

  //! Covariance matrix of the underlying random variable.
  virtual void   distributionVariance (M & covMatrix) const { queso_not_implemented(); };

  //! TODO: Computes the logarithm of the normalization factor.
  /*! \todo: implement me!*/
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

  //! Sets a value to be used in the normalization style of the prior density PDF (ie, protected attribute m_priorDensity).
  void   setNormalizationStyle    (unsigned int value) const;

  //! Returns the logarithm of the last computed Prior value. Access to protected attribute m_lastComputedLogPrior.
  double lastComputedLogPrior     () const;

  //! Returns the logarithm of the last computed likelihood value.  Access to protected attribute m_lastComputedLogLikelihood.
  double lastComputedLogLikelihood() const;

  //@}

  //! @name Using declarations
  //@{
  using BaseScalarFunction<V, M>::lnValue;
  //@}

protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;

  const BaseJointPdf      <V,M>& m_priorDensity;
  const BaseScalarFunction<V,M>& m_likelihoodFunction;
  double                                m_likelihoodExponent;
  mutable double                        m_lastComputedLogPrior;
  mutable double                        m_lastComputedLogLikelihood;

  mutable V  m_tmpVector1;
  mutable V  m_tmpVector2;
  mutable M* m_tmpMatrix;
};

}  // End namespace QUESO

#endif // UQ_BAYESIAN_JOINT_PROB_DENSITY_H
