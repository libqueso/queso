//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#ifndef UQ_JOINT_PROB_DENSITY_H
#define UQ_JOINT_PROB_DENSITY_H

#include <cmath>

#include <boost/math/special_functions.hpp> // for Boost isnan. Note parentheses are important in function call.

#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file JointPdf.h
 * \brief Classes to accommodate a probability density.
 *
 * \class BaseJointPdf
 * \brief A templated (base) class for handling joint PDFs.
 *
 * This class allows the mathematical definition of a Joint PDF, which is a scalar
 * function such as * \f$ \pi: B \subset R^n \rightarrow R \f$; ie a function of one
 * or more variables that has always one-dimensional range. QUESO currently supports
 * basic PDFs such as uniform and Gaussian and also more complex PDFs, such as the
 * ones coming from a Bayesian analysis. They are implemented in the derived classes
 * UniformJointPdf, GaussianJointPdf, and BayesianJointPdf,
 * respectively. The posterior PDF may be represented within QUESO by GenericJointPdf. */

template <class V = GslVector, class M = GslMatrix>
class BaseJointPdf : public BaseScalarFunction<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class, i.e. a scalar function, given a prefix and its domain.*/
  BaseJointPdf(const char*                  prefix,
		      const VectorSet<V,M>& domainSet);
  //! Destructor
  virtual ~BaseJointPdf();
  //@}

  //! @name Mathematical methods
  //@{
  //! Actual value of the PDF (scalar function).
  virtual double actualValue                    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;

  //! Logarithm of the value of the function.
  virtual double lnValue                        (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;

  //! Sets a value to be used in the normalization style (stored in the protected attribute m_normalizationStyle.)
  virtual void   setNormalizationStyle          (unsigned int value) const;

  //! Sets a logarithmic value to be used in the normalization factor (stored in the protected attribute m_normalizationStyle.)
  void   setLogOfNormalizationFactor    (double value) const;

  //! Computes the logarithm of the normalization factor. See template specialization.
  virtual double computeLogOfNormalizationFactor(unsigned int numSamples, bool m_logOfNormalizationFactor) const = 0;

  //const BaseScalarPdf<double>& component(unsigned int componentId) const;
  //@}
protected:
  //! Common method (to the derived classes) to compute the logarithm of the normalization factor.
  /*! The normalization factor is calculated by finding the max and min values of the domain set
   * and then drawing \c numSamples samples from a uniform distribution varying from \c min to
   * \c max. Such samples are averaged and the logarithmic value is assigned to protected attribute
   * m_logOfNormalizationFactor if the parameter \c m_logOfNormalizationFactor is true. */
  double commonComputeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;

  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;

  mutable unsigned int m_normalizationStyle;
  mutable double       m_logOfNormalizationFactor;

//std::vector<BaseScalarPdf<double>*> m_components; // FIXME: will need to be a parallel vector in case of a very large number of components
//BaseScalarPdf<double>               m_dummyComponent;
};

}  // End namespace QUESO

#endif // UQ_JOINT_PROB_DENSITY_H
