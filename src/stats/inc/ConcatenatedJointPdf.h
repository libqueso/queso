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

#ifndef UQ_CONCAT_JOINT_PROB_DENSITY_H
#define UQ_CONCAT_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Concatenated probability density class [PDF-11]
//*****************************************************
/*!
 * \class ConcatenatedJointPdf
 * \brief A class for handling concatenated PDFs.
 *
 * This class allows the user to defines concatenated probability density distributions,
 * i.e, two or more distinct PDFs can be concatenated into one single PDF.
 * This class used, for instance, to concatenate priors from two or more RVs, where one
 * of them has a uniform distribution whereas the other one(s) has a Gaussian distribution. */

template <class V = GslVector, class M = GslMatrix>
class ConcatenatedJointPdf : public BaseJointPdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Concatenates two PDFs: \c density1 and \c density2 into one vector PDF, given a prefix
   * and the concatenated domain of such PDFs.*/
  ConcatenatedJointPdf(const char*                     prefix,
                              const BaseJointPdf<V,M>& density1,
                              const BaseJointPdf<V,M>& density2,
                              const VectorSet   <V,M>& concatenatedDomain);

  //! Constructor
  /*! Concatenates a sequence of PDFs, given by: <c> std::vector<const BaseJointPdf<V,M>* >& densities </c>
   * into one single PDF, given a prefix and the concatenated domain of such PDFs.*/
  ConcatenatedJointPdf(const char*                                          prefix,
                              const std::vector<const BaseJointPdf<V,M>* >& densities,
                              const VectorSet<V,M>&                         concatenatedDomain);

  //! Destructor
  ~ConcatenatedJointPdf();
  //@}

  //! @name Math methods
  //@{
  //! Calculates the actual values of each density.
  /*! The final actual value is the multiplication of all values calculated.*/
  double actualValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Calculates the logarithm of the values of each density.
  /*! The final logarithm value is the addition of all values calculated.*/
  double lnValue              (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Mean value of the underlying random variable.
  virtual void   distributionMean (V & meanVector) const;

  //! Covariance matrix of the underlying random variable.
  virtual void   distributionVariance (M & covMatrix) const;

  //! Sets the normalization style of all densities to \c value.
  void   setNormalizationStyle(unsigned int value) const;

  //! Computes the logarithm of the normalization factor.
  /*! This method calls the computeLogOfNormalizationFactor() for each one of the densities that have
   * been concatenated.*/
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  //@}

protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;

  std::vector<const BaseJointPdf<V,M>* > m_densities;
};

}  // End namespace QUESO

#endif // UQ_CONCAT_JOINT_PROB_DENSITY_H
