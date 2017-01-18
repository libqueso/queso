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

#ifndef UQ_GAUSSIAN_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
#define UQ_GAUSSIAN_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H

#include <queso/VectorCdf.h>
#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/ScalarCdf.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Gaussian cumulative distribution function class
//*****************************************************
/*!
 * \class GaussianVectorCdf
 * \brief TODO: A class for handling Gaussian CDFs.
 *
 * This class \b will implement a Gaussian vector cumulative distribution function (CDF).
 * \todo: Implement me! */

template <class V = GslVector, class M = GslMatrix>
class GaussianVectorCdf : public BaseVectorCdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! TODO: Constructor.
  /*! \todo: implement me! This method calls commonConstructor() which is not yet implemented.
   *  Instantiates an object of the class given a prefix, the support of the related-PDF, and
   * the domain mean and expected values.
   */
  GaussianVectorCdf(const char*                  prefix,
                           const VectorSet<V,M>& pdfSupport,
                           const V&                     domainExpectedValues,
                           const V&                     domainVarianceValues);
  //! TODO: Constructor.
  /*! \todo: implement me! This method calls commonConstructor() which is not yet implemented.
   * Instantiates an object of the class given a prefix, the support of the related-PDF, and
   * the domain mean values and covariance matrix.*/
  GaussianVectorCdf(const char*                  prefix,
                           const VectorSet<V,M>& pdfSupport,
                           const V&                     domainExpectedValues,
                           const M&                     covMatrix);
  // Destructor
  ~GaussianVectorCdf();
  //@}

  //! @name Mathematical method
  //@{
  //! TODO: Returns the values of the vector CDF at each element of \c paramValues.
  /*! \todo: implement me!*/
  void values(const V& paramValues, V& cdfVec) const;
  //@}

  //! @name I/O method
  //@{
  //! TODO: Prints the vector CDF.
  /*! \todo: implement me!*/
  void print (std::ostream& os) const;
  //@}

protected:
  const M*                         m_covMatrix;

  using BaseVectorCdf<V,M>::m_env;
  using BaseVectorCdf<V,M>::m_prefix;
  using BaseVectorCdf<V,M>::m_pdfSupport;

  //! A common constructor to be used both class constructors.
  void commonConstructor();
};

}  // End namespace QUESO

#endif // UQ_GAUSSIAN_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
