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

#ifndef UQ_GAUSSIAN_VECTOR_MARGINAL_DENSITY_FUNCTION_H
#define UQ_GAUSSIAN_VECTOR_MARGINAL_DENSITY_FUNCTION_H

#include <queso/VectorMdf.h>
#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Gaussian marginal density function class
//*****************************************************
/*!
 * \class GaussianVectorMdf
 * \brief TODO: A class for handling Gaussian MDFs.
 *
 * This class \b will implement a Gaussian vector marginal density function function (MDF).
 * \todo: Implement me! */

template <class V = GslVector, class M = GslMatrix>
class GaussianVectorMdf : public BaseVectorMdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! TODO: Constructor.
  /*! \todo: implement me! This method calls commonConstructor() which is not yet implemented.
   * Instantiates an object of the class given a prefix, the domain set, and the domain mean
   * and expected values. */
  GaussianVectorMdf(const char*                    prefix,
                           const VectorSet<V,M>& domainSet,
                           const V&                       domainExpectedValues,
                           const V&                       domainVarianceValues);
  //! TODO: Constructor.
  /*! \todo: implement me! This method calls commonConstructor() which is not yet implemented.
   * Instantiates an object of the class given a prefix, the domain set, and the domain mean
   * and covariance matrix. */
  GaussianVectorMdf(const char*                    prefix,
                           const VectorSet<V,M>& domainSet,
                           const V&                       domainExpectedValues,
                           const M&                       covMatrix);
  //! Destructor
  ~GaussianVectorMdf();
  //@}

  //! @name Mathematical method
  //@{
  //! TODO: Returns the values of the vector CDF at each element of \c paramValues.
  /*! \todo: implement me!*/
  void values(const V& paramValues, V& mdfVec) const;
  //@}

  //! @name I/O method
  //@{
  //! TODO: Prints the vector CDF.
  /*! \todo: implement me!*/
  void print (std::ostream& os)                const;
  //@}

protected:
  const M*                         m_covMatrix;

  using BaseVectorMdf<V,M>::m_env;
  using BaseVectorMdf<V,M>::m_prefix;
  using BaseVectorMdf<V,M>::m_domainSet;

  void commonConstructor();
};

}  // End namespace QUESO

#endif // UQ_GAUSSIAN_VECTOR_MARGINAL_DENSITY_FUNCTION_H
