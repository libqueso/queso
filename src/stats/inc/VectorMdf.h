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

#ifndef UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H
#define UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H

#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorMdf.h
 * \brief Classes to accommodate a marginal density function of a vector RV.
 *
 * \class BaseVectorMdf
 * \brief A templated (base) class for handling MDFs of vector functions.
 *
 * To obtain the marginal distribution over a subset of multivariate (vector) RVs, one only
 * needs to drop the irrelevant variables (the variables that one wants to marginalize out).
 * If \b X is a Random Vector which contains the continuous random variables \f$ X_1, X_2, ..., X_n \f$.
 * Then each \f$ X_i \f$ is a continuous random variable with its own Probability Distribution called
 * the marginal distribution of \f$ X_i \f$.
 * This class handles MDFs of a vector RV.*/

template <class V = GslVector, class M = GslMatrix>
class BaseVectorMdf {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix and domain set of the MDF.*/
  BaseVectorMdf(const char*                  prefix,
                                const VectorSet<V,M>& domainSet);
  //! Virtual destructor.
  virtual ~BaseVectorMdf();
  //@}

  //! @name Mathematical methods
  //@{
  //! Returns the domain set; access to protected attribute m_domainSet.
  const   VectorSet<V,M>& domainSet() const;

  //! Finds the value of the vector MDF at each element of \c paramValue, and saves it in \c mdfVec. See template specialization.
  virtual void                   values   (const V& paramValues,
                                                 V& mdfVec)  const = 0;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the vector MDF. See template specialization.
  virtual void                   print    (std::ostream& os) const = 0;
  //@}

protected:

  const   BaseEnvironment& m_env;
          std::string             m_prefix;
  const   VectorSet<V,M>&  m_domainSet;
};

}  // End namespace QUESO

#endif // UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H
