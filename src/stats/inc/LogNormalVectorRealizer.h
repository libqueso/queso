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

#ifndef UQ_LOGNORMAL_REALIZER_H
#define UQ_LOGNORMAL_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// LogNormal class [R-10]
//*****************************************************
/*!
 * \class LogNormalVectorRealizer
 * \brief A class for handling sampling from a Log-Normal probability density distribution.
 *
 * This class handles sampling from a Log-Normal probability density distribution, of
 * mean and variance given.*/

template <class V = GslVector, class M = GslMatrix>
class LogNormalVectorRealizer : public BaseVectorRealizer<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object of the class, given a prefix and the image set, a  vector of
   * mean values, \c lawExpVector, and a lower triangular matrix resulting from Cholesky
   * decomposition of the covariance matrix, \c lowerCholLawCovMatrix.  */
  LogNormalVectorRealizer(const char*                  prefix,
                                 const VectorSet<V,M>& unifiedImageSet,
                                 const V&                     lawExpVector, // vector of mean values
                                 const M&                     lowerCholLawCovMatrix); // lower triangular matrix resulting from Cholesky decomposition of the covariance matrix

  //! Constructor
  /*! Constructs a new object of the class, given a prefix and the image set, a vector of
   * mean values, \c lawExpVector, and a set of two matrices and one vector resulting from the
   * Single Value Decomposition of the covariance matrix, \c matU, \c vecSsqrt and \c matVt.*/
  LogNormalVectorRealizer(const char*                  prefix,
                                 const VectorSet<V,M>& unifiedImageSet,
                                 const V&                     lawExpVector, // vector of mean values
                                 const M&                     matU,
                                 const V&                     vecSsqrt,
                                 const M&                     matVt);
  //! Destructor
  ~LogNormalVectorRealizer();
  //@}

  //! @name Realization-related methods
  //@{
  //! Access to the vector of mean values and private attribute:  m_unifiedLawExpVector.
  const V&   unifiedLawExpVector ()              const;

  //! Access to the vector of variance values and private attribute:  m_unifiedLawVarVector.
  const V&   unifiedLawVarVector ()              const;

  //! Draws a realization.
  /*! This function draws a realization of a LogNormal distribution and saves it in \c nextValues.*/
  void realization         (V& nextValues) const;
  //@}

private:
  V* m_unifiedLawExpVector;
  V* m_unifiedLawVarVector;
  M* m_lowerCholLawCovMatrix;
  M* m_matU;
  V* m_vecSsqrt;
  M* m_matVt;

  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;
};

}  // End namespace QUESO

#endif // UQ_LOGNORMAL_REALIZER_H
