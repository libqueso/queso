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

#ifndef UQ_CONCATENATED_REALIZER_H
#define UQ_CONCATENATED_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Concatenated class [R-11]
//*****************************************************
/*!
 * \class ConcatenatedVectorRealizer
 * \brief A class for handling sampling from concatenated probability density distributions.
 *
 * This class allows the user draw samples from concatenated probability density distributions (two
 * or more distinct probability distributions has(ve) been concatenated into one single vector RV).
 * This class used, for instance, to draw realization of concatenate priors from two or more RVs,
 * where one of them has a uniform distribution whereas the other one(s) has a Gaussian distribution. */

template <class V = GslVector, class M = GslMatrix>
class ConcatenatedVectorRealizer : public BaseVectorRealizer<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Concatenates two RVs: \c rv1 and \c rv2 into one vector RV, given a prefix and the image set of the vector RV.*/
  ConcatenatedVectorRealizer(const char*                           prefix,
                                    const BaseVectorRealizer<V,M>& realizer1,
                                    const BaseVectorRealizer<V,M>& realizer2,
                                    const VectorSet<V,M>&          unifiedImageSet);
  //! Constructor
  /*! Concatenates a sequence of RVs, given by: <c> std::vector<const BaseVectorRV<V,M>* >& rvs </c>
   * into one single vector RV, given a prefix and the image set of the resulting vector RV.*/
  ConcatenatedVectorRealizer(const char*                                                prefix,
                                    const std::vector<const BaseVectorRealizer<V,M>* >& realizers,
                                    unsigned int                                               minPeriod,
                                    const VectorSet<V,M>&                               unifiedImageSet);
  //! Destructor
  ~ConcatenatedVectorRealizer();
  //@}

  //! @name Realization-related methods
  //@{
  void realization(V& nextValues) const;
  //@}

private:
  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;

  std::vector<const BaseVectorRealizer<V,M>* > m_realizers;
};

}  // End namespace QUESO

#endif // UQ_CONCATENATED_REALIZER_H
