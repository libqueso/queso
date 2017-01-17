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

#ifndef UQ_GENERIC_REALIZER_H
#define UQ_GENERIC_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Generic class [R-01]
//*****************************************************
/*!
 * \class GenericVectorRealizer
 * \brief A class for handling sampling from generic probability density distributions.
 *
 * A realizer is an object that, simply put, contains a realization() operation
 * that returns a sample of a vector RV, or, particularly, a generic probability
 * density distribution. This is the class that handles generic sampling, used,
 * for example, to sample, posterior PDFs (the solution of a Bayesian problem).*/

template <class V = GslVector, class M = GslMatrix>
class GenericVectorRealizer : public BaseVectorRealizer<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer,
   * the sub period for the realizations and a pointer to a  generic routine. */
  GenericVectorRealizer(const char*                  prefix,
                               const VectorSet<V,M>& unifiedImageSet,
                               unsigned int                 subPeriod,
                               double (*routinePtr)(const void* routineDataPtr, V& nextParamValues),
                               const void* routineDataPtr);
 //! Destructor
  ~GenericVectorRealizer();
  //@}

  //! @name Realization-related methods
  //! Draws a realization.
  /*! This function draws a realization of \c this considering the generic routine \c m_routinePtr
   * and saves it in \c nextValues.*/
  void realization(V& nextValues) const;
  //@}

private:
  double (*m_routinePtr)(const void* routineDataPtr, V& nextParamValues);
  const void* m_routineDataPtr;

  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_REALIZER_H
