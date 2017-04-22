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

#ifndef UQ_UNIFORM_REALIZER_H
#define UQ_UNIFORM_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class UniformVectorRealizer
 * \brief A class for handling sampling from a Uniform probability density distribution.
 *
 * This class handles sampling from a uniform probability density distribution.*/

template <class V = GslVector, class M = GslMatrix>
class UniformVectorRealizer : public BaseVectorRealizer<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*!
   * Constructs a new object, given a prefix and the image set of the vector
   * realizer.
   *
   * Note: If \c unifiedImageSet is unbounded, it does not make sense to ask
   * for realizations.
   */
  UniformVectorRealizer(const char*                  prefix,
                               const VectorSet<V,M>& unifiedImageSet);
  //! Destructor
  ~UniformVectorRealizer();
  //@}

  //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*!
   * This function draws a realization of a uniform distribution and saves it
   * in \c nextValues.  It internally finds the minimum and the maximum values
   * of the distribution.
   *
   * Note: If the state space is unbounded, the distribution is improper and
   * realizations do not make sense.
   */
  void realization(V& nextValues) const;
  //@}

private:
  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;
};

}  // End namespace QUESO

#endif // UQ_UNIFORM_REALIZER_H
