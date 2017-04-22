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

#ifndef UQ_SEQUENTIAL_REALIZER_H
#define UQ_SEQUENTIAL_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class SequentialVectorRealizer
 * \brief A class for handling sequential draws (sampling) from probability density distributions.
 *
 * This class handles sequential sampling (it returns the next value of the chain) from a
 * probability density distribution.*/

template <class V = GslVector, class M = GslMatrix>
class SequentialVectorRealizer : public BaseVectorRealizer<V,M> {
public:

  //!@name Constructor/Destructor methods
  //! Default constructor.
  SequentialVectorRealizer(const char*                           prefix,
                                  const BaseVectorSequence<V,M>& chain);

  //! Destructor.
  ~SequentialVectorRealizer();
  //@}

  //!@name Sampling-related methods
  //! Returns the unified mean vector; access to private attribute m_unifiedSampleExpVector.
  const V&   unifiedSampleExpVector()              const;

  //! Returns the unified variance vector; access to private attribute m_unifiedSampleVarVector.
  const V&   unifiedSampleVarVector()              const;

  //! Draws the next value from this chain (\c m_chain) and saves it in \c nextValues
  void realization           (V& nextValues) const;
  //@}

private:
  const   BaseVectorSequence<V,M>& m_chain;
  mutable unsigned int                    m_currentChainPos;
          V*                              m_unifiedSampleExpVector;
          V*                              m_unifiedSampleVarVector;

  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;
};

}  // End namespace QUESO

#endif // UQ_SEQUENTIAL_REALIZER_H
