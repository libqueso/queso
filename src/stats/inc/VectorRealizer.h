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

#ifndef UQ_REALIZER_H
#define UQ_REALIZER_H

#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorRealizer.h
 * \brief A templated class for sampling from vector RVs (holding probability density distributions).
 *
 * \class BaseVectorRealizer
 * \brief A templated (base) class for handling sampling from vector RVs.
 *
 * A realizer is an object that, simply put, contains a realization() operation
 * that returns a sample of a vector RV. This is the base class. QUESO also support
 * uniform, Gaussian, Beta, Gamma, Inverse Gamma and LogNormal realizers, as described
 * and implemented in the derived classes. */

template <class V = GslVector, class M = GslMatrix>
class BaseVectorRealizer {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer.*/
  BaseVectorRealizer(const char*                  prefix,
			    const VectorSet<V,M>& unifiedImageSet,
			    unsigned int                 subPeriod);

  //! Virtual destructor
  virtual ~BaseVectorRealizer();
  //@}

  //! @name Realization-related methods
  //@{
  //! Image set where the realizations lie.  Access to protected attribute m_unifiedImageSet.
  const   VectorSet<V,M>& unifiedImageSet()              const;

  //! Sub-period of the realization. Access to protected attribute m_subPeriod.
  unsigned int           subPeriod      ()              const;

  //! Performs a realization (sample) from a probability density function. See template specialization.
  virtual void                   realization    (V& nextValues) const = 0;
  //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
  const VectorSet<V,M>&  m_unifiedImageSet;
        unsigned int            m_subPeriod;
};

}  // End namespace QUESO

#endif // UQ_REALIZER_H
