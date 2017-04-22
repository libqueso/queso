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

#ifndef UQ_WIGNER_REALIZER_H
#define UQ_WIGNER_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class WignerVectorRealizer
 * \brief A class for handling sampling from a Wigner probability density distribution.
 *
 * This class \b will handle sampling from an Wigner probability density distribution, with a
 * given center position and a radius.
 *
 * \todo: The method WignerVectorRealizer:realization() is not yet available,
 * thus this class does  nothing. */

template <class V = GslVector, class M = GslMatrix>
class WignerVectorRealizer : public BaseVectorRealizer<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, the
   * center position \c centerPos, and a radius \c radius.*/
  WignerVectorRealizer(const char*                  prefix,
                              const VectorSet<V,M>& unifiedImageSet,
                              const V&                     centerPos,
                              double                       radius);

  //! Destructor
 ~WignerVectorRealizer();
 //@}

  //! @name Realization-related methods
  //@{
  //! TODO: Draws a realization.
  /*! \todo: implement and explain me!*/
  void realization(V& nextValues) const;
  //@}

private:
  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;
  V*     m_centerPos;
  double m_radius;
};

}  // End namespace QUESO

#endif // UQ_WIGNER_REALIZER_H
