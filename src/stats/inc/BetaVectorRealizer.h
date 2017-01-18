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

#ifndef UQ_BETA_REALIZER_H
#define UQ_BETA_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Beta class [R-05]
//*****************************************************
/*!
 * \class BetaVectorRealizer
 * \brief A class for handling sampling from a Beta probability density distribution.
 *
 * This class handles sampling from a Beta probability density distribution, of
 * parameters \c alpha and \c beta.*/

template <class V = GslVector, class M = GslMatrix>
class BetaVectorRealizer : public BaseVectorRealizer<V,M> {
public:

   //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, and the
   * Beta distribution parameters \c alpha and \c beta, which are assigned to private attributes
   * m_alpha and m_beta.  */
  BetaVectorRealizer(const char*                  prefix,
                            const VectorSet<V,M>& unifiedImageSet,
                            const V&                     alpha,
                            const V&                     beta);

  //! Destructor
  ~BetaVectorRealizer();
  //@}

    //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of a Beta distribution and saves it in \c nextValues.
   * It internally checks whether the image set, where the realization should be drawn, belongs
   * to the interval (0, 1] - which is the range where Beta distribution is defined over. */
  void realization(V& nextValues) const;
  //@}
private:
  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;

  V m_alpha;
  V m_beta;
};

}  // End namespace QUESO

#endif // UQ_BETA_REALIZER_H
