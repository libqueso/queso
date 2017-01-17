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

/*
 * Brief description of this file:
 *
 * This file contains the code for the user defined qoi routine.
 */

#include <cmath>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <gravity_qoi.h>

template<class P_V, class P_M, class Q_V, class Q_M>
Qoi<P_V, P_M, Q_V, Q_M>::Qoi(const char * prefix,
    const QUESO::VectorSet<P_V, P_M> & domainSet,
    const QUESO::VectorSet<Q_V, Q_M> & imageSet)
  : QUESO::BaseVectorFunction<P_V, P_M, Q_V, Q_M>(prefix, domainSet, imageSet),
    m_angle(M_PI / 4.0),
    m_initialVelocity(5.0),
    m_initialHeight(0.0)
{
}

template<class P_V, class P_M, class Q_V, class Q_M>
Qoi<P_V, P_M, Q_V, Q_M>::~Qoi()
{
  // Deconstruct here
}

template<class P_V, class P_M, class Q_V, class Q_M>
void
Qoi<P_V, P_M, Q_V, Q_M>::compute(const P_V & domainVector,
    const P_V * domainDirection,
    Q_V & imageVector, QUESO::DistArray<P_V *> * gradVectors,
    QUESO::DistArray<P_M *> * hessianMatrices,
    QUESO::DistArray<P_V *> * hessianEffects) const
{
  if (domainVector.sizeLocal() != 1) {
    queso_error_msg("domainVector does not have size 1");
  }
  if (imageVector.sizeLocal() != 1) {
    queso_error_msg("imageVector does not have size 1");
  }

  double g = domainVector[0];  // Sample of the RV 'gravity acceleration'
  double distanceTraveled = 0.0;
  double aux = m_initialVelocity * std::sin(m_angle);
  distanceTraveled = (m_initialVelocity * std::cos(m_angle) / g) *
    (aux + std::sqrt(std::pow(aux, 2) + 2.0 * g * m_initialHeight));

  imageVector[0] = distanceTraveled;
}

template class Qoi<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector,
                   QUESO::GslMatrix>;
