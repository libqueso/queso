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
 * This is the header file from gravity_qoi.C.
 */

#ifndef QUESO_EXAMPLE_GRAVITY_QOI_H
#define QUESO_EXAMPLE_GRAVITY_QOI_H

#include <queso/VectorFunction.h>
#include <queso/DistArray.h>

template<class P_V = QUESO::GslVector, class P_M = QUESO::GslMatrix,
         class Q_V = QUESO::GslVector, class Q_M = QUESO::GslMatrix>
class Qoi : public QUESO::BaseVectorFunction<P_V, P_M, Q_V, Q_M>
{
public:
  Qoi(const char * prefix, const QUESO::VectorSet<P_V, P_M> & domainSet,
      const QUESO::VectorSet<Q_V, Q_M> & imageSet);
  virtual ~Qoi();
  virtual void compute(const P_V & domainVector, const P_V * domainDirection,
      Q_V & imageVector, QUESO::DistArray<P_V *> * gradVectors,
      QUESO::DistArray<P_M *> * hessianMatrices,
      QUESO::DistArray<P_V *> * hessianEffects) const;

  void setAngle(double angle);
  void setInitialVelocity(double velocity);
  void setInitialHeight(double height);

private:
  double m_angle;
  double m_initialVelocity;
  double m_initialHeight;
};

#endif
